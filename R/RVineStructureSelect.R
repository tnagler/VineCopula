RVineStructureSelect <- function(data, familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA, rotations = TRUE) {
    d <- ncol(data)
    n <- nrow(data)
    
    ## sanity checks
    if (type == 0)
        type <- "RVine" else if (type == 1)
            type <- "CVine"
    if (type != "RVine" & type != "CVine")
        stop("Vine model not implemented.")
    if (n < 2)
        stop("Number of observations has to be at least 2.")
    if (d < 3)
        stop("Dimension has to be at least 3.")
    if (any(data > 1) || any(data < 0))
        stop("Data has to be in the interval [0,1].")
    if (!is.na(familyset[1])) {
        for (i in 1:length(familyset)) {
            if (!(familyset[i] %in% c(0, 1:10, 13, 14, 16:20,
                                      23, 24, 26:30, 33, 34, 36:40,
                                      104, 114, 124, 134, 204, 214, 224, 234)))
                stop("Copula family not implemented.")
        }
    }
    if (selectioncrit != "AIC" && selectioncrit != "BIC")
        stop("Selection criterion not implemented.")
    if (level < 0 & level > 1)
        stop("Significance level has to be between 0 and 1.")
    
    ## set variable names and trunclevel if not provided
    if (is.null(colnames(data)))
        colnames(data) <- paste0("V", 1:d)
    if (is.na(trunclevel))
        trunclevel <- d
    
    ## adjust familyset
    if (trunclevel == 0)
        familyset <- 0
    if (rotations)
        familyset <- with_rotations(familyset)
    
    ## initialize object for results
    RVine <- list(Tree = NULL, Graph = NULL)
    
    ## estimation in first tree ----------------------------
    # find optimal tree
    g <- initializeFirstGraph2(data, weights)
    MST <- findMaximumTauTree2(g, mode = type)
    
    # estimate pair-copulas
    VineTree <- fit.FirstTreeCopulas2(MST,
                                      data,
                                      familyset,
                                      selectioncrit,
                                      indeptest,
                                      level,
                                      weights = weights)
    # store results
    RVine$Tree[[1]] <- VineTree
    RVine$Graph[[1]] <- g
    oldVineGraph <- VineTree
    
    ## estimation in higher trees --------------------------
    for (i in 2:(d - 1)) {
        # only estimate pair-copulas if not truncated
        if (trunclevel == i - 1)
            familyset <- 0
        # find optimal tree
        g <- buildNextGraph2(VineTree, weights)
        MST <- findMaximumTauTree2(g, mode = type)
        # estimate pair-copulas
        VineTree <- fit.TreeCopulas2(MST,
                                     VineTree,
                                     familyset,
                                     selectioncrit,
                                     indeptest,
                                     level,
                                     progress,
                                     weights = weights)
        # store results
        RVine$Tree[[i]] <- VineTree
        RVine$Graph[[i]] <- g
    }
    
    ## free memory and return results as 'RVineMatrix' object
    .RVine <- RVine
    rm(list = ls())
    as.RVM2(.RVine)
}

initializeFirstGraph2 <- function(data.univ, weights) {
    ## calculate Kendall's tau
    taus <- TauMatrix(data = data.univ, weights = weights)
    
    ## return full graph with tau as weights
    graphFromTauMatrix(taus)
}

findMaximumTauTree2 <- function(g, mode = "RVine") {
    ## construct adjency matrix
    A <- adjacencyMatrix(g)
    d <- ncol(A)
    
    if (mode == "RVine") {
        ## initialize
        tree <- NULL
        edges <- matrix(NA, d - 1, 2)
        w <- numeric(d - 1)
        i <- 1
        
        ## construct minimum spanning tree
        for (k in 1:(d - 1)) {
            # add selected edge to tree
            tree <- c(tree, i)
            
            # find edge with minimal weight
            m <- apply(as.matrix(A[, tree]), 2, min)
            a <- apply(as.matrix(A[, tree]), 2, 
                       function(x) order(rank(x)))[1, ]
            b <- order(rank(m))[1]
            j <- tree[b]
            i <- a[b]
            
            # store edge and weight
            edges[k, ] <- c(j, i)
            w[k] <- A[i, j]
            
            ## adjust adjecency matrix to prevent loops
            for (t in tree)
                A[i, t] <- A[t, i] <- Inf
        }
        
        ## reorder edges for backwads compatibility with igraph output
        edges <- t(apply(edges, 1, function(x) sort(x)))
        edges <- edges[order(edges[, 2], edges[, 1]), ]
        
        ## delete unused edges from graph
        E <- g$E$nums
        in.tree <- apply(matrix(edges, ncol = 2), 1,
                         function(x) which((x[1] == E[, 1]) & (x[2] == E[, 2])))
        if (is.list(in.tree))
            do.call()
        MST <- g
        g$E$todel <- rep(TRUE, nrow(E))
        if (any(g$E$todel)) {
            g$E$todel[in.tree] <- FALSE
            MST <- deleteEdges(g)
        }
    } else if (mode  == "CVine") {
        ## set root as vertex with minimal sum of weights
        A <- adjacencyMatrix(g)
        diag(A) <- 0
        sumtaus <- rowSums(A)
        root <- which.min(sumtaus)
        
        ## delete unused edges
        g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 1] == root))
        MST <- g
        if (any(g$E$todel ))
            MST <- deleteEdges(g)
    } else {
        stop("vine not implemented")
    }
    
    ## return result
    MST
}


# not required any longer? Use TauMatrix instead
fasttau <- function(x, y, weights = NA) {
    if (any(is.na(weights))) {
        m <- length(x)
        n <- length(y)
        if (m == 0 || n == 0)
            stop("both 'x' and 'y' must be non-empty")
        if (m != n)
            stop("'x' and 'y' must have the same length")
        out <- .C("ktau",
                  x = as.double(x),
                  y = as.double(y),
                  N = as.integer(n),
                  tau = as.double(0),
                  S = as.double(0),
                  D = as.double(0),
                  T = as.integer(0),
                  U = as.integer(0),
                  V = as.integer(0),
                  PACKAGE = "VineCopula")
        ktau <- out$tau
    } else {
        ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
    }
    return(ktau)
}

## fit pair-copulas for the first vine tree
fit.FirstTreeCopulas2 <- function(MST, data.univ, type, copulaSelectionBy, testForIndependence, testForIndependence.level, weights = NA) {
    
    ## initialize estimation results with empty list
    d <- nrow(MST$E$nums)
    parameterForACopula <- lapply(1:d, function(i) NULL)
    
    ## prepare for estimation and store names
    for (i in 1:d) {
        ## get edge and corresponding data
        a <- MST$E$nums[i, ]
        parameterForACopula[[i]]$zr1 <- data.univ[, a[1]]
        parameterForACopula[[i]]$zr2 <- data.univ[, a[2]]
        MST$E$Copula.Data.1[i] <- list(data.univ[, a[1]])
        MST$E$Copula.Data.2[i] <- list(data.univ[, a[2]])
        
        ## set names for this edge
        if (is.null(MST$V$names[a[1]])) {
            MST$E$Copula.CondName.1[i] <- a[1]
        } else {
            MST$E$Copula.CondName.1[i] <- MST$V$names[a[1]]
        }
        if (is.null(MST$V$names[a[2]])) {
            MST$E$Copula.CondName.2[i] <- a[2]
        } else {
            MST$E$Copula.CondName.2[i] <- MST$V$names[a[2]]
        }
        if (is.null(MST$V$names[a[1]]) || is.null(MST$V$names[a[2]])) {
            MST$E$Copula.Name[i] <- paste(a[1], a[2], sep = " , ")
        } else {
            MST$E$Copula.Name[i] <- paste(MST$V$names[a[1]],
                                          MST$V$names[a[2]],
                                          sep = " , ")
        }
    }
    
    ## estimate parameters and select family
    outForACopula <- lapply(X = parameterForACopula,
                            FUN = wrapper_fit.ACopula,
                            type,
                            copulaSelectionBy,
                            testForIndependence,
                            testForIndependence.level,
                            weights)
    
    ## store estimated model and pseudo-obversations for next tree
    for (i in 1:d) {
        MST$E$Copula.param[[i]] <- c(outForACopula[[i]]$par,
                                     outForACopula[[i]]$par2)
        MST$E$Copula.type[i] <- outForACopula[[i]]$family
        MST$E$Copula.out[i] <- list(outForACopula[[i]])
        
        MST$E$Copula.CondData.1[i] <- list(outForACopula[[i]]$CondOn.1)
        MST$E$Copula.CondData.2[i] <- list(outForACopula[[i]]$CondOn.2)
    }
    
    ## return results
    MST
}

## fit pair-copulas for vine trees 2,...
fit.TreeCopulas2 <- function(MST, oldVineGraph, type, copulaSelectionBy, testForIndependence, testForIndependence.level, progress, weights = NA) {
    
    ## initialize estimation results with empty list
    d <- nrow(MST$E$nums)
    parameterForACopula <- lapply(1:d, function(i) NULL)
    
    
    ## prepare for estimation
    for (i in 1:d) {
        ## get edge and corresponding data
        con <- MST$E$nums[i, ]
        temp <- oldVineGraph$E$nums[con, ]
        
        ## fetch corresponding data and names
        if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
            same <- temp[2, 1]
        } else {
            if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
                same <- temp[2, 2]
            }
        }
        
        if (temp[1, 1] == same) {
            zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
            n1  <- oldVineGraph$E$Copula.CondName.2[con[1]]
        } else {
            zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
            n1  <- oldVineGraph$E$Copula.CondName.1[con[1]]
        }
        if (temp[2, 1] == same) {
            zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
            n2  <- oldVineGraph$E$Copula.CondName.2[con[2]]
        } else {
            zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
            n2  <- oldVineGraph$E$Copula.CondName.1[con[2]]
        }
        
        if (is.list(zr1)) {
            zr1a <- as.vector(zr1[[1]])
            zr2a <- as.vector(zr2[[1]])
            n1a <- as.vector(n1[[1]])
            n2a <- as.vector(n2[[1]])
        } else {
            zr1a <- zr1
            zr2a <- zr2
            n1a <- n1
            n2a <- n2
        }
        
        if (progress == TRUE)
            message(n1a, " + ", n2a, " --> ", MST$E$names[i])
        
        parameterForACopula[[i]]$zr1 <- zr1a
        parameterForACopula[[i]]$zr2 <- zr2a
        
        MST$E$Copula.Data.1[i] <- list(zr1a)
        MST$E$Copula.Data.2[i] <- list(zr2a)
        
        MST$E$Copula.CondName.1[i] <- n1a
        MST$E$Copula.CondName.2[i] <- n2a
    }
    
    ## estimate parameters and select family
    outForACopula <- lapply(X = parameterForACopula,
                            FUN = wrapper_fit.ACopula,
                            type,
                            copulaSelectionBy,
                            testForIndependence,
                            testForIndependence.level,
                            weights)
    
    ## store estimated model and pseudo-obversations for next tree
    for (i in 1:d) {
        MST$E$Copula.param[[i]] <- c(outForACopula[[i]]$par,
                                     outForACopula[[i]]$par2)
        MST$E$Copula.type[i] <- outForACopula[[i]]$family
        MST$E$Copula.out[i] <- list(outForACopula[[i]])
        
        MST$E$Copula.CondData.1[i] <- list(outForACopula[[i]]$CondOn.1)
        MST$E$Copula.CondData.2[i] <- list(outForACopula[[i]]$CondOn.2)
    }
    
    ## return results
    MST
}

## initialize graph for next vine tree (possible edges)
buildNextGraph2 <- function(oldVineGraph, weights = NA) {
    
    d <- nrow(oldVineGraph$E$nums)
    
    ## initialize with full graph
    g <- makeFullGraph(d)
    g$V$names <- oldVineGraph$E$names
    g$V$conditionedSet <- oldVineGraph$E$conditionedSet
    g$V$conditioningSet <- oldVineGraph$E$conditioningSet
    
    ## get info for all edges
    out <- lapply(1:nrow(g$E$nums),
                  getEdgeInfo2,
                  g = g,
                  oldVineGraph = oldVineGraph,
                  weights = weights)
    
    ## annotate graph (same order as in old version of this function)
    g$E$weights         <- sapply(out, function(x) x$tau)
    g$E$names           <- sapply(out, function(x) x$name)
    g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
    g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
    g$E$todel           <- sapply(out, function(x) x$todel)
    
    ## delete edges that are prohibited by the proximity condition
    deleteEdges(g)
}

## function for obtaining edge information
getEdgeInfo2 <- function(i, g, oldVineGraph, weights) {
    
    ## get edge
    con <- g$E$nums[i, ]
    temp <- oldVineGraph$E$nums[con, ]
    
    ## check for proximity condition
    ok <- FALSE
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
        ok <- TRUE
        same <- temp[2, 1]
    } else {
        if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
            ok <- TRUE
            same <- temp[2, 2]
        }
    }
    
    ## dummy output
    tau <- nedSet <- ningSet <- name <- NA
    todel <- TRUE
    
    # info if proximity condition is fulfilled ...
    if (ok) {
        ## get data
        if (temp[1, 1] == same) {
            zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
        } else {
            zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
        }
        if (temp[2, 1] == same) {
            zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
        } else {
            zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
        }
        if (is.list(zr1)) {
            zr1a <- as.vector(zr1[[1]])
            zr2a <- as.vector(zr2[[1]])
        } else {
            zr1a <- zr1
            zr2a <- zr2
        }
        
        ## calculate Kendall's tau
        keine_nas <- !(is.na(zr1a) | is.na(zr2a))
        tau <- fasttau(zr1a[keine_nas], zr2a[keine_nas], weights)
        
        ## get names
        name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
        name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]
        
        ## infer conditioned set and conditioning set
        l1 <- c(g$V$conditionedSet[[con[1]]],
                g$V$conditioningSet[[con[1]]])
        l2 <- c(g$V$conditionedSet[[con[2]]],
                g$V$conditioningSet[[con[2]]])
        nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
        ningSet <- intersect(l1, l2)
        
        ## set edge name
        nmdiff <- c(setdiff(name.node1, name.node2),
                    setdiff(name.node2, name.node1))
        nmsect <- intersect(name.node1, name.node2)
        name <- paste(paste(nmdiff, collapse = ","),
                      paste(nmsect, collapse = ","),
                      sep = " ; ")
        
        ## mark as ok
        todel <- FALSE
    }
    
    ## return edge information
    list(tau = tau,
         nedSet = nedSet,
         ningSet = ningSet,
         name = name,
         todel = todel)
}


wrapper_fit.ACopula <- function(parameterForACopula, type, ...) {
    return(fit.ACopula(parameterForACopula$zr1,
                       parameterForACopula$zr2,
                       type,
                       ...))
}


## bivariate copula selection
fit.ACopula <- function(u1, u2, familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, weights = NA) {
    
    ## select family and estimate parameter(s) for the pair copula
    out <- BiCopSelect(u1, u2,
                       familyset,
                       selectioncrit,
                       indeptest,
                       level,
                       weights = weights,
                       rotations = FALSE)
    
    ## change rotation if family is not symmetric wrt the main diagonal
    if (out$family %in% c(23, 24, 26:30, 124, 224)) {
        out$family <- out$family + 10
    } else if (out$family %in% c(33, 34, 36:40, 134, 234)) {
        out$family <- out$family - 10
    }
    
    ## tawn copulas also have to change type
    if (out$family%/%100 == 1) {
        out$family <- out$family + 100
    } else if (out$family%/%200 == 1) {
        out$family <- out$family - 100
    }
    
    ## store pseudo-observations for estimation in next tree
    out$CondOn.1 <- .C("Hfunc1",
                       as.integer(out$family),
                       as.integer(length(u1)),
                       as.double(u1),
                       as.double(u2),
                       as.double(out$par),
                       as.double(out$par2),
                       as.double(rep(0, length(u1))),
                       PACKAGE = "VineCopula")[[7]]
    out$CondOn.2 <- .C("Hfunc2",
                       as.integer(out$family),
                       as.integer(length(u1)),
                       as.double(u2),
                       as.double(u1),
                       as.double(out$par),
                       as.double(out$par2),
                       as.double(rep(0, length(u1))),
                       PACKAGE = "VineCopula")[[7]]
    
    ## return results
    out
}

## build R-Vine matrix object based on nested set of trees
as.RVM2 <- function(RVine) {
    
    ## initialize objects
    n <- length(RVine$Tree) + 1
    con <- list()
    nam <- RVine$Tree[[1]]$V$names
    nedSets <- list()
    crspParams <- list()
    crspTypes <- list()
    
    ## get selected pairs, families and estimated parameters
    for (k in 1:(n - 2)) {
        nedSets[[k]]    <- RVine$Tree[[k]]$E$conditionedSet
        crspParams[[k]] <- as.list(RVine$Tree[[k]]$E$Copula.param)
        crspTypes[[k]]  <- as.list(RVine$Tree[[k]]$E$Copula.type)
    }
    crspParams[[n - 1]] <- as.list(RVine$Tree[[n - 1]]$E$Copula.param)
    crspTypes[[n - 1]]  <- as.list(RVine$Tree[[n - 1]]$E$Copula.type)
    if (is.list(RVine$Tree[[n - 1]]$E$conditionedSet)) {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet[[1]])
    } else {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet)
    }
    
    ## initialize matrices for RVineMatrix object
    Param <- array(dim = c(n, n))
    Params2 <- array(0, dim = c(n, n))
    Type <- array(dim = c(n, n))
    M <- matrix(NA, n, n)
    
    ## store structure, families and parameters in matrices
    for (k in 1:(n - 1)) {
        w <- nedSets[[n - k]][[1]][1]
        
        M[k, k] <- w
        M[(k + 1), k] <- nedSets[[n - k]][[1]][2]
        
        Param[(k + 1), k]   <- crspParams[[n - k]][[1]][1]
        Params2[(k + 1), k] <- crspParams[[n - k]][[1]][2]
        Type[(k + 1), k]    <- crspTypes[[n - k]][[1]]
        
        if (k == (n - 1)) {
            M[(k + 1), (k + 1)] <- nedSets[[n - k]][[1]][2]
        } else {
            for (i in (k + 2):n) {
                for (j in 1:length(nedSets[[n - i + 1]])) {
                    cs <- nedSets[[n - i + 1]][[j]]
                    cty <- crspTypes[[n - i + 1]][[j]]
                    if (cs[1] == w) {
                        M[i, k] <- cs[2]
                        Type[i, k] <- cty
                        break
                    } else if (cs[2] == w) {
                        # correct family for rotations
                        if (cty %in% c(23, 24, 26:30, 124, 224)) {
                            cty <- cty + 10
                        } else if (cty %in% c(33, 34, 36:40, 134, 234)) {
                            cty <- cty - 10
                        }
                        # change type for Tawn
                        if (cty%/%100 == 1) {
                            cty <- cty + 100
                        } else if (cty%/%200 == 1) {
                            cty <- cty - 100
                        }
                        M[i, k] <- cs[1]
                        Type[i, k] <- cty
                        break
                    }
                }
                Param[i, k]   <- crspParams[[n - i + 1]][[j]][1]
                Params2[i, k] <- crspParams[[n - i + 1]][[j]][2]
                nedSets[[n - i + 1]][[j]]    <- NULL
                crspParams[[n - i + 1]][[j]] <- NULL
                crspTypes[[n - i + 1]][[j]]  <- NULL
            }
        }
    }
    
    ## clean NAs
    M[is.na(M)] <- 0
    Type[is.na(Type)] <- 0
    
    ## return RVineMatrix object
    RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = nam)
}


## functions for handling the tree structure -------------------------
graphFromTauMatrix <- function(tau) {
    d <- ncol(tau)
    # get variable names
    nms <- colnames(tau)
    # construct edge set
    E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
               do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
    # add edge names
    E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
    # set weights
    w <- tau[upper.tri(tau)]
    
    ## return results
    list(V = list(names = nms,
                  conditionedSet = NULL,
                  conditioningSet = NULL),
         E = list(nums = E,
                  names = E.names,
                  weights = w,
                  conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                  conditioningSet = NULL))
}

makeFullGraph <- function(d) {
    ## create matrix of all combinations
    E <- cbind(do.call(c, lapply(1:(d-1), function(i) rep(i, d-i))),
               do.call(c, lapply(1:(d-1), function(i) (i+1):d)))
    E <- matrix(E, ncol = 2)
    
    ## output dummy list with edges set
    list(V = list(names = NULL,
                  conditionedSet = NULL,
                  conditioningSet = NULL),
         E = list(nums = E,
                  names = NULL,
                  weights = NULL,
                  conditionedSet = E,
                  conditioningSet = NULL))
}

adjacencyMatrix <- function(g) {
    ## create matrix of all combinations
    d <- length(g$V$names)
    v.all <- cbind(do.call(c, lapply(1:(d-1), function(i) seq.int(i))),
                   do.call(c, lapply(1:(d-1), function(i) rep(i+1, i))))
    
    ## fnd weight
    vals <- apply(v.all, 1, set_weight, E = g$E)
    
    ## create symmetric matrix of weights
    M <- matrix(0, d, d)
    M[upper.tri(M)] <- vals
    M <- M + t(M)
    diag(M) <- Inf
    
    ## return final matrix
    M
}

set_weight <- function(x, E) {
    is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
    if (!any(is.edge)) Inf else (1 - abs(E$weights[which(is.edge)]))
}


deleteEdges <- function(g) {
    ## reduce edge list
    keep <- which(!g$E$todel)
    E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
              names           = g$E$names[keep],
              weights         = g$E$weights[keep],
              conditionedSet  = g$E$conditionedSet[keep],
              conditioningSet = g$E$conditioningSet[keep])
    
    ## return reduced graph
    list(V = g$V, E = E)
}

