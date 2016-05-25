#' Sequential Specification of R- and C-Vine Copula Models
#'
#' This function fits either an R- or a C-vine copula model to a d-dimensional
#' copula data set. Tree structures are determined and appropriate pair-copula
#' families are selected using \code{\link{BiCopSelect}} and estimated
#' sequentially (forward selection of trees).
#'
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param familyset An integer vector of pair-copula families to select from.
#' The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases.  If \code{familyset = NA} (default), selection among all
#' possible families is performed.  Coding of pair-copula families is the same
#' as in \code{\link{BiCop}}.
#' @param type Type of the vine model to be specified:\cr
#' \code{0} or \code{"RVine"} = R-vine (default) \cr
#' \code{1} or \code{"CVine"} = C-vine \cr
#' C- and D-vine copula models with pre-specified order can be specified using
#' \code{CDVineCopSelect} of the package CDVine. Similarly, R-vine copula
#' models with pre-specified tree structure can be specified using
#' \code{\link{RVineCopSelect}}.
#' @param selectioncrit Character indicating the criterion for pair-copula
#' selection. Possible choices:\code{selectioncrit = "AIC"} (default),
#' \code{"BIC"}, or \code{"logLik"} (see \code{\link{BiCopSelect}}).
#' @param indeptest logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).  The
#' independence copula is chosen for a (conditional) pair if the null
#' hypothesis of independence cannot be rejected.
#' @param level numeri; significance level of the independence test
#' (default: \code{level = 0.05}).
#' @param trunclevel integer; level of truncation.
#' @param progress logical; whether the tree-wise specification progress is
#' printed (default: \code{progress = FALSE}).
#' @param weights numeric; weights for each observation (opitional).
#' @param treecrit edge weight for Dissman's structure selection algorithm, see
#' \emph{Details}.
#' @param rotations If \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#' @param cores integer; if \code{cores > 1}, estimation will be parallized
#' within each tree (using \code{\link[foreach]{foreach}}).
#'
#' @return An \code{\link{RVineMatrix}} object with the selected structure
#' (\code{RVM$Matrix}) and families (\code{RVM$family}) as well as sequentially
#' estimated parameters stored in \code{RVM$par} and \code{RVM$par2}. The object
#' is augmented by the following information about the fit:
#' \item{se, se2}{standard errors for the parameter estimates; note that these
#' are only approximate since they do not
#' account for the sequential nature of the estimation,}
#' \item{nobs}{number of observations,}
#' \item{logLik, pair.logLik}{log likelihood (overall and pairwise)}
#' \item{AIC, pair.AIC}{Aikaike's Informaton Criterion (overall and pairwise),}
#' \item{BIC, pair.BIC}{Bayesian's Informaton Criterion (overall and pairwise),}
#' \item{emptau}{matrix of empirical values of Kendall's tau,}
#' \item{p.value.indeptest}{matrix of p-values of the independence test.}
#'
#' @note For a comprehensive summary of the vine copula model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}.
#'
#' @details
#' R-vine trees are selected using maximum spanning trees w.r.t. some edge
#' weights. The most commonly used edge weigth is the absolute value of the
#' empirical Kendall's tau, say \eqn{\hat{\tau}_{ij}}. Then, the following o
#' ptimization problem is solved for each tree:
#' \deqn{\max \sum_{edges e_{ij} in spanning tree} |\hat{\tau}_{ij}|, }{
#' \max \sum_{edges e_{ij} in spanning tree} |\hat{\tau}_{ij}|, }
#' where a spanning tree is a tree on all nodes. The
#' setting of the first tree selection step is always a complete graph. For
#' subsequent trees, the setting depends on the R-vine construction principles,
#' in particular on the proximity condition.
#'
#' Some commonly used edge weights are implemented:
#' \tabular{ll}{
#' \code{"tau"} \tab absolute value of empirical Kendall's tau. \cr
#' \code{"rho"} \tab absolute value of empirical Spearman's rho. \cr
#' \code{"AIC"} \tab Akaike information (multiplied by -1). \cr
#' \code{"BIC"} \tab Bayesian information criterion (multiplied by -1). \cr
#' \code{"cAIC"}\tab corrected Akaike information criterion (multiplied by -1). \cr
#' }
#' The criteria \code{"AIC"}, \code{"BIC"}, and \code{"cAIC"} require estimation and
#' model selection for all possible pairs. This is computationally expensive and
#' much slower than \code{"tau"} or \code{"rho"}.
#' The user can also specify a custom function to calculate the edge weights.
#' The function has to be of type \code{function(u1, u2, weights) ...} and must
#' return a numeric value. The weigths argument must exist, but does not has to
#' be used. For example, \code{"tau"} (withouth using weights) can be implemented
#' as follows:\cr
#' \code{function(u1, u2, weights)
#'     abs(cor(u1, u2, method = "kendall", use = "complete.obs"))}
#'
#'
#' The root nodes of C-vine trees are determined similarly by identifying the
#' node with strongest dependencies to all other nodes. That is we take the
#' node with maximum column sum in the empirical Kendall's tau matrix.
#'
#' Note that a possible way to determine the order of the nodes in the D-vine
#' is to identify a shortest Hamiltonian path in terms of weights
#' \eqn{1-|\tau_{ij}|}. This can be established for example using the package
#' TSP. Example code is shown below.
#'
#' @author Jeffrey Dissmann, Eike Brechmann, Ulf Schepsmeier, Thomas Nagler
#'
#' @seealso
#' \code{\link{RVineMatrix}},
#' \code{\link{BiCop}},
#' \code{\link{RVineCopSelect}},
#' \code{\link{plot.RVineMatrix}},
#' \code{\link{contour.RVineMatrix}},
#' \code{\link[foreach]{foreach}}
#'
#' @references Brechmann, E. C., C. Czado, and K. Aas (2012). Truncated regular
#' vines in high dimensions with applications to financial data. Canadian
#' Journal of Statistics 40 (1), 68-85.
#'
#' Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
#' Selecting and estimating regular vine copulae and application to financial
#' returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' @examples
#'
#' # load data set
#' data(daxreturns)
#'
#' # select the R-vine structure, families and parameters
#' # using only the first 4 variables and the first 750 observations
#' # we allow for the copula families: Gauss, t, Clayton, Gumbel, Frank and Joe
#' RVM <- RVineStructureSelect(daxreturns[1:750,1:4], c(1:6), progress = TRUE)
#'
#' ## see the object's content or a summary
#' str(RVM)
#' summary(RVM)
#'
#' ## inspect the fitted model using plots
#' \donttest{plot(RVM)  # tree structure}
#' contour(RVM)  # contour plots of all pair-copulas
#'
#' \donttest{## estimate a C-vine copula model with only Clayton, Gumbel and Frank copulas
#' CVM <- RVineStructureSelect(daxreturns, c(3,4,5), "CVine")}
#'
#' \donttest{## determine the order of the nodes in a D-vine using the package TSP
#' library(TSP)
#' d <- dim(daxreturns)[2]
#' M <- 1 - abs(TauMatrix(daxreturns))
#' hamilton <- insert_dummy(TSP(M), label = "cut")
#' sol <- solve_TSP(hamilton, method = "repetitive_nn")
#' order <- cut_tour(sol, "cut")
#' DVM <- D2RVine(order, family = rep(0,d*(d-1)/2), par = rep(0, d*(d-1)/2))
#' RVineCopSelect(daxreturns, c(1:6), DVM$Matrix)}
#'
RVineStructureSelect <- function(data, familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE,
                                 level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA,
                                 treecrit = "tau", rotations = TRUE, cores = 1) {
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
    if (any(data[!is.na(data)] > 1) || any(data[!is.na(data)] < 0))
        stop("Data has to be in the interval [0,1].")
    if (!is.na(familyset[1])) {
        if (!all(abs(familyset) %in% allfams))
            stop("Copula family not implemented.")
        if (length(unique(sign(familyset))) != 1)
            stop("'familyset' must not contain positive AND negative numbers")
    } else {
        familyset <- allfams
    }
    if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
        stop("Selection criterion not implemented.")
    if (level < 0 & level > 1)
        stop("Significance level has to be between 0 and 1.")
    treecrit <- set_treecrit(treecrit, familyset)

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
    g <- initializeFirstGraph(data, treecrit, weights)
    MST <- findMaxTree(g, mode = type)

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    # estimate pair-copulas
    VineTree <- fit.FirstTreeCopulas(MST,
                                     data,
                                     familyset,
                                     selectioncrit,
                                     indeptest,
                                     level,
                                     weights = weights,
                                     cores = cores)
    # store results
    RVine$Tree[[1]] <- VineTree
    RVine$Graph[[1]] <- g
    oldVineGraph <- VineTree

    ## estimation in higher trees --------------------------
    for (tree in 2:(d - 1)) {
        ## old pseudo-observations are uneccessary in RVine object
        RVine$Tree[[tree - 1]]$E$Copula.CondData.1 <- NULL
        RVine$Tree[[tree - 1]]$E$Copula.CondData.2 <- NULL

        # only estimate pair-copulas if not truncated
        if (trunclevel == tree - 1)
            familyset <- 0
        # find optimal tree
        g <- buildNextGraph(VineTree, weights, treecrit = treecrit, cores > 1)
        MST <- findMaxTree(g, mode = type)
        # estimate pair-copulas
        VineTree <- fit.TreeCopulas(MST,
                                    VineTree,
                                    familyset,
                                    selectioncrit,
                                    indeptest,
                                    level,
                                    progress,
                                    weights = weights,
                                    cores = cores)
        # store results
        RVine$Tree[[tree]] <- VineTree
        RVine$Graph[[tree]] <- g
    }

    ## free memory and return results as 'RVineMatrix' object
    .RVine <- RVine
    .data <- data
    .callexp <- match.call()
    rm(list = ls())
    as.RVM2(.RVine, .data, .callexp)
}


set_treecrit <- function(treecrit, famset) {
    ## check if function is appropriate or type is implemented
    if (is.function(treecrit)) {
        if (!all(names(formals(treecrit)) == c("u1", "u2", "weights")))
            stop("treecrit must be of the form 'function(u1, u2, weights)'")
        if (!is.numeric(treecrit(runif(10), runif(10), rep(1, 10))))
            stop("treecrit does not return a numeric value")
    } else if (all(treecrit == "tau")) {
        treecrit <- function(u1, u2, weights) {
            complete.i <- which(!is.na(u1 + u2))
            complete.freq <- mean(!is.na(u1 + u2))
            tau <- abs(fasttau(u1[complete.i], u2[complete.i], weights))
            tau * sqrt(complete.freq)
        }
    } else if (all(treecrit == "rho")) {
        treecrit <- function(u1, u2, weights) {
            complete.freq <- mean(!is.na(u1 + u2))
            rho <- abs(cor(u1, u2, method = "spearman", use = "complete.obs"))
            rho * sqrt(complete.freq)
        }
    } else if (all(treecrit == "AIC")) {
        treecrit <- function(u1, u2, weights)
            - suppressWarnings(BiCopSelect(u1, u2,
                                           familyset = famset,
                                           weights = weights)$AIC)
    } else if (all(treecrit == "BIC")) {
        treecrit <- function(u1, u2, weights)
            - suppressWarnings(BiCopSelect(u1, u2,
                                           familyset = famset,
                                           weights = weights)$BIC)
    } else if (all(treecrit == "cAIC")) {
        treecrit <- function(u1, u2, weights) {
            fit <- suppressWarnings(BiCopSelect(u1, u2,
                                                familyset = famset,
                                                weights = weights))
            n <- fit$nobs
            p <- fit$npars
            - (fit$AIC + 2 * p * (p + 1) / (n - p - 1))
        }
    } else {
        txt1 <- 'treecrit must be one of "tau", "rho", "AIC", "BIC", "cAIC"'
        txt2 <- 'or a function like function(u1, u2, weights) ... returning a numeric value.'
        stop(paste(txt1, txt2))
    }

    ## return treecrit function
    treecrit
}

initializeFirstGraph <- function(data, treecrit, weights) {
    ## calculate edge weight for each possible pair
    all.pairs <- combn(1:ncol(data), 2)
    edge.ws <- apply(all.pairs, 2,
                     function(ind)
                         treecrit(data[, ind[1]], data[, ind[2]], weights))
    # number of pairwise complete observations / all observations
    rel.nobs <- apply(all.pairs, 2,
                      function(ind)
                          mean(!is.na(data[, ind[1]] + data[, ind[2]])))
    edge.ws <- edge.ws


    ## store in symmetric matrix with appropriate names
    W <- diag(ncol(data))
    W[lower.tri(W)] <- edge.ws
    W <- t(W)
    colnames(W) <- rownames(W) <- colnames(data)

    ## return full graph with edge weights
    graphFromWeightMatrix(W)
}

findMaxTree <- function(g, mode = "RVine") {
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
fit.FirstTreeCopulas <- function(MST, data.univ, type, copulaSelectionBy,
                                 testForIndependence, testForIndependence.level,
                                 weights = NA, cores = 1) {

    ## initialize estimation results with empty list
    d <- nrow(MST$E$nums)
    pc.data <- lapply(1:d, function(i) NULL)

    ## prepare for estimation and store names
    for (i in 1:d) {
        ## get edge and corresponding data
        a <- MST$E$nums[i, ]
        pc.data[[i]]$zr1 <- data.univ[, a[1]]
        pc.data[[i]]$zr2 <- data.univ[, a[2]]
        #         MST$E$Copula.Data.1[i] <- list(data.univ[, a[1]])
        #         MST$E$Copula.Data.2[i] <- list(data.univ[, a[2]])

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
    if (cores > 1) {
        pc.fits <- foreach(i = 1:length(pc.data),
                           .export = c("pcSelect",
                                       "fit.ACopula",
                                       "BiCopSelect")) %dopar%
            pcSelect(pc.data[[i]],
                     type,
                     copulaSelectionBy,
                     testForIndependence,
                     testForIndependence.level,
                     weights)
    } else {
        pc.fits <- lapply(X = pc.data,
                          FUN = pcSelect,
                          type,
                          copulaSelectionBy,
                          testForIndependence,
                          testForIndependence.level,
                          weights)
    }

    ## store estimated model and pseudo-obversations for next tree
    for (i in 1:d) {
        MST$E$Copula.param[[i]] <- c(pc.fits[[i]]$par,
                                     pc.fits[[i]]$par2)
        MST$E$Copula.type[i] <- pc.fits[[i]]$family
        MST$E$fits[[i]] <- pc.fits[[i]]

        MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
        MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)
    }

    ## return results
    MST
}

## fit pair-copulas for vine trees 2,...
fit.TreeCopulas <- function(MST, oldVineGraph, type, copulaSelectionBy,
                            testForIndependence, testForIndependence.level,
                            progress, weights = NA, cores = 1) {

    ## initialize estimation results with empty list
    d <- nrow(MST$E$nums)
    pc.data <- lapply(1:d, function(i) NULL)


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

        pc.data[[i]]$zr1 <- zr1a
        pc.data[[i]]$zr2 <- zr2a

        #         MST$E$Copula.Data.1[i] <- list(zr1a)
        #         MST$E$Copula.Data.2[i] <- list(zr2a)

        MST$E$Copula.CondName.1[i] <- n1a
        MST$E$Copula.CondName.2[i] <- n2a
    }

    ## estimate parameters and select family
    if (cores > 1) {
        pc.fits <- foreach(i = 1:length(pc.data),
                           .export = c("pcSelect", "fit.ACopula"),
                           .packages = "VineCopula") %dopar%
            pcSelect(pc.data[[i]],
                     type,
                     copulaSelectionBy,
                     testForIndependence,
                     testForIndependence.level,
                     weights)
    } else {
        pc.fits <- lapply(X = pc.data,
                          FUN = pcSelect,
                          type,
                          copulaSelectionBy,
                          testForIndependence,
                          testForIndependence.level,
                          weights)
    }

    ## store estimated model and pseudo-obversations for next tree

    for (i in 1:d) {
        MST$E$Copula.param[[i]] <- c(pc.fits[[i]]$par,
                                     pc.fits[[i]]$par2)
        MST$E$Copula.type[i] <- pc.fits[[i]]$family
        MST$E$fits[[i]] <- pc.fits[[i]]
        MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
        MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)
    }

    ## return results
    MST
}

## initialize graph for next vine tree (possible edges)
buildNextGraph <- function(oldVineGraph, treecrit, weights = NA, parallel) {

    d <- nrow(oldVineGraph$E$nums)

    ## initialize with full graph
    g <- makeFullGraph(d)
    g$V$names <- oldVineGraph$E$names
    g$V$conditionedSet <- oldVineGraph$E$conditionedSet
    g$V$conditioningSet <- oldVineGraph$E$conditioningSet

    ## get info for all edges
    if (parallel) {
        i <- NA  # dummy for CRAN checks
        out <- foreach(i = 1:nrow(g$E$nums)) %dopar%
            getEdgeInfo(i,
                        g = g,
                        oldVineGraph = oldVineGraph,
                        treecrit = treecrit,
                        weights = weights)
    } else {
        out <- lapply(1:nrow(g$E$nums),
                      getEdgeInfo,
                      g = g,
                      oldVineGraph = oldVineGraph,
                      treecrit = treecrit,
                      weights = weights)
    }

    ## annotate graph (same order as in old version of this function)
    g$E$weights         <- sapply(out, function(x) x$w)
    g$E$names           <- sapply(out, function(x) x$name)
    g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
    g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
    g$E$todel           <- sapply(out, function(x) x$todel)

    ## delete edges that are prohibited by the proximity condition
    deleteEdges(g)
}

## function for obtaining edge information
getEdgeInfo <- function(i, g, oldVineGraph, treecrit, weights) {

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
    w <- nedSet <- ningSet <- name <- NA
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
        w <- treecrit(zr1a[keine_nas], zr2a[keine_nas], weights)

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
    list(w = w,
         nedSet = nedSet,
         ningSet = ningSet,
         name = name,
         todel = todel)
}


pcSelect <- function(parameterForACopula, type, ...) {
    return(fit.ACopula(parameterForACopula$zr1,
                       parameterForACopula$zr2,
                       type,
                       ...))
}


## bivariate copula selection
fit.ACopula <- function(u1, u2, familyset = NA, selectioncrit = "AIC",
                        indeptest = FALSE, level = 0.05, weights = NA) {

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
as.RVM2 <- function(RVine, data, callexp) {

    ## initialize objects
    n <- length(RVine$Tree) + 1
    con <- list()
    nam <- RVine$Tree[[1]]$V$names
    nedSets <- list()
    crspParams <- list()
    crspTypes <- list()
    crspfits <- list()

    ## get selected pairs, families and estimated parameters
    for (k in 1:(n - 2)) {
        nedSets[[k]]    <- RVine$Tree[[k]]$E$conditionedSet
        crspParams[[k]] <- as.list(RVine$Tree[[k]]$E$Copula.param)
        crspTypes[[k]]  <- as.list(RVine$Tree[[k]]$E$Copula.type)
        crspfits[[k]]   <- as.list(RVine$Tree[[k]]$E$fits)

    }
    crspParams[[n - 1]] <- as.list(RVine$Tree[[n - 1]]$E$Copula.param)
    crspTypes[[n - 1]]  <- as.list(RVine$Tree[[n - 1]]$E$Copula.type)
    crspfits[[n - 1]]   <- as.list(RVine$Tree[[n - 1]]$E$fits)
    if (is.list(RVine$Tree[[1]]$E$conditionedSet)) {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet[[1]])
    } else {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet)
    }

    ## initialize matrices for RVineMatrix object
    Param <- array(dim = c(n, n))
    Params2 <- array(0, dim = c(n, n))
    Type <- array(dim = c(n, n))
    M <- matrix(NA, n, n)
    Ses     <- matrix(0, n, n)
    tmps    <- matrix(0, n, n)
    Se2s    <- matrix(0, n, n)
    emptaus <- matrix(0, n, n)
    pvals   <- matrix(0, n, n)

    ## store structure, families and parameters in matrices
    for (k in 1:(n - 1)) {
        w <- nedSets[[n - k]][[1]][1]

        M[k, k] <- w
        M[(k + 1), k] <- nedSets[[n - k]][[1]][2]

        Param[(k + 1), k]   <- crspParams[[n - k]][[1]][1]
        Params2[(k + 1), k] <- crspParams[[n - k]][[1]][2]
        Type[(k + 1), k]    <- crspTypes[[n - k]][[1]]
        Ses[(k + 1), k]     <- crspfits[[n - k]][[1]]$se
        tmpse2              <- crspfits[[n - k]][[1]]$se2
        Se2s[(k + 1), k]    <- ifelse(is.null(tmpse2), NA, tmpse2)
        emptaus[(k + 1), k] <- crspfits[[n - k]][[1]]$emptau
        pvals[(k + 1), k]   <- crspfits[[n - k]][[1]]$p.value.indeptest

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
                Ses[i, k]     <- crspfits[[n - i + 1]][[j]]$se
                tmpse2        <- crspfits[[n - i + 1]][[j]]$se2
                Se2s[i, k]    <- ifelse(is.null(tmpse2), NA, tmpse2)
                emptaus[i, k] <- crspfits[[n - i + 1]][[j]]$emptau
                pvals[i, k]   <- crspfits[[n - i + 1]][[j]]$p.value.indeptest
                nedSets[[n - i + 1]][[j]]    <- NULL
                crspParams[[n - i + 1]][[j]] <- NULL
                crspTypes[[n - i + 1]][[j]]  <- NULL
                crspfits[[n - i + 1]][[j]] <- NULL
            }
        }
    }

    ## clean NAs
    M[is.na(M)] <- 0
    Type[is.na(Type)] <- 0

    ## create RVineMatrix object
    RVM <- RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = nam)
    RVM$call <- callexp

    ## add information about pair-copulas
    RVM$se <- Ses
    RVM$se2 <- Se2s
    RVM$nobs <- crspfits[[1]][[1]]$nobs
    like <- RVineLogLik(data, RVM)
    RVM$logLik <- like$loglik
    RVM$pair.logLik <- like$V$value
    npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- (RVM$family %in% allfams[onepar]) +
        2 * (RVM$family %in% allfams[twopar])
    N <- nrow(data)
    RVM$AIC <- -2 * like$loglik + 2 * npar
    RVM$pair.AIC <- -2 * like$V$value + 2 * npar_pair
    RVM$BIC <- -2 * like$loglik + log(N) * npar
    RVM$pair.BIC <- -2 * like$V$value + log(N) * npar_pair
    RVM$emptau <- emptaus

    ## return final object
    RVM
}


## functions for handling the tree structure -------------------------
graphFromWeightMatrix <- function(W) {
    d <- ncol(W)
    # get variable names
    nms <- colnames(W)
    if (is.null(nms))
        nms <- paste0("V", 1:d)
    # construct edge set
    E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
               do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
    # add edge names
    E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
    # set weights
    w <- W[upper.tri(W)]

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
    ## convert weights so that minimum spanning tree algorithm can be applied
    is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
    if (!any(is.edge)) Inf else -E$weights[which(is.edge)]
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

