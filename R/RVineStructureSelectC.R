#' Simplified structure selection function
#'
#' This is a simplified version of \code{\link{RVineStructureSelect}}. It shall
#' serve as a blue-print for a pure C++ implementation. Note that a dot does not
#' have a special meaning in R. So \code{a.b} is just a regular name as is
#' \code{a_b} in C++.
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param familyset An integer vector of pair-copula families to select from.
#' The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence.  Coding of pair-copula families is the same
#' as in \code{\link{BiCop}}.
#' @param trunclevel integer; level of truncation.
#'
#' @return A list of trees, where each tree is represented in our internal data
#' structure. The final structure corresponding to the example code below is
#'
#' \preformatted{
#' List of 3
#'  $ tree_1:List of 2
#'   ..$ V:List of 3
#'   .. ..$ d              : int 4
#'   .. ..$ conditionedSet : NULL
#'   .. ..$ conditioningSet: NULL
#'   ..$ E:List of 5
#'   .. ..$ nums           : num [1:3, 1:2] 1 2 1 2 3 4
#'   .. ..$ weights        : num [1:3] 0.389 0.432 0.339
#'   .. ..$ conditionedSet :List of 3
#'   .. .. ..$ : num [1:2] 1 2
#'   .. .. ..$ : num [1:2] 2 3
#'   .. .. ..$ : num [1:2] 1 4
#'   .. ..$ conditioningSet: NULL
#'   .. ..$ fits           :List of 3
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'  $ tree_2:List of 2
#'   ..$ V:List of 2
#'   .. ..$ d             : int 3
#'   .. ..$ conditionedSet:List of 3
#'   .. .. ..$ : num [1:2] 1 2
#'   .. .. ..$ : num [1:2] 2 3
#'   .. .. ..$ : num [1:2] 1 4
#'   ..$ E:List of 5
#'   .. ..$ nums           : int [1:2, 1:2] 1 1 2 3
#'   .. ..$ weights        : num [1:2] 0.112 0.187
#'   .. ..$ conditionedSet :List of 2
#'   .. .. ..$ : num [1:2] 1 3
#'   .. .. ..$ : num [1:2] 2 4
#'   .. ..$ conditioningSet: num [1:2] 2 1
#'   .. ..$ fits           :List of 2
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'  $ tree_3:List of 2
#'   ..$ V:List of 3
#'   .. ..$ d              : int 2
#'   .. ..$ conditionedSet :List of 2
#'   .. .. ..$ : num [1:2] 1 3
#'   .. .. ..$ : num [1:2] 2 4
#'   .. ..$ conditioningSet: num [1:2] 2 1
#'   ..$ E:List of 7
#'   .. ..$ nums             : int [1, 1:2] 1 2
#'   .. ..$ weights          : num 0.126
#'   .. ..$ conditionedSet   :List of 1
#'   .. .. ..$ : num [1:2] 3 4
#'   .. ..$ conditioningSet  :List of 1
#'   .. .. ..$ : num [1:2] 1 2
#'   .. ..$ fits             :List of 1
#'   .. .. ..$ :List of 17
#'   .. .. .. ..- attr(*, "class")= chr "BiCop"
#'   .. ..$ Copula.CondData.1:List of 1
#'   .. .. ..$ : num [1:250] 0.0641 0.0319 0.425 0.3221 0.7609 ...
#'   .. ..$ Copula.CondData.2:List of 1
#'   .. .. ..$ : num [1:250] 0.933 0.876 0.563 0.396 0.425 ...
#' }
#'
#' @examples
#' data(daxreturns)
#'
#' daxreturns <- daxreturns[1:250, 1:4]
#' RVM <- RVineStructureSelect(daxreturns, 1:6)
#' RVM
#'
#' RVMC <- RVineStructureSelectC(daxreturns, 1:6)
#' sapply(RVMC, function(x) x$E$fits)
RVineStructureSelectC <- function(data, familyset = 0:6, trunclevel = Inf) {
    d <- ncol(data)  # number of variables
    n <- nrow(data)  # number of observations

    # Initialize object for final results.
    RVine <- vector("list", d - 1)
    names(RVine) <- paste0("tree_", 1:(d - 1))
    # List of 3
    # $ tree_1: NULL
    # $ tree_2: NULL
    # $ tree_3: NULL

    ## Estimation in the first tree ----------------------------

    # We start with a full graph, i.e. a graph with
    #   * node set V={1, ..., d},
    #   * edge set E containing all possible edges on V.
    #   * each edge is assigned a weight: the absolute value of Kendall's tau.
    g <- initializeFirstGraphC(data)
    # List of 2
    # $ V:List of 3
    # ..$ d              : int 4
    # ..$ conditionedSet : NULL
    # ..$ conditioningSet: NULL
    # $ E:List of 4
    # ..$ nums           : num [1:6, 1:2] 1 1 2 1 2 3 2 3 3 4 ...
    # ..$ weights        : num [1:6] 0.389 0.315 0.432 0.339 0.324 ...
    # ..$ conditionedSet :List of 6
    # .. ..$ : num [1:2] 1 2
    # .. ..$ : num [1:2] 1 3
    # .. ..$ : num [1:2] 2 3
    # .. ..$ : num [1:2] 1 4
    # .. ..$ : num [1:2] 2 4
    # .. ..$ : num [1:2] 3 4
    # ..$ conditioningSet: NULL

    # Find the maximum spanning tree and remove all other edges from g.
    # Note: The algorithm 'findMaxTreeC' actually finds a minimum spanning tree.
    # That's why all edge weights are multiplied with -1 before running it (see
    # 'set_weightsC' called from 'adjacencyMatrixC' at the beginniung of
    # 'findMaxTreeC').
    g <- findMaxTreeC(g)
    # List of 2
    # $ V:List of 3
    # ..$ d              : int 4
    # ..$ conditionedSet : NULL
    # ..$ conditioningSet: NULL
    # $ E:List of 4
    # ..$ nums           : num [1:3, 1:2] 1 2 1 2 3 4
    # ..$ weights        : num [1:3] 0.389 0.432 0.339
    # ..$ conditionedSet :List of 3
    # .. ..$ : num [1:2] 1 2
    # .. ..$ : num [1:2] 2 3
    # .. ..$ : num [1:2] 1 4
    # ..$ conditioningSet: NULL

    # Estimate a copula for each edge in the tree. The result will be the graph
    # object g from before, but with information about the copula fits added.
    RVine[[1]] <- fit.FirstTreeCopulasC(g, data, familyset)
    # List of 2
    # $ V:List of 3
    # ..$ d              : int 4
    # ..$ conditionedSet : NULL
    # ..$ conditioningSet: NULL
    # $ E:List of 7
    # ..$ nums             : num [1:3, 1:2] 1 2 1 2 3 4
    # ..$ weights          : num [1:3] 0.389 0.432 0.339
    # ..$ conditionedSet   :List of 3
    # .. ..$ : num [1:2] 1 2
    # .. ..$ : num [1:2] 2 3
    # .. ..$ : num [1:2] 1 4
    # ..$ conditioningSet  : NULL
    # ..$ fits             :List of 3
    # .. ..$ :List of 17
    # .. .. ..- attr(*, "class")= chr "BiCop"
    # .. ..$ :List of 17
    # .. .. ..- attr(*, "class")= chr "BiCop"
    # .. ..$ :List of 17
    # .. .. ..- attr(*, "class")= chr "BiCop"
    # ..$ Copula.CondData.1:List of 3
    # .. ..$ : num [1:250] 0.663 0.317 0.76 0.415 0.595 ...
    # .. ..$ : num [1:250] 0.695 0.722 0.467 0.729 0.166 ...
    # .. ..$ : num [1:250] 0.33 0.134 0.629 0.528 0.441 ...
    # ..$ Copula.CondData.2:List of 3
    # .. ..$ : num [1:250] 0.32 0.409 0.314 0.639 0.261 ...
    # .. ..$ : num [1:250] 0.1244 0.0406 0.511 0.2732 0.7695 ...
    # .. ..$ : num [1:250] 0.837 0.747 0.476 0.372 0.423 ...

    ## Estimation in higher trees --------------------------

    # Loop over remaining trees.
    for (t in 2:(d - 1)) {
        # If we have exceeded the truncation level, all copulas will be set to
        # the independence copula from here on. The simplest way to do this is
        # to reduce the familyset to the independence copula only.
        if (t > trunclevel)
            familyset <- 0

        # The edge set of tree_t becomes the node set of the tree_(t+1).
        # The new graph contains all possible edges on this node set - except
        # those prohibited by the proximity condition ("edges corresponding to
        # two connected nodes in tree_i share a common node in tree_(i-1)").
        g <- buildNextGraphC(RVine[[t - 1]])

        # Find the maximum spanning tree and remove all other edges from g.
        g <- findMaxTreeC(g)

        # Estimate a copula for each edge in the tree. The result will be the
        # graph object g from before, but with information about the copula fits
        # added.
        RVine[[t]] <- fit.TreeCopulasC(g, RVine[[t - 1]], familyset)

        # Pseudo-observations from previous trees are useless from here on; we
        # can delete them to save memory.
        RVine[[t - 1]]$E$Copula.CondData.1 <- NULL
        RVine[[t - 1]]$E$Copula.CondData.2 <- NULL
    }

    return(RVine)
}

initializeFirstGraphC <- function(data) {
    d <- ncol(data)  # number of variables

    # Create a matrix of all possible pairs.
    all.pairs <- t(combn(1:d, 2))
    #       [,1] [,2]
    # [1,]    1    2
    # [2,]    1    3
    # [3,]    1    4
    # [4,]    2    3
    # [5,]    2    4
    # [6,]    3    4

    # Compute the Kendall's tau for each pair of variables.
    edge.ws <- numeric(nrow(all.pairs))  # empty vector of length d choose 2.
    for (i in 1:nrow(all.pairs)) {
        e1 <- all.pairs[i, 1]  # index of first node contained in edge i
        e2 <- all.pairs[i, 2]  # index of second node contained in edge i
        edge.ws[i] <- fasttauC(data[, e1], data[, e2])
    }

    # We now store the information in an adjancency matrix W. The matrix W is
    # a symmetric dxd matrix where the (i, j) entry
    # of W contains the weight of the edge (i, j) connecting nodes i and j.
    W <- diag(d)                 # dxd identity matrix
    W[lower.tri(W)] <- edge.ws   # fill lower triangle with weights
    W <- W + t(W)                # add transpose to make it symmetric

    # Convert the proximity matrix into a graph object.
    graph <- graphFromWeightMatrixC(W)

    return(graph)
}

## Initialize graph for next vine tree (possible edges)
buildNextGraphC <- function(oldVineGraph) {
    d <- nrow(oldVineGraph$E$nums)  # number of nodes in the next graph

    # We start with a full graph on d variables.
    g <- makeFullGraphC(d)
    # We also store the conditioned and conditioning sets from the edges
    # in the previous tree. That's important to stay aware of the
    # interconnection between subsequent trees.
    # (Remember: Each node in Tree k has been an edge in Tree k-1.)
    g$V$conditionedSet <- oldVineGraph$E$conditionedSet
    g$V$conditioningSet <- oldVineGraph$E$conditioningSet

    # Calculate necessary info for each edge. This includes the weight
    # (Kendall's tau) and whether the edge is allowed by the proximity
    # condition (if allowed, todel = FALSE).
    for (i in 1:nrow(g$E$nums)) {
        edgeInfo <- getEdgeInfoC(i, g, oldVineGraph)
        g$E$weights[i]           <- edgeInfo$w
        g$E$conditionedSet[[i]]  <- edgeInfo$nedSet
        g$E$conditioningSet[[i]] <- edgeInfo$ningSet
        g$E$todel[i]             <- edgeInfo$todel
    }

    # Delete all edges that are prohibited by the proximity condition.
    g <- deleteEdgesC(g)

    return(g)
}

## This is Prim's algorithm for computing the minimum spanning tree based
## on the adjacency matrix. One can find several implementations in C or C++
## online. If we write our own implementation, I can make the code below read
## more like C++.
findMaxTreeC <- function(g) {
    # Construct adjency matrix from the graph object g.
    A <- adjacencyMatrixC(g)
    d <- ncol(A)  # number of nodes in the tree

    # Set up dummys for the algorithm and output.
    tree <- NULL
    edges <- matrix(NA, d - 1, 2)
    w <- numeric(d - 1)
    i <- 1

    # The actual algorithm:
    for (k in 1:(d - 1)) {
        # add selected edge to tree
        tree <- c(tree, i)

        # find the edge connected to one of the nodes in 'tree' with minimal
        # weight
        m <- apply(as.matrix(A[, tree]), 2, min)
        a <- apply(as.matrix(A[, tree]), 2,
                   function(x) order(rank(x)))[1, ]
        b <- order(rank(m))[1]
        j <- tree[b]
        i <- a[b]

        # store edge and weight
        edges[k, ] <- c(j, i)
        w[k] <- A[i, j]

        # adjust adjecency matrix to prevent loops
        for (t in tree)
            A[i, t] <- A[t, i] <- Inf
    }

    # reorder edges for backwards compatibility with igraph output
    edges <- t(apply(edges, 1, function(x) sort(x)))
    edges <- edges[order(edges[, 2], edges[, 1]), ]

    # delete unused edges from graph
    E <- g$E$nums
    in.tree <- apply(matrix(edges, ncol = 2), 1,
                     function(x) which((x[1] == E[, 1]) & (x[2] == E[, 2])))
    MST <- g
    g$E$todel <- rep(TRUE, nrow(E))
    if (any(g$E$todel)) {
        g$E$todel[in.tree] <- FALSE
        MST <- deleteEdgesC(g)
    }

    return(MST)
}

## Fit pair-copulas for the first vine tree
fit.FirstTreeCopulasC <- function(g, data, familyset) {
    # Loop through all edges.
    d <- nrow(g$E$nums)  # number of edges
    for (i in 1:d) {
        # Get ith edge.
        a <- g$E$nums[i, ]

        # Extract data belonging to the nodes that are connected by edge i.
        u1 <- data[, a[1]]
        u2 <- data[, a[2]]

        # Fit the copula for edge i. Note the reversed order of the arguments
        # u1 and u2.
        pc.fit <- fit.ACopulaC(u2, u1, familyset)

        # Store the fitted copula.
        g$E$fits[[i]] <- pc.fit

        # Store the conditional data (which is required for estimation in the
        # next tree).
        g$E$Copula.CondData.1[[i]] <- pc.fit$CondOn.1
        g$E$Copula.CondData.2[[i]] <- pc.fit$CondOn.2
    }

    return(g)
}

## Fit pair-copulas for vine trees 2, ..., d-1
# Note that there is a lot of overlap with the functino getEdgeInfoC. We do the
# work twice to keep the memory demand low.
fit.TreeCopulasC <- function(g, oldVineGraph, familyset) {
    d <- nrow(g$E$nums)  # number of edges
    # Loop through all edges.
    for (i in 1:d) {
        # Get ith edge.
        con <- g$E$nums[i, ]

        # Get edges from previous tree that corresponds to the nodes that are
        # connected by edge i.
        # (Remember: Each node in Tree k has been an edge in Tree k-1.)
        temp <- oldVineGraph$E$nums[con, ]

        # Check which variable has been conditioned on. It is the number that
        # shows up in both edges from the previous tree ('same').
        if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
            same <- temp[2, 1]
        } else {
            if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
                same <- temp[2, 2]
            }
        }

        # Extract the appropriate conditional data from the edge of the previous
        # tree.
        if (temp[1, 1] == same) {
            u1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
        } else {
            u1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
        }
        if (temp[2, 1] == same) {
            u2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
        } else {
            u2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
        }

        # Fit the copula for edge i. Note the reversed order of the arguments
        # u1 and u2.
        pc.fit <- fit.ACopulaC(u2, u1, familyset)

        # Store the fitted copula.
        g$E$fits[[i]] <- pc.fit

        # Store the conditional data (which is required for estimation in the
        # next tree).
        g$E$Copula.CondData.1[[i]] <- pc.fit$CondOn.1
        g$E$Copula.CondData.2[[i]] <- pc.fit$CondOn.2
    }

    return(g)
}

## Function for obtaining edge information
getEdgeInfoC <- function(i, g, oldVineGraph) {
    # Get ith edge.
    con <- g$E$nums[i, ]

    # Get edges from previous tree that corresponds to the nodes that are
    # connected by edge i.
    # (Remember: Each node in Tree k has been an edge in Tree k-1.)
    temp <- oldVineGraph$E$nums[con, ]

    # Check which variable has been conditioned on. It is the number that
    # shows up in both edges from the previous tree ('same'). If no number shows
    # about in both edges, it must be prohibited by the proximity condition. In
    # this case, no extra work has to be done and the edge will be marked as
    # "to delete" (todel = TRUE).
    todel <- TRUE
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
        todel <- FALSE
        same <- temp[2, 1]
    } else {
        if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
            todel <- FALSE
            same <- temp[2, 2]
        }
    }


    # Only do more work if proximity condition is fulfilled!
    if (!todel) {
        # Extract the appropriate conditional data from the edge of the previous
        # tree.
        if (temp[1, 1] == same) {
            u1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
        } else {
            u1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
        }
        if (temp[2, 1] == same) {
            u2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
        } else {
            u2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
        }

        # Calculate Kendall's tau between the conditional data corresponding to
        # edge i.
        w <- fasttau(u1, u2)

        # Infer conditioned set and conditioning set of edge i:
        #   * The conditioning set consists of all variables that appear in both
        #     edges in the previous tree.
        #   * The remaining two variables are the conditioned set.
        l1 <- c(g$V$conditionedSet[[con[1]]], g$V$conditioningSet[[con[1]]])
        l2 <- c(g$V$conditionedSet[[con[2]]], g$V$conditioningSet[[con[2]]])
        ningSet <- intersect(l1, l2)
        nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    } else {
        # Set dummy output.
        w <- nedSet <- ningSet <- NA
    }

    # Return edge information.
    return(list(w = w,
                nedSet = nedSet,
                ningSet = ningSet,
                todel = todel))
}


## Bivariate copula selection
fit.ACopulaC <- function(u1, u2, familyset) {
    # Fit the copula.
    out <- BiCopSelect(u1, u2, familyset)

    # Calculate conditional data (required for estimation in next tree).
    out$CondOn.1 <- BiCopHfunc1(u1, u2, out)
    out$CondOn.2 <- BiCopHfunc2(u1, u2, out)

    return(out)
}

## Functions for handling the tree structure -------------------------
graphFromWeightMatrixC <- function(W) {
    d <- ncol(W)  # number of nodes in the tree
    # Construct set of all possible edges.
    E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
               do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
    #       [,1] [,2]
    # [1,]    1    2
    # [2,]    1    3
    # [3,]    1    4
    # [4,]    1    5
    # [5,]    2    3
    # [6,]    2    4
    # [7,]    2    5
    # [8,]    3    4
    # [9,]    3    5
    # [10,]   4    5

    # Extract weights.
    w <- W[upper.tri(W)]

    return(list(V = list(d = d,
                         conditionedSet = NULL,
                         conditioningSet = NULL),
                E = list(nums = E,
                         weights = w,
                         # conditionedSet is a list where the ith entry contains
                         # the ith row of E
                         conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                         conditioningSet = NULL)))
}

makeFullGraphC <- function(d) {
    # Create matrix of all combinations of the set V={1, ..., d}.
    E <- cbind(do.call(c, lapply(1:(d-1), function(i) rep(i, d-i))),
               do.call(c, lapply(1:(d-1), function(i) (i+1):d)))
    E <- matrix(E, ncol = 2)
    #       [,1] [,2]
    # [1,]    1    2
    # [2,]    1    3
    # [3,]    1    4
    # [4,]    1    5
    # [5,]    2    3
    # [6,]    2    4
    # [7,]    2    5
    # [8,]    3    4
    # [9,]    3    5
    # [10,]   4    5

    # Output dummy list with edges set.
    return(list(V = list(d = d,
                         conditionedSet = NULL,
                         conditioningSet = NULL),
                E = list(nums = E,
                         weights = NULL,
                         conditionedSet = list(),
                         conditioningSet = list())))
}

adjacencyMatrixC <- function(g) {
    d <- g$V$d  # number of nodes in graph g

    # all.pairs contains all possible edges.
    all.pairs <- t(combn(1:d, 2))
    # They should be ordered like this:
    #       [,1] [,2]
    # [1,]    1    2
    # [2,]    1    3
    # [3,]    1    4
    # [4,]    2    3
    # [5,]    2    4
    # [6,]    3    4

    # Weights are set to minus the absolute value of Kendall's tau;
    # all edges in v.all that are not edges in E get weight Inf.
    edge.ws <- numeric(nrow(all.pairs))
    for (i in 1:length(edge.ws)) {
        edge.ws[i] <- set_weightC(all.pairs[i, ], E = g$E)
    }

    # We now store the information in an adjancency matrix W. The matrix W is
    # a symmetric dxd matrix where the (i, j) entry
    # of W contains the weight of the edge (i, j) connecting nodes i and j.
    W <- diag(d)                 # dxd identity matrix
    W[lower.tri(W)] <- edge.ws   # fill lower triangle with weights
    W <- W + t(W)                # add transpose to make it symmetric

    return(W)
}

set_weightC <- function(x, E) {
    # Check if x is actually and edge in E.
    is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])

    # Put minus before weights so that *minimum* spanning tree algorithm can be
    # applied; return Inf if x is not a member of E
    if (!any(is.edge)) return(Inf) else return(-E$weights[which(is.edge)])
}


deleteEdgesC <- function(g) {
    # Find all edges which have todel = FALSE.
    keep <- which(!g$E$todel)
    # The new edge set is the old edge set with all entries corresponding the
    # todel-edges removed.
    E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
              weights         = g$E$weights[keep],
              conditionedSet  = g$E$conditionedSet[keep],
              conditioningSet = g$E$conditioningSet[keep])

    # Return reduced graph.
    return(list(V = g$V, E = E))
}

## C implementation of Kendall's tau: function 'ktau' in file tools.c
fasttauC <- function(x, y) {
    out <- .C("ktau",
              x = as.double(x),
              y = as.double(y),
              N = as.integer(length(x)),
              tau = as.double(0),
              S = as.double(0),
              D = as.double(0),
              T = as.integer(0),
              U = as.integer(0),
              V = as.integer(0),
              PACKAGE = "VineCopula")

    return(out$tau)
}
