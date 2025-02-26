#' Simulation from an R-Vine Copula Model
#'
#' This function simulates from a given R-vine copula model.
#'
#'
#' @param N Number of d-dimensional observations to simulate.
#' @param RVM An [RVineMatrix()] object containing the information of
#' the R-vine copula model. Optionally, a length-`N` list of
#'   [RVineMatrix()]  objects sharing the same structure, but possibly
#'   different family/parameter can be supplied.
#' @param U If not [NULL()], an (N,d)-matrix of \eqn{U[0,1]} random
#' variates to be transformed to the copula sample.
#' @return An `N` x d matrix of data simulated from the given R-vine
#' copula model.
#' @author Jeffrey Dissmann
#' @seealso [RVineMatrix()], [BiCopSim()]
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#'   (2013). Selecting and estimating regular vine copulae and application to
#'   financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#' @examples
#'
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#'
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#'
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#'
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' # simulate a sample of size 300 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(300, RVM)
#'
RVineSim <- function(N, RVM, U = NULL) {
    ## vectorized call
    if (is.list(RVM) & !is(RVM, "RVineMatrix"))
        return(RVineSimVec(N, RVM, U))

    ## sanity checks
    stopifnot(N >= 1)
    if (!is(RVM, "RVineMatrix"))
        stop("'RVM' has to be an RVineMatrix object.")


    if (any(RVM$family == 1004)) {
        ix <- which(RVM$family == 1004)
        RVM$family[ix] <- 4 + 20 * (RVM$par[ix] < 0)
        RVM$par[ix] <- sign(RVM$par[ix] + 1e-100) * (1 + abs(RVM$par[ix]))
        print(RVM$par)
        print(RVM$family)
    }


    ## reorder matrix and U (if provided)
    n <- dim(RVM)
    o <- diag(RVM$Matrix)
    RVM <- normalizeRVineMatrix(RVM)
    takeU <- !is.null(U)
    if (takeU) {
        if (!is.matrix(U))
            U <- rbind(U, deparse.level = 0L)
        if ((d <- ncol(U)) < 2)
            stop("U should be at least bivariate")  # should be an (N, n) matrix
        U <- U[, rev(o)]
    }

    ## create objects for C-call
    matri <- as.vector(RVM$Matrix)
    w1 <- as.vector(RVM$family)
    th <- as.vector(RVM$par)
    th2 <- as.vector(RVM$par2)
    maxmat <- as.vector(RVM$MaxMat)
    conindirect <- as.vector(RVM$CondDistr$indirect)
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0
    tmp <- rep(0, n * N)

    ## simulate R-Vine
    tmp <- .C("SimulateRVine",
              as.integer(N),
              as.integer(n),
              as.integer(w1),
              as.integer(maxmat),
              as.integer(matri),
              as.integer(conindirect),
              as.double(th),
              as.double(th2),
              as.double(tmp),
              as.double(U),
              as.integer(takeU),
              PACKAGE = "VineCopula")[[9]]

    ## store results, bring back to initial order and return
    out <- matrix(tmp, ncol = n, byrow = TRUE)
    if (!is.null(RVM$names)) {
        colnames(out) <- RVM$names
    }
    out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    return(out)
}


RVineSimVec <- function(N, RVM, U = NULL) {
    ## sanity checks
    stopifnot(N >= 1)
    if (length(RVM) == 0)
        stop("'RVM' is an empty list")
    if (!all(vapply(RVM, is, F, "RVineMatrix")))
        stop("'RVM' has to be an RVineMatrix object or a list thereof.")

    dims <- dim(RVM[[1]])
    if (any(vapply(RVM, function(RVM) dim(RVM) != dims, FALSE)))
        stop("'RVM' is a list of RVineMatrix objects of different dimensions")

    mat <- RVM[[1]]$Matrix
    if (!all(vapply(RVM, function(RVM) identical(RVM$Matrix, mat), FALSE)))
        stop("'RVM' is a list of RVineMatrix objects with different structures")

    ## reorder matrix and U (if provided)
    n <- dim(RVM[[1]])
    o <- diag(RVM[[1]]$Matrix)
    RVM[[1]] <- normalizeRVineMatrix(RVM[[1]])
    takeU <- !is.null(U)
    if (takeU) {
        if (!is.matrix(U))
            U <- rbind(U, deparse.level = 0L)
        if ((d <- ncol(U)) < 2)
            stop("U should be at least bivariate")  # should be an (N, n) matrix
        U <- U[, rev(o)]
    }

    ## create objects for C-call
    matri <- RVM[[1]]$Matrix
    w1 <- extractMat(RVM, 'family')
    th <- extractMat(RVM, 'par')
    th2 <- extractMat(RVM, 'par2')
    maxmat <- RVM[[1]]$MaxMat
    conindirect <- RVM[[1]]$CondDistr$indirect
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0

    ## simulate R-Vine
    tmp <- SimulateRVineVec(N, n, w1, maxmat, matri, conindirect, th, th2,
                            U, takeU)

    ## store results, bring back to initial order and return
    out <- matrix(tmp, ncol = n, byrow = TRUE)
    if (!is.null(RVM[[1]]$names)) {
        colnames(out) <- RVM[[1]]$names
    }
    out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    return(out)
}

extractMat <- function(RVM, name) {
    prop <- vapply(RVM, function(RVM) RVM[[name]], RVM[[1]][[name]])
    array(prop, dim = c(dim(prop)[1:2], length(RVM)))
}

revert <- function(m) {
    if (length(dim(m)) == 2) {
        return(m[nrow(m):1, ncol(m):1, drop = F])
    } else {
        return(m[nrow(m):1, ncol(m):1, , drop = F])
    }
}

SimulateRVineVec <- function(N, d, family, maxmat, mat, cindirect, par, par2,
                             U, takeU) {
    vdirect <- array(dim = c(d, d, N))
    vindirect <- array(dim = c(d, d, N))

    family <- revert(family)
    par <- revert(par)
    par2 <- revert(par2)
    maxmat <- revert(maxmat)
    mat <- revert(mat)
    cindirect <- revert(cindirect)

    ## fill diagonal entries with independent uniforms
    id <- 1:d + d*(0:(d - 1)) + rep((0:(N - 1)) * d^2, each = d)
    vdirect[id] <- if (takeU) t(U) else runif(N * d)
    vindirect[1, 1,] <- vdirect[1, 1,]


    for (i in 2:d) {
        for (k in (i - 1):1) {
            m <- maxmat[k, i]
            u1 <- (if (mat[k, i] == m) vdirect else vindirect)[k, m,]
            vdirect[k, i,] <- BiCopHinv1(u1, vdirect[k + 1, i,],
                                         family[k, i,],
                                         par[k, i,],
                                         par2 = par2[k, i,],
                                         check.pars = FALSE)
            if (i < d) {
                if (cindirect[k + 1, i]) {
                    vindirect[k + 1, i,] <- BiCopHfunc2(u1, vdirect[k, i,],
                                                        family[k, i,],
                                                        par[k, i,],
                                                        par2 = par2[k, i,],
                                                        check.pars = FALSE)
                }
            }
        }
    }

    vdirect[1, , ]
}

