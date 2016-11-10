#' Simulation from an R-Vine Copula Model
#'
#' This function simulates from a given R-vine copula model.
#'
#'
#' @param N Number of d-dimensional observations to simulate.
#' @param RVM An \code{\link{RVineMatrix}} object containing the information of
#' the R-vine copula model.
#' @param U If not \code{\link{NULL}}, an (N,d)-matrix of U[0,1] random
#' variates to be transformed to the copula sample.
#' @return An \code{N} x d matrix of data simulated from the given R-vine
#' copula model.
#' @author Jeffrey Dissmann
#' @seealso \code{\link{RVineMatrix}}, \code{\link{BiCopSim}}
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
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
#' @export RVineSim
RVineSim <- function(N, RVM, U = NULL) {

    ## sanity checks
    stopifnot(N >= 1)
    if (!is(RVM, "RVineMatrix"))
        stop("'RVM' has to be an RVineMatrix object.")

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


transform <- function(M) {
    n <- dim(M)[1]

    M.new <- matrix(rep(0, n * n), n, n)
    for (i in 1:n) {
        for (j in 1:i) {
            M.new[(n - i + 1), (n - j + 1)] <- M[i, j]
        }
    }

    return(M.new)
}

