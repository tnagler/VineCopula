#' Sequential Estimation of an R-Vine Copula Model
#'
#' This function sequentially estimates the pair-copula parameters of a
#' d-dimensional R-vine copula model as specified by the corresponding
#' \code{\link{RVineMatrix}} object.
#'
#' The pair-copula parameter estimation is performed tree-wise, i.e., for each
#' R-vine tree the results from the previous tree(s) are used to calculate the
#' new copula parameters using \code{\link{BiCopEst}}.
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An \code{\link{RVineMatrix}} object including the structure, the
#' pair-copula families and the pair-copula parameters (if they are known).
#' @param method Character indicating the estimation method: either pairwise
#' maximum likelihood estimation (\code{method = "mle"}; default) or inversion
#' of Kendall's tau (\code{method = "itau"}; see \code{\link{BiCopEst}}.  For
#' \code{method = "itau"} only one parameter pair-copula families can be used
#' (\code{family = 1}, \code{3}, \code{4}, \code{5}, \code{6}, \code{13},
#' \code{14}, \code{16}, \code{23}, \code{24}, \code{26}, \code{33}, \code{34}
#' or \code{36}).
#' @param se Logical; whether standard errors are estimated (default: \code{se
#' = FALSE}).
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: \code{max.df = 30}; for more
#' details see \code{\link{BiCopEst}}).
#' @param max.BB List; upper bounds for the estimation of the two parameters
#' (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr (default:
#' \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).
#' @param progress Logical; whether the pairwise estimation progress is printed
#' (default: \code{progress = FALSE}).
#' @param weights Numerical; weights for each observation (opitional).
#' @param cores integer; if \code{cores > 1}, estimation will be parallized
#' within each tree (using \code{\link[foreach]{foreach}}).
#'
#' @return \item{RVM}{\code{\link{RVineMatrix}} object with the sequentially
#' estimated parameters stored in \code{RVM$par} and \code{RVM$par2}.}
#' \item{se}{Lower triangular d x d matrix with estimated standard errors of
#' the (first) pair-copula parameters for each (conditional) pair defined in
#' the \code{\link{RVineMatrix}} object (if \code{se = TRUE}).}
#' \item{se2}{Lower triangular d x d matrix with estimated standard errors of
#' the second parameters for pair-copula families with two parameters for each
#' (conditional) pair defined in the \code{\link{RVineMatrix}} object (if
#' \code{se = TRUE}).}
#'
#' @author Ulf Schepsmeier, Jeffrey Dissmann, Thomas Nagler
#'
#' @seealso \code{\link{BiCopEst}}, \code{\link{BiCopHfunc}},
#' \code{\link{RVineLogLik}}, \code{\link{RVineMLE}}, \code{\link{RVineMatrix}},
#' \code{\link[foreach]{foreach}}
#'
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
#' # sequential estimation
#' RVineSeqEst(simdata, RVM, method = "itau", se = TRUE)
#' RVineSeqEst(simdata, RVM, method = "mle", se = TRUE)
#'
#' @export RVineSeqEst
RVineSeqEst <- function(data, RVM, method = "mle", se = FALSE, max.df = 30,
                         max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)),
                         progress = FALSE, weights = NA, cores = 1) {
    data <- as.matrix(data)
    Matrix <- RVM$Matrix
    d <- n <- ncol(data)
    N <- nrow(data)

    ## sanity checks
    if (!is(RVM, "RVineMatrix"))
        stop("'RVM' has to be an RVineMatrix object.")
    if (nrow(Matrix) != ncol(Matrix))
        stop("Structure matrix has to be quadratic.")
    if (max(Matrix) > nrow(Matrix))
        stop("Error in the structure matrix.")
    if (N < 2)
        stop("Number of observations has to be at least 2.")
    if (d < 2)
        stop("Dimension has to be at least 2.")
    if (any(data > 1) || any(data < 0))
        stop("Data has be in the interval [0,1].")
    if (method != "mle" && method != "itau")
        stop("Estimation method has to be either 'mle' or 'itau'.")
    if (is.logical(se) == FALSE)
        stop("'se' has to be a logical variable (TRUE or FALSE).")
    if (max.df <= 1)
        stop("The upper bound for the degrees of freedom parameter has to be larger than 1.")
    if (!is.list(max.BB))
        stop("'max.BB' has to be a list.")
    if (max.BB$BB1[1] < 0.001)
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB1[2] < 1.001)
        stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[1] < 1.001)
        stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[2] < 1.001)
        stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[1] < 1.001)
        stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[2] < 0.001)
        stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[1] < 1.001)
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1)
        stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")

    ## set variable names and trunclevel if not provided
    if (is.null(colnames(data)))
        colnames(data) <- paste("V", 1:d, sep = "")
    varnames <- colnames(data)

    ## reorder matrix to natural order
    Matrix <- ToLowerTri(Matrix)
    M <- Matrix
    Mold <- M
    o <- diag(M)
    M <- reorderRVineMatrix(M)
    data <- data[, o[length(o):1]]

    ## create matrices required for selection of h-functions
    MaxMat <- createMaxMat(M)
    CondDistr <- neededCondDistr(M)

    ## create objects for results
    Types   <- matrix(0, d, d)
    Params  <- matrix(0, d, d)
    Params2 <- matrix(0, d, d)
    if (se)
        Se <- Se2 <- matrix(0, d, d)
    V <- list()
    V$direct <- array(NA, dim = c(d, N))
    V$indirect <- array(NA, dim = c(d, N))
    V$direct <- t(data[, d:1])

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

    ## loop over all trees and pair-copulas
    for (k in d:2) {
        doEst <- function(i) {
            if (k > i) {
                ## get pseudo-observaions
                m <- MaxMat[k, i]
                zr1 <- V$direct[i, ]

                zr2 <- if (m == M[k, i]) {
                    V$direct[(d - m + 1), ]
                } else {
                    V$indirect[(d - m + 1), ]
                }

                if (progress == TRUE) {
                    if (k == n) {
                        message(Mold[i, i],
                                ",",
                                Mold[k, i])
                    } else {
                        message(Mold[i, i],
                                ",",
                                Mold[k, i],
                                "|",
                                paste(Mold[(k + 1):n, i],
                                      collapse = ","))
                    }
                }

                ## select pair-copula
                cfit <- BiCopEst(zr2,
                                 zr1,
                                 RVM$family[k, i],
                                 method,
                                 se,
                                 max.df,
                                 max.BB,
                                 weights)

                ## transform data to pseudo-oberstavions in next tree
                direct <- indirect <- NULL
                if (CondDistr$direct[k - 1, i])
                    direct <- BiCopHfunc(zr1,
                                         zr2,
                                         cfit,
                                         check.pars = FALSE)$hfunc2
                if (CondDistr$indirect[k - 1, i])
                    indirect <- BiCopHfunc(zr1,
                                           zr2,
                                           cfit,
                                           check.pars = FALSE)$hfunc1

                ## return results
                list(direct = direct, indirect = indirect, cfit = cfit)
            } else {
                list(cfit = BiCop(0, 0))
            }
        }

        ## run pair-copula estimation for tree k
        res.k <- if (cores > 1) {
            foreach(i = 1:(k-1),
                    .packages = c("VineCopula"),
                    .export = ) %dopar% doEst(i)
        } else {
            lapply(1:(k-1), doEst)
        }

        for (i in 1:(k-1)) {
            ## store results in matrices
            Types[k, i] <- res.k[[i]]$cfit$family
            Params[k, i] <- res.k[[i]]$cfit$par
            Params2[k, i] <- res.k[[i]]$cfit$par2
            if (se == TRUE) {
                # se1 <- par.out$se
                Se[k, i] <- res.k[[i]]$cfit$se
                Se2[k, i] <- res.k[[i]]$cfit$se2
            }
            if (!is.null(res.k[[i]]$direct))
                V$direct[i, ] <- res.k[[i]]$direct
            if (!is.null(res.k[[i]]$indirect))
                V$indirect[i, ] <- res.k[[i]]$indirect
        } # end i = 1:(d-1)
    } # end k = d:2

    ## free memory and return results
    .RVM <- RVineMatrix(Mold,
                        family = Types,
                        par = Params,
                        par2 = Params2,
                        names = varnames)
    if (se) {
        .RVM$se <- Se
        .RVM$se2 <- Se2
    }
    rm(list = ls())
    .RVM
}
