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
#' @return \item{RVM}{\code{\link{RVineMatrix}} object with the sequentially
#' estimated parameters stored in \code{RVM$par} and \code{RVM$par2}.}
#' \item{se}{Lower triangular d x d matrix with estimated standard errors of
#' the (first) pair-copula parameters for each (conditional) pair defined in
#' the \code{\link{RVineMatrix}} object (if \code{se = TRUE}).}
#' \item{se2}{Lower triangular d x d matrix with estimated standard errors of
#' the second parameters for pair-copula families with two parameters for each
#' (conditional) pair defined in the \code{\link{RVineMatrix}} object (if
#' \code{se = TRUE}).}
#' @author Ulf Schepsmeier, Jeffrey Dissmann
#' @seealso \code{\link{BiCopEst}}, \code{\link{BiCopHfunc}},
#' \code{\link{RVineLogLik}}, \code{\link{RVineMLE}}, \code{\link{RVineMatrix}}
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
                        progress = FALSE, weights = NA) {
    data <- as.matrix(data)
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    n <- dim(RVM)
    N <- nrow(data)
    if (dim(data)[2] != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (N < 2) 
        stop("Number of observations has to be at least 2.")
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
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
    
    o <- diag(RVM$Matrix)
    
    oldRVM <- RVM
    
    if (any(o != length(o):1)) {
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }
    
    Params <- RVM$par
    Params2 <- RVM$par2
    
    if (se == TRUE) {
        seMat1 <- matrix(0, nrow = n, ncol = n)
        seMat2 <- matrix(0, nrow = n, ncol = n)
    }
    
    V <- list()
    V$direct <- array(NA, dim = c(n, n, N))
    V$indirect <- array(NA, dim = c(n, n, N))
    
    V$direct[n, , ] <- t(data[, n:1])
    
    for (i in (n - 1):1) {
        
        for (k in n:(i + 1)) {
            
            m <- RVM$MaxMat[k, i]
            zr1 <- V$direct[k, i, ]
            
            if (m == RVM$Matrix[k, i]) {
                zr2 <- V$direct[k, (n - m + 1), ]
            } else {
                zr2 <- V$indirect[k, (n - m + 1), ]
            }
            
            
            if (RVM$family[k, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20,
                                        27, 28, 29, 30, 37, 38, 39, 40)) {
                if (progress == TRUE) {
                    if (k == n) {
                        message(oldRVM$Matrix[i, i],
                                ",",
                                oldRVM$Matrix[k, i]) 
                    } else {
                        message(oldRVM$Matrix[i, i], 
                                ",", 
                                oldRVM$Matrix[k, i],
                                "|", 
                                paste(oldRVM$Matrix[(k + 1):n, i],
                                      collapse = ","))
                    }
                }
                par.out <- BiCopEst(zr2, 
                                    zr1, 
                                    RVM$family[k, i], 
                                    method,
                                    se,
                                    max.df, 
                                    max.BB,
                                    weights)
                # par1 <- out.par$par
                Params[k, i] <- par.out$par
                Params2[k, i] <- par.out$par2
                if (se == TRUE) {
                    # se1 <- par.out$se
                    seMat1[k, i] <- par.out$se
                    seMat2[k, i] <- par.out$se2
                }
            } else {
                if (progress == TRUE) {
                    if (k == n) {
                        message(oldRVM$Matrix[i, i],
                                ",",
                                oldRVM$Matrix[k, i])
                    } else {
                        message(oldRVM$Matrix[i, i], 
                                ",",
                                oldRVM$Matrix[k, i], 
                                "|", 
                                paste(oldRVM$Matrix[(k + 1):n, i],
                                      collapse = ","))
                    }
                }
                par.out <- BiCopEst(zr2,
                                    zr1,
                                    RVM$family[k, i], 
                                    method, 
                                    se,
                                    max.df,
                                    max.BB,
                                    weights)
                # par1 <- out.par$par
                Params[k, i] <- par.out$par
                Params2[k, i] <- par.out$par2
                if (se == TRUE) {
                    # se1 <- par.out$se
                    seMat1[k, i] <- par.out$se
                    seMat2[k, i] <- par.out$se2
                }
            }
            
            
            if (RVM$CondDistr$direct[k - 1, i]) {
                V$direct[k - 1, i, ] <- .C("Hfunc1",
                                           as.integer(RVM$family[k, i]),
                                           as.integer(length(zr1)), 
                                           as.double(zr1), 
                                           as.double(zr2), 
                                           as.double(Params[k, i]),
                                           as.double(Params2[k, i]), 
                                           as.double(rep(0, length(zr1))), 
                                           PACKAGE = "VineCopula")[[7]]
            }
            if (RVM$CondDistr$indirect[k - 1, i]) {
                V$indirect[k - 1, i, ] <- .C("Hfunc2", 
                                             as.integer(RVM$family[k, i]), 
                                             as.integer(length(zr2)),
                                             as.double(zr2),
                                             as.double(zr1),
                                             as.double(Params[k, i]), 
                                             as.double(Params2[k, i]), 
                                             as.double(rep(0, length(zr1))), 
                                             PACKAGE = "VineCopula")[[7]]
            }
            
        }
    }
    
    ## store results
    oldRVM$par <- Params
    oldRVM$par2 <- Params2
    if (se == FALSE) {
        .out <- list(RVM = oldRVM)
    } else {
        .out <- list(RVM = oldRVM, se = seMat1, se2 = seMat2)
    }
    
    ## free memory and return results
    rm(list = ls())
    .out
}
