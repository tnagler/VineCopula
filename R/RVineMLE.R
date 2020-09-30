#' Maximum Likelihood Estimation of an R-Vine Copula Model
#'
#' This function calculates the maximum likelihood estimate (MLE) of the
#' R-vine copula model parameters using sequential estimates as initial values
#' (if not provided).
#'
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An [RVineMatrix()] object including the structure and
#' the pair-copula families and parameters (if known).
#' @param start Lower triangular d x d matrix with zero diagonal entries with
#' starting values for the pair-copula parameters (optional; otherwise they are
#' calculated via \cr [RVineSeqEst()]; default: `start =
#' RVM$par`).
#' @param start2 Lower triangular d x d matrix with zero diagonal entries with
#' starting values for the second parameters of pair-copula families with two
#' parameters (optional; otherwise they are calculated via
#' [RVineSeqEst()]; default: `start2 = RVM$par2`).
#' @param maxit The maximum number of iteration steps (optional; default:
#' `maxit = 200`).
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: `max.df = 30`; for more
#' details see [BiCopEst()]).
#' @param max.BB List; upper bounds for the estimation of the two parameters
#' (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr (default:
#' `max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))`).
#' @param grad If RVM$family only contains one parameter copula families or the
#' t-copula the analytical gradient can be used for maximization of the
#' log-likelihood (see [RVineGrad()]; default: `grad = FALSE`).
#' @param hessian Logical; whether the Hessian matrix of parameter estimates is
#' estimated (default: `hessian = FALSE`). Note that this is not the
#' Hessian Matrix calculated via [RVineHessian()] but via finite
#' differences.
#' @param se Logical; whether standard errors of parameter estimates are
#' estimated on the basis of the Hessian matrix (see above; default: `se =
#' FALSE`).
#' @param ... Further arguments for `optim` (e.g. `factr` controls
#' the convergence of the "L-BFGS-B" method, or `trace`, a non-negative
#' integer, determines if tracing information on the progress of the
#' optimization is produced.) \cr For more details see the documentation of
#' [optim()].
#'
#' @return \item{RVM}{[RVineMatrix()] object with the calculated
#' parameters stored in `RVM$par` and `RVM$par2`. Additional
#' information about the fit is added (e.g., log-likelihood, AIC, BIC).}
#' \item{value}{Optimized log-likelihood value corresponding to the estimated
#' pair-copula parameters.} \item{convergence}{An integer code indicating
#' either successful convergence (`convergence = 0`) or an error:\cr
#' `1` = the iteration limit `maxit` has been reached \cr `51` =
#' a warning from the "L-BFGS-B" method; see component `message` for
#' further details \cr `52` = an error from the "L-BFGS-B" method; see
#' component `message` for further details} \item{message}{A character
#' string giving any additional information returned by [optim()], or
#' `NULL`.} \item{counts}{A two-element integer vector giving the number
#' of calls to `fn` and `gr` respectively.  This excludes those calls
#' needed to compute the Hessian, if requested, and any calls to `fn` to
#' compute a finite-difference approximation to the gradient.}
#' \item{hessian}{If `hessian = TRUE`, the Hessian matrix is returned. Its
#' calculation is on the basis of finite differences (output of `optim`).}
#'
#' @note `RVineMLE` uses the L-BFGS-B method for optimization. \cr If the
#' analytical gradient is used for maximization, computations may be up to 10
#' times faster than using finite differences.
#'
#' @author Ulf Schepsmeier, Jeffrey Dissmann
#'
#' @seealso [RVineSeqEst()],
#' [RVineStructureSelect()],
#' [RVineMatrix()],
#' [RVineGrad()],
#' [RVineHessian()]
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular
#' vine copula models. Computational Statistics, 1-29
#' <https://link.springer.com/article/10.1007/s00180-013-0423-8#>.
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
#' # compute the MLE
#' mle <- RVineMLE(simdata, RVM, grad = TRUE, trace = 0)
#'
#' # compare parameters
#' round(mle$RVM$par - RVM$par, 2)
#'
RVineMLE <- function(data, RVM, start = RVM$par, start2 = RVM$par2, maxit = 200, max.df = 30,
                     max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6),  BB8 = c(6, 1)),
                     grad = FALSE, hessian = FALSE, se = FALSE, ...) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    check_est_pars,
                    check_RVMs,
                    prep_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    if (all(RVM$family == 0)) {
        warning("RVM contains only Independence copulas, RVineMLE does nothing")
        return(RVM)
    }

    ## sanity checks
    if (maxit <= 0)
        stop("'maxit' has to be greater than zero.")

    ## sanity checks for start parameters
    Matrix <- RVM$Matrix
    sel <- lower.tri(Matrix)
    if (!all(start %in% c(0, NA))) {
        BiCopCheck(RVM$family[sel], start[sel], start2[sel])

        if (grad == TRUE) {
            if (any(RVM$family %in% c(7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234))) {
                message("The combination 'grad=TRUE' and a copula family ",
                        "of the vector (7:10,17:20,27:30,37:40,104,114,124,134,204,214,224,234) ",
                        "is not possible. The algorithm will continue with 'grad=FALSE'.")
                grad <- FALSE
            }
        }
    }

    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T

    ## normalization of R-vine matrix
    o <- diag(RVM$Matrix)
    oldRVM <- RVM
    RVM <- normalizeRVineMatrix(RVM)
    if (is.null(colnames(data)))
        colnames(data) <- paste0("V", 1:ncol(data))
    cnms <- colnames(data)
    data <- data[, o[length(o):1]]


    n <- dim(RVM)
    N <- dim(data)[1]

    ## sequential estimation of start parameters if not provided
    if (all(start == 0)) {
        est_start <- RVineSeqEst(data, RVM, max.df = max.df, max.BB = max.BB)
        start <- est_start$par
        start2 <- est_start$par2
    }

    ## Position of parameters in R-vine matrix
    posParams <- (RVM$family > 0)
    posParams2 <- (RVM$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20,
                                     27, 28, 29, 30, 37, 38, 39, 40,
                                     104, 114, 124, 134, 204, 214, 224, 234))
    posParams[is.na(posParams)] <- FALSE
    posParams2[is.na(posParams2)] <- FALSE

    ## number of parameters
    nParams <- sum(posParams, na.rm = TRUE)
    nParams2 <- sum(posParams2, na.rm = TRUE)

    ## vectors of start parameters and corresponding pair-copula families
    startpar <- double(nParams + nParams2)
    Copula.Types <- RVM$family[posParams]
    startpar[1:nParams] <- start[posParams]
    startpar[nParams + seq_len(nParams2)] <- start2[posParams2]

    ## lower and upper bounds
    lb <- double(nParams + nParams2)
    ub <- double(nParams + nParams2)
    for (i in seq_len(nParams)) {
        if (Copula.Types[i] %in% c(1, 2)) {
            # Normal
            lb[i] <- -0.99
            ub[i] <- 0.99
        } else if (Copula.Types[i] %in% c(3, 13)) {
            # clayton
            lb[i] <- 1e-04
            ub[i] <- 100
        } else if (Copula.Types[i] %in% c(23, 33)) {
            # rotated clayton
            lb[i] <- -100
            ub[i] <- -1e-04
        } else if (Copula.Types[i] %in% c(4, 14)) {
            # gumbel
            lb[i] <- 1.0001
            ub[i] <- 50
        } else if (Copula.Types[i] %in% c(24, 34)) {
            # rotated gumbel
            lb[i] <- -50
            ub[i] <- -1.0001
        } else if (Copula.Types[i] == 5) {
            # frank
            lb[i] <- -100
            ub[i] <- 100
        } else if (Copula.Types[i] %in% c(6, 16)) {
            # joe
            lb[i] <- 1.0001
            ub[i] <- 50
        } else if (Copula.Types[i] %in% c(26, 36)) {
            # rotated joe
            lb[i] <- -50
            ub[i] <- -1.0001
        } else if (Copula.Types[i] %in% c(7, 17)) {
            # bb1
            lb[i] <- 0.001
            ub[i] <- max.BB$BB1[1]
        } else if (Copula.Types[i] %in% c(8, 18)) {
            # bb6
            lb[i] <- 1.001
            ub[i] <- max.BB$BB6[1]
        } else if (Copula.Types[i] %in% c(9, 19)) {
            # bb7
            lb[i] <- 1.001
            ub[i] <- max.BB$BB7[1]
        } else if (Copula.Types[i] %in% c(10, 20)) {
            # bb8
            lb[i] <- 1.001
            ub[i] <- max.BB$BB8[1]
        } else if (Copula.Types[i] %in% c(27, 37)) {
            # rotated bb1
            lb[i] <- -max.BB$BB1[1]
            ub[i] <- -0.001
        } else if (Copula.Types[i] %in% c(28, 38)) {
            # rotated bb6
            lb[i] <- -max.BB$BB6[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(29, 39)) {
            # rotated bb7
            lb[i] <- -max.BB$BB7[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(30, 40)) {
            # rotated bb8
            lb[i] <- -max.BB$BB8[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(43, 44)) {
            # double Clayton, Gumbel
            lb[i] <- -0.9999
            ub[i] <- 0.9999
        } else if (Copula.Types[i] %in% c(104, 114, 204, 214)) {
            # Tawn
            lb[i] <- 1.001
            ub[i] <- 20
        } else if (Copula.Types[i] %in% c(124, 134, 224, 234)) {
            # Tawn
            lb[i] <- -20
            ub[i] <- -1.001
        }
    }

    todo <- which(Copula.Types %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234))

    for (i in seq_len(nParams2)) {
        if (Copula.Types[todo[i]] == 2) {
            # t
            lb[nParams + i] <- 2.001
            ub[nParams + i] <- max.df
        } else if (Copula.Types[todo[i]] %in% c(7, 17)) {
            # bb1
            lb[nParams + i] <- 1.001
            ub[nParams + i] <- max.BB$BB1[2]
        } else if (Copula.Types[todo[i]] %in% c(8, 18)) {
            # bb6
            lb[nParams + i] <- 1.001
            ub[nParams + i] <- max.BB$BB6[2]
        } else if (Copula.Types[todo[i]] %in% c(9, 19)) {
            # bb7
            lb[nParams + i] <- 0.001
            ub[nParams + i] <- max.BB$BB7[2]
        } else if (Copula.Types[todo[i]] %in% c(10, 20)) {
            # bb8
            lb[nParams + i] <- 0.001
            ub[nParams + i] <- max.BB$BB8[2]
        } else if (Copula.Types[todo[i]] %in% c(27, 37)) {
            # rotated bb1
            lb[nParams + i] <- -max.BB$BB1[2]
            ub[nParams + i] <- -1.001
        } else if (Copula.Types[todo[i]] %in% c(28, 38)) {
            # rotated bb6
            lb[nParams + i] <- -max.BB$BB6[2]
            ub[nParams + i] <- -1.001
        } else if (Copula.Types[todo[i]] %in% c(29, 39)) {
            # rotated bb7
            lb[nParams + i] <- -max.BB$BB7[2]
            ub[nParams + i] <- -0.001
        } else if (Copula.Types[todo[i]] %in% c(30, 40)) {
            # rotated bb8
            lb[nParams + i] <- -max.BB$BB8[2]
            ub[nParams + i] <- -0.001
        } else if (Copula.Types[todo[i]] %in% c(104, 114, 124, 134, 204, 214, 224, 234)) {
            # Tawn
            lb[nParams + i] <- 0.001
            ub[nParams + i] <- 0.99
        }
    }

    ## log-likelihood function to be maximized
    optim_LL <- function(parm, data, posParams, posParams2, Copula.Types, start, start2, RVM, calcupdate = NA) {

        nParams <- sum(posParams, na.rm = TRUE)
        nParams2 <- sum(posParams2, na.rm = TRUE)

        matrixParams <- start
        matrixParams2 <- start2

        matrixParams[posParams] <- parm[1:nParams]
        matrixParams2[posParams2] <- parm[nParams + seq_len(nParams2)]


        ll <- RVineLogLik(data, RVM, par = matrixParams, par2 = matrixParams2)


        if (is.finite(ll$loglik)) {
            return(ll$loglik)
        } else {
            if (is.na(ll$loglik)) {
                message(parm)
            }
            message(ll$loglik)
            return(-10^305)
        }
    }

    ## gradient
    gradient <- function(parm, data, posParams, posParams2, Copula.Types, start, start2, RVM, calcupdate) {
        nParams <- sum(posParams, na.rm = TRUE)
        nParams2 <- sum(posParams2, na.rm = TRUE)

        matrixParams <- start
        matrixParams2 <- start2

        matrixParams[posParams] <- parm[1:nParams]
        if (nParams2 > 0) {
            matrixParams2[posParams2] <- parm[nParams + seq_len(nParams2)]
        }

        grad <- RVineGrad(data = data,
                          RVM = RVM,
                          par = matrixParams,
                          par2 = matrixParams2,
                          posParams = posParams)$gradient
        return(grad)
    }

    ## default values for parscale (see optim)
    pscale <- numeric()
    for (i in seq_len(nParams)) {
        pscale[i] <- ifelse(Copula.Types[i] %in% c(1, 2, 43, 44), 0.01, 1)
    }
    pscale2 <- numeric()
    for (i in seq_len(nParams2)) {
        pscale2[i] <- ifelse(Copula.Types[i] %in% c(104, 114, 124, 134, 204, 214, 224, 234), 0.05, 1)
    }
    pscale <- c(pscale, pscale2)

    ## (default) values for control parameters of optim
    ctrl <- list(fnscale = -1,
                 maxit = maxit,
                 trace = 1,
                 parscale = pscale,
                 factr = 1e+08)
    ctrl <- modifyList(ctrl, list(...))

    ## optimization
    if (all(Copula.Types %in% c(0, 1, 2, 3:6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 43, 44)) && grad == TRUE) {

        if (hessian == TRUE || se == TRUE) {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          gr = gradient,
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2, Copula.Types = Copula.Types,
                          start = start,
                          start2 = start2,
                          RVM = RVM,
                          method = "L-BFGS-B",
                          control = ctrl,
                          lower = lb,
                          upper = ub,
                          hessian = TRUE)
        } else {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          gr = gradient,
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2,
                          Copula.Types = Copula.Types,
                          start = start,
                          start2 = start2,
                          RVM = RVM,
                          method = "L-BFGS-B",
                          control = ctrl,
                          lower = lb,
                          upper = ub)
        }
    } else {
        if (hessian == TRUE || se == TRUE) {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2,
                          Copula.Types = Copula.Types,
                          start = start,
                          start2 = start2,
                          RVM = RVM,
                          calcupdate = NA,
                          method = "L-BFGS-B",
                          control = ctrl,
                          lower = lb,
                          upper = ub,
                          hessian = TRUE)
        } else {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2,
                          Copula.Types = Copula.Types,
                          start = start,
                          start2 = start2,
                          RVM = RVM,
                          calcupdate = NA,
                          method = "L-BFGS-B",
                          control = ctrl,
                          lower = lb,
                          upper = ub)
        }
    }

    # list for final output
    out <- list()
    # store results in out
    out$value <- out1$value
    out$convergence <- out1$convergence
    out$message <- out1$message
    out$counts <- out1$counts
    if (hessian == TRUE)
        out$hessian <- out1$hessian
    if (se == TRUE)
        out1$se <- sqrt((diag(solve(-out1$hessian))))
    # create parameter matrices
    kk <- 1
    for (ll in seq_len(nParams)) {
        out1$par[ll] <- out1$par[ll]
        if (Copula.Types[ll] %in% c(2, 7:10, 17:20, 27:30, 37:40,
                                    104, 114, 124, 134, 204, 214, 224, 234)) {
            out1$par[nParams + kk] <- out1$par[nParams + kk]
            kk <- kk + 1
        }
    }
    newpar <- newpar2 <- matrix(0, d, d)
    newpar[posParams]  <- out1$par[1:nParams]
    newpar2[posParams2] <- out1$par[nParams + seq_len(nParams2)]
    # create RVineMatrix object
    out$RVM <- RVineMatrix(Matrix = oldRVM$Matrix,
                           family = oldRVM$family,
                           par = newpar,
                           par2 = newpar2,
                           names = oldRVM$names)
    # add standad errors
    if (se == TRUE) {
        out$RVM$se <- matrix(0, d, d)
        out$RVM$se2 <- matrix(0, d, d)
        out$RVM$se[posParams] <- out1$se[1:nParams]
        out$RVM$se2[posParams2] <- out1$se[nParams + seq_len(nParams2)]
    }
    # add summary statistics
    like <- RVineLogLik(data[, cnms], out$RVM)
    out$RVM$logLik <- like$loglik
    out$RVM$pair.logLik <- like$V$value
    npar <- sum(out$RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(out$RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- (out$RVM$family %in% allfams[onepar]) +
        2 * (out$RVM$family %in% allfams[twopar])
    out$RVM$AIC <- -2 * like$loglik + 2 * npar
    out$RVM$pair.AIC <- -2 * like$V$value + 2 * npar_pair
    out$RVM$BIC <- -2 * like$loglik + log(N) * npar
    out$RVM$pair.BIC <- -2 * like$V$value + log(N) * npar_pair
    out$RVM$emptau <- oldRVM$emptau
    out$RVM$p.value.indeptest <- oldRVM$p.value.indeptest

    ## return results
    out
}
