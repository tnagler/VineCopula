#' Sequential Estimation of an R-Vine Copula Model
#'
#' This function sequentially estimates the pair-copula parameters of a
#' d-dimensional R-vine copula model as specified by the corresponding
#' [RVineMatrix()] object.
#'
#' The pair-copula parameter estimation is performed tree-wise, i.e., for each
#' R-vine tree the results from the previous tree(s) are used to calculate the
#' new copula parameters using [BiCopEst()].
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An [RVineMatrix()] object including the structure, the
#' pair-copula families and the pair-copula parameters (if they are known).
#' @param method indicates the estimation method: either maximum
#' likelihood estimation (`method = "mle"`; default) or inversion of
#' Kendall's tau (`method = "itau"`). For `method = "itau"` only
#' one parameter families and the Student t copula can be used (`family =
#' 1,2,3,4,5,6,13,14,16,23,24,26,33,34` or `36`). For the t-copula,
#' `par2` is found by a crude profile likelihood optimization over the
#' interval (2, 10].
#' @param se Logical; whether standard errors are estimated (default: `se
#' = FALSE`).
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: `max.df = 30`; for more
#' details see [BiCopEst()]).
#' @param max.BB List; upper bounds for the estimation of the two parameters
#' (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr (default:
#' `max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))`).
#' @param progress Logical; whether the pairwise estimation progress is printed
#' (default: `progress = FALSE`).
#' @param weights Numerical; weights for each observation (optional).
#' @param cores integer; if `cores > 1`, estimation will be parallelized
#' within each tree (using [foreach::foreach()]). However, the
#' overhead caused by parallelization is likely to make the function run slower
#' unless sample size is really large and `method = "itau"`.
#'
#' @return An [RVineMatrix()] object with the sequentially
#' estimated parameters stored in `RVM$par` and `RVM$par2`. The object
#' is augmented by the following information about the fit:
#' \item{se, se2}{standard errors for the parameter estimates (if
#' `se = TRUE`); note that these are only approximate since they do not
#' account for the sequential nature of the estimation,}
#' \item{nobs}{number of observations,}
#' \item{logLik, pair.logLik}{log likelihood (overall and pairwise)}
#' \item{AIC, pair.AIC}{Aikaike's Informaton Criterion (overall and pairwise),}
#' \item{BIC, pair.BIC}{Bayesian's Informaton Criterion (overall and pairwise),}
#' \item{emptau}{matrix of empirical values of Kendall's tau,}
#' \item{p.value.indeptest}{matrix of p-values of the independence test.}
#'
#' @note For a comprehensive summary of the fitted model, use
#' `summary(object)`; to see all its contents, use `str(object)`.
#'
#' @author Ulf Schepsmeier, Jeffrey Dissmann, Thomas Nagler
#'
#' @seealso
#' [RVineMatrix()],
#' [BiCop()],
#' [BiCopEst()],
#' [plot.RVineMatrix()],
#' [contour.RVineMatrix()]
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
#' summary(RVineSeqEst(simdata, RVM, method = "itau", se = TRUE))
#' summary(RVineSeqEst(simdata, RVM, method = "mle", se = TRUE))
#'
RVineSeqEst <- function(data, RVM, method = "mle", se = FALSE, max.df = 30,
                        max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)),
                        progress = FALSE, weights = NA, cores = 1) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    check_if_01,
                    check_est_pars,
                    check_RVMs,
                    prep_RVMs)
    list2env(args, environment())
    if (any(is.na(data)))
        warning(" In ", args$call[1], ": ",
                "Some of the data are NA. ",
                "Only pairwise complete observations are used.",
                call. = FALSE)

    d <- n <- ncol(data)
    N <- nrow(data)
    varnames <- colnames(data)

    ## reorder matrix to natural order
    Matrix <- RVM$Matrix
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
    emptaus <- matrix(0, d, d)
    logLiks <- matrix(0, d, d)
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
            cl <- makePSOCKcluster(cores)
            setDefaultCluster(cl)
            on.exit(try(stopCluster(cl), silent = TRUE))
        }
    }

    ## loop over all trees and pair-copulas up to a truncation level
    warn <- NULL
    trc_lvl <- tail(which(apply(RVM$family, 1, function(x) all(x == 0))), 1)
    for(k in d:(trc_lvl+1)) {
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

                na.ind <- which(is.na(zr1 + zr2))
                if (length(na.ind) >= length(zr1) - 10) {
                    cfit <- BiCop(0)
                    ## add more information about the fit
                    cfit$se  <- NA
                    cfit$se2 <- NA
                    cfit$nobs   <- 0
                    cfit$logLik <- 0
                    cfit$AIC    <- 0
                    cfit$BIC    <- 0
                    cfit$emptau <- NA
                    cfit$p.value.indeptest <- NA
                    warn <- paste("Insufficient data for at least one pair;",
                                  "family has been set to 0.")
                } else {
                    ## select pair-copula
                    cfit <- suppressWarnings(BiCopEst(zr2,
                                                      zr1,
                                                      RVM$family[k, i],
                                                      method,
                                                      se,
                                                      max.df,
                                                      max.BB,
                                                      weights))
                    warn <- NULL
                }

                ## transform data to pseudo-oberstavions in next tree
                direct <- indirect <- NULL
                if (CondDistr$direct[k - 1, i])
                    direct <- suppressWarnings(BiCopHfunc1(zr2,
                                                           zr1,
                                                           cfit,
                                                           check.pars = FALSE))
                if (CondDistr$indirect[k - 1, i])
                    indirect <- suppressWarnings(BiCopHfunc2(zr2,
                                                             zr1,
                                                             cfit,
                                                             check.pars = FALSE))

                ## return results
                list(direct = direct,
                     indirect = indirect,
                     cfit = cfit,
                     warn = warn)
            } else {
                list(cfit = BiCop(0, 0))
            }
        }

        ## run pair-copula estimation for tree k
        if (cores > 1)
            lapply <- function(...) parallel::parLapply(getDefaultCluster(), ...)
        res.k <- lapply(seq_len(k - 1), doEst)

        for (i in 1:(k-1)) {
            ## store results in matrices
            Types[k, i]   <- res.k[[i]]$cfit$family
            Params[k, i]  <- res.k[[i]]$cfit$par
            Params2[k, i] <- res.k[[i]]$cfit$par2
            emptaus[k, i] <- res.k[[i]]$cfit$emptau
            logLiks[k, i] <- res.k[[i]]$cfit$logLik
            if (!is.null(res.k[[i]]$warn))
                warn <- res.k[[i]]$warn

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
    if (!is.null(warn))
        warning(" In ", args$call[1], ": ", warn, call. = FALSE)

    ## store results in RVineMatrix object
    .RVM <- RVineMatrix(Mold,
                        family = Types,
                        par = Params,
                        par2 = Params2,
                        names = varnames)
    if (se) {
        .RVM$se <- Se
        .RVM$se2 <- Se2
    }
    .RVM$nobs <- N
    revo <- sapply(1:d, function(i) which(o[length(o):1] == i))
    like <- suppressWarnings(RVineLogLik(data[, revo], .RVM, calculate.V = TRUE))
    .RVM$logLik <- like$loglik
    .RVM$pair.logLik <- logLiks
    npar <- sum(.RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(.RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- (.RVM$family %in% allfams[onepar]) +
        2 * (.RVM$family %in% allfams[twopar])
    .RVM$AIC <- -2 * like$loglik + 2 * npar
    .RVM$pair.AIC <- -2 * like$V$value + 2 * npar_pair
    .RVM$BIC <- -2 * like$loglik + log(N) * npar
    .RVM$pair.BIC <- -2 * like$V$value + log(N) * npar_pair
    .RVM$emptau <- emptaus

    ## free memory and return results
    rm(list = ls())
    .RVM
}
