#' Sequential Pair-Copula Selection and Estimation for R-Vine Copula Models
#'
#' This function fits a R-vine copula model to a d-dimensional copula data set.
#' Pair-copula families are selected using \code{\link{BiCopSelect}} and
#' estimated sequentially.
#'
#' R-vine copula models with unknown structure can be specified using
#' \code{\link{RVineStructureSelect}}.
#'
#' @param data N x d data matrix (with uniform margins).
#' @param familyset integer vector of pair-copula families to select from.
#' The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases.  If \code{familyset = NA} (default), selection among all
#' possible families is performed. If a vector of negative numbers is provided,
#' selection among all but \code{abs(familyset)} is performed. Coding of
#' pair copula families is the same as in \code{\link{BiCop}}.
#' @param Matrix lower or upper triangular d x d matrix that defines the R-vine
#' tree structure.
#' @param selectioncrit Character indicating the criterion for pair-copula
#' selection. Possible choices: \code{selectioncrit = "AIC"} (default),
#' \code{"BIC"}, or \code{"logLik"} (see \code{\link{BiCopSelect}}).
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).  The
#' independence copula is chosen for a (conditional) pair if the null
#' hypothesis of independence cannot be rejected.
#' @param level numeric; significance level of the independence test (default:
#' \code{level = 0.05}).
#' @param trunclevel integer; level of truncation.
#' @param rotations logical; if \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#' @param cores integer; if \code{cores > 1}, estimation will be parallized
#' within each tree (using \code{\link[foreach]{foreach}}).
#'
#' @return An \code{\link{RVineMatrix}} object with the selected families
#' (\code{RVM$family}) as well as sequentially
#' estimated parameters stored in \code{RVM$par} and \code{RVM$par2}. The object
#' is augmented by the following information about the fit:
#' \item{se, se2}{standard errors for the parameter estimates  (if
#' \code{se = TRUE}; note that these are only approximate since they do not
#' account for the sequential nature of the estimation,}
#' \item{nobs}{number of observations,}
#' \item{logLik, pair.logLik}{log likelihood (overall and pairwise)}
#' \item{AIC, pair.AIC}{Aikaike's Informaton Criterion (overall and pairwise),}
#' \item{BIC, pair.BIC}{Bayesian's Informaton Criterion (overall and pairwise),}
#' \item{emptau}{matrix of empirical values of Kendall's tau,}
#' \item{p.value.indeptest}{matrix of p-values of the independence test.}#'
#'
#' @note For a comprehensive summary of the vine copula model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}.
#'
#' @author Eike Brechmann, Thomas Nagler
#'
#' @seealso
#' \code{\link{RVineMatrix}},
#' \code{\link{BiCop}},
#' \code{\link{BiCopSelect}},
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
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' ## define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' ## simulate a sample of size 1000 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(1000, RVM)
#'
#' ## determine the pair-copula families and parameters
#' RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4, 5 ,6), Matrix)
#'
#' ## see the object's content or a summary
#' str(RVM1)
#' summary(RVM1)
#'
#' ## inspect the fitted model using plots
#' \dontrun{
#' plot(RVM1)  # tree structure
#' }
#' contour(RVM1)  # contour plots of all pair-copulas
#'
RVineCopSelect <- function(data, familyset = NA, Matrix, selectioncrit = "AIC", indeptest = FALSE,
                           level = 0.05, trunclevel = NA, rotations = TRUE, cores = 1) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    check_if_01,
                    prep_familyset,
                    check_matrix)
    list2env(args, environment())
    warning(" In ", args$call[1], ": ",
            "Some of the data are NA. ",
            "Only pairwise complete observations are used.",
            call. = FALSE)

    ## sanity checks
    if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
        stop("Selection criterion not implemented.")
    if (level < 0 & level > 1)
        stop("Significance level has to be between 0 and 1.")

    d <- n <- ncol(data)
    N <- nrow(data)
    ## set variable names and trunclevel
    varnames <- colnames(data)
    if (is.na(trunclevel))
        trunclevel <- d

    ## adjust familyset
    types <- familyset
    if (trunclevel == 0)
        types <- 0

    ## reorder matrix to natural order
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
    Ses     <- matrix(0, d, d)
    Se2s    <- matrix(0, d, d)
    emptaus <- matrix(0, d, d)
    pvals   <- matrix(0, d, d)
    nobs    <- matrix(0, d, d)
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
    warn <- NULL
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

                ## select pair-copula
                if (trunclevel <= (d-k))
                    familyset <- 0

                na.ind <- which(is.na(zr1 + zr2))
                if (length(na.ind) >= length(zr1) - 1) {
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
                    warn <- paste("Insufficient data for at least one pair.",
                                  "Independence has been selected automatically.")
                } else {
                    cfit <- suppressWarnings(BiCopSelect(zr2,
                                                         zr1,
                                                         familyset,
                                                         selectioncrit,
                                                         indeptest,
                                                         level,
                                                         weights = NA,
                                                         rotations,
                                                         se = TRUE))
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

        ## run pair-copula selection for tree k
        res.k <- if (cores > 1) {
            foreach(i = 1:(k-1),
                    .packages = c("VineCopula"),
                    .export = "familyset") %dopar% doEst(i)
        } else {
            lapply(1:(k-1), doEst)
        }

        for (i in 1:(k-1)) {
            ## store info about selected pair-copula in matrices
            Types[k, i]   <- res.k[[i]]$cfit$family
            Params[k, i]  <- res.k[[i]]$cfit$par
            Params2[k, i] <- res.k[[i]]$cfit$par2
            Ses[k, i]     <- res.k[[i]]$cfit$se
            tmpse2        <- res.k[[i]]$cfit$se2
            Se2s[k, i]    <- ifelse(is.null(tmpse2), NA, tmpse2)
            emptaus[k, i] <- res.k[[i]]$cfit$emptau
            pvals[k, i]   <- res.k[[i]]$cfit$p.value.indeptest
            if (!is.null(res.k[[i]]$warn))
                warn <- res.k[[i]]$warn
            ## replace pseudo observations for estimation of next tree
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
    .RVM$se <- Ses
    .RVM$se2 <- Se2s
    .RVM$nobs <- N
    revo <- sapply(1:d, function(i) which(o[length(o):1] == i))
    like <- suppressWarnings(RVineLogLik(data[, revo], .RVM))
    .RVM$logLik <- like$loglik
    .RVM$pair.logLik <- like$V$value
    npar <- sum(.RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(.RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- (.RVM$family %in% allfams[onepar]) +
        2 * (.RVM$family %in% allfams[twopar])
    .RVM$AIC <- -2 * like$loglik + 2 * npar
    .RVM$pair.AIC <- -2 * like$V$value + 2 * npar_pair
    .RVM$BIC <- -2 * like$loglik + log(N) * npar
    .RVM$pair.BIC <- -2 * like$V$value + log(N) * npar_pair
    .RVM$emptau <- emptaus
    .RVM$p.value.indeptest <- pvals

    ## free memory and return final object
    rm(list = ls())
    .RVM
}
