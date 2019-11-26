#' Gradient of the Log-Likelihood of an R-Vine Copula Model
#'
#' This function calculates the gradient of the log-likelihood of a
#' d-dimensional R-vine copula model with respect to the copula parameter and
#' evaluates it on a given copula data set.
#'
#' The ordering of the gradient is due to the ordering of the R-vine matrix.
#' The gradient starts at the lower right corner of the R-vine matrix and goes
#' column by column to the left and up, i.e. the first entry of the gradient is
#' the last entry of the second last column of the `par`-matrix followed
#' by the last entry of the third last column and the second last entry of this
#' column. If there is a copula family with two parameters, i.e. the t-copula,
#' the derivative with respect to the second parameter is at the end of the
#' gradient vector in order of their occurrence.
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An [RVineMatrix()] object including the structure and
#' the pair-copula families and parameters. \cr
#' Only the following copula
#' families are allowed in `RVM$family` \cr
#' `0` = independence copula \cr
#' `1` = Gaussian copula \cr
#' `2` = Student t copula (t-copula)\cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `13` = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' `14` = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' `16` = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' `23` = rotated Clayton copula (90 degrees) \cr
#' `24` = rotated Gumbel copula (90 degrees) \cr
#' `26` = rotated Joe copula (90 degrees) \cr
#' `33` = rotated Clayton copula (270 degrees) \cr
#' `34` = rotated Gumbel copula (270 degrees) \cr
#' `36` = rotated Joe copula (270 degrees) \cr
#' @param par A d x d matrix with the pair-copula parameters (optional;
#' default: `par = RVM$par`).
#' @param par2 A d x d matrix with the second parameters of pair-copula
#' families with two parameters (optional; default: `par2 = RVM$par2`).
#' @param start.V Transformations (h-functions and log-likelihoods of each
#' pair-copula) of previous calculations (see output; default: `start.V =
#' NA`).
#' @param posParams A d x d matrix indicating which copula has to be considered
#' in the gradient (default: `posParams = (RVM$family > 0)`).
#'
#' @return gradient The calculated gradient of the log-likelihood value
#' of the R-vine copula model. (three matrices: `direct`, `indirect`
#' and `value`).
#'
#' @note The gradient for R-vine copula models with two parameter Archimedean
#' copulas, i.e. BB1, BB6, BB7, BB8 and their rotated versions can not yet be calculated.
#' The derivatives of these bivariate copulas are more complicated.
#'
#' @author Ulf Schepsmeier, Jakob Stoeber
#'
#' @seealso [BiCopDeriv()],
#' [BiCopDeriv2()],
#' [BiCopHfuncDeriv()],
#' [BiCopHfuncDeriv2()], \cr
#' [RVineMatrix()],
#' [RVineMLE()],
#' [RVineHessian()]
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' Schepsmeier, U. and J. Stoeber (2014)
#' Derivatives and Fisher information of bivariate copulas.
#' Statistical Papers, 55(2), 525-542.
#' online first: <http://link.springer.com/article/10.1007/s00362-013-0498-x>.
#'
#' Web supplement: Derivatives and Fisher Information of bivariate copulas.
#' <http://mediatum.ub.tum.de/node?id=1119201>
#'
#' Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular
#' vine copula models. Computational Statistics, 28 (6), 2679-2707
#' <http://link.springer.com/article/10.1007/s00180-013-0423-8#>.
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
#' # compute the gradient of the first row of the data
#' out2 <- RVineGrad(simdata[1,], RVM)
#' out2$gradient
#'
RVineGrad <- function(data, RVM, par = RVM$par, par2 = RVM$par2, start.V = NA, posParams = (RVM$family > 0)) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36))))
        stop("Copula family not implemented.")

    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T

    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }

    if (any(is.na(start.V))) {
        loglik <- RVineLogLik(data, RVM, par = par, par2 = par2, separate = TRUE)
        V <- loglik$V

    } else {
        V <- start.V
        V$value[V$value %in% c(NA, NaN, -Inf)] <- -1e+10
        if (any(is.na(V$value)))
            message("NA in LogL call")
    }


    ll <- as.vector(V$value)
    vv <- as.vector(V$direct)
    vv2 <- as.vector(V$indirect)

    w1 <- as.vector(RVM$family)
    w1[is.na(w1)] <- 0
    th <- as.vector(par)
    th[is.na(th)] <- 0
    th2 <- as.vector(par2)
    th2[is.na(th2)] <- 0
    condirect <- as.vector(as.numeric(RVM$CondDistr$direct))
    conindirect <- as.vector(as.numeric(RVM$CondDistr$indirect))
    maxmat <- as.vector(RVM$MaxMat)
    matri <- as.vector(RVM$Matrix)
    matri[is.na(matri)] <- 0
    maxmat[is.na(maxmat)] <- 0
    condirect[is.na(condirect)] <- 0
    conindirect[is.na(conindirect)] <- 0


    out <- rep(0, sum(posParams[lower.tri(posParams, diag = FALSE)]) + sum(w1 == 2))

    out <- .C("VineLogLikRvineGradient",
              as.integer(T),
              as.integer(d),
              as.integer(w1),
              as.integer(maxmat),
              as.integer(matri),
              as.integer(condirect),
              as.integer(conindirect),
              as.double(th),
              as.double(th2),
              as.double(data),
              as.double(out),
              as.double(ll),
              as.double(vv),
              as.double(vv2),
              as.integer(as.vector(posParams)),
              PACKAGE = 'VineCopula')



    gradient2 <- out[[11]]
    gradient2[gradient2 %in% c(NA, NaN, -Inf)] <- -1e+10

    dd <- sum(RVM$family > 0)
    tt <- sum(w1 == 2)
    grad1 <- gradient2[1:dd]
    gradient <- grad1[dd:1]
    if (tt > 0) {
        grad2 <- gradient2[(dd + 1):(dd + tt)]
        gradient <- c(gradient, grad2[tt:1])
    }

    out2 <- list(gradient = gradient)
    return(out2)
}
