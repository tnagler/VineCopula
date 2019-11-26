#' Hessian Matrix of the Log-Likelihood of an R-Vine Copula Model
#'
#' This function calculates the Hessian matrix of the log-likelihood of a
#' d-dimensional R-vine copula model with respect to the copula parameter and
#' evaluates it on a given copula data set.
#'
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An [RVineMatrix()] object including the structure, the
#' pair-copula families, and the parameters. \cr
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
#'
#' @return \item{hessian}{The calculated Hessian matrix of the log-likelihood
#' value of the R-vine copula model.}
#' \item{der}{The product of the gradient
#' vector with its transposed version.}
#'
#' @note The Hessian matrix is not available for R-vine copula models with two
#' parameter Archimedean copulas, i.e. BB1, BB6, BB7, BB8 and their rotated
#' versions.
#'
#' @author Ulf Schepsmeier, Jakob Stoeber
#'
#' @seealso
#' [BiCopDeriv()],
#' [BiCopDeriv2()],
#' [BiCopHfuncDeriv()],
#' [BiCopHfuncDeriv2()], \cr
#' [RVineMatrix()],
#' [RVineMLE()],
#' [RVineGrad()]
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
#' # compute the Hessian matrix of the first row of the data
#' out2 <- RVineHessian(simdata[1,], RVM)
#' out2$hessian
#'
RVineHessian <- function(data, RVM) {
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

    n <- d <- args$d
    N <- T

    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- getFromNamespace("normalizeRVineMatrix", "VineCopula")(RVM)
        data <- data[, o[length(o):1]]
    }

    dd <- d * (d - 1)/2
    tt <- sum(RVM$family == 2)
    hessian <- matrix(0, dd + tt, dd + tt)
    subhess <- matrix(0, dd + tt, dd + tt)
    der <- matrix(0, dd + tt, dd + tt)
    subder <- matrix(0, dd + tt, dd + tt)

    out <- .C("hesse",
              as.integer(T),
              as.integer(d),
              as.integer(as.vector(RVM$family)),
              as.integer(as.vector(RVM$MaxMat)),
              as.integer(as.vector(RVM$Matrix)),
              as.integer(as.vector(RVM$CondDistr$direct)),
              as.integer(as.vector(RVM$CondDistr$indirect)),
              as.double(as.vector(RVM$par)),
              as.double(as.vector(RVM$par2)),
              as.double(as.vector(data)),
              as.double(as.vector(hessian)),
              as.double(as.vector(subhess)),
              as.double(as.vector(der)),
              as.double(as.vector(subder)),
              PACKAGE = 'VineCopula')

    hessian <- matrix(out[[11]], dd + tt, dd + tt)
    subhess <- matrix(out[[12]], dd + tt, dd + tt)
    der <- matrix(out[[13]], dd + tt, dd + tt)
    subder <- matrix(out[[14]], dd + tt, dd + tt)


    test <- apply(hessian, 2, function(x) max(abs(x)))
    hessian <- hessian[test > 0, test > 0]
    subhess <- subhess[test > 0, test > 0]
    der <- der[test > 0, test > 0]
    subder <- subder[test > 0, test > 0]


    out <- list(hessian = hessian, der = der)

    return(out)
}
