#' Standard Errors of an R-Vine Copula Model
#'
#' This function calculates the standard errors of a d-dimensional R-vine
#' copula model given the Hessian matrix.
#'
#'
#' @param hessian The Hessian matrix of the given R-vine.
#' @param RVM An \code{\link{RVineMatrix}} object including the structure, the
#' pair-copula families, and the parameters.
#'
#' @return \item{se}{The calculated standard errors for the first parameter
#' matrix. The entries are ordered with respect to the ordering of the
#' \code{RVM$par} matrix.} \item{se2}{The calculated standard errors for the
#' second parameter matrix.}
#'
#' @note The negative Hessian matrix should be positive semidefinite. Otherwise
#' NAs will be returned in some entries and the non-NA entries may be wrong. If
#' the negaive Hessian matrix is negative definite, then one could try a near
#' positive matrix. The package \code{Matrix} provides a function called
#' \code{nearPD} to estimate a matrix which is positive definite and close to
#' the given matrix.
#'
#' @author Ulf Schepsmeier, Jakob Stoeber
#'
#' @seealso
#' \code{\link{BiCopDeriv}},
#' \code{\link{BiCopDeriv2}},
#' \code{\link{BiCopHfuncDeriv}},
#' \code{\link{BiCopHfuncDeriv2}}, \cr
#' \code{\link{RVineMatrix}},
#' \code{\link{RVineHessian}},
#' \code{\link{RVineGrad}}
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' Schepsmeier, U. and J. Stoeber (2014)
#' Derivatives and Fisher information of bivariate copulas.
#' Statistical Papers, 55(2), 525-542.
#' online first: \url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.
#'
#' Web supplement: Derivatives and Fisher Information of bivariate copulas.
#' \url{http://mediatum.ub.tum.de/node?id=1119201}
#'
#' Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular
#' vine copula models. Computational Statistics, 28 (6), 2679-2707
#' \url{http://link.springer.com/article/10.1007/s00180-013-0423-8#}.
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
#' out2 <- RVineHessian(simdata,RVM)
#'
#' # get the standard errors
#' RVineStdError(out2$hessian, RVM)
#'
RVineStdError <- function(hessian, RVM) {
    # Test auf pos. semidef.
    se3 <- numeric()
    a <- eigen(-hessian, only.values = TRUE)
    if (any(a$values < 0)) {
        warning("The negative Hessian matrix is not positive definite. Thus NAs will be returned in some entries.")
    }
    se <- sqrt((diag(solve(-hessian))))

    d <- dim(RVM$family)[1]
    posParams <- (RVM$family > 0)
    posParams2 <- (RVM$family == 2)

    posParams[is.na(posParams)] <- FALSE
    posParams2[is.na(posParams2)] <- FALSE

    nParams <- sum(posParams, na.rm = TRUE)
    nParams2 <- sum(posParams2, na.rm = TRUE)

    d2 <- dim(hessian)[1]

    if (nParams < d2) {
        # t-copula involved
        se2 <- se[(nParams + 1):d2]
        se <- se[1:nParams]
    }

    SE <- matrix(0, d, d)
    t <- 1
    for (i in (d - 1):1) {
        for (j in d:(i + 1)) {
            if (posParams[j, i]) {
                SE[j, i] <- se[t]
                t <- t + 1
            }
        }
    }

    t <- 1
    if (nParams < d2) {
        SE2 <- matrix(0, d, d)
        for (i in (d - 1):1) {
            for (j in d:(i + 1)) {
                if (RVM$family[j, i] == 2) {
                    SE2[j, i] <- se2[t]
                    t <- t + 1
                }
            }
        }
    }

    out <- list()
    out$se <- SE
    if (nParams < d2) out$se2 <- SE2

    return(out)
}
