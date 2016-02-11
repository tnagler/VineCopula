#' Inverse Conditional Distribution Function of a Bivariate Copula
#'
#' This function evaluates the inverse conditional distribution function
#' (inverse h-function) of a given parametric bivariate copula.
#'
#' The h-function is defined as the conditional distribution function of a
#' bivariate copula, i.e.,
#' \deqn{h(u|v,\boldsymbol{\theta}) := F(u|v) = \frac{\partial C(u,v)}{\partial v}, }{
#' h(u|v,\theta) := F(u|v) = \partial C(u,v) / \partial v, }
#' where \eqn{C} is a bivariate copula distribution function with parameter(s)
#' \eqn{\boldsymbol{\theta}}{\theta}. For more details see Aas et al. (2009). \cr \cr
#'
#' If the family and parameter specification is stored in a \code{\link{BiCop}}
#' object \code{obj}, the alternative version
#' \preformatted{BiCopHinv(u1, u2, obj)}
#' can be used.
#'
#' @param u1,u2 numeric vectors of equal length with values in [0,1].
#' @param family integer; single number or vector of size \code{length(u1)};
#' defines the bivariate copula family: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par numeric; single number or vector of size \code{length(u1)};
#' copula parameter.
#' @param par2 numeric; single number or vector of size \code{length(u1)};
#' second parameter for bivariate copulas with two parameters (t, BB1, BB6,
#' BB7, BB8, Tawn type 1 and type 2; default: \code{par2 = 0}). \code{par2}
#' should be an positive integer for the Students's t copula \code{family = 2}.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#'
#' @return \item{hinv1}{Numeric vector of the inverse conditional distribution
#' function (inverse h-function) of the copula \code{family} with parameter(s)
#' \code{par}, \code{par2} evaluated at \code{u2} given \code{u1}, i.e.,
#' \eqn{h^{-1}(\code{u2}|\code{u1},\boldsymbol{\theta})}{h^{-1}(u2|u1,\theta)}.}
#' \item{hinv2}{Numeric vector of the inverse conditional distribution function
#' (inverse h-function) of the copula \code{family} with parameter(s)
#' \code{par}, \code{par2}evaluated at \code{u1} given \code{u2}, i.e.,
#' \eqn{h^{-1}(\code{u1}|\code{u2},\boldsymbol{\theta})}{h^{-1}(u1|u2,\theta)}.}
#'
#' @author Thomas Nagler
#'
#' @seealso \code{\link{BiCopPDF}}, \code{\link{BiCopCDF}},
#' \code{\link{RVineLogLik}}, \code{\link{RVineSeqEst}}, \code{\link{BiCop}}
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
#' Pair-copula constructions of multiple dependence. Insurance: Mathematics and
#' Economics 44 (2), 182-198.
#'
#' @examples
#'
#' # inverse h-functions of the Gaussian copula
#' cop <- BiCop(1, 0.5)
#' h1 <- BiCopHinv(0.1, 0.2, cop)
#' \dontshow{
#' h1
#' }
#' ## check if it is actually the inverse
#' cop <- BiCop(3, 3)
#' all.equal(0.2, BiCopHfunc(0.1, BiCopHinv(0.1, 0.2, cop)$hinv1, cop)$hfunc1)
#' all.equal(0.1, BiCopHfunc(BiCopHinv(0.1, 0.2, cop)$hinv2, 0.2, cop)$hfunc2)
#'
#'
#' @export BiCopHinv
BiCopHinv <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE)
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2))
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(c(u1, u2) > 1) || any(c(u1, u2) < 0))
        stop("Data has be in the interval [0,1].")
    n <- length(u1)

    ## extract family and parameters if BiCop object is provided
    if (missing(family))
        family <- NA
    if (missing(par))
        par <- NA
    # for short hand usage extract obj from family
    if (class(family) == "BiCop")
        obj <- family
    if (!is.null(obj)) {
        stopifnot(class(obj) == "BiCop")
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }

    ## adjust length for parameter vectors; stop if not matching
    if (missing(par) & (all(family == 0)))
        par <- 0
    if (any(c(length(family), length(par), length(par2)) == n)) {
        if (length(family) == 1)
            family <- rep(family, n)
        if (length(par) == 1)
            par <- rep(par, n)
        if (length(par2) == 1)
            par2 <- rep(par2, n)
    }
    if (!(length(family) %in% c(1, n)))
        stop("'family' has to be a single number or a size n vector")
    if (!(length(par) %in% c(1, n)))
        stop("'par' has to be a single number or a size n vector")
    if (!(length(par2) %in% c(1, n)))
        stop("'par2' has to be a single number or a size n vector")

    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }

    ## calculate inverse h-functions
    if (length(par) == 1) {
        # call for single parameters
        hinv1 <- .C("Hinv1",                      # h(u2|u1)
                    as.integer(family),
                    as.integer(n),
                    as.double(u2),
                    as.double(u1),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2",                      # h(u1|u2)
                    as.integer(family),
                    as.integer(n),
                    as.double(u1),
                    as.double(u2),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        hinv1 <- .C("Hinv1_vec",                      # h(u2|u1)
                    as.integer(family),
                    as.integer(n),
                    as.double(u2),
                    as.double(u1),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2_vec",                      # h(u1|u2)
                    as.integer(family),
                    as.integer(n),
                    as.double(u1),
                    as.double(u2),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
    }

    ## return results
    list(hinv1 = hinv1, hinv2 = hinv2)
}

