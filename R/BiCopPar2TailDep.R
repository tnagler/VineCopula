#' Tail Dependence Coefficients of a Bivariate Copula
#' 
#' This function computes the theoretical tail dependence coefficients of a
#' bivariate copula for given parameter values.
#' 
#' If the family and parameter specification is stored in a \code{BiCop} object
#' \code{obj}, the alternative version \cr \preformatted{BiCopPar2TailDep(obj)}
#' can be used.
#' 
#' @param family integer; single number or vector of size \code{m}; defines the
#' bivariate copula family: \cr \code{0} = independence copula \cr \code{1} =
#' Gaussian copula \cr \code{2} = Student t copula (t-copula) \cr \code{3} =
#' Clayton copula \cr \code{4} = Gumbel copula \cr \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr \code{7} = BB1 copula \cr \code{8} = BB6 copula
#' \cr \code{9} = BB7 copula \cr \code{10} = BB8 copula \cr \code{13} = rotated
#' Clayton copula (180 degrees; ``survival Clayton'') \cr \code{14} = rotated
#' Gumbel copula (180 degrees; ``survival Gumbel'') \cr \code{16} = rotated Joe
#' copula (180 degrees; ``survival Joe'') \cr \code{17} = rotated BB1 copula
#' (180 degrees; ``survival BB1'')\cr \code{18} = rotated BB6 copula (180
#' degrees; ``survival BB6'')\cr \code{19} = rotated BB7 copula (180 degrees;
#' ``survival BB7'')\cr \code{20} = rotated BB8 copula (180 degrees; ``survival
#' BB8'')\cr \code{23} = rotated Clayton copula (90 degrees) \cr \code{24} =
#' rotated Gumbel copula (90 degrees) \cr \code{26} = rotated Joe copula (90
#' degrees) \cr \code{27} = rotated BB1 copula (90 degrees) \cr \code{28} =
#' rotated BB6 copula (90 degrees) \cr \code{29} = rotated BB7 copula (90
#' degrees) \cr \code{30} = rotated BB8 copula (90 degrees) \cr \code{33} =
#' rotated Clayton copula (270 degrees) \cr \code{34} = rotated Gumbel copula
#' (270 degrees) \cr \code{36} = rotated Joe copula (270 degrees) \cr \code{37}
#' = rotated BB1 copula (270 degrees) \cr \code{38} = rotated BB6 copula (270
#' degrees) \cr \code{39} = rotated BB7 copula (270 degrees) \cr \code{40} =
#' rotated BB8 copula (270 degrees) \cr \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr \code{124} =
#' rotated Tawn type 1 copula (90 degrees) \cr \code{134} = rotated Tawn type 1
#' copula (270 degrees) \cr \code{204} = Tawn type 2 copula \cr \code{214} =
#' rotated Tawn type 2 copula (180 degrees) \cr \code{224} = rotated Tawn type
#' 2 copula (90 degrees) \cr \code{234} = rotated Tawn type 2 copula (270
#' degrees) \cr
#' @param par numeric; single number or vector of size \code{m}; copula
#' parameter.
#' @param par2 numeric; single number or vector of size \code{m}; second
#' parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8,
#' Tawn type 1 and type 2; default: \code{par2 = 0}). \code{par2} should be an
#' positive integer for the Students's t copula \code{family = 2}.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#' @return \item{lower}{Lower tail dependence coefficient for the given
#' bivariate copula \code{family} and parameter(s) \code{par}, \code{par2}:
#' \deqn{ }{ \lambda_L = lim_{u->0} C(u,u)/u }\deqn{ \lambda_L =
#' \lim_{u\searrow 0}\frac{C(u,u)}{u} }{ \lambda_L = lim_{u->0} C(u,u)/u
#' }\deqn{ }{ \lambda_L = lim_{u->0} C(u,u)/u } } \item{upper}{Upper tail
#' dependence coefficient for the given bivariate copula family \code{family}
#' and parameter(s) \code{par}, \code{par2}: \deqn{ }{ \lambda_U =
#' lim_{u->1}(1-2u+C(u,u))/(1-u) }\deqn{ \lambda_U = \lim_{u\nearrow
#' 1}\frac{1-2u+C(u,u)}{1-u} }{ \lambda_U = lim_{u->1}(1-2u+C(u,u))/(1-u)
#' }\deqn{ }{ \lambda_U = lim_{u->1}(1-2u+C(u,u))/(1-u) } } Lower and upper
#' tail dependence coefficients for bivariate copula families and parameters
#' (\eqn{\theta} for one parameter families and the first parameter of the
#' t-copula with \eqn{\nu} degrees of freedom, \eqn{\theta} and \eqn{\delta}
#' for the two parameter BB1, BB6, BB7 and BB8 copulas) are given in the
#' following table.  \tabular{lll}{ No. \tab Lower tail dependence \tab Upper
#' tail dependence \cr \code{1} \tab - \tab - \cr \code{2} \tab
#' \eqn{2t_{\nu+1}\left(-\sqrt{\nu+1}\sqrt{\frac{1-\theta}{1+\theta}}\right)}{2t_{\nu+1}(-\sqrt{\nu+1}\sqrt{(1-\theta)/(1+\theta)})}
#' \tab
#' \eqn{2t_{\nu+1}\left(-\sqrt{\nu+1}\sqrt{\frac{1-\theta}{1+\theta}}\right)}{2t_{\nu+1}(-\sqrt{\nu+1}\sqrt{(1-\theta)/(1+\theta)})}
#' \cr \code{3} \tab \eqn{2^{-1/\theta}} \tab - \cr \code{4} \tab - \tab
#' \eqn{2-2^{1/\theta}} \cr \code{5} \tab - \tab - \cr \code{6} \tab - \tab
#' \eqn{2-2^{1/\theta}} \cr \code{7} \tab \eqn{2^{-1/(\theta\delta)}} \tab
#' \eqn{2-2^{1/\delta}} \cr \code{8} \tab - \tab \eqn{2-2^{1/(\theta\delta)}}
#' \cr \code{9} \tab \eqn{2^{-1/\delta}} \tab \eqn{2-2^{1/\theta}} \cr
#' \code{10} \tab - \tab \eqn{2-2^{1/\theta}} if \eqn{\delta=1} otherwise 0 \cr
#' \code{13} \tab - \tab \eqn{2^{-1/\theta}} \cr \code{14} \tab
#' \eqn{2-2^{1/\theta}} \tab - \cr \code{16} \tab \eqn{2-2^{1/\theta}} \tab -
#' \cr \code{17} \tab \eqn{2-2^{1/\delta}} \tab \eqn{2^{-1/(\theta\delta)}} \cr
#' \code{18} \tab \eqn{2-2^{1/(\theta\delta)}} \tab - \cr \code{19} \tab
#' \eqn{2-2^{1/\theta}} \tab \eqn{2^{-1/\delta}} \cr \code{20} \tab
#' \eqn{2-2^{1/\theta}} if \eqn{\delta=1} otherwise 0 \tab - \cr \code{23, 33}
#' \tab - \tab - \cr \code{24, 34} \tab - \tab - \cr \code{26, 36} \tab - \tab
#' - \cr \code{27, 37} \tab - \tab - \cr \code{28, 38} \tab - \tab - \cr
#' \code{29, 39} \tab - \tab - \cr \code{30, 40} \tab - \tab - \cr
#' \code{104,204} \tab - \tab \eqn{\delta+1-(\delta^{\theta}+1)^{1/\theta}} \cr
#' \code{114, 214} \tab \eqn{1+\delta-(\delta^{\theta}+1)^{1/\theta}} \tab -
#' \cr \code{124, 224} \tab - \tab - \cr \code{134, 234} \tab - \tab - \cr }
#' @note The number \code{m} can be chosen arbitrarily.
#' @author Eike Brechmann
#' @seealso \code{\link{BiCopPar2Tau}}
#' @references Joe, H. (1997). Multivariate Models and Dependence Concepts.
#' Chapman and Hall, London.
#' @examples
#' 
#' ## Example 1: Gaussian copula
#' BiCopPar2TailDep(1, 0.7)
#' 
#' ## Example 2: t copula
#' BiCopPar2TailDep(2, 0.7, 4)
#' 
#' @export BiCopPar2TailDep
BiCopPar2TailDep <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
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
    n <- max(length(family), length(par), length(par2))
    if (length(family) == 1) 
        family <- rep(family, n)
    if (length(par) == 1) 
        par <- rep(par, n)
    if (length(par2) == 1)
        par2 <- rep(par2, n)
    if (!all(c(length(family), length(par), length(par2)) %in% c(1, n)))
        stop("Input lenghts don't match")
    
    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }
    
    ## calculate tail dependence coefficient
    if (length(par) == 1) {
        # call for single parameters
        out <- matrix(calcTD(family, par, par2), ncol = 2)
    } else {
        # vectorized call
        out <- t(vapply(1:length(par),
                        function(i) calcTD(family[i], par[i], par2[i]),
                        numeric(2)))
    }
    
    ## return result
    list(lower = out[, 1], upper = out[, 2])
}

calcTD <- function(family, par, par2) {
    if (family == 0 | family == 1 | family == 5 | family %in% c(23, 24, 26, 27, 28, 29,
                                                                30, 33, 34, 36, 37, 38, 39,
                                                                40, 124, 134, 224, 234)) {
        lower <- 0
        upper <- 0
    } else if (family == 2) {
        lower <- 2 * pt((-sqrt(par2 + 1) * sqrt((1 - par)/(1 + par))), df = par2 + 
                            1)
        upper <- lower
    } else if (family == 3) {
        lower <- 2^(-1/par)
        upper <- 0
    } else if (family == 4 | family == 6) {
        lower <- 0
        upper <- 2 - 2^(1/par)
    } else if (family == 7) {
        lower <- 2^(-1/(par * par2))
        upper <- 2 - 2^(1/par2)
    } else if (family == 8) {
        lower <- 0
        upper <- 2 - 2^(1/(par * par2))
    } else if (family == 9) {
        lower <- 2^(-1/par2)
        upper <- 2 - 2^(1/par)
    } else if (family == 10) {
        lower <- 0
        if (par2 == 1) 
            upper <- 2 - 2^(1/par) else upper <- 0
    } else if (family == 13) {
        lower <- 0
        upper <- 2^(-1/par)
    } else if (family == 14 | family == 16) {
        lower <- 2 - 2^(1/par)
        upper <- 0
    } else if (family == 17) {
        lower <- 2 - 2^(1/par2)
        upper <- 2^(-1/par * par2)
    } else if (family == 18) {
        lower <- 2 - 2^(1/(par * par2))
        upper <- 0
    } else if (family == 19) {
        lower <- 2 - 2^(1/par)
        upper <- 2^(-1/par2)
    } else if (family == 20) {
        if (par2 == 1) 
            lower <- 2 - 2^(1/par) else lower <- 0
            upper <- 0
    } else if (family == 104) {
        par3 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 114) {
        par3 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    } else if (family == 204) {
        par3 <- par2
        par2 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 214) {
        par3 <- par2
        par2 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    }
    
    ## return result
    c(upper, lower)
}
