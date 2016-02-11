#' Kendall's Tau Value of a Bivariate Copula
#'
#' This function computes the theoretical Kendall's tau value of a bivariate
#' copula for given parameter values.
#'
#' If the family and parameter specification is stored in a \code{\link{BiCop}}
#' object \code{obj}, the alternative version \cr
#' \preformatted{BiCopPar2Tau(obj)} can be used.
#'
#' @param family integer; single number or vector of size \code{m}; defines the
#' bivariate copula family: \cr
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
#' @param par numeric; single number or vector of size \code{m}; copula
#' parameter.
#' @param par2 numeric; single number or vector of size \code{m}; second
#' parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8,
#' Tawn type 1 and type 2; default: \code{par2 = 0}).  Note that the degrees of
#' freedom parameter of the t-copula does not need to be set, because the
#' theoretical Kendall's tau value of the t-copula is independent of this
#' choice.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#' @return Theoretical value of Kendall's tau (vector) corresponding to the
#' bivariate copula \code{family} and parameter vector \eqn{(\theta, \delta) =}
#' \code{(par, par2)}.
#' \tabular{ll}{ No. (\code{family}) \tab Kendall's tau (\code{tau}) \cr
#' \code{1, 2} \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
#' \code{3, 13} \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
#' \code{4, 14} \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
#' \code{5} \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta +
#' 4 D_1(\theta)/\theta} \cr
#' \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{x/\theta}{\exp(x)-1}dx}{D_1(\theta)=
#' \int_0^\theta (x/\theta)/(exp(x)-1)dx} (Debye function) \cr
#' \code{6, 16} \tab \eqn{1+\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx}{1+4/\theta^2\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx} \cr
#' \code{7, 17} \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
#' \code{8, 18} \tab \eqn{1+4\int_0^1 -\log(-(1-t)^\theta+1)
#' (1-t-(1-t)^{-\theta}+(1-t)^{-\theta}t)/(\delta\theta) dt} \cr
#' \code{9, 19} \tab \eqn{1+4\int_0^1 ( (1-(1-t)^{\theta})^{-\delta} - 1)
#' /( -\theta\delta(1-t)^{\theta-1}(1-(1-t)^{\theta})^{-\delta-1} ) dt} \cr
#' \code{10, 20} \tab \eqn{1+4\int_0^1
#' -\log \left(((1-t\delta)^\theta-1)/((1-\delta)^\theta-1) \right) } \cr
#' \tab \eqn{* (1-t\delta-(1-t\delta)^{-\theta}+(1-t\delta)^{-\theta}t\delta)/(\theta\delta) dt} \cr
#' \code{23, 33} \tab \eqn{\frac{\theta}{2-\theta}}{\theta/(2-\theta)} \cr
#' \code{24, 34} \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
#' \code{26, 36} \tab \eqn{-1-\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{-2(1+\theta)/\theta}dx}{-1-4/\theta^2
#' \int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx} \cr
#' \code{27, 37} \tab \eqn{-1-\frac{2}{\delta(2-\theta)}}{1-2/(\delta(\theta+2))} \cr
#' \code{28, 38} \tab \eqn{-1-4\int_0^1 -\log(-(1-t)^{-\theta}+1)
#' (1-t-(1-t)^{\theta}+(1-t)^{\theta}t)/(\delta\theta) dt} \cr
#' \code{29, 39} \tab \eqn{-1-4\int_0^1 ( (1-(1-t)^{-\theta})^{\delta} - 1)
#' /( -\theta\delta(1-t)^{-\theta-1}(1-(1-t)^{-\theta})^{\delta-1} ) dt} \cr
#' \code{30, 40} \tab \eqn{-1-4\int_0^1 -\log
#' \left( ((1+t\delta)^{-\theta}-1)/((1+\delta)^{-\theta}-1) \right)} \cr
#' \tab \eqn{* (1+t\delta-(1+t\delta)^{\theta}-(1+t\delta)^{\theta}t\delta)/(\theta\delta) dt} \cr
#' \code{104,114} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{\theta}+t^{\theta}]^{1/\theta}} \cr
#' \code{204,214} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt}  \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#' \code{124,134} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{-\theta}+t^{-\theta}]^{-1/\theta}} \cr
#' \code{224,234} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#'
#' }
#' @note The number \code{m} can be chosen arbitrarily.
#' @author Ulf Schepsmeier, Tobias Erhardt
#' @seealso \code{\link{BiCopTau2Par}}, \code{\link{BiCop}}
#' @references Joe, H. (1997). Multivariate Models and Dependence Concepts.
#' Chapman and Hall, London.
#'
#' Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation
#' of mixed C-vines with application to exchange rates. Statistical Modelling,
#' 12(3), 229-255.
#' @examples
#'
#' ## Example 1: Gaussian copula
#' tau0 <- 0.5
#' rho <- BiCopTau2Par(family = 1, tau = tau0)
#'
#' # transform back
#' tau <- BiCopPar2Tau(family = 1, par = rho)
#' tau - 2/pi*asin(rho)
#'
#'
#' ## Example 2: Student-t copula
#' rho <- BiCopTau2Par(family = 2, tau = c(0.4, 0.5, 0.6))
#' obj <- BiCop(family = 2, par = rho, par2 = 6:8)
#' BiCopPar2Tau(obj)
#'
#'
#' ## Example 3:
#' vpar <- seq(from = 1.1, to = 10, length.out = 100)
#' tauC <- BiCopPar2Tau(family = 3, par = vpar)
#' tauG <- BiCopPar2Tau(family = 4, par = vpar)
#' tauF <- BiCopPar2Tau(family = 5, par = vpar)
#' tauJ <- BiCopPar2Tau(family = 6, par = vpar)
#' plot(tauC ~ vpar, type = "l", ylim = c(0,1))
#' lines(tauG ~ vpar, col = 2)
#' lines(tauF ~ vpar, col = 3)
#' lines(tauJ ~ vpar, col = 4)
#'
#'
#' ## Example 4: different copula families
#' theta <- BiCopTau2Par(family = c(3,4,6), tau = c(0.4, 0.5, 0.6))
#' BiCopPar2Tau(family = c(3,4,6), par = theta)
#'
#' \dontshow{
#' # Test BiCopPar2Tau (one parametric families)
#' theta <- BiCopTau2Par(family = 0, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 0, par = theta)
#' theta <- BiCopTau2Par(family = 1, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 1, par = theta)
#' theta <- BiCopTau2Par(family = 3, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 3, par = theta)
#' theta <- BiCopTau2Par(family = 4, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 4, par = theta)
#' theta <- BiCopTau2Par(family = 5, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 5, par = theta)
#' theta <- BiCopTau2Par(family = 6, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 6, par = theta)
#' theta <- BiCopTau2Par(family = 13, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 13, par = theta)
#' theta <- BiCopTau2Par(family = 14, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 14, par = theta)
#' theta <- BiCopTau2Par(family = 16, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 16, par = theta)
#' theta <- BiCopTau2Par(family = 23, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 23, par = theta)
#' theta <- BiCopTau2Par(family = 24, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 24, par = theta)
#' theta <- BiCopTau2Par(family = 26, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 26, par = theta)
#' theta <- BiCopTau2Par(family = 33, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 33, par = theta)
#' theta <- BiCopTau2Par(family = 34, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 34, par = theta)
#' theta <- BiCopTau2Par(family = 36, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 36, par = theta)
#' theta <- BiCopTau2Par(family = 41, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 41, par = theta)
#' theta <- BiCopTau2Par(family = 51, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 51, par = theta)
#' theta <- BiCopTau2Par(family = 61, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 61, par = theta)
#' theta <- BiCopTau2Par(family = 71, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 71, par = theta)
#' theta <- BiCopTau2Par(family = 41, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 41, par = theta)
#' theta <- BiCopTau2Par(family = 51, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 51, par = theta)
#' theta <- BiCopTau2Par(family = 61, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 61, par = theta)
#' theta <- BiCopTau2Par(family = 71, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 71, par = theta)
#'
#' # Test BiCopPar2Tau (two parametric families)
#' theta <- BiCopTau2Par(family = 2, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 2, par = theta)
#' theta <- 1:3
#' delta <- 1:3
#' BiCopPar2Tau(family = 7, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 17, par = theta, par2 = delta)
#' theta <- -(1:3)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 27, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 37, par = theta, par2 = delta)
#' theta <- 2:4
#' delta <- 1:3
#' BiCopPar2Tau(family = 8, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 18, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 28, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 38, par = theta, par2 = delta)
#' theta <- 1:3
#' delta <- 1:3
#' BiCopPar2Tau(family = 9, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 19, par = theta, par2 = delta)
#' theta <- -(1:3)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 29, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 39, par = theta, par2 = delta)
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 10, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 20, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- -c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 30, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 40, par = theta, par2 = delta)
#'
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 104, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 114, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 124, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 134, par = theta, par2 = delta)
#'
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 204, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 214, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 224, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 234, par = theta, par2 = delta)
#' }
#'
#' @export BiCopPar2Tau
BiCopPar2Tau <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
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
    n <- max(length(family), length(par), length(par2))
    if (length(family) == 1)
        family <- rep(family, n)
    if (length(par) == 1)
        par <- rep(par, n)
    if (length(par2) == 1)
        par2 <- rep(par2, n)
    if (!all(c(length(family), length(par), length(par2)) %in% c(1, n)))
        stop("Input lenghts don't match")

    # set arbitrary par2 for t-copula
    par2[family == 2] <- par2[family == 2] + 4

    ## check for family/parameter consistency
    if (check.pars)
        BiCopCheck(family, par, par2)

    ## calculate Kendall's tau
    if (length(par) == 1) {
        # call for single parameters
        out <- calcTau(family, par, par2)
    } else {
        # vectorized call
        out <- vapply(1:length(par),
                      function(i) calcTau(family[i], par[i], par2[i]),
                      numeric(1))
    }

    ## return result
    out
}

calcTau <- function(family, par, par2) {
    ## calculation of tau(s) depending on pair-copula family
    if (family == 0) {
        tau <- rep(0, times = length(par))
    } else if (family == 1 | family == 2) {
        tau <- 2/pi * asin(par)
    } else if (family == 3 || family == 13) {
        tau <- par/(par + 2)
    } else if (family == 4 || family == 14) {
        tau <- 1 - 1/par
    } else if (family == 5) {
        tau <- if (par == 0) 0 else 1 - 4/par + 4/par * debye1(par)
    } else if (family == 6 || family == 16) {
        # tau = 1 + 4/par^2 * integrate(function(x) log(x)*x*(1-x)^(2*(1-par)/par), 0,
        # 1)$value
        param1 <- 2/par + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - par)
        tau[par == 2] <- 1 - trigamma(2)
    } else if (family == 7 || family == 17) {
        theta <- par
        delta <- par2
        tau <- 1 - 2/(delta * (theta + 2))
    } else if (family == 8 || family == 18) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 9 || family == 19) {
        theta <- par
        delta <- par2

        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^-de - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t)  kt(t, th = theta, de = delta), 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 10 || family == 20) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 23 || family == 33) {
        tau <- par/(-par + 2)
    } else if (family == 24 || family == 34) {
        tau <- -1 - 1/par
    } else if (family == 26 || family == 36) {
        theta <- -par
        param1 <- 2/theta + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - theta)
        tau[theta == 2] <- 1 - trigamma(2)
        tau <- -tau
    } else if (family == 27 || family == 37) {
        theta <- -par
        delta <- -par2
        tau <- 1 - 2/(delta * (theta + 2))
        tau <- -tau
    } else if (family == 28 || family == 38) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 29 || family == 39) {
        theta <- -par
        delta <- -par2

        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^(-de) - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 30 || family == 40) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 41 || family == 51) {
        de <- par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
    } else if (family == 61 || family == 71) {
        de <- -par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
        tau <- -tau
    } else if (family == 42) {
        tau <- (75 * par2 - par2^2 + par * (25 - par2))/450
    } else if (family == 104 || family == 114 || family == 204 || family == 214) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2",
                        as.double(t),
                        as.integer(length(t)),
                        as.double(th),
                        as.double(de),
                        as.double(1),
                        as.double(rep(0, length(t))),
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn",
                              as.double(t),
                              as.integer(length(t)),
                              as.double(th),
                              as.double(de),
                              as.double(1),
                              as.double(rep(0, length(t))),
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- try(mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 124 || family == 134 || family == 224 || family == 234) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2",
                        as.double(t),
                        as.integer(length(t)),
                        as.double(-th),
                        as.double(de),
                        as.double(1),
                        as.double(rep(0, length(t))),
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn",
                              as.double(t),
                              as.integer(length(t)),
                              as.double(-th),
                              as.double(de),
                              as.double(1),
                              as.double(rep(0, length(t))),
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- try(mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    }

    ## return result
    tau
}
