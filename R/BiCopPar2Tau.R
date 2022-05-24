#' Kendall's Tau Value of a Bivariate Copula
#'
#' This function computes the theoretical Kendall's tau value of a bivariate
#' copula for given parameter values.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version \cr
#' \preformatted{BiCopPar2Tau(obj)} can be used.
#'
#' @param family integer; single number or vector of size `m`; defines the
#' bivariate copula family: \cr
#' `0` = independence copula \cr
#' `1` = Gaussian copula \cr
#' `2` = Student t copula (t-copula) \cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `7` = BB1 copula \cr
#' `8` = BB6 copula \cr
#' `9` = BB7 copula \cr
#' `10` = BB8 copula \cr
#' `13` = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' `14` = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' `16` = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' `17` = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' `18` = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' `19` = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' `20` = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' `23` = rotated Clayton copula (90 degrees) \cr
#' `24` = rotated Gumbel copula (90 degrees) \cr
#' `26` = rotated Joe copula (90 degrees) \cr
#' `27` = rotated BB1 copula (90 degrees) \cr
#' `28` = rotated BB6 copula (90 degrees) \cr
#' `29` = rotated BB7 copula (90 degrees) \cr
#' `30` = rotated BB8 copula (90 degrees) \cr
#' `33` = rotated Clayton copula (270 degrees) \cr
#' `34` = rotated Gumbel copula (270 degrees) \cr
#' `36` = rotated Joe copula (270 degrees) \cr
#' `37` = rotated BB1 copula (270 degrees) \cr
#' `38` = rotated BB6 copula (270 degrees) \cr
#' `39` = rotated BB7 copula (270 degrees) \cr
#' `40` = rotated BB8 copula (270 degrees) \cr
#' `104` = Tawn type 1 copula \cr
#' `114` = rotated Tawn type 1 copula (180 degrees) \cr
#' `124` = rotated Tawn type 1 copula (90 degrees) \cr
#' `134` = rotated Tawn type 1 copula (270 degrees) \cr
#' `204` = Tawn type 2 copula \cr
#' `214` = rotated Tawn type 2 copula (180 degrees) \cr
#' `224` = rotated Tawn type 2 copula (90 degrees) \cr
#' `234` = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par numeric; single number or vector of size `n`; copula
#' parameter.
#' @param par2 numeric; single number or vector of size `n`; second
#' parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8,
#' Tawn type 1 and type 2; default: `par2 = 0`).  Note that the degrees of
#' freedom parameter of the t-copula does not need to be set, because the
#' theoretical Kendall's tau value of the t-copula is independent of this
#' choice.
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#' @return Theoretical value of Kendall's tau (vector) corresponding to the
#' bivariate copula `family` and parameter vector \eqn{(\theta, \delta) =}
#' `(par, par2)`.
#' \tabular{ll}{ No. (`family`) \tab Kendall's tau (`tau`) \cr
#' `1, 2` \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
#' `3, 13` \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
#' `4, 14` \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
#' `5` \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta +
#' 4 D_1(\theta)/\theta} \cr
#' \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{x/\theta}{\exp(x)-1}dx}{D_1(\theta)=
#' \int_0^\theta (x/\theta)/(exp(x)-1)dx} (Debye function) \cr
#' `6, 16` \tab \eqn{1+\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx}{1+4/\theta^2\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx} \cr
#' `7, 17` \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
#' `8, 18` \tab \eqn{1+4\int_0^1 -\log(-(1-t)^\theta+1)
#' (1-t-(1-t)^{-\theta}+(1-t)^{-\theta}t)/(\delta\theta) dt} \cr
#' `9, 19` \tab \eqn{1+4\int_0^1 ( (1-(1-t)^{\theta})^{-\delta} - 1)
#' /( -\theta\delta(1-t)^{\theta-1}(1-(1-t)^{\theta})^{-\delta-1} ) dt} \cr
#' `10, 20` \tab \eqn{1+4\int_0^1
#' -\log \left(((1-t\delta)^\theta-1)/((1-\delta)^\theta-1) \right) } \cr
#' \tab \eqn{* (1-t\delta-(1-t\delta)^{-\theta}+(1-t\delta)^{-\theta}t\delta)/(\theta\delta) dt} \cr
#' `23, 33` \tab \eqn{\frac{\theta}{2-\theta}}{\theta/(2-\theta)} \cr
#' `24, 34` \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
#' `26, 36` \tab \eqn{-1-\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{-2(1+\theta)/\theta}dx}{-1-4/\theta^2
#' \int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx} \cr
#' `27, 37` \tab \eqn{-1-\frac{2}{\delta(2-\theta)}}{1-2/(\delta(\theta+2))} \cr
#' `28, 38` \tab \eqn{-1-4\int_0^1 -\log(-(1-t)^{-\theta}+1)
#' (1-t-(1-t)^{\theta}+(1-t)^{\theta}t)/(\delta\theta) dt} \cr
#' `29, 39` \tab \eqn{-1-4\int_0^1 ( (1-(1-t)^{-\theta})^{\delta} - 1)
#' /( -\theta\delta(1-t)^{-\theta-1}(1-(1-t)^{-\theta})^{\delta-1} ) dt} \cr
#' `30, 40` \tab \eqn{-1-4\int_0^1 -\log
#' \left( ((1+t\delta)^{-\theta}-1)/((1+\delta)^{-\theta}-1) \right)} \cr
#' \tab \eqn{* (1+t\delta-(1+t\delta)^{\theta}-(1+t\delta)^{\theta}t\delta)/(\theta\delta) dt} \cr
#' `104,114` \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{\theta}+t^{\theta}]^{1/\theta}} \cr
#' `204,214` \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt}  \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#' `124,134` \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{-\theta}+t^{-\theta}]^{-1/\theta}} \cr
#' `224,234` \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#'
#' }
#'
#' @note The number `n` can be chosen arbitrarily, but must agree across
#' arguments.
#'
#' @author Ulf Schepsmeier, Tobias Erhardt
#'
#' @seealso [BiCopTau2Par()], [BiCop()]
#'
#' @references Joe, H. (1997). Multivariate Models and Dependence Concepts.
#' Chapman and Hall, London.
#'
#' Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation
#' of mixed C-vines with application to exchange rates. Statistical Modelling,
#' 12(3), 229-255.
#'
#' @examples
#' ## Example 1: Gaussian copula
#' tau0 <- 0.5
#' rho <- BiCopTau2Par(family = 1, tau = tau0)
#' # transform back
#' tau <- BiCopPar2Tau(family = 1, par = rho)
#' tau - 2/pi*asin(rho)
#'
#' ## Example 2:
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
#' ## Example 3: different copula families
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
BiCopPar2Tau <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    # fix for SemiParBIVProbit package
    dims <- set_dims(family, par, par2)
    # set arbitrary par2 for t-copula
    if (!inherits(family, "BiCop"))
        par2[family == 2] <- par2[family == 2] + 4
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate Kendall's tau
    out <- vapply(1:length(par),
                  function(i) calcTau(family[i], par[i], par2[i]),
                  numeric(1))

    ## return result
    if (length(dims) > 1)
        out <- array(out, dim = dims)
    out
}


frankParGrid <- c(-10^(5:2), seq(-36, 36, l = 100), 10^(2:5))
# frankTauVals <- 1 - 4/frankParGrid + 4/frankParGrid * copula::debye1(frankParGrid)
frankTauVals <- c(-0.99996000, -0.99960007, -0.99600658, -0.96065797, -0.89396585, -0.89188641, -0.88972402, -0.88747361,
                  -0.88512973, -0.88268644, -0.88013730, -0.87747533, -0.87469288, -0.87178163, -0.86873246, -0.86553538,
                  -0.86217942, -0.85865252, -0.85494133, -0.85103114, -0.84690560, -0.84254657, -0.83793380, -0.83304468,
                  -0.82785385, -0.82233280, -0.81644934, -0.81016705, -0.80344454, -0.79623459, -0.78848313, -0.78012801,
                  -0.77109744, -0.76130821, -0.75066334, -0.73904940, -0.72633308, -0.71235706, -0.69693500, -0.67984550,
                  -0.66082501, -0.63955977, -0.61567712, -0.58873731, -0.55822774, -0.52356346, -0.48410024, -0.43916961,
                  -0.38814802, -0.33057147, -0.26629772, -0.19569472, -0.11979813, -0.04035073,  0.04035073,  0.11979813,
                  0.19569472,  0.26629772,  0.33057147,  0.38814802,  0.43916961,  0.48410024, 0.52356346,  0.55822774,
                  0.58873731,  0.61567712,  0.63955977,  0.66082501,  0.67984550,  0.69693500,  0.71235706,  0.72633308,
                  0.73904940,  0.75066334,  0.76130821,  0.77109744,  0.78012801,  0.78848313,  0.79623459,  0.80344454,
                  0.81016705,  0.81644934,  0.82233280,  0.82785385,  0.83304468,  0.83793380,  0.84254657,  0.84690560,
                  0.85103114,  0.85494133,  0.85865252,  0.86217942,  0.86553538,  0.86873246,  0.87178163,  0.87469288,
                  0.87747533,  0.88013730,  0.88268644,  0.88512973,  0.88747361,  0.88972402,  0.89188641,  0.89396585,
                  0.96065797,  0.99600658,  0.99960007,  0.99996000)
frankTau <- function(par) {
    approx(x = frankParGrid, y = frankTauVals, xout = par)$y
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
        tau <- frankTau(par)
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
