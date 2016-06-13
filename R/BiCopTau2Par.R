#' Parameter of a Bivariate Copula for a given Kendall's Tau Value
#'
#' This function computes the parameter of a (one parameter) bivariate copula
#' for a given value of Kendall's tau.
#'
#'
#' @param family integer; single number or vector of size \code{n}; defines the
#' bivariate copula family: \cr \code{0} = independence copula \cr \code{1} =
#' Gaussian copula \cr \code{2} = Student t copula (Here only the first
#' parameter can be computed) \cr \code{3} = Clayton copula \cr \code{4} =
#' Gumbel copula \cr \code{5} = Frank copula \cr \code{6} = Joe copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr \code{23}
#' = rotated Clayton copula (90 degrees) \cr \code{24} = rotated Gumbel copula
#' (90 degrees) \cr \code{26} = rotated Joe copula (90 degrees) \cr \code{33} =
#' rotated Clayton copula (270 degrees) \cr \code{34} = rotated Gumbel copula
#' (270 degrees) \cr \code{36} = rotated Joe copula (270 degrees)\cr Note that
#' (with exception of the t-copula) two parameter bivariate copula families
#' cannot be used.
#' @param tau numeric; single number or vector of size \code{n}; Kendall's tau
#' value (vector with elements in [-1,1]).
#' @param check.taus logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/tau-consistency are ommited (should only be used with care).
#'
#' @return Parameter (vector) corresponding to the bivariate copula family and
#' the value(s) of Kendall's tau (\eqn{\tau}). \tabular{ll}{ No.
#' (\code{family}) \tab Parameter (\code{par}) \cr \code{1, 2} \tab
#' \eqn{\sin(\tau \frac{\pi}{2})}{sin(\tau \pi/2)} \cr \code{3, 13} \tab
#' \eqn{2\frac{\tau}{1-\tau}}{2\tau/(1-\tau)} \cr \code{4, 14} \tab
#' \eqn{\frac{1}{1-\tau}}{1/(1-\tau)} \cr \code{5} \tab no closed form
#' expression (numerical inversion) \cr \code{6, 16} \tab no closed form
#' expression (numerical inversion) \cr \code{23, 33} \tab
#' \eqn{2\frac{\tau}{1+\tau}}{2\tau/(1+\tau)} \cr \code{24, 34} \tab
#' \eqn{-\frac{1}{1+\tau}}{-1/(1+\tau)} \cr \code{26, 36} \tab no closed form
#' expression (numerical inversion) }
#'
#' @note The number \code{n} can be chosen arbitrarily, but must agree across
#' arguments.
#'
#' @author Jakob Stoeber, Eike Brechmann, Tobias Erhardt
#'
#' @seealso \code{\link{BiCopPar2Tau}}
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
#' BiCop(1, tau = tau0)$par  # alternative
#'
#' ## Example 2:
#' vtau <- seq(from = 0.1, to = 0.8, length.out = 100)
#' thetaC <- BiCopTau2Par(family = 3, tau = vtau)
#' thetaG <- BiCopTau2Par(family = 4, tau = vtau)
#' thetaF <- BiCopTau2Par(family = 5, tau = vtau)
#' thetaJ <- BiCopTau2Par(family = 6, tau = vtau)
#' plot(thetaC ~ vtau, type = "l", ylim = range(thetaF))
#' lines(thetaG ~ vtau, col = 2)
#' lines(thetaF ~ vtau, col = 3)
#' lines(thetaJ ~ vtau, col = 4)
#'
#' ## Example 3: different copula families
#' theta <- BiCopTau2Par(family = c(3,4,6), tau = c(0.4, 0.5, 0.6))
#' BiCopPar2Tau(family = c(3,4,6), par = theta)
#'
#' \dontshow{
#' # Test BiCopTau2Par
#' BiCopTau2Par(family = 0, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 1, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 2, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 3, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 4, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 5, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 6, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 13, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 14, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 16, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 23, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 24, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 26, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 33, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 34, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 36, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 41, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 51, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 61, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 71, tau = c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 41, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 51, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 61, tau = -c(0.4,0.5,0.6))
#' BiCopTau2Par(family = 71, tau = -c(0.4,0.5,0.6))
#' }
#'
#' @export BiCopTau2Par
BiCopTau2Par <- function(family, tau, check.taus = TRUE) {
    ## sanity check
    if (any(family %in% setdiff(allfams[twopar], 2)))
        stop("For two parameter copulas (except t) Kendall's tau cannot be inverted.")
    if (any(abs(tau) > 0.99999))
        stop("some tau is too close to -1 or 1")

    # fix for SemiParBIVProbit package
    dims <- set_dims(family, tau = tau)

    ## adjust length for input vectors; stop if not matching
    family <- c(family)
    tau <- c(tau)
    n <- max(length(family), length(tau))
    if (length(family) == 1)
        family <- rep(family, n)
    if (length(tau) == 1)
        par <- rep(tau, n)
    if (!all(c(length(family), length(tau)) %in% c(1, n)))
        stop("Input lenghts don't match")

    ## check for family/tau consistency
    if (check.taus)
        BiCopCheckTaus(family, tau)

    ## calculate the parameter
    out <- vapply(1:length(tau),
                  function(i) calcPar(family[i], tau[i]),
                  numeric(1))

    ## return result
    if (length(dims) > 1)
        out <- array(out, dim = dims)
    out
}

calcPar <- function(family, tau) {
    ## calculation of parameter(s) depending on pair-copula family
    if (family == 0) {
        par <- rep(0, times = length(tau))
    } else if (family %in% 1:2) {
        par <- sin(pi * tau/2)
    } else if (family %in% c(3, 13)) {
        par <- 2 * tau/(1 - tau)
    } else if (family %in% c(4, 14)) {
        par <- 1/(1 - tau)
    } else if (family == 5) {
        par <- if (tau == 0) 0 else Frank.itau.JJ(tau)
    } else if (family %in% c(6, 16)) {
        par <- Joe.itau.JJ(tau)
    } else if (family %in% c(23, 33)) {
        par <- 2 * tau/(1 + tau)
    } else if (family %in% c(24, 34)) {
        par <- -(1/(1 + tau))
    } else if (family %in% c(26, 36)) {
        par <- -Joe.itau.JJ(-tau)
    } else if (family %in% c(41, 51)) {
        par <- ipsA.tau2cpar(tau)
    } else if (family %in% c(61, 71)) {
        par <- -ipsA.tau2cpar(-tau)
    }

    ## return result
    par
}

Frank.itau.JJ <- function(tau) {
    a <- 1
    if (tau < 0) {
        a <- -1
        tau <- -tau
    }
    v <- uniroot(function(x) tau - (1 - 4/x + 4/x * debye1(x)),
                 lower = 0 + .Machine$double.eps^0.5, upper = 5e5,
                 tol = .Machine$double.eps^0.5)$root
    return(a*v)
}


Joe.itau.JJ <- function(tau) {
    if (tau < 0) {
        return(1.000001)
    } else {
        tauF <- function(par) {
            param1 <- 2/par + 1
            tem <- digamma(2) - digamma(param1)
            tau <- 1 + tem * 2/(2 - par)
            tau[par == 2] <- 1 - trigamma(2)
            tau
        }

        v <- uniroot(function(x) tau - tauF(x),
                     lower = 1,
                     upper = 5e5,
                     tol = .Machine$double.eps^0.5)$root
        return(v)
    }
}

ipsA.tau2cpar <- function(tau, mxiter = 20, eps = 1e-06, dstart = 0, iprint = FALSE) {
    con <- log((1 - tau) * sqrt(pi)/2)
    de <- dstart
    if (dstart <= 0)
        de <- tau + 1
    iter <- 0
    diff <- 1
    while (iter < mxiter & max(abs(diff)) > eps) {
        g <- con + lgamma(1 + de) - lgamma(de + 0.5)
        gp <- digamma(1 + de) - digamma(de + 0.5)
        iter <- iter + 1
        diff <- g/gp
        de <- de - diff
        while (min(de) <= 0) {
            diff <- diff/2
            de <- de + diff
        }
        if (iprint)
            cat(iter, " ", de, " ", diff, "\n")
    }
    if (iter >= mxiter)
        cat("did not converge\n")
    de
}

