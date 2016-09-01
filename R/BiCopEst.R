#' Parameter Estimation for Bivariate Copula Data
#'
#' This function estimates the parameter(s) of a bivariate copula using either
#' inversion of empirical Kendall's tau (for one parameter copula families only) or
#' maximum likelihood estimation for implemented copula families.
#'
#' If \code{method = "itau"}, the function computes the empirical Kendall's tau
#' of the given copula data and exploits the one-to-one relationship of copula
#' parameter and Kendall's tau which is available for many one parameter
#' bivariate copula families (see \code{\link{BiCopPar2Tau}} and
#' \code{\link{BiCopTau2Par}}). The inversion of Kendall's tau is however not
#' available for all bivariate copula families (see above). If a two parameter
#' copula family is chosen and \code{method = "itau"}, a warning message is
#' returned and the MLE is calculated.
#'
#' For \code{method = "mle"} copula parameters are estimated by maximum
#' likelihood using starting values obtained by \code{method = "itau"}.  If no
#' starting values are available by inversion of Kendall's tau, starting values
#' have to be provided given expert knowledge and the boundaries \code{max.df}
#' and \code{max.BB} respectively. Note: The MLE is performed via numerical
#' maximazation using the L_BFGS-B method. For the Gaussian, the t- and the
#' one-parametric Archimedean copulas we can use the gradients, but for the BB
#' copulas we have to use finite differences for the L_BFGS-B method.
#'
#' A warning message is returned if the estimate of the degrees of freedom
#' parameter of the t-copula is larger than \code{max.df}. For high degrees of
#' freedom the t-copula is almost indistinguishable from the Gaussian and it is
#' advised to use the Gaussian copula in this case. As a rule of thumb
#' \code{max.df = 30} typically is a good choice. Moreover, standard errors of
#' the degrees of freedom parameter estimate cannot be estimated in this case.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param family An integer defining the bivariate copula family: \cr
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
#' @param method Character indicating the estimation method: either maximum
#' likelihood estimation (\code{method = "mle"}; default) or inversion of
#' Kendall's tau (\code{method = "itau"}).\cr For \code{method = "itau"} only
#' one parameter bivariate copula families can be used (\code{family =
#' 1,3,4,5,6,13,14,16,23,24,26,33,34} or \code{36}).
#' @param se Logical; whether standard error(s) of parameter estimates is/are
#' estimated (default: \code{se = FALSE}).
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: \code{max.df = 30}).
#' @param max.BB List; upper bounds for the estimation of the two parameters
#' (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr (default:
#' \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).
#' @param weights Numerical; weights for each observation (opitional).
#'
#' @return An object of class \code{\link{BiCop}}, augmented with the following
#' entries:
#' \item{se, se2}{standard errors for the parameter estimates (if
#' \code{se = TRUE},}
#' \item{nobs}{number of observations,}
#' \item{logLik}{log likelihood}
#' \item{AIC}{Aikaike's Informaton Criterion,}
#' \item{BIC}{Bayesian's Informaton Criterion,}
#' \item{emptau}{empirical value of Kendall's tau,}
#' \item{p.value.indeptest}{p-value of the independence test.}
#'
#' @note For a comprehensive summary of the fitted model, use \code{summary(object)};
#' to see all its contents, use \code{str(object)}.
#'
#' @author Ulf Schepsmeier, Eike Brechmann, Jakob Stoeber, Carlos Almeida
#'
#' @seealso
#' \code{\link{BiCop}},
#' \code{\link{BiCopPar2Tau}},
#' \code{\link{BiCopTau2Par}},
#' \code{\link{RVineSeqEst}},
#' \code{\link{BiCopSelect}},
#'
#' @references Joe, H. (1997). Multivariate Models and Dependence Concepts.
#' Chapman and Hall, London.
#'
#' @examples
#'
#' ## Example 1: bivariate Gaussian copula
#' dat <- BiCopSim(500, 1, 0.7)
#' u1 <- dat[, 1]
#' v1 <- dat[, 2]
#'
#' # estimate parameters of Gaussian copula by inversion of Kendall's tau
#' est1.tau <- BiCopEst(u1, v1, family = 1, method = "itau")
#' est1.tau  # short overview
#' summary(est1.tau)  # comprehensive overview
#' str(est1.tau)  # see all contents of the object
#'
#' # check if parameter actually coincides with inversion of Kendall's tau
#' tau1 <- cor(u1, v1, method = "kendall")
#' all.equal(BiCopTau2Par(1, tau1), est1.tau$par)
#'
#' # maximum likelihood estimate for comparison
#' est1.mle <- BiCopEst(u1, v1, family = 1, method = "mle")
#' summary(est1.mle)
#'
#'
#' ## Example 2: bivariate Clayton and survival Gumbel copulas
#' # simulate from a Clayton copula
#' dat <- BiCopSim(500, 3, 2.5)
#' u2 <- dat[, 1]
#' v2 <- dat[, 2]
#'
#' # empirical Kendall's tau
#' tau2 <- cor(u2, v2, method = "kendall")
#'
#' # inversion of empirical Kendall's tau for the Clayton copula
#' BiCopTau2Par(3, tau2)
#' BiCopEst(u2, v2, family = 3, method = "itau")
#'
#' # inversion of empirical Kendall's tau for the survival Gumbel copula
#' BiCopTau2Par(14, tau2)
#' BiCopEst(u2, v2, family = 14, method = "itau")
#'
#' # maximum likelihood estimates for comparison
#' BiCopEst(u2, v2, family = 3, method = "mle")
#' BiCopEst(u2, v2, family = 14, method = "mle")
#'
#'
BiCopEst <- function(u1, u2, family, method = "mle", se = FALSE, max.df = 30,
                     max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)),
                     weights = NA) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    check_est_pars,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    ## calculate empirical Kendall's tau and invert for initial estimate
    tau <- fasttau(u1, u2, weights)
    if (family %in% c(0, 2, allfams[onepar]))
        theta <- BiCopTau2Par(family, tau)

    ## inversion of kendall's tau -----------------------------
    if (method == "itau") {

        ## standard errors for method itau
        se1 <- 0
        if (se == TRUE) {
            p <- 2
            n <- length(u1)
            ec <- numeric(n)
            u <- cbind(u1, u2)
            v <- matrix(0, n, p * (p - 1)/2)

            if (family == 1)
                tauder <- function(x) {
                    2/(pi * sqrt(1 - x^2))
                } else if (family %in% c(3, 13, 23, 33)) {
                    tauder <- function(x) 2 * (2 + x)^(-2)
                } else if (family %in% c(4, 14, 24, 34)) {
                    tauder <- function(x) x^(-2)
                } else if (family == 5) {
                    f <- function(x) x/(exp(x) - 1)
                    tauder <- function(x) {
                        lwr <- 0 + .Machine$double.eps^0.5
                        intgrl <- integrate(f,
                                            lower = lwr,
                                            upper = x)$value
                        4/x^2 - 8/x^3 * intgrl + 4/(x * (exp(x) - 1))
                    }
                } else if (family %in% c(6, 16, 26, 36)) {
                    tauder <- function(x) {
                        euler <- 0.577215664901533
                        -((-2 + 2 * euler + 2 * log(2) + digamma(1/x) +
                               digamma(1/2 * (2 + x)/x) + x)/(-2 + x)^2) +
                            ((-trigamma(1/x)/x^2 + trigamma(1/2 * (2 + x)/x) *
                                  (1/(2 + x) - (2 + x)/(2 * x^2)) + 1)/(-2 + x))
                    }
                } else if (family %in% c(41, 51, 61, 71)) {
                    tauder <- function(x) {
                        2 * sqrt(pi) * gamma(0.5 + x) *
                            (digamma(1 + x) - digamma(0.5 + x))/gamma(1 + x)
                    }
                }

            l <- 1
            for (j in 1:(p - 1)) {
                for (i in (j + 1):p) {
                    for (k in 1:n)
                        ec[k] <- sum(u[, i] <= u[k, i] & u[, j] <= u[k, j])/n
                    v[, l] <- 2 * ec - u[, i] - u[, j]
                    l <- l + 1
                }
            }

            if (family == 0) {
                D <- 0
            } else if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 41, 51)) {
                D <- 1/tauder(theta)
            } else if (family %in% c(23, 33, 24, 34, 26, 36, 61, 71)) {
                D <- 1/tauder(-theta)
            }


            se1 <- as.numeric(sqrt(16/n * var(v %*% D)))
        }  # end if (se == TRUE)
    }  # end if (method == "itau")

    ## MLE ------------------------------------------
    if (method == "mle") {
        ## set starting parameters for maximum likelihood estimation
        theta1 <- 0
        delta <- 0

        if (!(family %in% c(2, 6, 7, 8, 9, 10,
                            17, 18, 19, 20,
                            27, 28, 29, 30,
                            37, 38, 39, 40,
                            104, 114, 124, 134,
                            204, 214, 224, 234))) {
            theta1 <- theta
        }
        if (family == 2) {
            ## t
            theta1 <- sin(tau * pi/2)
            delta <- 8
        } else if (family == 7 || family == 17) {
            ## BB1
            if (tau < 0) {
                print("The BB1 or survival BB1 copula cannot be used for negatively dependent data.")
                delta <- 1.001
                theta1 <- 0.001
            } else {
                delta <- min(1.5, max((max.BB$BB1[2] + 1.001)/2, 1.001))
                theta1 <- min(0.5, max((max.BB$BB1[1] + 0.001)/2, 0.001))
            }
        } else if (family == 27 || family == 37) {
            ## BB1
            if (tau > 0) {
                print("The rotated BB1 copulas cannot be used for positively dependent data.")
                delta <- -1.001
                theta1 <- -0.001
            } else {
                delta <- max(-1.5, -max((max.BB$BB1[2] + 1.001)/2, 1.001))
                theta1 <- max(-0.5, -max((max.BB$BB1[1] + 0.001)/2, 0.001))
            }
        } else if (family == 8 || family == 18) {
            ## BB6
            if (tau < 0) {
                print("The BB6 or survival BB6 copula cannot be used for negatively dependent data.")
                delta <- 1.001
                theta1 <- 1.001
            } else {
                delta <- min(1.5, max((max.BB$BB6[2] + 1.001)/2, 1.001))
                theta1 <- min(1.5, max((max.BB$BB6[1] + 1.001)/2, 1.001))
            }
        } else if (family == 28 || family == 38) {
            ## BB6
            if (tau > 0) {
                print("The rotated BB6 copulas cannot be used for positively dependent data.")
                delta <- -1.001
                theta1 <- -1.001
            } else {
                delta <- max(-1.5, -max((max.BB$BB6[2] + 1.001)/2, 1.001))
                theta1 <- max(-1.5, -max((max.BB$BB6[1] + 1.001)/2, 1.001))
            }
        } else if (family == 9 || family == 19) {
            ## BB7
            if (tau < 0) {
                print("The BB7 or survival BB7 copula cannot be used for negatively dependent data.")
                delta <- 0.001
                theta <- 1.001
            } else {
                delta <- min(0.5, max((max.BB$BB7[2] + 0.001)/2, 0.001))
                theta1 <- min(1.5, max((max.BB$BB7[1] + 1.001)/2, 1.001))
            }
        } else if (family == 29 || family == 39) {
            ## BB7
            if (tau > 0) {
                print("The rotated BB7 copulas cannot be used for positively dependent data.")
                delta <- -0.001
                theta1 <- -1.001
            } else {
                delta <- max(-0.5, -max((max.BB$BB7[2] + 0.001)/2, 0.001))
                theta1 <- max(-1.5, -max((max.BB$BB7[1] + 1.001)/2, 1.001))
            }
        } else if (family == 10 || family == 20) {
            ## BB8
            if (tau < 0) {
                print("The BB8 or survival BB8 copula cannot be used for negatively dependent data.")
                delta <- 0.001
                theta <- 1.001
            } else {
                delta <- min(0.5, max((max.BB$BB8[2] + 0.001)/2, 0.001))
                theta1 <- min(1.5, max((max.BB$BB8[1] + 1.001)/2, 1.001))
            }
        } else if (family == 30 || family == 40) {
            ## BB8
            if (tau > 0) {
                print("The rotated BB8 copulas cannot be used for positively dependent data.")
                delta <- -0.001
                theta1 <- -1.001
            } else {
                delta <- max(-0.5, -max((max.BB$BB8[2] + 0.001)/2, 0.001))
                theta1 <- max(-1.5, -max((max.BB$BB8[1] + 1.001)/2, 1.001))
            }
        } else if (family %in% allfams[tawns]) {
            ## Tawn
            # the folllowing gives a theoretical kendall's tau close to the empirical one
            delta <- min(abs(tau) + 0.1, 0.999)
            theta1 <- 1 + 6 * abs(tau)

            # check if data can be modeled by selected family
            if (family %in% c(104, 114)) {
                if (tau < 0) {
                    print("The Tawn or survival Tawn copula cannot be used for negatively dependent data.")
                    delta <- 1
                    theta1 <- 1.001
                }
            } else if (family %in% c(124, 134)) {
                if (tau > 0) {
                    print("The rotated Tawn copula cannot be used for positively dependent data.")
                    delta <- 1
                    theta1 <- -1.001
                } else theta1 <- -theta1

            } else if (family %in% c(204, 214)) {
                if (tau < 0) {
                    print("The Tawn2 or survival Tawn2 copula cannot be used for negatively dependent data.")
                    delta <- 1
                    theta1 <- 1.001
                }
            } else if (family %in% c(224, 234)) {
                if (tau > 0) {
                    print("The rotated Tawn2 copula cannot be used for positively dependent data.")
                    delta <- 1
                    theta1 <- -1.001
                } else theta1 <- -theta1
            }
        }

        ## maximum likelihood optimization
        if (family == 0) {
            theta <- 0
            se1 <- 0
            out <- list(value = 0)
        } else if (family < 100) {
            out <- MLE_intern(cbind(u1, u2),
                              c(theta1, delta),
                              family = family,
                              se,
                              max.df,
                              max.BB,
                              weights)
            theta <- out$par
            if (se == TRUE)
                se1 <- out$se
        } else if (family > 100) {
            # New
            out <- MLE_intern_Tawn(cbind(u1, u2),
                                   c(theta1, delta),
                                   family = family,
                                   se)
            theta <- out$par
            if (se == TRUE)
                se1 <- out$se
        }
    }

    ## store estimated parameters
    if (length(theta) == 1)
        theta <- c(theta, 0)
    obj <- BiCop(family, theta[1], theta[2])

    ## store standard errors (if asked for)
    if (se == TRUE) {
        if (length(se1) == 1)
            se1 <- c(se1, 0)
        obj$se <- se1[1]
        obj$se2 <- se1[2]
    }

    ## add more information about the fit
    obj$nobs   <- length(u1)
    # for method "itau" the log-likelihood hasn't been calculated yet
    obj$logLik <- switch(method,
                         "itau" = sum(log(BiCopPDF(u1, u2,
                                                   obj$family,
                                                   obj$par,
                                                   obj$par2,
                                                   check.pars = FALSE))),
                         "mle"  = out$value)
    obj$AIC    <- - 2 * obj$logLik + 2 * obj$npars
    obj$BIC    <- - 2 * obj$logLik + log(obj$nobs) * obj$npars
    obj$emptau <- tau
    obj$p.value.indeptest <- BiCopIndTest(u1, u2)$p.value

    ## store the call that created the BiCop object
    obj$call <- match.call()

    ## return results
    obj
}


## internal version without checking and option for reduced outpout
BiCopEst.intern <- function(u1, u2, family, method = "mle", se = TRUE, max.df = 30,
                            max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)),
                            weights = NA, as.BiCop = TRUE) {

    ## calculate empirical Kendall's tau and invert for initial estimate
    tau <- fasttau(u1, u2, weights)
    if (family %in% c(0, 2, allfams[onepar]))
        theta <- BiCopTau2Par(family, tau, check.taus = FALSE)

    ## inversion of kendall's tau -----------------------------
    if (method == "itau") {

        ## standard errors for method itau
        se1 <- 0
        if (se == TRUE) {
            p <- 2
            n <- length(u1)
            ec <- numeric(n)
            u <- cbind(u1, u2)
            v <- matrix(0, n, p * (p - 1)/2)

            if (family == 1)
                tauder <- function(x) {
                    2/(pi * sqrt(1 - x^2))
                } else if (family %in% c(3, 13, 23, 33)) {
                    tauder <- function(x) 2 * (2 + x)^(-2)
                } else if (family %in% c(4, 14, 24, 34)) {
                    tauder <- function(x) x^(-2)
                } else if (family == 5) {
                    f <- function(x) x/(exp(x) - 1)
                    tauder <- function(x) {
                        lwr <- 0 + .Machine$double.eps^0.5
                        intgrl <- integrate(f,
                                            lower = lwr,
                                            upper = x)$value
                        4/x^2 - 8/x^3 * intgrl + 4/(x * (exp(x) - 1))
                    }
                } else if (family %in% c(6, 16, 26, 36)) {
                    tauder <- function(x) {
                        euler <- 0.577215664901533
                        -((-2 + 2 * euler + 2 * log(2) + digamma(1/x) +
                               digamma(1/2 * (2 + x)/x) + x)/(-2 + x)^2) +
                            ((-trigamma(1/x)/x^2 + trigamma(1/2 * (2 + x)/x) *
                                  (1/(2 + x) - (2 + x)/(2 * x^2)) + 1)/(-2 + x))
                    }
                } else if (family %in% c(41, 51, 61, 71)) {
                    tauder <- function(x) {
                        2 * sqrt(pi) * gamma(0.5 + x) *
                            (digamma(1 + x) - digamma(0.5 + x))/gamma(1 + x)
                    }
                }

            l <- 1
            for (j in 1:(p - 1)) {
                for (i in (j + 1):p) {
                    for (k in 1:n)
                        ec[k] <- sum(u[, i] <= u[k, i] & u[, j] <= u[k, j])/n
                    v[, l] <- 2 * ec - u[, i] - u[, j]
                    l <- l + 1
                }
            }

            if (family == 0) {
                D <- 0
            } else if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 41, 51)) {
                D <- 1/tauder(theta)
            } else if (family %in% c(23, 33, 24, 34, 26, 36, 61, 71)) {
                D <- 1/tauder(-theta)
            }


            se1 <- as.numeric(sqrt(16/n * var(v %*% D)))
        }  # end if (se == TRUE)
    }  # end if (method == "itau")

    ## MLE ------------------------------------------
    if (method == "mle") {
        ## set starting parameters for maximum likelihood estimation
        theta1 <- 0
        delta <- 0

        if (!(family %in% c(2, 6, 7, 8, 9, 10,
                            17, 18, 19, 20,
                            27, 28, 29, 30,
                            37, 38, 39, 40,
                            104, 114, 124, 134,
                            204, 214, 224, 234))) {
            theta1 <- theta
        }
        if (family == 2) {
            ## t
            theta1 <- sin(tau * pi/2)
            delta <- 8
        } else if (family == 7 || family == 17) {
            ## BB1
            delta <- 1.5
            theta1 <- 0.5
        } else if (family == 27 || family == 37) {
            ## BB1
            delta <- -1.5
            theta1 <- -0.5
        } else if (family == 8 || family == 18) {
            ## BB6
            delta <- 1.5
            theta1 <- 1.5
        } else if (family == 28 || family == 38) {
            ## BB6
            delta <- -1.5
            theta1 <- -1.5
        } else if (family == 9 || family == 19) {
            ## BB7
            delta <- 0.5
            theta1 <- 1.5
        } else if (family == 29 || family == 39) {
            ## BB7
            delta <- max(-0.5, -max((max.BB$BB7[2] + 0.001)/2, 0.001))
            theta1 <- max(-1.5, -max((max.BB$BB7[1] + 1.001)/2, 1.001))
        } else if (family == 10 || family == 20) {
            ## BB8
            delta <- 0.5
            theta1 <- 1.5
        } else if (family == 30 || family == 40) {
            ## BB8
            delta <- -0.5
            theta1 <- -1.5
        } else if (family %in% allfams[tawns]) {
            ## Tawn
            # the folllowing gives a theoretical kendall's tau close
            #  to the empirical one
            delta <- min(abs(tau) + 0.1, 0.999)
            theta1 <- 1 + 6 * abs(tau)
            if (family %in% negfams)
                theta1 <- - theta1
        }

        ## maximum likelihood optimization
        if (family == 0) {
            out <- list(value = 0)
            theta <- 0
            se1 <- 0
        } else if (family < 100) {
            out <- MLE_intern(cbind(u1, u2),
                              c(theta1, delta),
                              family = family,
                              se,
                              max.df,
                              max.BB,
                              weights)
            theta <- out$par
            if (se == TRUE)
                se1 <- out$se
        } else if (family > 100) {
            # New
            out <- MLE_intern_Tawn(cbind(u1, u2),
                                   c(theta1, delta),
                                   family = family,
                                   se)
            theta <- out$par
            if (se == TRUE)
                se1 <- out$se
        }
    }

    ## store estimated parameters
    if (length(theta) == 1)
        theta <- c(theta, 0)
    if (!as.BiCop) {
        obj <- list(family = family, par = theta[1], par2 = theta[2])
        ## store standard errors (if asked for)
        if (se == TRUE) {
            if (length(se1) == 1)
                se1 <- c(se1, 0)
            obj$se <- se1[1]
            obj$se2 <- se1[2]
        }
    } else {
        obj <- BiCop(family, theta[1], theta[2])

        ## store standard errors (if asked for)
        if (se == TRUE) {
            if (length(se1) == 1)
                se1 <- c(se1, 0)
            obj$se <- se1[1]
            obj$se2 <- se1[2]
        }

        ## add more information about the fit
        obj$nobs   <- length(u1)
        # for method "itau" the log-likelihood hasn't been calculated yet
        obj$logLik <- switch(method,
                             "itau" = sum(log(BiCopPDF(u1, u2,
                                                       obj$family,
                                                       obj$par,
                                                       obj$par2,
                                                       check.pars = FALSE))),
                             "mle"  = out$value)
        obj$AIC    <- - 2 * obj$logLik + 2 * obj$npars
        obj$BIC    <- - 2 * obj$logLik + log(obj$nobs) * obj$npars
        obj$emptau <- tau
        obj$p.value.indeptest <- BiCopIndTest(u1, u2)$p.value

        ## store the call that created the BiCop object
        obj$call <- match.call()
    }

    ## return results
    obj
}




#############################################################
# bivariate MLE function
#
#------------------------------------------------------------
# INPUT:
#   data    Data for which to estimate parameter
#   start.parm	Start parameter for the MLE
#   Maxiter	max number of iterations
#   se		TRUE or FALSE
# OUTPUT:
#   out     Estimated Parameters and standard error (if se==TRUE)
#--------------------------------------------------------------
# Author: Ulf Schepsmeier
# Date: 2011-02-04
# Version: 1.1
#---------------------------------------------------------------

MLE_intern <- function(data, start.parm, family, se = FALSE, max.df = 30,
                       max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6),
                                     BB7 = c(5, 6), BB8 = c(6, 1)),
                       weights = NULL, cor.fixed = FALSE) {

    n <- dim(data)[1]
    if (any(is.na(weights)))
        weights <- NULL

    if (family %in% c(7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        t_LL <- function(param) {

            if (is.null(weights)) {
                ll <- .C("LL_mod2",
                         as.integer(family),
                         as.integer(n),
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param[1]),
                         as.double(param[2]),
                         as.double(0),
                         PACKAGE = "VineCopula")[[7]]
            } else {
                ll <- .C("LL_mod_seperate",
                         as.integer(family),
                         as.integer(n),
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param[1]),
                         as.double(param[2]),
                         as.double(rep(0, n)),
                         PACKAGE = "VineCopula")[[7]] %*% weights
            }

            if (is.infinite(ll) || is.na(ll) || ll < -10^250)
                ll <- -10^250

            return(ll)
        }

        if (family == 7 || family == 17) {
            low <- c(0.001, 1.001)
            up <- pmin(max.BB$BB1, c(7, 7))
        } else if (family == 8 || family == 18) {
            low <- c(1.001, 1.001)
            up <- pmin(max.BB$BB6, c(6, 8))
        } else if (family == 9 | family == 19) {
            low <- c(1.001, 0.001)
            up <- pmin(max.BB$BB7, c(6, 75))
        } else if (family == 10 | family == 20) {
            low <- c(1.001, 0.001)
            up <- pmin(max.BB$BB8, c(8, 1))
        } else if (family == 27 | family == 37) {
            up <- c(-0.001, -1.001)
            low <- -pmin(max.BB$BB1, c(7, 7))
        } else if (family == 28 | family == 38) {
            up <- c(-1.001, -1.001)
            low <- -pmin(max.BB$BB6, c(6, 8))
        } else if (family == 29 | family == 39) {
            up <- c(-1.001, -0.001)
            low <- -pmin(max.BB$BB7, c(6, 75))
        } else if (family == 30 | family == 40) {
            up <- c(-1.001, -0.001)
            low <- -pmin(max.BB$BB8, c(8, 1))
        }

        if (se == TRUE) {
            optimout <- optim(par = start.parm,
                              fn = t_LL,
                              method = "L-BFGS-B",
                              lower = low,
                              upper = up,
                              control = list(fnscale = -1, maxit = 500),
                              hessian = TRUE)
        } else {
            optimout <- optim(par = start.parm,
                              fn = t_LL,
                              method = "L-BFGS-B",
                              lower = low,
                              upper = up,
                              control = list(fnscale = -1, maxit = 500))
        }

    } else if (family == 2) {

        if (cor.fixed == FALSE) {

            t_LL <- function(param) {
                if (param[1] < -0.9999 | param[1] > 0.9999 | param[2] < 2.0001 | param[2] > max.df) {
                    ll <- -10^10
                } else {
                    if (is.null(weights)) {
                        ll <- .C("LL_mod2",
                                 as.integer(family),
                                 as.integer(n),
                                 as.double(data[, 1]),
                                 as.double(data[, 2]),
                                 as.double(param[1]),
                                 as.double(param[2]),
                                 as.double(0),
                                 PACKAGE = "VineCopula")[[7]]
                    } else {
                        ll <- .C("LL_mod_seperate",
                                 as.integer(family),
                                 as.integer(n),
                                 as.double(data[, 1]),
                                 as.double(data[, 2]),
                                 as.double(param[1]),
                                 as.double(param[2]),
                                 as.double(rep(0, n)),
                                 PACKAGE = "VineCopula")[[7]] %*% weights
                    }

                    if (is.infinite(ll) || is.na(ll) || ll < -10^10)
                        ll <- -10^10
                }
                return(ll)
            }

            gr_LL <- function(param) {
                gr <- rep(0, 2)
                gr[1] <- sum(BiCopDeriv(data[, 1],
                                        data[, 2],
                                        family = 2,
                                        par = param[1],
                                        par2 = param[2],
                                        deriv = "par",
                                        log = TRUE,
                                        check.pars = FALSE))
                gr[2] <- sum(BiCopDeriv(data[, 1],
                                        data[, 2],
                                        family = 2,
                                        par = param[1],
                                        par2 = param[2],
                                        deriv = "par2",
                                        log = TRUE,
                                        check.pars = FALSE))
                return(gr)
            }

            if (is.null(weights)) {
                if (se == TRUE) {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL,
                                      gr = gr_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE,
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                } else {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL,
                                      gr = gr_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                }
            } else {
                if (se == TRUE) {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE,
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                } else {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                }
            }

        } else {
            t_LL <- function(param) {
                if (is.null(weights)) {
                    ll <- .C("LL_mod2",
                             as.integer(family),
                             as.integer(n),
                             as.double(data[, 1]),
                             as.double(data[, 2]),
                             as.double(start.parm[1]),
                             as.double(param[1]),
                             as.double(0),
                             PACKAGE = "VineCopula")[[7]]
                } else {
                    ll <- .C("LL_mod_seperate",
                             as.integer(family),
                             as.integer(n),
                             as.double(data[, 1]),
                             as.double(data[, 2]),
                             as.double(start.parm[1]),
                             as.double(param[1]),
                             as.double(rep(0, n)),
                             PACKAGE = "VineCopula")[[7]] %*% weights
                }

                if (is.infinite(ll) || is.na(ll) || ll < -10^250)
                    ll <- -10^250

                return(ll)
            }

            gr_LL <- function(param) {
                gr <- sum(BiCopDeriv(data[, 1],
                                     data[, 2],
                                     family = 2,
                                     par = start.parm[1],
                                     par2 = param[1],
                                     deriv = "par2",
                                     log = TRUE,
                                     check.pars = FALSE))
                return(gr)
            }

            if (se == TRUE) {
                if (is.null(weights)) {
                    optimout <- optim(par = start.parm[2],
                                      fn = t_LL,
                                      gr = gr_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE,
                                      lower = 2.0001,
                                      upper = max.df)
                } else {
                    optimout <- optim(par = start.parm[2],
                                      fn = t_LL,
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE,
                                      lower = 2.0001,
                                      upper = max.df)
                }
            } else {
                optimout <- optimize(f = t_LL,
                                     maximum = TRUE,
                                     interval = c(2.0001, max.df))
                optimout$par <- optimout$maximum
                optimout$value <- optimout$objective
            }
            optimout$par <- c(0, optimout$par)

        }

    } else {

        t_LL <- function(param) {
            if (is.null(weights)) {
                ll <- .C("LL_mod2", as.integer(family),
                         as.integer(n),
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param),
                         as.double(0), as.double(0),
                         PACKAGE = "VineCopula")[[7]]
            } else {
                ll <- .C("LL_mod_seperate",
                         as.integer(family),
                         as.integer(n),
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param[1]),
                         as.double(0),
                         as.double(rep(0, n)),
                         PACKAGE = "VineCopula")[[7]] %*% weights
            }
            if (is.infinite(ll) || is.na(ll) || ll < -10^250)
                ll <- -10^250

            return(ll)
        }

        gr_LL <- function(param) {
            gr <- sum(BiCopDeriv(data[, 1],
                                 data[, 2],
                                 family = family,
                                 par = param,
                                 deriv = "par",
                                 log = TRUE,
                                 check.pars = FALSE))
            return(gr)
        }

        low <- -Inf
        up <- Inf

        if (family == 1) {
            low <- -0.9999
            up <- 0.9999
        } else if (family %in% c(3, 13)) {
            low <- 1e-04
            up <- 100
        } else if (family %in% c(4, 14)) {
            low <- 1.0001
            up <- 100
        } else if (family %in% c(5)) {
            low <- -100
            up <- 100
        } else if (family %in% c(6, 16)) {
            low <- 1.0001
            up <- 50
        } else if (family %in% c(23, 33)) {
            up <- -1e-04
            low <- -100
        } else if (family %in% c(24, 34)) {
            up <- -1.0001
            low <- -100
        } else if (family %in% c(26, 36)) {
            up <- -1.0001
            low <- -50
        } else if (family %in% c(41, 51)) {
            low <- 1e-04
            up <- BiCopTau2Par(family, 0.85)
            # if(t_LL(up)==-10^250) up=BiCopTau2Par(family,0.95) if(t_LL(up)==-10^250) up=BiCopTau2Par(family,0.9)
        } else if (family %in% c(61, 71)) {
            up <- -1e-04
            low <- BiCopTau2Par(family, -0.85)
            if (t_LL(low) == -10^250)
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250)
                low <- BiCopTau2Par(family, -0.9)
        }

        ## ensure that starting parameters are withing bounds
        start.parm[1] <- min(max(start.parm[1], low), up)
        optimout <- optimize(f = t_LL,
                             maximum = TRUE,
                             interval = c(low, up))
        optimout$par <- c(optimout$maximum, 0)
        optimout$value <- optimout$objective
        if (se == TRUE) {
            d2 <- BiCopDeriv2(data[, 1], data[, 1],
                              family,
                              optimout$par[1],
                              check.pars = FALSE)
            optimout$hessian <- sum(-d2)
        }

    }

    out <- list()

    if (se == TRUE) {
        if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
            out$par <- optimout$par

            if (!is.finite(det(optimout$hessian))) {
                var <- matrix(NA)
            } else if (det(optimout$hessian) == 0) {
                var <- diag(1, dim(optimout$hessian)[1])
            } else {
                var <- try((-solve(optimout$hessian)), silent = TRUE)
                if (inherits(var, 'try-error'))
                    var <- c(NA, NA)
            }
            out$se <- suppressWarnings(sqrt(diag(var)))

            if ((family == 2) && (out$par[2] >= (max.df - 1e-04)))
                out$se[2] <- NA

        } else {
            out$par <- optimout$par[1]

            if (!is.finite(optimout$hessian)) {
                var <- NA
            } else if (optimout$hessian == 0) {
                var <- 1
            } else {
                if (optimout$hessian < 0)
                    optimout$hessian <- -optimout$hessian
                var <- 1/optimout$hessian
            }

            out$se <- as.numeric(sqrt(var))
        }
    } else {
        if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
            out$par <- optimout$par
        } else {
            out$par[1] <- optimout$par[1]
        }
    }

    out$value <- optimout$value
    out
}


# New for Tawn

MLE_intern_Tawn <- function(data, start.parm, family, se = FALSE) {

    n <- dim(data)[1]

    ## set bounds for optimization
    tau <- fasttau(data[, 1], data[, 2])
    if (family == 104 || family == 114 || family == 204 || family == 214) {
        parlower <- c(1.001, max(tau - 0.1, 1e-04))
        parupper <- c(20, min(tau + 0.2, 0.99))
    } else if (family == 124 || family == 134 || family == 224 || family == 234) {
        parlower <- c(-20, max(-tau - 0.1, 1e-04))
        parupper <- c(-1.001, min(-tau + 0.2, 0.99))
    }

    ## log-liklihood function
    loglikfunc <- function(param) {
        ll <- .C("LL_mod2",
                 as.integer(family),
                 as.integer(n),
                 as.double(data[, 1]),
                 as.double(data[, 2]),
                 as.double(param[1]),
                 as.double(param[2]),
                 as.double(0),
                 PACKAGE = "VineCopula")[[7]]
        if (is.infinite(ll) || is.na(ll) || ll < -10^250)
            ll <- -10^250
        # print(param) print(ll)
        return(ll)
    }

    ## optimize log-likelihood
    out <- list()
    # print(start.parm)
    if (se == TRUE) {
        optimout <- optim(par = start.parm,
                          fn = loglikfunc,
                          method = c("L-BFGS-B"),
                          lower = parlower,
                          upper = parupper,
                          control = list(fnscale = -1, maxit = 500),
                          hessian = TRUE)
        if (!is.finite(det(optimout$hessian))) {
            var <- matrix(NA, 2, 2)
        } else if (det(optimout$hessian) == 0) {
            var <- diag(NA, dim(optimout$hessian)[1])
        } else {
            var <- try((-solve(optimout$hessian)), silent = TRUE)
            if (inherits(var, 'try-error'))
                var <- c(NA, NA)
        }
        out$se <- suppressWarnings(sqrt(diag(var)))
    } else {
        optimout <- optim(par = start.parm,
                          fn = loglikfunc,
                          method = c("L-BFGS-B"),
                          lower = parlower,
                          upper = parupper,
                          control = list(fnscale = -1, maxit = 500))
    }

    ## return results
    out$par <- optimout$par
    out$value <- optimout$value
    return(out)
}



fasttau <- function(x, y, weights = NA) {
    if (any(is.na(weights))) {
        m <- length(x)
        n <- length(y)
        if (m == 0 || n == 0)
            stop("both 'x' and 'y' must be non-empty")
        if (m != n)
            stop("'x' and 'y' must have the same length")
        out <- .C("ktau",
                  x = as.double(x),
                  y = as.double(y),
                  N = as.integer(n),
                  tau = as.double(0),
                  S = as.double(0),
                  D = as.double(0),
                  T = as.integer(0),
                  U = as.integer(0),
                  V = as.integer(0),
                  PACKAGE = "VineCopula")
        ktau <- out$tau
    } else {
        ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
    }
    return(ktau)
}
