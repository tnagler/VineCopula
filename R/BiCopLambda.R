############ plot of the theoretical and empirical lambda-function


#' Lambda-Function (Plot) for Bivariate Copula Data
#'
#' This function plots/returns the lambda-function of given bivariate copula
#' data.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative versions \cr
#' \preformatted{BiCopLambda(obj, PLOT = TRUE, ...)}
#' and
#' \preformatted{BiCopLambda((u1, u2, obj, PLOT = TRUE, ...)}
#' can be used.
#'
#' @param u1,u2 Data vectors of equal length with values in \eqn{[0,1]} (default:
#' `u1` and `u2 = NULL`).
#' @param family An integer defining the bivariate copula family or indicating
#' the empirical lambda-function: \cr
#' `"emp"` = empirical lambda-function
#' (default) \cr
#' `1` = Gaussian copula; the theoretical lambda-function is simulated
#' (no closed formula available) \cr
#' `2` = Student-t copula; the theoretical lambda-function is simulated
#' (no closed formula available) \cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `7` = BB1 copula \cr
#' `8` = BB6 copula \cr
#' `9` = BB7 copula \cr
#' `10` = BB8
#' copula
#' @param par Copula parameter; if the empirical lambda-function is chosen,
#' `par = NULL` or `0` (default).
#' @param par2 Second copula parameter for t-, BB1, BB6, BB7 and BB8 copulas
#' (default: `par2 = 0`).
#' @param PLOT Logical; whether the results are plotted. If `PLOT =
#' FALSE`, the values \cr `empLambda` and/or `theoLambda` are
#' returned (see below; default: `PLOT = TRUE`).
#' @param ... Additional plot arguments.
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#'
#' @return \item{empLambda}{If the empirical lambda-function is chosen and
#' `PLOT = FALSE`, a vector of the empirical lambda's is returned.}
#' \item{theoLambda}{If the theoretical lambda-function is chosen and
#' `PLOT = FALSE`, a vector of the theoretical lambda's is returned.}
#' @note The \eqn{\lambda}-function is characteristic for each bivariate copula
#' family and defined by Kendall's distribution function \eqn{K}:
#' \deqn{ \lambda(v,\theta) := v - K(v,\theta) } with
#' \deqn{K(v,\theta) := P(C_{\theta}(U_1,U_2) \leq v),\ \ v\in [0,1]. }{
#' K(v,\theta) := P(C_{\theta}(U_1,U_2) <= v), v \in [0,1]. }
#' For Archimedean copulas one has the following closed form expression
#' in terms of the generator function \eqn{\varphi}{\phi} of the copula
#' \eqn{C_{\theta}}:
#' \deqn{\lambda(v,\theta) = \frac{\varphi(v)}{\varphi '(v)}, }{
#' \lambda(v,\theta) = \phi(v) / \phi'(v), }
#' where \eqn{\varphi '}{\phi'} is the derivative of \eqn{\varphi}{\phi}. For
#' more details see Genest and Rivest (1993) or Schepsmeier (2010).
#'
#' For the bivariate Gaussian and Student-t copula no closed form expression for
#' the theoretical \eqn{\lambda}-function exists. Therefore it is simulated based
#' on samples of size 1000. For all other implemented copula families there are
#' closed form expressions available.
#'
#' The plot of the theoretical \eqn{\lambda}-function also shows the limits of
#' the \eqn{\lambda}-function corresponding to Kendall's tau \eqn{=0} and
#' Kendall's tau \eqn{=1} (\eqn{\lambda=0}).
#'
#' For rotated bivariate copulas one has to transform the input arguments
#' `u1` and/or `u2`. In particular, for copulas rotated by 90 degrees
#' `u1` has to be set to `1-u1`, for 270 degrees `u2` to
#' `1-u2` and for survival copulas `u1` and `u2` to `1-u1`
#' and `1-u2`, respectively. Then \eqn{\lambda}-functions for the
#' corresponding non-rotated copula families can be considered.
#'
#' @author Ulf Schepsmeier
#'
#' @seealso [BiCopMetaContour()], [BiCopKPlot()],
#' [BiCopChiPlot()], [BiCop()]
#'
#' @references Genest, C. and L.-P. Rivest (1993). Statistical inference
#' procedures for bivariate Archimedean copulas. Journal of the American
#' Statistical Association, 88 (423), 1034-1043.
#'
#' Schepsmeier, U. (2010). Maximum likelihood estimation of C-vine pair-copula
#' constructions based on bivariate copulas from different families. Diploma
#' thesis, Technische Universitaet Muenchen.\cr
#' <https://mediatum.ub.tum.de/?id=1079296>.
#'
#' @examples
#' \dontshow{set.seed(123)}
#' # simulate from Clayton copula
#' cop <- BiCop(3, tau = 0.5)
#' dat <- BiCopSim(1000, cop)
#'
#' # create lambda-function plots
#' op <- par(mfrow = c(1, 3))
#' BiCopLambda(dat[, 1], dat[, 2])  # empirical lambda-function
#' BiCopLambda(cop)	# theoretical lambda-function
#' BiCopLambda(dat[, 1], dat[, 2], cop)	# both
#' par(op)
#'
#'
BiCopLambda <- function(u1 = NULL, u2 = NULL, family = "emp", par = 0, par2 = 0,
                        PLOT = TRUE, obj = NULL, ...) {
    ## extract family and parameters if BiCop object is provided
    if (!is.null(obj)) {
        stopifnot(inherits(obj, "BiCop"))
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }
    if (inherits(u1, "BiCop")) {
        # for short hand usage extract from u1
        if (inherits(u2, "logical"))
            PLOT <- u2
        obj <- u1
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
        u1 <- NULL
    }
    if (inherits(family, "BiCop")) {
        # for short hand usage extract from family
        if (inherits(par, "logical"))
            PLOT <- par
        obj <- family
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }

    if ((is.null(u1) == TRUE) && (is.null(u2) == TRUE) && (family == 0 || par == 0))
        stop("Either 'u1' and 'u2' have to be set for the emp.
             lambda-function or 'family' and 'par' for the theo. lambda-function.")
    if (length(u1) != length(u2))
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (!(family %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "emp")))
        stop("Copula family not implemented.")
    if ((family %in% c(2, 7, 8, 9, 10)) && (par2 == 0))
        stop("For t-, BB1 and BB7 copulas, 'par2' must be set.")
    if (any(is.na(u1 + u2))) {
        # send warning message
        warning("Some of the data are NA. ",
                "Only complete observations are used.")
        # remove NAs
        na.ind <- which(is.na(u1 + u2))
        u1 <- u1[-na.ind]
        u2 <- u2[-na.ind]
    }
    if ((is.null(u1) == FALSE) && (any(u1 > 1) || any(u1 < 0)))
        stop("Data has be in the interval [0,1].")
    if ((is.null(u2) == FALSE) && (any(u2 > 1) || any(u2 < 0)))
        stop("Data has be in the interval [0,1].")

    if (PLOT != TRUE && PLOT != FALSE)
        stop("The parameter 'PLOT' has to be set to 'TRUE' or 'FALSE'.")



    ## check for parameter consistency
    if (family != "emp")
        BiCopCheck(family, par, par2)

    ## create grids
    if (!is.null(u1))
        v <- seq(0.001, 1, length.out = length(u1)) else v <- seq(0.001, 1, 0.001)
    v1 <- v
    theoLambda <- rep(0, length(v))
    lambdaFull <- rep(0, length(v))

    ## calculate theoretical lambda
    main <- ""
    for (i in 2:length(v)) {
        if (family == 3) {
            theoLambda[i] <- -v1[i] * (1 - v1[i]^(par))/par
            main <- "Clayton copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 4) {
            theoLambda[i] <- v1[i] * log(v1[i])/(par)
            main <- "Gumbel copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 5) {
            theoLambda[i] <- -log((1 - exp(-par))/(1 - exp(-par * v1[i]))) * (1 - exp(-par * v1[i]))/(par * exp(-par * v1[i]))
            main <- "Frank copula"
            lambdaFull[i] <- (-v1[i] * log(1/v1[i], exp(1)))
        } else if (family == 6) {
            theoLambda[i] <- (log(1 - (1 - v[i])^par) * (1 - (1 - v[i])^par))/(par * (1 - v[i])^(par - 1))
            main <- "Joe copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 7) {
            theta <- par
            delta <- par2
            theoLambda[i] <- -1/(theta * delta) * (v[i]^(-theta) - 1)/(v[i]^(-1 - theta))
            main <- "BB1 copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 8) {
            theta <- par
            delta <- par2
            theoLambda[i] <- -log(-(1 - v[i])^theta + 1) * (1 - v[i] - (1 - v[i])^(-theta) + (1 - v[i])^(-theta) * v[i])/(delta * theta)
            main <- "BB6 copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 9) {
            theta <- par
            delta <- par2
            theoLambda[i] <- -1/(theta * delta) * ((1 - (1 - v[i])^theta)^(-delta) - 1)/((1 - v[i])^(theta - 1) * (1 - (1 - v[i])^theta)^(-delta - 1))
            main <- "BB7 copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        } else if (family == 10) {
            theta <- par
            delta <- par2
            theoLambda[i] <- -log(((1 - v[i] * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - v[i] * delta - (1 - v[i] * delta)^(-theta) + (1 - v[i] * delta)^(-theta) * v[i] * delta)/(theta * delta)
            main <- "BB8 copula"
            lambdaFull[i] <- (v1[i] * log(v1[i], exp(1)))
        }

    }

    if (all(is.null(u1) == FALSE))
        len <- length(u1) else len <- 1000

    if (family == 1) {
        theoLambda <- gtLambda(1, par, len = len)
        main <- "Gaussian copula"
        lambdaFull <- gtLambda(1, 0, len = len)
    } else if (family == 2) {
        theoLambda <- gtLambda(2, c(par, par2), len = len)
        main <- "t-copula"
        lambdaFull <- gtLambda(2, c(0, par2), len = len)
    }

    ## calculate empircal lambda
    if (family == "emp" || is.null(u1) == FALSE) {
        n <- length(v)
        nn <- length(u1)
        empLambda <- rep(0, n)
        K <- rep(0, n)
        V <- rep(0, nn)

        # Berechnung der konkordanten Paare

        V[1:nn] <- mapply(function(x, y) length(which(x > u1 & y > u2)), u1, u2)


        # Berechnung des empirischen K's

        V1 <- V/(n - 1)
        V1 <- V1

        K <- sapply(v1, function(x) (1/nn) * length(which(V1[1:nn] <= x)))

        # Berechnung der emp. lambdas

        empLambda <- v1 - K

    }

    ## create plot
    if (PLOT) {
        if (family == "emp") {
            # only empirical one
            if (min(empLambda[2:(nn - 1)]) < (-0.4)) {
                plot(v1,
                     empLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.6, 0),
                     ...)
            } else {
                plot(v1,
                     empLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.4, 0),
                     ...)
            }

        } else if (family != "emp" && is.null(u1) == TRUE) {
            # only theoretical one
            if (min(theoLambda[2:999]) < (-0.4)) {
                plot(v1,
                     theoLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.6, 0),
                     main = main,
                     ...)
            } else {
                plot(v1,
                     theoLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.4, 0),
                     main = main,
                     ...)
            }

            lines(v1, lambdaFull, type = "l", lty = 2)
            lines(c(0, 1), c(0, 0), type = "l", lty = 2)
        } else {
            # both
            minLambda <- min(min(theoLambda[2:(nn - 1)]), min(empLambda))
            if (minLambda < (-0.4)) {
                plot(v1,
                     empLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.6, 0),
                     ...)
            } else {
                plot(v1,
                     empLambda,
                     type = "l",
                     ylab = expression(lambda(v)),
                     xlab = "v",
                     ylim = c(-0.4, 0),
                     ...)
            }

            lines(v1, theoLambda, xlab = "v", type = "l", lwd = 2, col = "grey")
            lines(v1, lambdaFull, type = "l", lty = 2)
            lines(c(0, 1), c(0, 0), type = "l", lty = 2)

        }
    } else {
        out <- list()
        if (family == "emp")
            out$empLambda <- empLambda else if (family != "emp" && is.null(u1) == TRUE) {
                out$theoLambda <- theoLambda
            } else {
                out$empLambda <- empLambda
                out$theoLambda <- theoLambda
            }
        return(out)
    }

}


###### lambda-function for Gaussian- and t-copula #
# Input:
# copula Copula family (1='N',2='t')
# param Parameter
# Output:
# lambda lambda-function #

gtLambda <- function(copula, param, len = 10000) {
    v <- seq(0.001, 1, length.out = len)
    v1 <- v
    n <- length(v)
    nn <- len
    lambda <- rep(0, n)
    K <- rep(0, n)
    V <- rep(0, nn)
    mu <- c(0, 0)

    if (param[1] == 0) {
        for (i in 1:n) {
            lambda[i] <- (v1[i] * log(v1[i], exp(1)))
        }
    } else {
        if (copula == 1) {
            rho <- param
            # sigma=matrix(c(1,rho,rho,1), c(2,2)) Z=rmvnorm(nn, mu, sigma) u1=pnorm(Z[,1])
            # u2=pnorm(Z[,2])

            uu1 <- runif(nn)
            vv2 <- runif(nn)
            uu2 <- .C("Hinv1",
                      as.integer(1),
                      as.integer(nn),
                      as.double(vv2),
                      as.double(uu1),
                      as.double(rho),
                      as.double(0),
                      as.double(rep(0, nn)),
                      PACKAGE = "VineCopula")[[7]]


            # x1=qnorm(U1) x2=qnorm(U2) C=dmvnorm(cbind(x1,x2), mu, sigma)
        } else if (copula == 2) {
            rho <- param[1]
            nu <- param[2]

            # sigma=matrix(c(1,rho,rho,1), c(2,2)) Z=rmvt(nn, sigma, nu) u1=pt(Z[,1], df=nu)
            # u2=pt(Z[,2], df=nu)

            uu1 <- runif(nn)
            vv2 <- runif(nn)
            uu2 <- .C("Hinv1",
                      as.integer(2),
                      as.integer(nn),
                      as.double(vv2),
                      as.double(uu1),
                      as.double(rho),
                      as.double(nu),
                      as.double(rep(0, nn)),
                      PACKAGE = "VineCopula")[[7]]

        }


        # Berechnung der konkordanten Paare
        V[1:nn] <- mapply(function(x, y) length(which(x > uu1 & y > uu2)), uu1, uu2)

        # Berechnung des empirischen K's

        V1 <- V/(n - 1)
        V1 <- V1

        K <- sapply(v1, function(x) (1/nn) * length(which(V1[1:nn] <= x)))

        # Berechnung der emp. lambdas

        lambda <- v1 - K

        # lambda=smooth(lambda)
    }
    return(lambda)
}
