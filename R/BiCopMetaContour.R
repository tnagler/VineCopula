#' Contour Plot of Bivariate Meta Distribution
#'
#' Note: This function is deprecated and only available for backwards
#' compatibility. See \code{\link{contour.BiCop}} for contour plots of
#' parametric copulas, and \code{\link{BiCopKDE}} for kernel estimates.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1] (default:
#' \code{u1} and \code{u2 = NULL}).
#' @param bw Bandwidth (smoothing factor; default: \code{bw = 1}).
#' @param size Number of grid points; default: \code{size = 100}.
#' @param levels Vector of contour levels. For Gaussian, Student-t or
#' exponential margins the default value (\code{levels = c(0.01, 0.05, 0.1,
#' 0.15, 0.2)}) typically is a good choice. For uniform margins we
#' recommend\cr \code{levels = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)}\cr
#' and for Gamma margins\cr \code{levels = c(0.005, 0.01, 0.03, 0.05, 0.07,
#' 0.09)}.
#' @param family An integer defining the bivariate copula family or indicating
#' an empirical contour plot: \cr
#' \code{"emp"} = empirical contour plot
#' (default; margins can be specified by \code{margins}) \cr
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
#' @param par Copula parameter; if empirical contour plot, \code{par = NULL} or
#' \code{0} (default).
#' @param par2 Second copula parameter for t-, BB1, BB6, BB7, BB8, Tawn type 1
#' and type 2 copulas (default: \code{par2 = 0}).
#' @param PLOT Logical; whether the results are plotted.  If \code{PLOT =
#' FALSE}, the values \code{x}, \code{y} and \code{z} are returned (see below;
#' default: \code{PLOT = TRUE}).
#' @param margins Character; margins for the bivariate copula contour plot.
#' Possible margins are:\cr
#' \code{"norm"} = standard normal margins (default)\cr
#' \code{"t"} = Student t margins with degrees of freedom as
#' specified by \code{margins.par}\cr
#' \code{"gamma"} = Gamma margins with shape and scale as
#' specified by \code{margins.par}\cr
#' \code{"exp"} = Exponential margins with rate as
#' specified by \code{margins.par}\cr
#' \code{"unif"} = uniform margins
#' @param margins.par Parameter(s) of the distribution of the margins if
#' necessary (default: \code{margins.par = 0}), i.e.,
#' \itemize{
#' \item a positive real number for the degrees of freedom of
#' Student t margins (see \code{\link{dt}}),
#' \item a 2-dimensional vector of positive real numbers for
#' the shape and scale parameters of Gamma margins (see \code{\link{dgamma}}),
#' \item a positive real number for the rate parameter of
#' exponential margins (see \code{\link{dexp}}).
#' }
#' @param xylim A 2-dimensional vector of the x- and y-limits.  By default
#' (\code{xylim = NA}) standard limits for the selected margins are used.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param ... Additional plot arguments.
#'
#' @return \item{x}{A vector of length \code{size} with the x-values of the
#' kernel density estimator with Gaussian kernel if the empirical contour plot
#' is chosen and a sequence of values in \code{xylim} if the theoretical
#' contour plot is chosen.}
#' \item{y}{A vector of length \code{size} with the
#' y-values of the kernel density estimator with Gaussian kernel if the
#' empirical contour plot is chosen and a sequence of values in \code{xylim} if
#' the theoretical contour plot is chosen.}
#' \item{z}{A matrix of dimension
#' \code{size} with the values of the density of the meta distribution with
#' chosen margins (see \code{margins} and \code{margins.par}) evaluated at the
#' grid points given by \code{x} and \code{y}.}
#'
#' @note The combination \code{family = 0} (independence copula) and
#' \code{margins = "unif"} (uniform margins) is not possible because all
#' \code{z}-values are equal.
#'
#' @author Ulf Schepsmeier, Alexander Bauer
#'
#' @seealso \code{\link{BiCopChiPlot}}, \code{\link{BiCopKPlot}},
#' \code{\link{BiCopLambda}}
#'
#' @examples
#' ## meta Clayton distribution  with Gaussian margins
#' cop <- BiCop(family = 1, tau = 0.5)
#' BiCopMetaContour(obj = cop, main = "Clayton - normal margins")
#' # better:
#' contour(cop, main = "Clayton - normal margins")
#'
#' ## empirical contour plot with standard normal margins
#' dat <- BiCopSim(1000, cop)
#' BiCopMetaContour(dat[, 1], dat[, 2], bw = 2, family = "emp",
#'                  main = "empirical - normal margins")
#' # better:
#' BiCopKDE(dat[, 1], dat[, 2],
#'         main = "empirical - normal margins")
#'
#' ## empirical contour plot with exponential margins
#' BiCopMetaContour(dat[, 1], dat[, 2], bw = 2,
#'                  main = "empirical - exponential margins",
#'                  margins = "exp", margins.par = 1)
#' # better:
#' BiCopKDE(dat[, 1], dat[, 2],
#'          main = "empirical - exponential margins",
#'          margins = "exp")
#'
BiCopMetaContour <- function(u1 = NULL, u2 = NULL, bw = 1, size = 100,
                             levels = c(0.01, 0.05, 0.1, 0.15, 0.2), family = "emp",
                             par = 0, par2 = 0, PLOT = TRUE, margins = "norm",
                             margins.par = 0, xylim = NA, obj = NULL,...) {
    warning("This function is deprecated. ",
            "See ?contour.BiCop for contour plots of parametric copulas\n",
            "  and ?BiCopKDE for kernel estimates.")

    ## preprocessing of arguments
    if ((family == "emp") & is.null(obj)) {
        args <- preproc(c(as.list(environment()), call = match.call()),
                        check_u,
                        remove_nas,
                        check_if_01,
                        na.txt = " Only complete observations are used.")
    } else {
        args <- preproc(c(as.list(environment()),
                          call = match.call(),
                          check.pars = TRUE),
                        extract_from_BiCop,
                        check_fam_par)
    }
    list2env(args, environment())

    # check plot option
    if (PLOT != TRUE && PLOT != FALSE)
        stop("The parameter 'PLOT' has to be set to 'TRUE' or 'FALSE'.")
    # limits for size parameter
    if (size > 1000)
        stop("Size parameter should not be greater than 1000. Otherwise computational time and memory space are too large.")
    if (size < 50)
        stop("Size parameter should not be smaller than 50.")
    # limits bandwidth parameter
    if (bw < 1)
        stop("The bandwidth parameter 'bw' should be greater or equal to 1.")
    if (bw > 5)
        stop("The bandwidth parameter 'bw' should not be greater than 5.")

    ## check for appropriate call w.r.t. margins
    if (margins != "norm" && margins != "t" && margins != "exp" && margins != "gamma" &&
        margins != "unif")
        stop("The function only supports Gaussian ('norm'), Student t ('t'), exponential ('exp'), Gamma ('gamma') and uniform ('unif') margins.")
    if (margins == "t" && margins.par <= 0)
        stop("The degrees of freedom parameter for the Student t margins has to positive.")
    if (margins == "Gamma" && length(margins.par) != 2)
        stop("For Gamma margins two parameters are required in 'margins.par'.")
    if (margins == "exp" && margins.par == 0)
        stop("Exponential margins require one parameter in 'margins.par'.")
    if (margins == "unif" && family == 0)
        stop("The combination independence copula and uniform margins is not possible because all z-values are equal.")

    ## set margins for theoretical contour plot
    if (is.null(u1) && is.null(u2) && family != "emp") {
        u1 <- runif(1000)
        u2 <- runif(1000)
    }
    if (!is.na(xylim) && length(xylim) != 2)
        stop("'xylim' has to be a vector of length 2.")

    ## transform grid marginally
    if (margins == "norm") {
        x1 <- qnorm(p = u1)
        x2 <- qnorm(p = u2)
        if (any(is.na(xylim)))
            xylim <- c(-3, 3)
    } else if (margins == "t") {
        x1 <- qt(p = u1, df = margins.par)
        x2 <- qt(p = u2, df = margins.par)
        if (any(is.na(xylim)))
            xylim <- c(-3, 3)
    } else if (margins == "exp") {
        x1 <- qexp(p = u1, rate = margins.par)
        x2 <- qexp(p = u2, rate = margins.par)
        if (any(is.na(xylim)))
            xylim <- c(0, 5)
    } else if (margins == "gamma") {
        x1 <- qgamma(p = u1, shape = margins.par[1], scale = margins.par[2])
        x2 <- qgamma(p = u2, shape = margins.par[1], scale = margins.par[2])
        if (any(is.na(xylim)))
            xylim <- c(0, 5)
    } else if (margins == "unif") {
        x1 <- u1
        x2 <- u2
        if (any(is.na(xylim)))
            xylim <- c(0, 1)
    }

    x <- y <- seq(from = xylim[1], to = xylim[2], length.out = size)

    if (family != "emp") {
        ## calculate theoretical contours
        if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38,
                          39, 40, 42, 52, 62, 72, 104, 114, 124, 134, 204, 214, 224, 234)) {
            z <- matrix(data = meta.dens(x1 = rep(x = x, each = size),
                                         x2 = rep(x = y, times = size),
                                         param = c(par, par2),
                                         copula = family,
                                         margins = margins,
                                         margins.par = margins.par),
                        nrow = size,
                        byrow = TRUE)
        } else {
            z <- matrix(data = meta.dens(x1 = rep(x = x, each = size),
                                         x2 = rep(x = y, times = size),
                                         param = par,
                                         copula = family,
                                         margins = margins,
                                         margins.par = margins.par),
                        nrow = size,
                        byrow = TRUE)
        }
    } else {
        ## calculate empirical contours
        bw1 <- bw * bandwidth.nrd(x1)
        bw2 <- bw * bandwidth.nrd(x2)

        kd.est <- kde2d(x = x1, y = x2, h = c(bw1, bw2), n = size)

        x <- kd.est$x
        y <- kd.est$y
        z <- kd.est$z
    }

    if (PLOT) {
        ## plot contour lines
        contour(x = x,
                y = y,
                z = z,
                levels = levels,
                ylim = xylim,
                xlim = xylim,
                ...)
    } else {
        ## output bivarate meta density z(x,y)
        out <- list()
        out$x <- x
        out$y <- y
        out$z <- z

        return(out)
    }
}


### gen
# Input:
# u data vector
# param copula parameters
# copula copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)
# Output:
# out generator

gen <- function(u, param, copula) {
    # param == c(theta, delta)
    out <- numeric(length(u))

    if (copula == 7) {
        out <- (u^(-param[1]) - 1)^param[2]
    } else if (copula == 8) {
        out <- (-log(-(1 - u)^param[1] + 1))^param[2]
    } else if (copula == 9) {
        out <- (1 - (1 - u)^param[1])^(-param[2]) - 1
    } else if (copula == 10) {
        out <- -log((1 - (1 - param[2] * u)^param[1])/(1 - (1 - param[2])^param[1]))
    }

    return(out)
}


### gen.inv
# Input:
# u data vector
# param copula parameters
# copula copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)
# Output:
# out inverse generator

gen.inv <- function(u, param, copula) {
    out <- numeric(length(u))

    if (copula == 7) {
        out <- (1 + u^(1/param[2]))^(-1/param[1])
    } else if (copula == 8) {
        out <- 1 - (1 - exp(-u^(1/param[2])))^(1/param[1])
    } else if (copula == 9) {
        out <- 1 - (1 - (1 + u)^(-1/param[2]))^(1/param[1])
    } else if (copula == 10) {
        out <- 1/param[2] * (1 - (1 - (1 - (1 - param[2])^param[1]) * exp(-u))^(1/param[1]))
    }

    return(out)
}


### gen.drv
# Input:
# u data vector
# param copula parameters
# copula copula  family (7=BB1, 8=BB6, 9=BB7, 10=BB8)
# Output:
# out First derivative of the generator

gen.drv <- function(u, param, copula) {
    out <- numeric(length(u))

    if (copula == 7) {
        out <- -prod(param) * (u^-(param[1]) - 1)^(param[2] - 1) * u^(-1 - param[1])
    } else if (copula == 8) {
        out <- (-log(-(1 - u)^param[1] + 1))^(param[2] - 1) * param[2] * (1 - u)^(param[1] - 1) * param[1]/((1 - u)^param[1] - 1)
    } else if (copula == 9) {
        out <- -prod(param) * (1 - u)^(param[1] - 1) * (1 - (1 - u)^param[1])^(-1 - param[2])
    } else if (copula == 10) {
        out <- -prod(param) * ((1 - param[2] * u)^(param[1] - 1))/(1 - (1 - param[2] * u)^param[1])
    }

    return(out)
}


#### gen.drv2
# Input:
# u data vector
# param copula parameters
# copula copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)
# Output:
# out Second derivative of the generator

gen.drv2 <- function(u, param, copula) {
    out <- numeric(length(u))

    if (copula == 7) {
        out <- prod(param) * u^(-2 - param[1]) * (u^(-param[1]) - 1)^(param[2] - 2) * ((1 + prod(param)) * u^(-param[1]) - param[1] - 1)
    } else if (copula == 8) {
        out <- (prod(param) * ((-log(-(1 - u)^param[1] + 1))^(param[2] - 2) * (1 - u)^(2 * param[1] - 2) * prod(param) - (-log(-(1 - u)^param[1] + 1))^(param[2] - 2) * (1 - u)^(2 * param[1] - 2) * param[1] - (-log(-(1 - u)^param[1] +  1))^(param[2] - 1) * (1 - u)^(param[1] - 2) * param[1] - (-log(-(1 - u)^param[1] + 1))^(param[2] - 1) * (1 - u)^(2 * param[1] - 2) + (-log(-(1 - u)^param[1] + 1))^(param[2] - 1) * (1 - u)^(param[1] - 2)))/((1 - u)^param[1] - 1)^2
    } else if (copula == 9) {
        out <- prod(param) * (1 - u)^(param[1] - 2) * (1 - (1 - u)^param[1])^(-2 - param[2]) * ((1 + prod(param)) * (1 - u)^param[1] + param[1] - 1)
    } else if (copula == 10) {
        out <- (param[2]^2 * param[1] * ((1 - u * param[2])^(param[1] - 2) * param[1] + (1 - u * param[2])^(2 * param[1] - 2) - (1 - u * param[2])^(param[1] - 2)))/(((1 - u * param[2])^param[1] - 1)^2)
    }

    return(out)
}


### cop.cdf
# Input:
# u1,u2 data vectors
# param copula parameters
# copula copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)
# Output:
# out copula

cop.cdf <- function(u1, u2, param, copula) {
    return(gen.inv(u = gen(u = u1, param = param, copula = copula) + gen(u = u2, param = param, copula = copula),
                   param = param,
                   copula = copula))
}


######################### pdf for BB6 #

bb6pdf <- function(u, v, th, de) {
    t1 <- 1 - u
    t2 <- t1^th
    t3 <- 1 - t2
    t4 <- log(t3)
    t5 <- (-t4)^de
    t12 <- 1/de
    t16 <- 1/th
    t32 <- de - 1
    t38 <- 2 * de
    t39 <- -1 + t38
    t40 <- (-t4)^t39
    t47 <- 3 * de - 1
    t50 <- (-t4)^t32
    t61 <- (-t4)^t47
    t90 <- (-t4)^t38
    # above depend on u and parameters only loop below for solving for conditional
    # quantile
    t6 <- 1 - v
    t7 <- t6^th
    t8 <- 1 - t7
    t9 <- log(t8)
    t10 <- (-t9)^de
    t11 <- t5 + t10
    t13 <- t11^t12
    t14 <- exp(-t13)
    t15 <- 1 - t14
    t17 <- t15^t16
    t35 <- t11^(-2 * t32 * t12)
    t36 <- t35 * th
    t37 <- exp(t13)
    t42 <- (-t9)^t39
    t48 <- (-t9)^t47
    t53 <- t13 * de
    t56 <- (-t9)^t32
    t57 <- t37 * t50 * t56
    t59 <- t13 * th
    t78 <- t37 - 1
    t80 <- (t78 * t14)^t16
    t87 <- t78 * t78
    t93 <- (-t9)^t38
    # c21 = -t17*t13*t5*t2/t1/t3/t4/t11*t14/t15;
    pdf <- (2 * t36 * t37 * t40 * t42 + t36 * t37 * t48 * t50 + t53 * th * t57 - t59 * t57 + t36 * t37 * t61 * t56 - 2 * t35 * t40 * t42 - t35 * t61 * t56 -  t53 * th * t50 * t56 + t59 * t50 * t56 - t35 * t48 * t50) * t80 * t7 * t2/t3/t8/t87/(t90 + 2 * t5 * t10 + t93)/t1/t6
    pdf
}



#### cop.pdf
# Input:
# u1,u2 data vectors
# param copula parameters
# copula copula
# family (1,2,3,...14)
# Output:
# out copula density #

cop.pdf <- function(u1, u2, param, copula) {
    if (copula == 7 | copula == 9 | copula == 10) {
        return(-gen.drv2(u = cop.cdf(u1 = u1, u2 = u2, param = param, copula = copula), param = param, copula = copula) *
                   gen.drv(u = u1, param = param, copula = copula) *
                   gen.drv(u = u2, param = param, copula = copula)/
                   gen.drv(u = cop.cdf(u1 = u1, u2 = u2, param = param, copula = copula), param = param, copula = copula)^3)
    } else if (copula == 8) {
        return(bb6pdf(u1, u2, param[1], param[2]))
    } else if (copula == 17 | copula == 19 | copula == 20) {
        d1 <- 1 - u1
        d2 <- 1 - u2
        return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 10), param = param, copula = copula - 10) *
                   gen.drv(u = d1, param = param, copula = copula - 10) *
                   gen.drv(u = d2, param = param, copula = copula - 10)/
                   gen.drv(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 10), param = param, copula = copula - 10)^3)
    } else if (copula == 18) {
        d1 <- 1 - u1
        d2 <- 1 - u2
        return(bb6pdf(d1, d2, param[1], param[2]))
    } else if (copula == 27 | copula == 29 | copula == 30) {
        d1 <- 1 - u1
        d2 <- u2
        param <- -param
        return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 20), param = param,  copula = copula - 20) *
                   gen.drv(u = d1, param = param, copula = copula - 20) *
                   gen.drv(u = d2, param = param, copula = copula - 20)/
                   gen.drv(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 20), param = param, copula = copula - 20)^3)
    } else if (copula == 28) {
        d1 <- 1 - u1
        d2 <- u2
        return(bb6pdf(d1, d2, -param[1], -param[2]))
    } else if (copula == 37 | copula == 39 | copula == 40) {
        d1 <- u1
        d2 <- 1 - u2
        param <- -param
        return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 30), param = param, copula = copula - 30) *
                   gen.drv(u = d1, param = param, copula = copula - 30) *
                   gen.drv(u = d2, param = param, copula = copula - 30)/
                   gen.drv(u = cop.cdf(u1 = d1,  u2 = d2, param = param, copula = copula - 30), param = param, copula = copula - 30)^3)
    } else if (copula == 38) {
        d1 <- u1
        d2 <- 1 - u2
        return(bb6pdf(d1, d2, -param[1], -param[2]))
    } else if (copula == 0) {
        # independent
        return(rep(1, length(u1)))
    } else if (copula == 1) {
        # Gaussian
        t1 <- qnorm(p = u1)
        t2 <- qnorm(p = u2)
        rho <- param
        return(1/sqrt(1 - rho^2) * exp(-(rho^2 * (t1^2 + t2^2) - 2 * rho * t1 * t2)/(2 * (1 - rho^2))))
    } else if (copula == 2) {
        # t-copula
        rho <- param[1]
        nu <- param[2]

        t1 <- qt(u1, nu)
        t2 <- qt(u2, nu)
        return(1/(2 * pi * sqrt(1 - rho^2) * dt(t1, nu) * dt(t2, nu)) * (1 + (t1^2 + t2^2 - 2 * rho * t1 * t2)/(nu * (1 - rho^2)))^(-(nu + 2)/2))
    } else if (copula == 3) {
        # Clayton
        theta <- param
        return((1 + theta) * (u1 * u2)^(-1 - theta) * (u1^(-theta) + u2^(-theta) - 1)^(-2 - 1/theta))
    } else if (copula == 4) {
        # Gumbel
        theta <- param
        t1 <- (-log(u1))^(theta) + (-log(u2))^(theta)
        t2 <- exp(-t1^(1/theta))
        return(t2/(u1 * u2) * t1^(-2 + 2/theta) * (log(u1) * log(u2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
    } else if (copula == 5) {
        # Frank
        theta <- param
        return((theta * (exp(theta) - 1) * exp(theta * u2 + theta * u1 + theta))/(exp(theta * u2 + theta * u1) - exp(theta * u2 + theta) - exp(theta * u1 + theta) + exp(theta))^2)
    } else if (copula == 6) {
        # Joe
        theta <- param
        return(((1 - u1)^(theta) + (1 - u2)^(theta) - (1 - u1)^(theta) * (1 - u2)^(theta))^(1/(theta) - 2) * (1 - u1)^(theta - 1) * (1 - u2)^(theta - 1) * (theta - 1 + (1 -
                                                                                                                                                                             u1)^(theta) + (1 - u2)^(theta) - (1 - u1)^(theta) * (1 - u2)^(theta)))
    } else if (copula == 13) {
        # rotated Clayton (180)
        theta <- param
        d1 <- 1 - u1
        d2 <- 1 - u2
        return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
    } else if (copula == 14) {
        # rotated Gumbel (180)
        theta <- param
        d1 <- 1 - u1
        d2 <- 1 - u2
        t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
        t2 <- exp(-t1^(1/theta))
        return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
    } else if (copula == 16) {
        # rotated Joe (180)
        theta <- param
        d1 <- 1 - u1
        d2 <- 1 - u2
        return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 -
                                                                                                                                                                             d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta)))
    } else if (copula == 23) {
        # rotated Clayton (90)
        theta <- -param
        d1 <- 1 - u1
        d2 <- u2
        return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
    } else if (copula == 24) {
        # rotated Gumbel (90)
        theta <- -param
        d1 <- 1 - u1
        d2 <- u2
        t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
        t2 <- exp(-t1^(1/theta))
        return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
    } else if (copula == 26) {
        # rotated Joe (90)
        theta <- -param
        d1 <- 1 - u1
        d2 <- u2
        return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 -
                                                                                                                                                                             d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta)))
    } else if (copula == 33) {
        # rotaed Clayton (270)
        theta <- -param
        d1 <- u1
        d2 <- 1 - u2
        return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
    } else if (copula == 34) {
        # rotated Gumbel (270)
        theta <- -param
        d1 <- u1
        d2 <- 1 - u2
        t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
        t2 <- exp(-t1^(1/theta))
        return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
    } else if (copula == 36) {
        # rotated Joe (270)
        theta <- -param
        d1 <- u1
        d2 <- 1 - u2
        return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 -
                                                                                                                                                                             d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta)))
    } else if (copula == 41) {
        # New: Archimedean copula based on integrated positive stable LT; reflection
        # asymmetric copula (from Harry Joe)
        de <- param
        tem1 <- qgamma(1 - u1, de)
        tem2 <- qgamma(1 - u2, de)
        con <- gamma(1 + de)/de
        sm <- tem1^de + tem2^de
        tem <- sm^(1/de)
        pdf <- con * tem * exp(-tem + tem1 + tem2)/sm
        return(pdf)
    } else if (copula == 51) {
        # rotated reflection asymmetric copula (180)
        de <- param
        d1 <- 1 - u1
        d2 <- 1 - u2
        tem1 <- qgamma(1 - d1, de)
        tem2 <- qgamma(1 - d2, de)
        con <- gamma(1 + de)/de
        sm <- tem1^de + tem2^de
        tem <- sm^(1/de)
        pdf <- con * tem * exp(-tem + tem1 + tem2)/sm
        return(pdf)
    } else if (copula == 61) {
        # rotated reflection asymmetric copula (90)
        de <- -param
        d1 <- 1 - u1
        d2 <- u2
        tem1 <- qgamma(1 - d1, de)
        tem2 <- qgamma(1 - d2, de)
        con <- gamma(1 + de)/de
        sm <- tem1^de + tem2^de
        tem <- sm^(1/de)
        pdf <- con * tem * exp(-tem + tem1 + tem2)/sm
        return(pdf)
    } else if (copula == 71) {
        # rotated reflection asymmetric copula (270)
        de <- -param
        d1 <- u1
        d2 <- 1 - u2
        tem1 <- qgamma(1 - d1, de)
        tem2 <- qgamma(1 - d2, de)
        con <- gamma(1 + de)/de
        sm <- tem1^de + tem2^de
        tem <- sm^(1/de)
        pdf <- con * tem * exp(-tem + tem1 + tem2)/sm
        return(pdf)
    } else if (copula == 42) {
        # 2-parametric asymmetric copula (thanks to Benedikt Graeler)
        a <- param[1]
        b <- param[2]
        return(pmax(a * u2 * (((12 - 9 * u1) * u1 - 3) * u2 + u1 * (6 * u1 - 8) + 2) + b * (u2 * ((u1 * (9 * u1 - 12) + 3) * u2 + (12 - 6 * u1) * u1 -  4) - 2 * u1 + 1) + 1, 0))
    } else if (copula == 52) {
        # rotated 2-parametric asymmetric copula (180)
        a <- param[1]
        b <- param[2]
        d1 <- 1 - u1
        d2 <- 1 - u2
        return(pmax(a * d2 * (((12 - 9 * d1) * d1 - 3) * d2 + d1 * (6 * d1 - 8) + 2) + b * (d2 * ((d1 * (9 * d1 - 12) + 3) * d2 + (12 - 6 * d1) * d1 - 4) - 2 * d1 + 1) + 1, 0))
    } else if (copula == 62) {
        # rotated 2-parametric asymmetric copula (180)
        a <- -param[1]
        b <- -param[2]
        d1 <- 1 - u1
        d2 <- u2
        return(pmax(a * d2 * (((12 - 9 * d1) * d1 - 3) * d2 + d1 * (6 * d1 - 8) + 2) + b * (d2 * ((d1 * (9 * d1 - 12) + 3) * d2 + (12 - 6 * d1) * d1 - 4) - 2 * d1 + 1) + 1, 0))
    } else if (copula == 52) {
        # rotated 2-parametric asymmetric copula (180)
        a <- -param[1]
        b <- -param[2]
        d1 <- u1
        d2 <- 1 - u2
        return(pmax(a * d2 * (((12 - 9 * d1) * d1 - 3) * d2 + d1 * (6 * d1 - 8) + 2) + b * (d2 * ((d1 * (9 * d1 - 12) + 3) * d2 + (12 - 6 * d1) * d1 - 4) - 2 * d1 + 1) + 1, 0))
    } else if (copula == 104) {
        # Tawn copula (psi2 fix)
        par <- param[1]
        par2 <- param[2]
        fam <- 104
        return(BiCopPDF(u1, u2, fam, par, par2))
    } else if (copula == 114) {
        # Tawn copula (psi2 fix)
        par <- param[1]
        par2 <- param[2]
        fam <- 104
        return(BiCopPDF(1 - u1, 1 - u2, fam, par, par2))
    } else if (copula == 124) {
        # Tawn copula (psi2 fix)
        par <- -param[1]
        par2 <- param[2]
        fam <- 104
        return(BiCopPDF(1 - u1, u2, fam, par, par2))
    } else if (copula == 134) {
        # Tawn copula (psi2 fix)
        par <- -param[1]
        par2 <- param[2]
        fam <- 104
        return(BiCopPDF(u1, 1 - u2, fam, par, par2))
    } else if (copula == 204) {
        # Tawn copula (psi1 fix)
        par <- param[1]
        par2 <- param[2]
        fam <- 204
        return(BiCopPDF(u1, u2, fam, par, par2))
    } else if (copula == 214) {
        # Tawn copula (psi1 fix)
        par <- param[1]
        par2 <- param[2]
        fam <- 204
        return(BiCopPDF(1 - u1, 1 - u2, fam, par, par2))
    } else if (copula == 224) {
        # Tawn copula (psi1 fix)
        par <- -param[1]
        par2 <- param[2]
        fam <- 204
        return(BiCopPDF(1 - u1, u2, fam, par, par2))
    } else if (copula == 234) {
        # Tawn copula (psi1 fix)
        par <- -param[1]
        par2 <- param[2]
        fam <- 204
        return(BiCopPDF(u1, 1 - u2, fam, par, par2))
    }
}



### Density of a meta-distribution with normal margins
# Input:
# x1,x2 vectors
# param Copula parameter(s)
# copula copula family
# Output:
# density


meta.dens <- function(x1, x2, param, copula, margins, margins.par) {
    if (margins == "norm") {
        return(cop.pdf(u1 = pnorm(x1),
                       u2 = pnorm(x2),
                       param = param,
                       copula = copula)  * dnorm(x1) * dnorm(x2))
    } else if (margins == "t") {
        return(cop.pdf(u1 = pt(x1, df = margins.par),
                       u2 = pt(x2, df = margins.par),
                       param = param,
                       copula = copula) * dt(x1, df = margins.par) * dt(x2, df = margins.par))
    } else if (margins == "unif") {
        return(cop.pdf(u1 = x1,
                       u2 = x2,
                       param = param,
                       copula = copula))
    } else if (margins == "gamma") {
        return(cop.pdf(u1 = pgamma(x1, shape = margins.par[1], scale = margins.par[2]),
                       u2 = pgamma(x2, shape = margins.par[1], scale = margins.par[2]),
                       param = param,
                       copula = copula) *
                   dgamma(x1, shape = margins.par[1], scale = margins.par[2]) *
                   dgamma(x2, shape = margins.par[1], scale = margins.par[2]))
    } else if (margins == "exp") {
        return(cop.pdf(u1 = pexp(x1, rate = margins.par),
                       u2 = pexp(x2, rate = margins.par),
                       param = param, copula = copula) * dexp(x1, rate = margins.par) * dexp(x2,  rate = margins.par))
    }
}

