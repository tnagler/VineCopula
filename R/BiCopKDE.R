#' Kernel estimate of  a Bivariate Copula Density
#'
#' A kernel density estimate of the copula density is visualized. The function
#' provides the same options as [plot.BiCop()]. Further arguments can
#' be passed to [kdecopula::kdecop()] to modify the estimate. The
#' [kdecopula::kdecopula-package()] must be installed to use
#' this function.
#'
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param type plot type; either `"contour"` or `"surface"` (partial
#' matching is activated) for a contour or perspective/surface plot
#' respectively.
#' @param margins only relevant for types `"contour"` and
#' `"surface"`; options are: `"unif"` for the original copula density,
#' `"norm"` for the transformed density with standard normal margins,
#' `"exp"` with standard exponential margins, and  `"flexp"` with
#' flipped exponential margins. Default is `"norm"` for `type =
#' "contour"`, and `"unif"` for `type = "surface"`.
#' `"norm"` for the transformed density with standard normal margins
#' (partial matching is activated). Default is `"norm"` for `type =
#' "contour"`, and `"unif"` for `type = "surface"`.
#' @param size integer; the plot is based on values on a `size x size`
#' grid; default is 100 for `type = "contour"`, and 25 for `type =
#' "surface"`.
#' @param kde.pars list of arguments passed to
#'  [kdecopula::kdecop()].
#' @param \dots optional arguments passed to [contour()] or
#' [lattice::wireframe()].
#'
#' @details
#' For further details on estimation see [kdecopula::kdecop()].
#'
#' @author Thomas Nagler
#'
#' @examples
#' # simulate data from Joe copula
#' cop <- BiCop(3, tau = 0.3)
#' u <- BiCopSim(1000, cop)
#' contour(cop)  # true contours
#'
#' # kernel contours with standard normal margins
#' BiCopKDE(u[, 1], u[, 2])
#' BiCopKDE(u[, 1], u[, 2], kde.pars = list(mult = 0.5))  # undersmooth
#' BiCopKDE(u[, 1], u[, 2], kde.pars = list(mult = 2))  # oversmooth
#'
#' # kernel density with uniform margins
#' BiCopKDE(u[, 1], u[, 2], type = "surface", zlim = c(0, 4))
#' plot(cop, zlim = c(0, 4))  # true density
#'
#' # kernel contours are also used in pairs.copuladata
#' \donttest{data(daxreturns)
#' data <- as.copuladata(daxreturns)
#' pairs(data[c(4, 5, 14, 15)])}
#'
BiCopKDE <- function(u1, u2, type = "contour", margins, size,
                     kde.pars = list(), ...) {
    if (!requireNamespace("kdecopula", quietly = TRUE))
        stop("The 'kdecopula' package must be installed.")

    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_if_01,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    ## prepare the data for usage with plot.kdecopula function
    args <- list(udata = cbind(u1, u2))
    if (all(colnames(args$udata) == c("u1", "u2")))
        args$udata <- unname(args$udata)

    ## estimate copula density with kde.pars
    args <- modifyList(args, kde.pars)
    est <- do.call(kdecopula::kdecop, args)

    ## choose margins if missing
    if (missing(margins)) {
        margins <- switch(type,
                          "contour" = "norm",
                          "surface" = "unif")
    }

    # plot
    return(plot(est,
                type = type,
                margins = margins,
                size = size,
                ...))
}
