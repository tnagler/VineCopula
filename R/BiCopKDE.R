#' Kernel estimate of  a Bivariate Copula Density
#'
#' A kernel density estimate of the copula density is visualized. The function
#' provides the same options as \code{\link{plot.BiCop}}. Further arguments can
#' be passed to \code{\link[kdecopula::kdecop]{kdecop}} to modify the estimate.
#'
#' @param u1, u2 data vecotrs.
#' @param type plot type; either \code{"contour"} or \code{"surface"} (partial
#' matching is activated) for a contour or perspective/surface plot
#' respectively.
#' @param margins only relevant for types \code{"contour"} and
#' \code{"surface"}; options are: \code{"unif"} for the original copula density,
#' \code{"norm"} for the transformed density with standard normal margins,
#' \code{"exp"} with standard exponential margins, and  \code{"flexp"} with
#' flipped exponential margins. Default is \code{"norm"} for \code{type =
#' "contour"}, and \code{"unif"} for \code{type = "surface"}.
#' \code{"norm"} for the transformed density with standard normal margins
#' (partial matching is activated). Default is \code{"norm"} for \code{type =
#' "contour"}, and \code{"unif"} for \code{type = "surface"}.
#' @param size integer; the plot is based on values on a \eqn{size x size}
#' grid; default is 100 for \code{type = "contour"}, and 25 for \code{type =
#' "surface"}.
#' @param kde.pars list of arguments passed to
#'  \code{\link[kdecopula::kdecop]{kdecop}}.
#' @param \dots optional arguments passed to \code{\link{contour}} or
#' \code{\link{wireframe}}.
#'
#' @details
#' For further details on estimation see \code{\link[kdecopula::kdecop]{kdecop}}.
#'
#' @author Thomas Nagler
#'
#' @examples
#'
#'
BiCopKDE <- function(u1, u2, type = "contour", margins, size,
                     kde.pars = list(), ...) {
    ## prepare the data for usage with plot.kdecopula function
    udata <- list(udata = cbind(u1, u2))
    if (all(colnames(udata$udata) == c("u1", "u2")))
        udata$udata <- unname(udata$udata)

    ## estimate copula density with kde.pars
    # args <- modifyList(udata, kde.pars)
    # est <- do.call(kdecopula::kdecop, args)
    #
    # ## choose margins if missing
    # if (missing(margins)) {
    #     margins <- switch(type,
    #                       "contour" = "norm",
    #                       "surface" = "unif")
    # }
    #
    # # plot
    # return(plot(est,
    #             type = type,
    #             margins = margins,
    #             size = size,
    #             ...))
}
