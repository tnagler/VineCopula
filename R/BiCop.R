#' Constructing BiCop-objects
#'
#' This function creates an object of class `BiCop` and checks for
#' family/parameter consistency.
#'
#' @param family An integer defining the bivariate copula family: \cr
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
#' @param par Copula parameter.
#' @param par2 Second parameter for bivariate copulas with two parameters (t,
#' BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default is `par2 = 0`).
#' `par2` should be an positive integer for the Students's t copula
#' `family = 2`.
#' @param tau numeric; value of Kendall's tau; has to lie in the interval
#' (-1, 1). Can only be used with one-parameter families and the t copula.
#' If `tau` is provided, `par` will be ignored.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#'
#' @return An object of class [BiCop()]. It is a list containing
#' information about the bivariate copula. Its components are:
#' \item{family, par, par2}{copula family number and parameter(s),}
#' \item{npars}{number of parameters,}
#' \item{familyname}{name of the copula family,}
#' \item{tau}{Kendall's tau,}
#' \item{beta}{Blomqvist's beta,}
#' \item{taildep}{lower and upper tail dependence coefficients,}
#' \item{call}{the call that created the object.}
#' Objects of this class are also returned by the [BiCopEst()] and
#' [BiCopSelect()] functions. In this case, further information about
#' the fit is added.
#'
#' @note For a comprehensive summary of the model, use `summary(object)`;
#' to see all its contents, use `str(object)`.
#'
#' @author Thomas Nagler
#'
#' @seealso
#' [BiCopPDF()],
#' [BiCopHfunc()],
#' [BiCopSim()],
#' [BiCopEst()],
#' [BiCopSelect()],
#' [plot.BiCop()],
#' [contour.BiCop()]
#'
#' @examples
#'
#' ## create BiCop object for bivariate t-copula
#' obj <- BiCop(family = 2, par = 0.4, par2 = 6)
#' obj
#'
#' ## see the object's content or a summary
#' str(obj)
#' summary(obj)
#'
#' ## a selection of functions that can be used with BiCop objects
#' simdata <- BiCopSim(300, obj)  # simulate data
#' BiCopPDF(0.5, 0.5, obj) # evaluate density in (0.5,0.5)
#' plot(obj)  # surface plot of copula density
#' contour(obj)  # contour plot with standard normal margins
#' print(obj)  # brief overview of BiCop object
#' summary(obj)  # comprehensive overview of BiCop object
#'
BiCop <- function(family, par, par2 = 0, tau = NULL, check.pars = TRUE) {
    ## use tau to construct object (if provided)
    if (!is.null(tau))
        par <- BiCopTau2Par(family, tau)
    if (missing(par) & (family == 0))
        par <- 0
    stopifnot(is.logical(check.pars))
    if (length(c(family, par, par2)) != 3)
        stop("family, par, and par2 have to be a single number.")

    ## family/parameter consistency checks
    if (check.pars) {
        # check for consistency
        BiCopCheck(family, par, par2, call = match.call()[1])
        # warn if par2 is unused
        if ((family %in% allfams[onepar]) && (par2 != 0)) {
            warning("The ",
                    BiCopName(family, short = FALSE),
                    " copula has only one parameter; 'par2' is useless.")
        }
    }

    # calculate dependence measures
    if (family == 0) {
        tau <- 0
        beta <- 0
        taildep <- list()
        taildep$upper <- 0
        taildep$lower <- 0
    } else {
        tau <- BiCopPar2Tau(family, par, par2, check.pars = FALSE)
        taildep <- BiCopPar2TailDep(family, par, par2, check.pars = FALSE)
        beta <- NA  # beta does not work for t copula (cdf disabled)
        if (family != 2)
            beta <- BiCopPar2Beta(family, par, par2, check.pars = FALSE)
    }

    ## get full family name and calculate number of parameters
    familyname <- BiCopName(family, short = FALSE)
    npars <- if (family == 0) 0 else ifelse(family %in% allfams[onepar], 1, 2)

    ## return BiCop object
    out <- list(family     = family,
                par        = par,
                par2       = par2,
                npars      = npars,
                familyname = familyname,
                tau        = tau,
                beta       = beta,
                taildep    = taildep,
                call       = match.call())
    class(out) <- "BiCop"
    out
}

## sets of families
allfams <- c(0:10,
             13, 14, 16:20,
             23, 24, 26:30, 33, 34, 36:40,
             104, 114, 124, 134, 204, 214, 224, 234,
             1004)
tawns <- which(allfams > 100 & allfams < 1000)
onepar <- setdiff(which(allfams %% 10 %in% c(1, 3, 4, 5, 6)), tawns)
twopar <- seq_along(allfams)[-c(1, onepar)]
negfams <- c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40, 124, 134, 224, 234)
posfams <- c(1:10, 13, 14, 16:20, 104, 114, 204, 214)
# families with more dependence near (u, v) = (1, 1)
fams11 <- c(13, 4, 6, 7, 17, 8, 9, 19, 10, 104, 204)
# families with more dependence near (u, v) = (0, 0)
fams00 <- c(3, 14, 16, 7, 17, 18, 9, 19, 20, 114, 214)
# families with more dependence near (u, v) = (1, 0)
fams10 <- c(23, 34, 36, 27, 37, 38, 29, 39, 40, 134, 234)
# families with more dependence near (u, v) = (0, 1)
fams01 <- c(33, 24, 26, 27, 37, 28, 29, 39, 30, 124, 224)

#' @noRd
#' @export
print.BiCop <- function(x, ...) {
    cat("Bivariate copula: ")
    cat(x$familyname, " (par = ", round(x$par, 2), sep = "")
    if (x$family %in% allfams[twopar])
        cat(", par2 = ", round(x$par2, 2), sep = "")
    cat(", tau = ", round(x$tau, 2), sep = "")
    cat(") \n")

    ## return BiCop object invsibly
    invisible(x)
}

#' @noRd
#' @export
summary.BiCop <- function(object, ...) {
    ## print family name
    cat("Family\n")
    cat("------ \n")
    cat("No:   ", object$family)
    cat("\n")
    cat("Name: ", object$familyname)
    cat("\n")
    cat("\n")

    ## print parameters and standard errors
    cat("Parameter(s)\n")
    cat("------------\n")
    cat("par: ", as.character(round(object$par, 2)))
    if (!is.null(object$se))
        cat("  (SE = ", as.character(round(object$se[1], 2)), ")", sep = "")
    cat("\n")
    if (object$family %in% allfams[twopar]) {
        cat("par2:", as.character(round(object$par2, 2)))
        if (!is.null(object$se))
            cat("  (SE = ", as.character(round(object$se2, 2)), ")", sep = "")
    }
    cat("\n")

    ## show dependence measures
    #     object$rho <- BiCopPar2Rho(object)
    cat("Dependence measures\n")
    cat("-------------------\n")
    cat("Kendall's tau:   ", as.character(round(object$tau, 2)))
    if (!is.null(object$emptau)) {
        p <- object$p.value.indeptest
        cat(" (empirical = ",
            as.character(round(object$emptau, 2)),
            ", ",
            "p value ",
            ifelse(p < 0.01,
                   "< 0.01",
                   paste0("= ", as.character(round(p, 2)))),
            ")",
            sep = "")
    }
    cat("\n")
    cat("Upper TD:        ", as.character(round(object$taildep$upper, 2)), "\n")
    cat("Lower TD:        ", as.character(round(object$taildep$lower, 2)), "\n")
#     cat("Blomqvist's beta:", as.character(round(object$beta, 2)), "\n")
    cat("\n")

    ## print fit statistics if available
    if (!is.null(object$nobs)) {
        cat("Fit statistics\n")
        cat("--------------\n")
        cat("logLik: ", as.character(round(object$logLik, 2)), "\n")
        cat("AIC:   ", as.character(round(object$AIC, 2)), "\n")
        cat("BIC:   ", as.character(round(object$BIC, 2)), "\n")
        cat("\n")

    }

    ## return BiCop object invsibly
    invisible(object)
}
