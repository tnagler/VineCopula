#' Constructing BiCop-objects
#'
#' This function creates an object of class \code{BiCop} and checks for
#' family/parameter consistency.
#'
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
#' \code{114} = rotated Tawn type 1 copula
#' (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par Copula parameter.
#' @param par2 Second parameter for bivariate copulas with two parameters (t,
#' BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default is \code{par2 = 0}).
#' \code{par2} should be an positive integer for the Students's t copula
#' \code{family = 2}.
#' @param tau numeric; value of Kendall's tau; has to lie in the interval
#' (-1, 1). If \code{tau} is provided, \code{par} will be ignored.
#'
#' @return An object of class \code{\link{BiCop}}. Objects of this class are
#' also returned by the \code{\link{BiCopEst}} and \code{\link{BiCopSelect}}
#' functions.
#'
#' @author Thomas Nagler
#'
#' @seealso
#' \code{\link{BiCopPDF}},
#' \code{\link{BiCopHfunc}},
#' \code{\link{BiCopSim}},
#' \code{\link{BiCopEst}},
#' \code{\link{BiCopSelect}},
#' \code{\link{plot.BiCop}}
#'
#' @examples
#'
#' ## create BiCop object for bivariate t-copula
#' obj <- BiCop(family = 2, par = 0.4, par2 = 6)
#'
#' ## a selection of function that can be used with BiCop objects
#' simdata <- BiCopSim(300, obj)  # simulate data
#' BiCopPDF(0.5, 0.5, obj) # evaluate density in (0.5,0.5)
#' plot(obj)  # normal contour plot
#'
BiCop <- function(family, par, par2 = 0, tau = NULL) {
    ## use tau to construct object (if provided)
    if (!is.null(tau))
        par <- BiCopTau2Par(family, tau)

    ## family/parameter consistency checks
    BiCopCheck(family, par, par2)

    # calculate dependence measures
    tau <- BiCopPar2Tau(family, par, par2, check.pars = FALSE)
    taildep <- BiCopPar2TailDep(family, par, par2, check.pars = FALSE)
    beta <- BiCopPar2Beta(family, par, par2, check.pars = FALSE)

    ## get full family name
    familyname <- BiCopName(family, short = FALSE)

    ## return BiCop object
    out <- list(family     = family,
                par        = par,
                par2       = par2,
                familyname = familyname,
                tau     = tau,
                beta    = beta,
                taildep = taildep)
    class(out) <- "BiCop"
    out
}

print.BiCop <- function(x, ...) {
    cat("Bivariate copula: \n")
    cat(x$familyname, " (par = ", round(x$par, 2), sep = "")
    if (x$par2 != 0)
        cat(", par2 = ", round(x$par2, 2), sep = "")
    cat(") \n")

    ## return BiCop object invsibly
    invisible(x)
}

summary.BiCop <- function(object, ...) {
    ## print object
    print.BiCop(object)
    cat("\n")

    ## show dependence measures as table
    ms <- c(tau = object$tau,
            bet = object$beta,
            utd = object$taildep$upper,
            ltd = object$taildep$lower)
    names(ms) <- c("Kendall's tau",
                   "Blomqvist's beta",
                   "Upper TD",
                   "Lower TD")
    cat("Dependence measures: \n")
    print(ms, digits = 2)
    cat("\n")

    ## return BiCop object invsibly
    invisible(object)
}
