BiCopPar2Beta <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
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
    ## check for reasonable input
    if (any(is.na(family)) | any(is.na(par)))
        stop("Provide either 'family' and 'par' or 'obj'")
    n <- max(length(family), length(par), length(par2))
    if (!all(c(length(family), length(par), length(par2)) %in% c(1, n)))
        stop("Input lengths don't match")
    
    ## calculate beta
    Cuv <- BiCopCDF(rep(0.5, n),
                    rep(0.5, n),
                    family,
                    par,
                    par2,
                    check.pars = check.pars)
    4 * Cuv - 1
}