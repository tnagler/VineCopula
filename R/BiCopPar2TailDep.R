BiCopPar2TailDep <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
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
    
    ## adjust length for parameter vectors; stop if not matching
    n <- max(length(family), length(par), length(par2))
    if (length(family) == 1) 
        family <- rep(family, n)
    if (length(par) == 1) 
        par <- rep(par, n)
    if (length(par2) == 1)
        par2 <- rep(par2, n)
    if (!all(c(length(family), length(par), length(par2)) %in% c(1, n)))
        stop("Input lenghts don't match")
    
    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }
    
    ## calculate tail dependence coefficient
    if (length(par) == 1) {
        # call for single parameters
        out <- matrix(calcTD(family, par, par2), ncol = 2)
    } else {
        # vectorized call
        out <- t(vapply(1:length(par),
                        function(i) calcTD(family[i], par[i], par2[i]),
                        numeric(2)))
    }
    
    ## return result
    list(lower = out[, 1], upper = out[, 2])
}

calcTD <- function(family, par, par2) {
    if (family == 0 | family == 1 | family == 5 | family %in% c(23, 24, 26, 27, 28, 29,
                                                                30, 33, 34, 36, 37, 38, 39,
                                                                40, 124, 134, 224, 234)) {
        lower <- 0
        upper <- 0
    } else if (family == 2) {
        lower <- 2 * pt((-sqrt(par2 + 1) * sqrt((1 - par)/(1 + par))), df = par2 + 
                            1)
        upper <- lower
    } else if (family == 3) {
        lower <- 2^(-1/par)
        upper <- 0
    } else if (family == 4 | family == 6) {
        lower <- 0
        upper <- 2 - 2^(1/par)
    } else if (family == 7) {
        lower <- 2^(-1/(par * par2))
        upper <- 2 - 2^(1/par2)
    } else if (family == 8) {
        lower <- 0
        upper <- 2 - 2^(1/(par * par2))
    } else if (family == 9) {
        lower <- 2^(-1/par2)
        upper <- 2 - 2^(1/par)
    } else if (family == 10) {
        lower <- 0
        if (par2 == 1) 
            upper <- 2 - 2^(1/par) else upper <- 0
    } else if (family == 13) {
        lower <- 0
        upper <- 2^(-1/par)
    } else if (family == 14 | family == 16) {
        lower <- 2 - 2^(1/par)
        upper <- 0
    } else if (family == 17) {
        lower <- 2 - 2^(1/par2)
        upper <- 2^(-1/par * par2)
    } else if (family == 18) {
        lower <- 2 - 2^(1/(par * par2))
        upper <- 0
    } else if (family == 19) {
        lower <- 2 - 2^(1/par)
        upper <- 2^(-1/par2)
    } else if (family == 20) {
        if (par2 == 1) 
            lower <- 2 - 2^(1/par) else lower <- 0
            upper <- 0
    } else if (family == 104) {
        par3 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 114) {
        par3 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    } else if (family == 204) {
        par3 <- par2
        par2 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 214) {
        par3 <- par2
        par2 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    }
    
    ## return result
    c(upper, lower)
}