BiCopHinv <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(c(u1, u2) > 1) || any(c(u1, u2) < 0)) 
        stop("Data has be in the interval [0,1].")
    n <- length(u1)
    
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
    if (any(c(length(family), length(par), length(par2)) == n)) {
        if (length(family) == 1) 
            family <- rep(family, n)
        if (length(par) == 1) 
            par <- rep(par, n)
        if (length(par2) == 1)
            par2 <- rep(par2, n)
    }
    if (!(length(family) %in% c(1, n)))
        stop("'family' has to be a single number or a size n vector")
    if (!(length(par) %in% c(1, n)))
        stop("'par' has to be a single number or a size n vector")
    if (!(length(par2) %in% c(1, n)))
        stop("'par2' has to be a single number or a size n vector")
    
    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }
    
    ## calculate inverse h-functions
    if (length(par) == 1) {
        # call for single parameters
        hinv1 <- .C("Hinv1",                      # h(u2|u1)
                    as.integer(family),
                    as.integer(n), 
                    as.double(u2), 
                    as.double(u1), 
                    as.double(par),
                    as.double(par2), 
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2",                      # h(u1|u2)
                    as.integer(family),
                    as.integer(n), 
                    as.double(u1),
                    as.double(u2), 
                    as.double(par),
                    as.double(par2), 
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        hinv1 <- .C("Hinv1_vec",                      # h(u2|u1) 
                    as.integer(family),
                    as.integer(n), 
                    as.double(u2), 
                    as.double(u1), 
                    as.double(par),
                    as.double(par2), 
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2_vec",                      # h(u1|u2) 
                    as.integer(family),
                    as.integer(n), 
                    as.double(u1),
                    as.double(u2), 
                    as.double(par),
                    as.double(par2), 
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
    }
    
    ## return results
    list(hinv1 = hinv1, hinv2 = hinv2)
}

