BiCopSim <- function(N, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
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
    if (any(c(length(family), length(par), length(par2)) == N)) {
        if (length(family) == 1) 
            family <- rep(family, N)
        if (length(par) == 1) 
            par <- rep(par, N)
        if (length(par2) == 1)
            par2 <- rep(par2, N)
    }
    if (!(length(family) %in% c(1, N)))
        stop("'family' has to be a single number or a size N vector")
    if (!(length(par) %in% c(1, N)))
        stop("'par' has to be a single number or a size N vector")
    if (!(length(par2) %in% c(1, N)))
        stop("'par2' has to be a single number or a size N vector")
    
    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }
    
    
    ## start with independent uniforms (byrow for backwards compatibility)
    w <- matrix(runif(2*N), ncol = 2, byrow = TRUE)
    
    ## simulate from copula by inverse rosenblatt transform
    if (length(par) == 1) {
        # call for single parameters
        tmp <- .C("Hinv1", 
                  as.integer(family),
                  as.integer(N), 
                  as.double(w[, 2]), 
                  as.double(w[, 1]), 
                  as.double(par),
                  as.double(par2), 
                  as.double(rep(0, N)),
                  PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        tmp <- .C("Hinv1_vec", 
                  as.integer(family),
                  as.integer(N), 
                  as.double(w[, 2]), 
                  as.double(w[, 1]), 
                  as.double(par),
                  as.double(par2), 
                  as.double(rep(0, N)),
                  PACKAGE = "VineCopula")[[7]]
    }
    
    ## return results
    U <- matrix(c(w[, 1], tmp), ncol = 2)
    U
}
