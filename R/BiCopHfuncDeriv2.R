BiCopHfuncDeriv2 <- function(u1, u2, family, par, par2 = 0, deriv = "par", obj = NULL, check.pars = TRUE) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
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
    
    ## check for reasonable input
    if (any(is.na(family)) | any(is.na(par))) 
        stop("Provide either 'family' and 'par' or 'obj'")
    if (!all(family %in% c(0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36))) 
        stop("Copula family not implemented.")
    if (any((family == 2) & (par2 == 0))) 
        stop("For t-copulas, 'par2' must be set.")
    if ((deriv %in% c("par2", "par1par2", "par2u2")) && any(family != 2)) 
        stop("The derivative with respect to the second parameter can only be derived for the t-copula.")
    
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
    
    ## check for family/parameter consistency
    if (check.pars)
        BiCopCheck(family, par, par2)
    
    ## call C routines for specified 'deriv' case 
    if (length(par) == 1) {
        ## call for single parameters
        if (deriv == "par") {
            if (family == 2) {
                out <- .C("diff2hfunc_rho_tCopula_new", 
                          as.double(u1),
                          as.double(u2), 
                          as.integer(n), as.double(c(par, par2)),
                          as.integer(2), 
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diff2hfunc_mod",
                          as.double(u1), 
                          as.double(u2),
                          as.integer(n), 
                          as.double(par),
                          as.integer(family), 
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2") {
            out <- .C("diff2hfunc_nu_tCopula_new", 
                      as.double(u1), 
                      as.double(u2),
                      as.integer(n), 
                      as.double(c(par, par2)), 
                      as.integer(2), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "u2") {
            out <- .C("diff2hfunc_v_mod", 
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(c(par, par2)), 
                      as.integer(family),
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "par1par2") {
            out <- .C("diff2hfunc_rho_nu_tCopula_new",
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(c(par, par2)),
                      as.integer(2), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "par1u2") {
            if (family == 2) {
                out <- .C("diff2hfunc_rho_v_tCopula_new",
                          as.double(u1),
                          as.double(u2), 
                          as.integer(n), 
                          as.double(c(par, par2)),
                          as.integer(2), 
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diff2hfunc_par_v_mod",
                          as.double(u1),
                          as.double(u2), 
                          as.integer(n), 
                          as.double(c(par, par2)),
                          as.integer(family), 
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2u2") {
            out <- .C("diff2hfunc_nu_v_tCopula_new",
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(c(par, par2)), 
                      as.integer(2), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[6]]
        } else {
            stop("This kind of derivative is not implemented")
        }
    } else {
        # vectorized calls
        if (deriv == "par") {
            out <- .C("diff2hfunc_mod_vec",
                      as.double(u1), 
                      as.double(u2),
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family), 
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2") {
            out <- .C("diff2hfunc_nu_tCopula_new_vec", 
                      as.double(u1), 
                      as.double(u2),
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "u2") {
            out <- .C("diff2hfunc_v_mod_vec", 
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par1par2") {
            out <- .C("diff2hfunc_rho_nu_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par1u2") {
            out <- .C("diff2hfunc_par_v_mod_vec",
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family), 
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2u2") {
            out <- .C("diff2hfunc_nu_v_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2), 
                      as.integer(n), 
                      as.double(par),
                      as.double(par2),
                      as.integer(family), 
                      as.double(rep(0, n)), 
                      PACKAGE = "VineCopula")[[7]]
        } else {
            stop("This kind of derivative is not implemented")
        }
    }
    
    ## return results
    out
}
