BiCopTau2Par <- function(family, tau, check.taus = TRUE) {
    ## sanity check
    if (any(abs(tau) > 0.99999))
        stop("some tau is too close to -1 or 1")
    
    ## adjust length for input vectors; stop if not matching
    n <- max(length(family), length(tau))
    if (length(family) == 1) 
        family <- rep(family, n)
    if (length(tau) == 1) 
        par <- rep(tau, n)
    if (!all(c(length(family), length(tau)) %in% c(1, n)))
        stop("Input lenghts don't match")
    
    ## check for family/tau consistency
    if (check.taus)
        BiCopCheckTaus(family, tau)
    
    
    ## calculate the parameter
    if (length(tau) == 1) {
        # call for single parameters
        out <- calcPar(family, tau)
    } else {
        # vectorized call
        out <- vapply(1:length(tau),
                      function(i) calcPar(family[i], tau[i]),
                      numeric(1))
    }
    
    ## return result
    out
}

calcPar <- function(family, tau) {
    ## calculation of parameter(s) depending on pair-copula family
    if (family == 0) {
        par <- rep(0, times = length(tau))
    } else if (family %in% 1:2) {
        par <- sin(pi * tau/2)
    } else if (family %in% c(3, 13)) {
        par <- 2 * tau/(1 - tau)
    } else if (family %in% c(4, 14)) {
        par <- 1/(1 - tau)
    } else if (family == 5) {
        par <- if (tau == 0) 0 else Frank.itau.JJ(tau)
    } else if (family %in% c(6, 16)) {
        par <- Joe.itau.JJ(tau)
    } else if (family %in% c(23, 33)) {
        par <- 2 * tau/(1 + tau)
    } else if (family %in% c(24, 34)) {
        par <- -(1/(1 + tau))
    } else if (family %in% c(26, 36)) {
        par <- -Joe.itau.JJ(-tau)
    } else if (family %in% c(41, 51)) {
        par <- ipsA.tau2cpar(tau)
    } else if (family %in% c(61, 71)) {
        par <- -ipsA.tau2cpar(-tau)
    }
    
    ## return result
    par
}

