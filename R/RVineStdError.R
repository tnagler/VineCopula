RVineStdError <- function(hessian, RVM) {
    # Test auf pos. semidef.
    se3 <- numeric()
    a <- eigen(-hessian, only.values = TRUE)
    if (any(a$values < 0)) {
        warning("The negative Hessian matrix is not positive definite. Thus NAs will be returned in some entries.")
    }
    se <- sqrt((diag(solve(-hessian))))
    
    d <- dim(RVM$family)[1]
    posParams <- (RVM$family > 0)
    posParams2 <- (RVM$family == 2)
    
    posParams[is.na(posParams)] <- FALSE
    posParams2[is.na(posParams2)] <- FALSE
    
    nParams <- sum(posParams, na.rm = TRUE)
    nParams2 <- sum(posParams2, na.rm = TRUE)
    
    d2 <- dim(hessian)[1]
    
    if (nParams < d2) {
        # t-copula involved
        se2 <- se[(nParams + 1):d2]
        se <- se[1:nParams]
    }
    
    SE <- matrix(0, d, d)
    t <- 1
    for (i in (d - 1):1) {
        for (j in d:(i + 1)) {
            if (posParams[j, i]) {
                SE[j, i] <- se[t]
                t <- t + 1
            }
        }
    }
    
    t <- 1
    if (nParams < d2) {
        SE2 <- matrix(0, d, d)
        for (i in (d - 1):1) {
            for (j in d:(i + 1)) {
                if (RVM$family[j, i] == 2) {
                    SE2[j, i] <- se2[t]
                    t <- t + 1
                }
            }
        }
    }
    
    out <- list()
    out$se <- SE
    if (nParams < d2) out$se2 <- SE2
    
    return(out)
}
