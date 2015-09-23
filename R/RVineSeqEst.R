RVineSeqEst <- function(data, RVM, method = "mle", se = FALSE, max.df = 30,
                        max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), 
                        progress = FALSE, weights = NA) {
    data <- as.matrix(data)
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    n <- dim(RVM)
    N <- nrow(data)
    if (dim(data)[2] != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (N < 2) 
        stop("Number of observations has to be at least 2.")
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    if (method != "mle" && method != "itau") 
        stop("Estimation method has to be either 'mle' or 'itau'.")
    if (is.logical(se) == FALSE) 
        stop("'se' has to be a logical variable (TRUE or FALSE).")
    
    if (max.df <= 1) 
        stop("The upper bound for the degrees of freedom parameter has to be larger than 1.")
    if (!is.list(max.BB)) 
        stop("'max.BB' has to be a list.")
    if (max.BB$BB1[1] < 0.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB1[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[2] < 0.001) 
        stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
        stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
    
    o <- diag(RVM$Matrix)
    
    oldRVM <- RVM
    
    if (any(o != length(o):1)) {
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }
    
    Params <- RVM$par
    Params2 <- RVM$par2
    
    if (se == TRUE) {
        seMat1 <- matrix(0, nrow = n, ncol = n)
        seMat2 <- matrix(0, nrow = n, ncol = n)
    }
    
    V <- list()
    V$direct <- array(NA, dim = c(n, n, N))
    V$indirect <- array(NA, dim = c(n, n, N))
    
    V$direct[n, , ] <- t(data[, n:1])
    
    for (i in (n - 1):1) {
        
        for (k in n:(i + 1)) {
            
            m <- RVM$MaxMat[k, i]
            zr1 <- V$direct[k, i, ]
            
            if (m == RVM$Matrix[k, i]) {
                zr2 <- V$direct[k, (n - m + 1), ]
            } else {
                zr2 <- V$indirect[k, (n - m + 1), ]
            }
            
            
            if (RVM$family[k, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20,
                                        27, 28, 29, 30, 37, 38, 39, 40)) {
                if (progress == TRUE) {
                    if (k == n) {
                        message(oldRVM$Matrix[i, i],
                                ",",
                                oldRVM$Matrix[k, i]) 
                    } else {
                        message(oldRVM$Matrix[i, i], 
                                ",", 
                                oldRVM$Matrix[k, i],
                                "|", 
                                paste(oldRVM$Matrix[(k + 1):n, i],
                                      collapse = ","))
                    }
                }
                par.out <- BiCopEst(zr2, 
                                    zr1, 
                                    RVM$family[k, i], 
                                    method,
                                    se,
                                    max.df, 
                                    max.BB,
                                    weights)
                # par1 <- out.par$par
                Params[k, i] <- par.out$par
                Params2[k, i] <- par.out$par2
                if (se == TRUE) {
                    # se1 <- par.out$se
                    seMat1[k, i] <- par.out$se
                    seMat2[k, i] <- par.out$se2
                }
            } else {
                if (progress == TRUE) {
                    if (k == n) {
                        message(oldRVM$Matrix[i, i],
                                ",",
                                oldRVM$Matrix[k, i])
                    } else {
                        message(oldRVM$Matrix[i, i], 
                                ",",
                                oldRVM$Matrix[k, i], 
                                "|", 
                                paste(oldRVM$Matrix[(k + 1):n, i],
                                      collapse = ","))
                    }
                }
                par.out <- BiCopEst(zr2,
                                    zr1,
                                    RVM$family[k, i], 
                                    method, 
                                    se,
                                    max.df,
                                    max.BB,
                                    weights)
                # par1 <- out.par$par
                Params[k, i] <- par.out$par
                Params2[k, i] <- par.out$par2
                if (se == TRUE) {
                    # se1 <- par.out$se
                    seMat1[k, i] <- par.out$se
                    seMat2[k, i] <- par.out$se2
                }
            }
            
            
            if (RVM$CondDistr$direct[k - 1, i]) {
                V$direct[k - 1, i, ] <- .C("Hfunc1",
                                           as.integer(RVM$family[k, i]),
                                           as.integer(length(zr1)), 
                                           as.double(zr1), 
                                           as.double(zr2), 
                                           as.double(Params[k, i]),
                                           as.double(Params2[k, i]), 
                                           as.double(rep(0, length(zr1))), 
                                           PACKAGE = "VineCopula")[[7]]
            }
            if (RVM$CondDistr$indirect[k - 1, i]) {
                V$indirect[k - 1, i, ] <- .C("Hfunc2", 
                                             as.integer(RVM$family[k, i]), 
                                             as.integer(length(zr2)),
                                             as.double(zr2),
                                             as.double(zr1),
                                             as.double(Params[k, i]), 
                                             as.double(Params2[k, i]), 
                                             as.double(rep(0, length(zr1))), 
                                             PACKAGE = "VineCopula")[[7]]
            }
            
        }
    }
    
    ## store results
    oldRVM$par <- Params
    oldRVM$par2 <- Params2
    if (se == FALSE) {
        .out <- list(RVM = oldRVM)
    } else {
        .out <- list(RVM = oldRVM, se = seMat1, se2 = seMat2)
    }
    
    ## free memory and return results
    rm(list = ls())
    .out
}
