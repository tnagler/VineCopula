gof_White <- function(data, RVM, B = 200) {
    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36)))) 
        stop("Copula family not implemented.")
    
    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    T <- dim(data)[1]
    d <- dim(data)[2]
    
    if (d != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (!("RVineMatrix" %in% is(RVM))) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    dd <- sum(RVM$family != 0)
    tt <- sum(RVM$family == 2)
    D <- rep(0, (dd + tt) * (dd + tt + 1)/2)
    V0 <- matrix(0, (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
    
    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        # RVM = normalizeRVineMatrix(RVM)
        RVM <- getFromNamespace("normalizeRVineMatrix", "VineCopula")(RVM)
        data <- data[, o[length(o):1]]
    }
    
    out <- .C("White", 
              as.integer(T), 
              as.integer(d), 
              as.integer(as.vector(RVM$family)), 
              as.integer(as.vector(RVM$MaxMat)),
              as.integer(as.vector(RVM$Matrix)), 
              as.integer(as.vector(RVM$CondDistr$direct)),
              as.integer(as.vector(RVM$CondDistr$indirect)), 
              as.double(as.vector(RVM$par)), 
              as.double(as.vector(RVM$par2)), 
              as.double(as.vector(data)), 
              as.double(as.vector(D)),
              as.double(as.vector(V0)), 
              PACKAGE = "VineCopula")
    
    D <- out[[11]]
    # print(round(D,1))
    V0 <- matrix(out[[12]], (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
    handle2 <- try(solve(V0, D), TRUE)
    if (is.null(dim(handle2))) 
        handle2 <- ginv(V0) %*% D
    test <- T * t(D) %*% handle2
    # test=T*t(D)%*%solve(V0,D)
    test <- as.numeric(test)
    
    if (B == 0) {
        pvalue <- 1 - pchisq(test, df = length(D))
    } else {
        # bootstrap
        pvalue <- 0
        for (i in 1:B) {
            D <- rep(0, (dd + tt) * (dd + tt + 1)/2)
            V0 <- matrix(0, (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
            f <- sample(T)
            out <- .C("White", 
                      as.integer(T),
                      as.integer(d),
                      as.integer(as.vector(RVM$family)), 
                      as.integer(as.vector(RVM$MaxMat)),
                      as.integer(as.vector(RVM$Matrix)), 
                      as.integer(as.vector(RVM$CondDistr$direct)),
                      as.integer(as.vector(RVM$CondDistr$indirect)), 
                      as.double(as.vector(RVM$par)), 
                      as.double(as.vector(RVM$par2)), 
                      as.double(as.vector(data[f, ])),
                      as.double(as.vector(D)), 
                      as.double(as.vector(V0)), 
                      PACKAGE = "VineCopula")
            
            D <- out[[11]]
            V0 <- matrix(out[[12]], (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
            handle2 <- try(solve(V0, D), TRUE)
            if (is.null(dim(handle2))) 
                handle2 <- ginv(V0) %*% D
            test_b <- T * t(D) %*% handle2
            # test_b=T*t(D)%*%solve(V0,D)
            
            if (test_b >= test) 
                pvalue <- pvalue + 1/B
        }
    }
    
    out <- list(White = test, p.value = pvalue)
    return(out)
}


gof_White2 <- function(data, RVM, B = 200) {
    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36)))) 
        stop("Copula family not implemented.")
    
    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    T <- dim(data)[1]
    d <- dim(data)[2]
    
    if (d != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (!("RVineMatrix" %in% is(RVM))) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    dd <- sum(RVM$family != 0)
    tt <- sum(RVM$family == 2)
    D <- rep(0, (dd + tt) * (dd + tt + 1)/2)
    V0 <- matrix(0, (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
    
    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        # RVM = normalizeRVineMatrix(RVM)
        RVM <- getFromNamespace("normalizeRVineMatrix", "VineCopula")(RVM)
        data <- data[, o[length(o):1]]
    }
    
    dd <- sum(RVM$family != 0)
    tt <- sum(RVM$family == 2)
    Dprime <- rep(0, (dd + tt) * (dd + tt + 1)/2)
    Vt <- matrix(0, (dd + tt) * (dd + tt + 1)/2, (dd + tt) * (dd + tt + 1)/2)
    Dprime2 <- matrix(0, T, (dd + tt) * (dd + tt + 1)/2)  # => d_t
    gradD <- matrix(0, (dd + tt) * (dd + tt + 1)/2, (dd + tt))
    
    for (t2 in 1:T) {
        out <- RVineHessian(data[t2, ], RVM)
        handle <- out$hessian
        Hprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
        handle <- out$der
        Cprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
        Dprime2[t2, ] <- Hprime + Cprime
        Dprime <- Dprime + Dprime2[t2, ]
        
        # Finite Differenzen f?r gradD
        for (i in 1:(dd + tt)) {
            # laufe ?ber alle Parameter
            if (i <= dd) {
                par_plus <- RVM$par
                par_plus[i] <- par_plus[i] + 1e-05  # + epsilon
                par2_plus <- RVM$par2
                par_minus <- RVM$par
                par_minus[i] <- par_minus[i] - 1e-05  # - epsilon
                par2_minus <- RVM$par2
            } else {
                par_plus <- RVM$par
                par2_plus <- RVM$par2
                par2_plus[i] <- par2_plus[i] + 1e-05  # + epsilon
                par_minus <- RVM$par
                par2_minus <- RVM$par2
                par2_minus[i] <- par2_minus[i] - 1e-05  # - epsilon
            }
            RVM_plus <- RVineMatrix(Matrix = RVM$Matrix,
                                    family = RVM$family, 
                                    par = matrix(par_plus, d, d), 
                                    par2 = matrix(par2_plus, d, d))
            RVM_minus <- RVineMatrix(Matrix = RVM$Matrix, 
                                     family = RVM$family, 
                                     par = matrix(par_minus, d, d),
                                     par2 = matrix(par2_minus, d, d))
            
            out <- RVineHessian(data[t2, ], RVM_plus)
            handle <- out$hessian
            Hprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
            handle <- out$der
            Cprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
            Dplus <- Hprime + Cprime
            out <- RVineHessian(data[t2, ], RVM_minus)
            handle <- out$hessian
            Hprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
            handle <- out$der
            Cprime <- as.vector(handle[lower.tri(handle, diag = TRUE)])
            Dminus <- Hprime + Cprime
            
            gradD[, i] <- gradD[, i] + (Dplus - Dminus)/2/1e-05  # finite Differenzen
        }
        
    }
    
    D <- Dprime/T  # => d
    gradD <- gradD/T  # => Erwartungswert der Ableitung von d
    H <- RVineHessian(data, RVM)$hessian  # H bzw. -B
    
    # Berechne V_0
    for (t2 in 1:T) {
        handle1 <- solve(H, RVineGrad(data[t2, ], RVM)$gradient)
        handle <- Dprime2[t2, ] - gradD %*% handle1
        Vt <- Vt + (handle %*% t(handle))
    }
    
    V0 <- Vt/T
    test <- T * t(D) %*% solve(V0, D)
    
    pvalue <- 1 - pchisq(test, df = length(D))
    
    out <- list(White = test, p.value = pvalue)
    return(out)
}