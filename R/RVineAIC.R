RVineAIC <- function(data, RVM, par = RVM$par, par2 = RVM$par2) {
    
    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T
    if (n != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    par[is.na(par)] <- 0
    par[upper.tri(par, diag = T)] <- 0
    par2[is.na(par2)] <- 0
    par2[upper.tri(par2, diag = T)] <- 0
    
    if (any(par != NA) & dim(par)[1] != dim(par)[2]) 
        stop("Parameter matrix has to be quadratic.")
    if (any(par2 != NA) & dim(par2)[1] != dim(par2)[2]) 
        stop("Second parameter matrix has to be quadratic.")
    
    family <- RVM$family
    
    if (!all(par %in% c(0, NA))) {
        for (i in 2:dim(RVM$Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if ((family[i, j] == 1 || family[i, j] == 2) && abs(par[i, j]) >= 1) 
                    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
                if (family[i, j] == 2 && par2[i, j] <= 2) 
                    stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
                if ((family[i, j] == 3 || family[i, j] == 13) && par[i, j] <= 0) 
                    stop("The parameter of the Clayton copula has to be positive.")
                if ((family[i, j] == 4 || family[i, j] == 14) && par[i, j] < 1) 
                    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
                if ((family[i, j] == 6 || family[i, j] == 16) && par[i, j] <= 1) 
                    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
                if (family[i, j] == 5 && par[i, j] == 0) 
                    stop("The parameter of the Frank copula has to be unequal to 0.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par[i, j] <= 0) 
                    stop("The first parameter of the BB1 copula has to be positive.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par2[i,  j] < 1) 
                    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par[i,  j] <= 0) 
                    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par2[i, j] < 1) 
                    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par[i, j] < 1) 
                    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par2[i, j] <= 0) 
                    stop("The second parameter of the BB7 copula has to be positive.")
                if ((family[i, j] == 10 || family[i, j] == 20) && par[i, j] < 1) 
                    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 10 || family[i, j] == 20) && (par2[i, j] <= 0 || par2[i, j] > 1)) 
                    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
                if ((family[i, j] == 23 || family[i, j] == 33) && par[i, j] >= 0) 
                    stop("The parameter of the rotated Clayton copula has to be negative.")
                if ((family[i, j] == 24 || family[i, j] == 34) && par[i, j] > -1) 
                    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 26 || family[i, j] == 36) && par[i, j] >= -1) 
                    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
                if ((family[i, j] == 27 || family[i, j] == 37) && par[i, j] >= 0) 
                    stop("The first parameter of the rotated BB1 copula has to be negative.")
                if ((family[i, j] == 27 || family[i, j] == 37) && par2[i, j] > -1) 
                    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par[i, j] >= 0) 
                    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par2[i, j] > -1) 
                    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par[i, j] > -1) 
                    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par2[i, j] >= 0) 
                    stop("The second parameter of the rotated BB7 copula has to be negative.")
                if ((family[i, j] == 30 || family[i, j] == 40) && par[i, j] > -1) 
                    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 30 || family[i, j] == 40) && (par2[i, j] >= 0 || par2[i, j] < (-1))) 
                    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && par[i, j] < 1) 
                    stop("Please choose 'par' of the Tawn copula in [1,oo).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && (par2[i, j] < 0 ||  par2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && par[i, j] > -1) 
                    stop("Please choose 'par' of the Tawn copula in (-oo,-1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && (par2[i, j] < 0 ||  par2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
            }
        }
    }
    
    npar <- sum(RVM$family >= 1, na.rm = TRUE) + sum(RVM$family %in% c(2,  7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), 
                                                     na.rm = TRUE)
    npar_pair <- (RVM$family >= 1) + (RVM$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234))
    
    RVM2 <- RVM
    RVM2$par <- par
    RVM2$par2 <- par2
    
    like <- RVineLogLik(data, RVM2)
    
    AIC <- -2 * like$loglik + 2 * npar
    pair.AIC <- -2 * like$V$value + 2 * npar_pair
    
    return(list(AIC = AIC, pair.AIC = pair.AIC))
}

RVineBIC <- function(data, RVM, par = RVM$par, par2 = RVM$par2) {
    
    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T
    if (n != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (is(RVM)[1] != "RVineMatrix") 
        stop("'RVM' has to be an RVineMatrix object.")
    
    par[is.na(par)] <- 0
    par[upper.tri(par, diag = T)] <- 0
    par2[is.na(par2)] <- 0
    par2[upper.tri(par2, diag = T)] <- 0
    
    if (any(par != NA) & dim(par)[1] != dim(par)[2]) 
        stop("Parameter matrix has to be quadratic.")
    if (any(par2 != NA) & dim(par2)[1] != dim(par2)[2]) 
        stop("Second parameter matrix has to be quadratic.")
    
    family <- RVM$family
    
    if (!all(par %in% c(0, NA))) {
        for (i in 2:dim(RVM$Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if ((family[i, j] == 1 || family[i, j] == 2) && abs(par[i, j]) >= 1) 
                    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
                if (family[i, j] == 2 && par2[i, j] <= 2) 
                    stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
                if ((family[i, j] == 3 || family[i, j] == 13) && par[i, j] <= 0) 
                    stop("The parameter of the Clayton copula has to be positive.")
                if ((family[i, j] == 4 || family[i, j] == 14) && par[i, j] < 1) 
                    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
                if ((family[i, j] == 6 || family[i, j] == 16) && par[i, j] <= 1) 
                    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
                if (family[i, j] == 5 && par[i, j] == 0) 
                    stop("The parameter of the Frank copula has to be unequal to 0.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par[i, j] <= 0) 
                    stop("The first parameter of the BB1 copula has to be positive.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par2[i, j] < 1) 
                    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par[i, j] <= 0) 
                    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par2[i, j] < 1) 
                    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par[i, j] < 1) 
                    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par2[i, j] <= 0) 
                    stop("The second parameter of the BB7 copula has to be positive.")
                if ((family[i, j] == 10 || family[i, j] == 20) && par[i, j] < 1) 
                    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 10 || family[i, j] == 20) && (par2[i, j] <= 0 || par2[i, j] > 1)) 
                    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
                if ((family[i, j] == 23 || family[i, j] == 33) && par[i, j] >= 0) 
                    stop("The parameter of the rotated Clayton copula has to be negative.")
                if ((family[i, j] == 24 || family[i, j] == 34) && par[i, j] > -1) 
                    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 26 || family[i, j] == 36) && par[i, j] >= -1) 
                    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
                if ((family[i, j] == 27 || family[i, j] == 37) && par[i, j] >= 0) 
                    stop("The first parameter of the rotated BB1 copula has to be negative.")
                if ((family[i, j] == 27 || family[i, j] == 37) && par2[i,j] > -1) 
                    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par[i, j] >= 0) 
                    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par2[i, j] > -1) 
                    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par[i, j] > -1) 
                    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par2[i, j] >= 0) 
                    stop("The second parameter of the rotated BB7 copula has to be negative.")
                if ((family[i, j] == 30 || family[i, j] == 40) && par[i, j] > -1) 
                    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 30 || family[i, j] == 40) && (par2[i, j] >= 0 || par2[i, j] < (-1))) 
                    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i,  j] == 204 || family[i, j] == 214) && par[i, j] < 1) 
                    stop("Please choose 'par' of the Tawn copula in [1,oo).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && (par2[i, j] < 0 || par2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && par[i, j] > -1) 
                    stop("Please choose 'par' of the Tawn copula in (-oo,-1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i,  j] == 224 || family[i, j] == 234) && (par2[i, j] < 0 || par2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
            }
        }
    }
    
    npar <- sum(RVM$family >= 1, na.rm = TRUE) + sum(RVM$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), 
                                                     na.rm = TRUE)
    npar_pair <- (RVM$family >= 1) + (RVM$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234))
    
    RVM2 <- RVM
    RVM2$par <- par
    RVM2$par2 <- par2
    
    like <- RVineLogLik(data, RVM2)
    
    BIC <- -2 * like$loglik + log(T) * npar
    pair.BIC <- -2 * like$V$value + log(T) * npar_pair
    
    return(list(BIC = BIC, pair.BIC = pair.BIC))
}
