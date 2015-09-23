gof_PIT <- function(data, RVM, weight = "Breymann", B = 200, statisticName = "AD", alpha = 2) {
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
    
    if (weight == "Breymann") {
        method <- 1 
    } else if (weight == "Berg") {
        method <- 2 
    } else if (weight == "Berg2") {
        method <- 3
    }
    
    if (statisticName == "Cramer-von Mises" || statisticName == "CvM") {
        statisticName <- 3 
    } else if (statisticName == "Kolmogorov-Smirnov" || statisticName == "KS") {
        statisticName <- 2 
    } else if (statisticName == "Anderson-Darling" || statisticName == "AD") {
        statisticName <- 1
    }
    
    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }
    
    N <- T
    n <- d
    V <- list()
    V$direct <- array(0, dim = c(n, n, N))
    V$indirect <- array(0, dim = c(n, n, N))
    if (is.vector(data)) {
        V$direct[n, , ] <- data[n:1]
    } else {
        V$direct[n, , ] <- t(data[, n:1])
    }
    
    vv <- as.vector(V$direct)
    vv2 <- as.vector(V$indirect)
    calcup <- as.vector(matrix(1, dim(RVM), dim(RVM)))
    
    w1 <- as.vector(RVM$family)
    w1[is.na(w1)] <- 0
    th <- as.vector(RVM$par)
    th[is.na(th)] <- 0
    th2 <- as.vector(RVM$par2)
    th2[is.na(th2)] <- 0
    condirect <- as.vector(as.numeric(RVM$CondDistr$direct))
    conindirect <- as.vector(as.numeric(RVM$CondDistr$indirect))
    maxmat <- as.vector(RVM$MaxMat)
    matri <- as.vector(RVM$Matrix)
    matri[is.na(matri)] <- 0
    maxmat[is.na(maxmat)] <- 0
    condirect[is.na(condirect)] <- 0
    conindirect[is.na(conindirect)] <- 0
    
    # Anderson-Darling test statistic (mittlerweile auch KS und CvM)
    tmp <- .C("gofPIT_AD",
              as.integer(T),
              as.integer(d),
              as.integer(w1), 
              as.integer(maxmat),
              as.integer(matri),
              as.integer(condirect),
              as.integer(conindirect), 
              as.double(th),
              as.double(th2), 
              as.double(data), 
              as.double(0), 
              as.double(vv), 
              as.double(vv2), 
              as.integer(calcup),
              as.integer(method),
              as.integer(alpha), 
              as.integer(B),
              as.integer(statisticName), 
              PACKAGE = "VineCopula")[[11]]
    
    statistic <- tmp
    
    if (B == 0) {
        dat <- RVinePIT(data, RVM)
        
        tmp <- matrix(0, T, d)
        for (i in 1:d) {
            if (method == 1) {
                tmp[, i] <- qnorm(dat[, i])^2
            } else if (method == 2) {
                tmp[, i] <- abs(dat[, i] - 0.5) 
            } else if (method == 3) {
                tmp[, i] <- (dat[, i] - 0.5)^alpha
            }
        }
        
        S <- rep(0, T)
        for (i in 1:T) {
            S[i] <- sum(tmp[i, ])
        }
        if (statisticName == 1 && method == 1) {
            # Macht nur Sinn bei Breymann
            if (requireNamespace("ADGofTest", quietly = TRUE)) {
                pvalue <- ADGofTest::ad.test(S, pchisq, df = d)$p.value
            } else {
                stop("For Anderson-Darling 'statistic' the package 'ADGofTest' has to be installed")
            }
        } else if (statisticName == 2 && method == 1) {
            pvalue <- ks.test(S, "pchisq", df = d)$p.value 
        } else { 
            # stop('This statistic is not implemented or for this statistic no
            # analytic form of the test statistic and p-value is available')
            pvalue <- NULL
        }
    } else {
        # p-value according to the statistic
        tmp <- .C("gofPIT_AD_pvalue",
                  as.integer(T),
                  as.integer(d), 
                  as.integer(w1), 
                  as.integer(maxmat),
                  as.integer(matri), 
                  as.integer(condirect), 
                  as.integer(conindirect),
                  as.double(th), 
                  as.double(th2), 
                  as.double(data), 
                  as.double(statistic), 
                  as.double(vv), 
                  as.double(vv2), 
                  as.integer(calcup), 
                  as.integer(method), 
                  as.integer(alpha), 
                  as.integer(B), 
                  as.double(0), 
                  as.integer(statisticName), 
                  PACKAGE = "VineCopula")[[18]]
        
        pvalue <- tmp
    }
    
    if (statisticName == 1) {
        out <- list(AD = statistic, p.value = pvalue) 
    } else if (statisticName == 2) {
        out <- list(KS = statistic, p.value = pvalue) 
    } else {
        out <- list(CvM = statistic, p.value = pvalue) 
    }
    
    return(out)
}
