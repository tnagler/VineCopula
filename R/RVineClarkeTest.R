RVineClarkeTest <- function(data, RVM1, RVM2) {
    
    N <- dim(data)[1]
    
    if (dim(data)[2] < 2) 
        stop("Dimension has to be at least 2.")
    if (N < 2) 
        stop("Number of observations has to be at least 2.")
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    if (is(RVM1)[1] != "RVineMatrix") 
        stop("'RVM1' has to be an RVineMatrix object.")
    if (is(RVM2)[1] != "RVineMatrix") 
        stop("'RVM2' has to be an RVineMatrix object.")
    
    Model1.ll <- RVineLogLik(data, RVM1, separate = TRUE)$loglik
    Model2.ll <- RVineLogLik(data, RVM2, separate = TRUE)$loglik
    
    anz.1 <- sum(RVM1$family >= 1, na.rm = TRUE) + sum(RVM1$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    anz.2 <- sum(RVM2$family >= 1, na.rm = TRUE) + sum(RVM2$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    
    B <- sum(Model1.ll - Model2.ll > 0)
    B.Schwarz <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2) * log(N)/(2 * N) > 0)
    B.Akaike <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2)/N > 0)
    
    if (B == 0 | B == N/2) {
        p <- 1 
    } else { 
        p <- 2 * min(pbinom(B, N, 0.5), 1 - pbinom(B - 1, N, 0.5))
    }
    
    if (B.Schwarz == 0 | B.Schwarz == N/2) {
        p.Schwarz <- 1
    } else {
        p.Schwarz <- 2 * min(pbinom(B.Schwarz, N, 0.5), 1 - pbinom(B.Schwarz - 1, N, 0.5))
    }
    
    if (B.Akaike == 0 | B.Akaike == N/2) {
        p.Akaike <- 1 
    }   else {
        p.Akaike <- 2 * min(pbinom(B.Akaike, N, 0.5), 1 - pbinom(B.Akaike - 1, N, 0.5))
    }
    
    return(list(statistic = B, 
                statistic.Akaike = B.Akaike,
                statistic.Schwarz = B.Schwarz,
                p.value = p, 
                p.value.Akaike = p.Akaike,
                p.value.Schwarz = p.Schwarz))
}