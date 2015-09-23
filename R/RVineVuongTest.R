RVineVuongTest <- function(data, RVM1, RVM2) {
    
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
    
    anz.1 <- sum(RVM1$family >= 1, na.rm = TRUE) + sum(RVM1$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    anz.2 <- sum(RVM2$family >= 1, na.rm = TRUE) + sum(RVM2$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    
    if (all(Model1.ll - Model2.ll == 0)) {
        # models are the same
        V <- 0
        V.Schwarz <- 0
        V.Akaike <- 0
        
        p <- 1
        p.Schwarz <- 1
        p.Akaike <- 1
    } else {
        w <- 1/N * sum((Model1.ll - Model2.ll)^2) + (1/N * sum(Model1.ll - Model2.ll))^2
        w <- sqrt(w)
        
        LR <- sum(Model1.ll) - sum(Model2.ll)
        LR.Schwarz <- LR - ((anz.1/2 * log(N) - anz.2/2 * log(N)))
        LR.Akaike <- LR - (anz.1 - anz.2)
        
        V <- LR/(sqrt(N) * w)
        V.Schwarz <- LR.Schwarz/(sqrt(N) * w)
        V.Akaike <- LR.Akaike/(sqrt(N) * w)
        
        p <- 2 * min(pnorm(V), 1 - pnorm(V))
        p.Schwarz <- 2 * min(pnorm(V.Schwarz), 1 - pnorm(V.Schwarz))
        p.Akaike <- 2 * min(pnorm(V.Akaike), 1 - pnorm(V.Akaike))
    }
    
    return(list(statistic = V, 
                statistic.Akaike = V.Akaike,
                statistic.Schwarz = V.Schwarz, 
                p.value = p,
                p.value.Akaike = p.Akaike,
                p.value.Schwarz = p.Schwarz))
}