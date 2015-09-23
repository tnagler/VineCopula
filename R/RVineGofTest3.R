RVineGofTest <- function(data, RVM, method = "White", statistic = "CvM", B = 200, alpha = 2) {
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
    if (is(RVM)[1] != "RVineMatrix") 
        stop("'RVM' has to be an RVineMatrix object.")
    
    if (statistic == "Cramer-von Mises") {
        statistic <- "CvM" 
    } else if (statistic == "Kolmogorov-Smirnov") {
        statistic <- "KS"  
    } else if (statistic == "Anderson-Darling") { 
        statistic <- "AD"
    }
    
    if (method == "White") {
        out <- gof_White(data, RVM, B)
    } else if (method == "Breymann" || method == "Berg" || method == "Berg2") {
        out <- gof_PIT(data, RVM, method, B, statistic, alpha)
    } else if (method == "ECP" || method == "ECP2") {
        if (statistic == "AD") 
            stop("The Anderson-Darling statistic is not available for the empirical copula process based goodness-of-fit tests.")
        out <- gof_ECP(data, RVM, B, method, statistic)
    } else if (method == "IR") {
        out2 <- RVineHessian(data, RVM)
        C <- out2$der
        H <- out2$hessian
        p <- dim(C)[1]
        Z <- solve(-C, H)
        IR <- sum(diag(Z))/p
        
        # TODO: p-value via bootstrap (?)
        out <- list(IR = IR, pvalue = NULL)
    }
    
    return(out)
}
