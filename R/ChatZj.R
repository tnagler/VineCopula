"ChatZj" <- function(data, u = NULL) {
    # Function that computes the empirical copula cdf 
    # $ \hat C ( u ) = \frac{1}{n+1} \sum_{j=1}^n I( data_j <= u ) , u \in [0,1]^d $
    #----------------------------------------------------------------
    # INPUT:
    #   data        Data to compare u with
    #   u           Data for which to compute empirical copula cdf
    # OUTPUT:
    #   Chat        Empirical copula cdf
    #----------------------------------------------------------------
    # Author: Daniel Berg <daniel at danielberg.no>
    # Date: 26 September 2006
    # Version: 1.0.2
    #----------------------------------------------------------------
    if (missing(u)) 
        u <- data
    n <- dim(u)[1]
    d <- dim(u)[2]
    m <- dim(data)[1]
    Chat <- vector("numeric", n)
    tdata <- t(data)
    tu <- t(u)
    for (j in 1:n) Chat[j] <- sum(colSums(tdata <= tu[, j]) == d)/(m + 1)
    Chat
}


ChatZj2 <- function(data, u = NULL) {
    if (missing(u)) 
        u <- data
    n <- dim(u)[1]
    d <- dim(u)[2]
    m <- dim(data)[1]
    tmp <- .C("ChatZj", 
              as.double(data),
              as.double(u), 
              as.integer(n), 
              as.integer(d), 
              as.integer(m),
              as.double(rep(0, n)),
              PACKAGE = "VineCopula")[[6]]
    return(tmp)
}
