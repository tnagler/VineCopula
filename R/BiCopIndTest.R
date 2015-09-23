BiCopIndTest <- function(u1, u2) {
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (length(u1) < 2) 
        stop("Number of observations has to be at least 2.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
        stop("Data has be in the interval [0,1].")
    
    # tau = cor(u1,u2,method='kendall')
    tau <- fasttau(u1, u2)
    
    N <- length(u1)
    f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)
    
    return(list(statistic = f, p.value = 2 * (1 - pnorm(f))))
}
