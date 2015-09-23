BetaMatrix <- function(data) {
    d <- dim(data)[2]
    
    betahat <- matrix(1, d, d)
    for (i in 1:(d - 1)) {
        u1 <- data[, i]
        for (j in (i + 1):d) {
            u2 <- data[, j]
            betahat[i, j] <- betaFunc(u1, u2, 1/2, 1/2)
            betahat[j, i] <- betahat[i, j]
        }
    }
    
    return(betahat)
}


# empirical copula
empcop <- function(u1, u2, u, v) {
    n <- length(u1)
    a <- which(u1 < u)
    b <- which(u2 < v)
    sc <- intersect(a, b)
    return(1/n * length(sc))
}

# survival copula
survivalcop <- function(u1, u2, u, v) {
    n <- length(u1)
    a <- which(u1 > u)
    b <- which(u2 > v)
    sc <- intersect(a, b)
    return(1/n * length(sc))
}

# h_d
h <- function(u, v) (min(u, v) + min(1 - u) - u * v - (1 - u) * (1 - v))^-1

# g_d
g <- function(u, v) (u * v) + (1 - u) * (1 - v)

# beta
betaFunc <- function(u1, u2, u, v) {
    h(u, v) * (empcop(u1, u2, u, v) + survivalcop(u1,  u2, u, v) - g(u, v))
}
