RVinePIT <- function(data, RVM) {
    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 
                                104, 114, 124, 134, 204, 214, 224, 234)))) 
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
    
    tmp <- .C("RvinePIT",
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
              as.double(rep(0,T*d)),
              as.double(vv),
              as.double(vv2),
              as.integer(calcup),
              PACKAGE = 'VineCopula')[[11]]
    U <- matrix(tmp, ncol = d)
    U <- U[, sort(o[length(o):1], index.return = TRUE)$ix]
    
    return(U)
}
