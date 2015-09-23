RVineLogLik <- function(data, RVM, par = RVM$par, par2 = RVM$par2, separate = FALSE, verbose = TRUE) {
    
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
    
    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }
    
    V <- list()
    V$direct <- array(0, dim = c(n, n, N))
    V$indirect <- array(0, dim = c(n, n, N))
    if (is.vector(data)) {
        V$direct[n, , ] <- data[n:1]
    } else {
        V$direct[n, , ] <- t(data[, n:1])
    }
    
    
    if (separate) {
        V$value <- array(0, c(n, n, N))
    } else {
        V$value <- matrix(0, n, n)
    }
    
    ll <- as.vector(V$value)
    vv <- as.vector(V$direct)
    vv2 <- as.vector(V$indirect)
    calcup <- as.vector(matrix(1, dim(RVM), dim(RVM)))
    
    w1 <- as.vector(RVM$family)
    w1[is.na(w1)] <- 0
    th <- as.vector(par)
    th[is.na(th)] <- 0
    th2 <- as.vector(par2)
    th2[is.na(th2)] <- 0
    condirect <- as.vector(as.numeric(RVM$CondDistr$direct))
    conindirect <- as.vector(as.numeric(RVM$CondDistr$indirect))
    maxmat <- as.vector(RVM$MaxMat)
    matri <- as.vector(RVM$Matrix)
    matri[is.na(matri)] <- 0
    maxmat[is.na(maxmat)] <- 0
    condirect[is.na(condirect)] <- 0
    conindirect[is.na(conindirect)] <- 0
    
    if (separate) {
        out <- rep(0, N)
    } else {
        out <- 0
    }
    
    out <- .C("VineLogLikRvine",
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
              as.double(out),
              as.double(ll),
              as.double(vv),
              as.double(vv2),
              as.integer(calcup),
              as.integer(separate),
              PACKAGE = 'VineCopula'
    )
    
    ll <- out[[12]]
    loglik <- out[[11]]
    loglik[loglik %in% c(NA, NaN, -Inf, Inf)] <- -1e+10
    vv <- out[[13]]
    vv2 <- out[[14]]
    V$direct <- array(vv, dim = c(n, n, N))
    V$indirect <- array(vv2, dim = c(n, n, N))
    if (separate) {
        V$value <- array(ll, c(n, n, N))
    } else {
        V$value <- matrix(ll, n, n)
    }
    if (any(V$value %in% c(NA, NaN, -Inf, Inf)) & verbose) {
        print(V$value[V$value %in% c(NA, NaN, -Inf, Inf)])
        print(th)
        print(th2)
    }
    V$value[V$value %in% c(NA, NaN, -Inf, Inf)] <- -1e+10
    # if(any(V$value %in% c(NA,NaN,-Inf,Inf))) print('Fehler in V$value')
    
    return(list(loglik = loglik, V = V))
}


RVinePDF <- function(newdata, RVM) {
    exp(RVineLogLik(newdata, RVM, separate = TRUE)$loglik)
}