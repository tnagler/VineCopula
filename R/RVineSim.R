RVineSim <- function(N, RVM, U = NULL) {
    
    ## sanity checks
    stopifnot(N >= 1)
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    ## reorder matrix and U (if provided)
    n <- dim(RVM)
    o <- diag(RVM$Matrix)
    RVM <- normalizeRVineMatrix(RVM)
    takeU <- !is.null(U)
    if (takeU) {
        if (!is.matrix(U)) 
            U <- rbind(U, deparse.level = 0L)
        if ((d <- ncol(U)) < 2) 
            stop("U should be at least bivariate")  # should be an (N, n) matrix
        U <- U[, rev(o)]
    }
    
    ## create objects for C-call
    matri <- as.vector(RVM$Matrix)
    w1 <- as.vector(RVM$family)
    th <- as.vector(RVM$par)
    th2 <- as.vector(RVM$par2)
    maxmat <- as.vector(RVM$MaxMat)
    conindirect <- as.vector(RVM$CondDistr$indirect)
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0
    tmp <- rep(0, n * N)
    
    ## simulate R-Vine
    tmp <- .C("SimulateRVine", 
              as.integer(N),
              as.integer(n),
              as.integer(w1), 
              as.integer(maxmat),
              as.integer(matri),
              as.integer(conindirect), 
              as.double(th), 
              as.double(th2),
              as.double(tmp), 
              as.double(U),
              as.integer(takeU), 
              PACKAGE = "VineCopula")[[9]]
    
    ## store results, bring back to initial order and return 
    out <- matrix(tmp, ncol = n)
    if (!is.null(RVM$names)) {
        colnames(out) <- RVM$names
    }
    out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    return(out)
}


transform <- function(M) {
    n <- dim(M)[1]
    
    M.new <- matrix(rep(0, n * n), n, n)
    for (i in 1:n) {
        for (j in 1:i) {
            M.new[(n - i + 1), (n - j + 1)] <- M[i, j]
        }
    }
    
    return(M.new)
}
