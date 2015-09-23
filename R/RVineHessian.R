RVineHessian <- function(data, RVM) {
    
    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36)))) 
        stop("Copula family not implemented.")
    
    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    if (is.null(dim(data))) {
        d <- length(data)
        T <- 1
    } else {
        d <- dim(data)[2]
        T <- dim(data)[1]
    }
    n <- d
    N <- T
    if (n != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    
    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- getFromNamespace("normalizeRVineMatrix", "VineCopula")(RVM)
        data <- data[, o[length(o):1]]
    }
    
    dd <- d * (d - 1)/2
    tt <- sum(RVM$family == 2)
    hessian <- matrix(0, dd + tt, dd + tt)
    subhess <- matrix(0, dd + tt, dd + tt)
    der <- matrix(0, dd + tt, dd + tt)
    subder <- matrix(0, dd + tt, dd + tt)
    
    out <- .C("hesse",
              as.integer(T),
              as.integer(d),
              as.integer(as.vector(RVM$family)),
              as.integer(as.vector(RVM$MaxMat)),
              as.integer(as.vector(RVM$Matrix)),
              as.integer(as.vector(RVM$CondDistr$direct)),
              as.integer(as.vector(RVM$CondDistr$indirect)),
              as.double(as.vector(RVM$par)),
              as.double(as.vector(RVM$par2)),
              as.double(as.vector(data)),
              as.double(as.vector(hessian)),
              as.double(as.vector(subhess)),
              as.double(as.vector(der)),
              as.double(as.vector(subder)),
              PACKAGE = 'VineCopula')
    
    hessian <- matrix(out[[11]], dd + tt, dd + tt)
    subhess <- matrix(out[[12]], dd + tt, dd + tt)
    der <- matrix(out[[13]], dd + tt, dd + tt)
    subder <- matrix(out[[14]], dd + tt, dd + tt)
    
    
    test <- apply(hessian, 2, function(x) max(abs(x)))
    hessian <- hessian[test > 0, test > 0]
    subhess <- subhess[test > 0, test > 0]
    der <- der[test > 0, test > 0]
    subder <- subder[test > 0, test > 0]
    
    
    out <- list(hessian = hessian, der = der)
    
    return(out)
}
