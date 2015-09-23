#################################################################
#    							#
# RVineGrad							#
#								#
# Function to calculate the derivative of one			#
# pair-copula in an R-vine					#
#								#
# Input:							#
# data		data set					#
# RVM		R-Vine matrix object				#
# calcupdate	array of Update-Matrices (output of RVineMatrixUpdate)	#
# par, par2	Copula parameter stored in an RVM-Matrix	#
# start.V	log-liklihoods (output of RVineLogLik)		#
#								#
# Output:							#
# gradient	gradient of the R-vine				#
#################################################################

RVineGrad <- function(data, RVM, par = RVM$par, par2 = RVM$par2, start.V = NA, posParams = (RVM$family > 0)) {
    
    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36)))) 
        stop("Copula family not implemented.")
    
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
    
    if (any(is.na(start.V))) {
        loglik <- RVineLogLik(data, RVM, par = par, par2 = par2, separate = TRUE)
        V <- loglik$V
        
    } else {
        V <- start.V
        V$value[V$value %in% c(NA, NaN, -Inf)] <- -1e+10
        if (any(is.na(V$value))) 
            message("NA in LogL call")
    }
    
    
    ll <- as.vector(V$value)
    vv <- as.vector(V$direct)
    vv2 <- as.vector(V$indirect)
    
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
    
    
    out <- rep(0, sum(posParams[lower.tri(posParams, diag = FALSE)]) + sum(w1 == 2))
    
    out <- .C("VineLogLikRvineGradient",
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
              as.integer(as.vector(posParams)),
              PACKAGE = 'VineCopula')
    
    
    
    gradient2 <- out[[11]]
    gradient2[gradient2 %in% c(NA, NaN, -Inf)] <- -1e+10
    
    dd <- sum(RVM$family > 0)
    tt <- sum(w1 == 2)
    grad1 <- gradient2[1:dd]
    gradient <- grad1[dd:1]
    if (tt > 0) {
        grad2 <- gradient2[(dd + 1):(dd + tt)]
        gradient <- c(gradient, grad2[tt:1])
    }
    
    
    out2 <- list(gradient = gradient)
    return(out2)
}
