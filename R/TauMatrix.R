###########################################
# TauMatrix
#
# Input:
# data    	data matrix
# weights   vector of weights
#
# Output:
# ktauMarix	Kendall's tau matrix
#
#
# Algorithm for weighted tau adapted from matlab code by http://www.mathworks.com/matlabcentral/fileexchange/27361-weighted-kendall-rank-correlation-matrix/content/kendalltau.m
############################################

TauMatrix <- function(data, weights = NA) {
    data <- as.matrix(data)
    if (any(is.na(weights))) {
        # if (any(data > 1) || any(data < 0)) 
        #     stop("Data has be in the interval [0,1].")
        d <- dim(data)[2]
        N <- dim(data)[1]
        
        ktau <- rep(0, d * (d - 1)/2)
        
        out <- .C("ktau_matrix", 
                  as.double(data), 
                  as.integer(d), 
                  as.integer(N), 
                  as.double(ktau), 
                  PACKAGE = "VineCopula")
        
        ktau <- out[[4]]
        
        ktauMatrix <- matrix(1, d, d)
        k <- 1
        for (i in 1:(d - 1)) {
            for (j in (i + 1):d) {
                ktauMatrix[i, j] <- ktau[k]
                ktauMatrix[j, i] <- ktau[k]
                k <- k + 1
            }
        }
        if (!is.null(colnames(data))) {
            rownames(ktauMatrix) <- colnames(ktauMatrix) <- colnames(data)
        }
    } else {
        T <- dim(data)[1]
        A <- data
        out <- combn(1:T, 2)
        i1 <- out[1, ]
        i2 <- out[2, ]
        w <- as.numeric(weights)/sqrt(sum(as.numeric(weights[i1] * weights[i2])))
        tau <- sign(A[i1, ] - A[i2, ])
        tau <- t(tau) %*% (tau * w[i1] * w[i2])
        temp <- diag(tau)
        ktauMatrix <- tau/sqrt(temp %*% t(temp))
    }
    return(ktauMatrix)
}
