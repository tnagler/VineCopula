###########################################
# Algorithm for weighted tau adapted from matlab code by
# http://www.mathworks.com/matlabcentral/fileexchange/27361-weighted-kendall-rank-correlation-matrix/content/kendalltau.m
############################################

#' Matrix of Empirical Kendall's Tau Values
#'
#' This function computes the empirical Kendall's tau using the algorithm by
#' Knight (1966).
#'
#'
#' @param data An N x d data matrix.
#' @param weights Numerical; weights for each observation (optional).
#' @return Matrix of the empirical Kendall's taus.
#' @author Ulf Schepsmeier
#' @seealso [BiCopTau2Par()], [BiCopPar2Tau()],
#' [BiCopEst()]
#' @references Knight, W. R. (1966). A computer method for calculating
#' Kendall's tau with ungrouped data. Journal of the American Statistical
#' Association 61 (314), 436-439.
#' @examples
#'
#' data(daxreturns)
#' Data <- as.matrix(daxreturns)
#'
#' # compute the empirical Kendall's taus
#' TauMatrix(Data)
#'
TauMatrix <- function(data, weights = NA) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_nobs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

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
        if (!is.null(args$na.ind))
            weights <- weights[-args$na.ind]
        weights <- weights / sum(weights)
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
