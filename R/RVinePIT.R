#' Probability Integral Transformation for R-Vine Copula Models
#'
#' This function applies the probability integral transformation (PIT) for
#' R-vine copula models to given copula data.
#'
#' The multivariate probability integral transformation (PIT) of Rosenblatt
#' (1952) transforms the copula data \eqn{u = (u_1,\ldots,u_d)} with a given
#' multivariate copula C into independent data in \eqn{[0,1]^d}, where d is the
#' dimension of the data set. \cr
#'
#' Let \eqn{u = (u_1,\ldots,u_d)} denote copula data of dimension d. Further
#' let C be the joint cdf of \eqn{u = (u_1,\ldots,u_d)}. Then Rosenblatt's
#' transformation of u, denoted as \eqn{y = (y_1,\ldots,y_d)}, is defined as
#' \deqn{ y_1 := u_1,\ \ y_2 := C(u_2|u_1), \ldots\ y_d :=
#' C(u_d|u_1,\ldots,u_{d-1}), } where \eqn{C(u_k|u_1,\ldots,u_{k-1})} is the
#' conditional copula of \eqn{U_k} given \eqn{U_1 = u_1,\ldots, U_{k-1} =
#' u_{k-1}, k = 2,\ldots,d}. The data vector \eqn{y = (y_1,\ldots,y_d)} is now
#' i.i.d. with \eqn{y_i \sim U[0, 1]}. The algorithm for the R-vine PIT is
#' given in the appendix of Schepsmeier (2015).
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM \code{\link{RVineMatrix}} objects of the R-vine model.
#' @return An \code{N} x d matrix of PIT data from the given R-vine copula
#' model.
#'
#' @author Ulf Schepsmeier
#'
#' @seealso \code{\link{RVineGofTest}}
#'
#' @references Rosenblatt, M. (1952).  Remarks on a Multivariate
#' Transformation. The Annals of Mathematical Statistics 23 (3), 470-472.
#'
#' Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for
#' vine copula models with fixed margins. Journal of Multivariate Analysis 138,
#' 34-52.
#'
#' @examples
#' # load data set
#' data(daxreturns)
#'
#' # select the R-vine structure, families and parameters
#' RVM <- RVineStructureSelect(daxreturns[,1:3], c(1:6))
#'
#' # PIT data
#' pit <- RVinePIT(daxreturns[,1:3], RVM)
#'
#' par(mfrow = c(1,2))
#' plot(daxreturns[,1], daxreturns[,2])	# correlated data
#' plot(pit[,1], pit[,2])	# i.i.d. data
#'
#' cor(pit, method = "kendall")
#'
RVinePIT <- function(data, RVM) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    fix_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs)
    list2env(args, environment())

    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36,
                                104, 114, 124, 134, 204, 214, 224, 234))))
        stop("Copula family not implemented.")

    T <- dim(data)[1]
    d <- dim(data)[2]
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
    U <- reset_nas(U, args)
    U <- U[, sort(o[length(o):1], index.return = TRUE)$ix]

    return(U)
}
