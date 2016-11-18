#' Log-Likelihood of an R-Vine Copula Model
#'
#' This function calculates the log-likelihood of a d-dimensional R-vine copula
#' model for a given copula data set.
#'
#' For observations
#' \eqn{\boldsymbol{u}=(\boldsymbol{u}_1^\prime,...,\boldsymbol{u}_N^\prime)^\prime}{u=(u_1',...,u_N')'}
#' the log-likelihood of a \eqn{d}-dimensional R-vine copula with \eqn{d-1}
#' trees and corresponding edge sets \eqn{E_1,...,E_{d-1}} is given by
#' \deqn{\texttt{loglik}:=l_{RVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right)
#' }{loglik:=ll_{RVine}(\theta|u)}
#' \deqn{=\sum_{i=1}^N \sum_{\ell=1}^{d-1} \sum_{e\in E_\ell}
#'  \ln[c_{j(e),k(e)|D(e)}(F(u_{i,j(e)}|u_{i,D(e)}),F(u_{i,k(e)}|u_{i,D(e)})|\theta_{j(e),k(e)|D(e)})] }{
#'  =\sum_{i=1}^N \sum_{k=1}^{d-1} \sum_{e\in E_k}
#' \ln[c_{j(e),k(e)|D(e)}(F(u_{i,j(e)}|u_{i,D(e)}),F(u_{i,k(e)}|u_{i,D(e)})|\theta_{j(e),k(e)|D(e)})],
#' }
#'  where \eqn{\boldsymbol{u}_i=(u_{i,1},...,u_{i,d})^\prime\in[0,1]^d,\
#' i=1,...,N}{u_i=(u_{i,1},...,u_{i,d})'\in[0,1]^d, i=1,...,N}. Further
#' \eqn{c_{j(e),k(e)|D(e)}} denotes a bivariate copula density associated to an
#' edge \eqn{e} and with parameter(s)
#' \eqn{\boldsymbol{\theta}_{j(e),k(e)|D(e)}}{\theta_{j(e),k(e)|D(e)}}.
#' Conditional distribution functions such as
#' \eqn{F(u_{i,j(e)}|\boldsymbol{u}_{i,D(e)})}{F(u_{i,j(e)}|u_{i,D(e)})} are
#' obtained recursively using the relationship
#' \deqn{h(u|\boldsymbol{v},\boldsymbol{\theta}) := F(u|\boldsymbol{v}) =
#' d C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j}))/d
#' F(v_j|v_{-j}), }{
#' h(u|v,\theta) := F(u|v) = d C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j}))/d
#' F(v_j|v_{-j}), }
#' where
#' \eqn{C_{uv_j|\boldsymbol{v}_{-j}}}{C_{uv_j|v_{-j}}} is a bivariate copula
#' distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}
#' and \eqn{\boldsymbol{v}_{-j}}{v_{-j}} denotes a vector with the \eqn{j}-th
#' component \eqn{v_j} removed. The notation of h-functions is introduced for
#' convenience. For more details see Dissmann et al. (2013).
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An \code{\link{RVineMatrix}} object including the structure and
#' the pair-copula families and parameters.
#' @param par A d x d matrix with the pair-copula parameters (optional;
#' default: \code{par = RVM$par}).
#' @param par2 A d x d matrix with the second parameters of pair-copula
#' families with two parameters (optional; default: \code{par2 = RVM$par2}).
#' @param separate Logical; whether log-likelihoods are returned point wisely
#' (default: \code{separate = FALSE}).
#' @param verbose In case something goes wrong, additional output will be
#' plotted.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#' @param calculate.V logical; whether V matrices should be calculated. Default
#' is \code{TRUE}, but requires a lot of memory when dimension is large.
#' Use \code{FALSE} for a memory efficient version.
#'
#' @return \item{loglik}{The calculated log-likelihood value of the R-vine
#' copula model.} \item{V}{The stored transformations (h-functions and
#' log-likelihoods of each pair-copula) which may be used for posterior updates
#' (three matrices: \code{direct}, \code{indirect} and \code{value}).}
#'
#' @author Ulf Schepsmeier, Jeffrey Dissmann, Jakob Stoeber
#'
#' @seealso \code{\link{BiCopHfunc}}, \code{\link{RVineMatrix}},
#' \code{\link{RVineMLE}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}}
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' @examples
#'
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#'
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#'
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#'
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' # simulate a sample of size 300 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(300, RVM)
#'
#' # compute the log-likelihood
#' ll <- RVineLogLik(simdata, RVM, separate = FALSE)
#' ll$loglik
#'
#' # compute the pointwise log-likelihoods
#' ll <- RVineLogLik(simdata, RVM, separate = TRUE)
#' ll$loglik
#'
RVineLogLik <- function(data, RVM, par = RVM$par, par2 = RVM$par2,
                        separate = FALSE, verbose = TRUE, check.pars = TRUE,
                        calculate.V = TRUE) {
    ## preprocessing of arguments
    RVM$par <- par
    RVM$par2 <- par2
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    fix_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs)
    list2env(args, environment())
    par <- RVM$par
    par2 <- RVM$par2

    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T

    o <- diag(RVM$Matrix)
    if (any(o != length(o):1)) {
        oldRVM <- RVM
        RVM <- normalizeRVineMatrix(RVM)
        data <- data[, o[length(o):1]]
    }

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

    if (calculate.V) {
        V <- list()
        V$direct <- array(0, dim = c(n, n, N))
        V$indirect <- array(0, dim = c(n, n, N))
        if (is.vector(data)) {
            V$direct[n, , ] <- data[n:1]
        } else {
            V$direct[n, , ] <- t(data[, n:1])
        }
        V$value <- array(0, c(n, n, N))

        ll <- as.vector(V$value)
        vv <- as.vector(V$direct)
        vv2 <- as.vector(V$indirect)
        calcup <- as.vector(matrix(1, dim(RVM), dim(RVM)))

        out <- rep(0, N)
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
                  as.integer(TRUE),
                  PACKAGE = 'VineCopula'
        )

        ll <- out[[12]]
        loglik <- out[[11]]
        loglik[loglik %in% c(NaN, -Inf, Inf)] <- -1e+10
        vv <- out[[13]]
        vv2 <- out[[14]]
        V$direct <- array(vv, dim = c(n, n, N))
        V$indirect <- array(vv2, dim = c(n, n, N))
        V$value <- array(ll, dim = c(n, n, N))
        V <- suppressWarnings(lapply(V, reset_nas, args = args))

        if (separate) {
            loglik <- reset_nas(loglik, args)
        } else {
            if (!is.null(args$msg))
                args$msg <- paste(args$msg, "Only complete observations are used.")
            loglik <- sum(reset_nas(loglik, args), na.rm = TRUE)
            V$value <- apply(V$value, 1:2, sum, na.rm = TRUE)
        }
        if (any(V$value %in% c(NaN, -Inf, Inf)) & verbose) {
            print(V$value[V$value %in% c(NaN, -Inf, Inf)])
            print(th)
            print(th2)
        }
        V$value[V$value %in% c(NaN, -Inf, Inf)] <- -1e+10

        return(list(loglik = loglik, V = V))

    } else {
        # memory efficient version (not calculating V matrices)
        out <- rep(0, N)
        out <- .C("VineLogLikRvine2",
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
                  PACKAGE = 'VineCopula'
        )

        loglik <- out[[11]]
        if (separate) {
            loglik <- reset_nas(loglik, args)
        } else {
            if (!is.null(args$msg))
                args$msg <- paste(args$msg, "Only complete observations are used.")
            loglik <- sum(reset_nas(loglik, args), na.rm = TRUE)
        }

        return(list(loglik = loglik))
    }
}

#' PDF of an R-Vine Copula Model
#'
#' This function calculates the probability density function of a d-dimensional
#' R-vine copula.
#'
#' The density of a \eqn{d}-dimensional R-vine copula with \eqn{d-1} trees and
#' corresponding edge sets \eqn{E_1,...,E_{d-1}} is given by
#' \deqn{\prod_{\ell=1}^{d-1} \prod_{e\in E_\ell }
#' c_{j(e),k(e)|D(e)}(F(u_{j(e)}|u_{D(e)}),F(u_{k(e)}|u_{D(e)})|\theta_{j(e),k(e)|D(e)}), }{
#'  =\prod_{m=1}^{d-1} \prod_{e\in E_m}
#' c_{j(e),k(e)|D(e)}(F(u_{j(e)}|u_{D(e)}),F(u_{k(e)}|u_{D(e)})|\theta_{j(e),k(e)|D(e)}),
#' }
#' where
#' \eqn{\boldsymbol{u}=(u_{1},...,u_{d})^\prime\in[0,1]^d}{u=(u_{1},...,u_{d})'\in[0,1]^d}.
#' Further \eqn{c_{j(e),k(e)|D(e)}} denotes a bivariate copula density
#' associated to an edge \eqn{e} and with parameter(s)
#' \eqn{\boldsymbol{\theta}_{j(e),k(e)|D(e)}}{\theta_{j(e),k(e)|D(e)}}.
#' Conditional distribution functions such as
#' \eqn{F(u_{j(e)}|\boldsymbol{u}_{D(e)})}{F(u_{j(e)}|u_{D(e)})} are obtained
#' recursively using the relationship
#' \deqn{h(u|\boldsymbol{v},\boldsymbol{\theta}) := F(u|\boldsymbol{v}) =
#' d C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j}))/d F(v_j|v_{-j}),}{
#' h(u|v,\theta) := F(u|v) = d C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j}))/d
#' F(v_j|v_{-j}), }
#' where
#' \eqn{C_{uv_j|\boldsymbol{v}_{-j}}}{C_{uv_j|v_{-j}}} is a bivariate copula
#' distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}
#' and \eqn{\boldsymbol{v}_{-j}}{v_{-j}} denotes a vector with the \eqn{j}-th
#' component \eqn{v_j} removed. The notation of h-functions is introduced for
#' convenience. For more details see Dissmann et al. (2013).
#'
#' The function is actually just a wrapper to \code{\link{RVineLogLik}}.
#'
#' @param newdata An N x d data matrix that specifies where the density shall
#' be evaluated.
#' @param RVM An \code{\link{RVineMatrix}} object including the structure and
#' the pair-copula families and parameters.
#'
#' @author Thomas Nagler
#'
#' @seealso \code{\link{BiCopHfunc}}, \code{\link{RVineMatrix}},
#' \code{\link{RVineMLE}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}}
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' @examples
#'
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#'
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#'
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#'
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' # compute the density at (0.1, 0.2, 0.3, 0.4, 0.5)
#' RVinePDF(c(0.1, 0.2, 0.3, 0.4, 0.5), RVM)
#'
RVinePDF <- function(newdata, RVM) {
    exp(RVineLogLik(newdata, RVM, separate = TRUE, calculate.V = FALSE)$loglik)
}
