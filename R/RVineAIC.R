#' AIC and BIC of an R-Vine Copula Model
#'
#' These functions calculate the Akaike and Bayesian Information criteria of a
#' d-dimensional R-vine copula model for a given copula data set.
#'
#' If \eqn{k} denotes the number of parameters of an R-vine copula model with
#' log-likelihood \eqn{l_{RVine}} and parameter set
#' \eqn{\boldsymbol{\theta}}{\theta}, then the Akaike Information Criterion (AIC)
#' by Akaike (1973) is defined as
#' \deqn{AIC := -2 l_{RVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right) +
#' 2 k,}{AIC := -2 l_{RVine}(\theta|u) + 2 k,}  for observations
#' \eqn{\boldsymbol{u}=(\boldsymbol{u}_1^\prime,...,
#' \boldsymbol{u}_N^\prime)^\prime}{u=(u_1',...,u_N')'}.
#'
#' Similarly, the Bayesian Information Criterion (BIC) by Schwarz (1978) is
#' given by \deqn{ BIC := -2
#' l_{RVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right) + \log(N) k. }
#'
#' @aliases RVineAIC RVineBIC
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM An \code{\link{RVineMatrix}} object including the structure and
#' the pair-copula families and parameters.
#' @param par A d x d matrix with the pair-copula parameters (optional;
#' default: \code{par = RVM$par}).
#' @param par2 A d x d matrix with the second parameters of pair-copula
#' families with two parameters (optional; default: \code{par2 = RVM$par2}).
#' @return \item{AIC, BIC}{The computed AIC or BIC value, respectively.}
#' \item{pair.AIC, pair.BIC}{A d x d matrix of individual contributions to the
#' AIC or BIC value for each pair-copula, respectively. Note: \code{AIC =
#' sum(pair.AIC)} and similarly \code{BIC = sum(pair.BIC)}.}
#' @author Eike Brechmann
#' @seealso \code{\link{RVineLogLik}}, \code{\link{RVineVuongTest}},
#' \code{\link{RVineClarkeTest}}
#' @references Akaike, H. (1973). Information theory and an extension of the
#' maximum likelihood principle. In B. N. Petrov and F. Csaki (Eds.),
#' Proceedings of the Second International Symposium on Information Theory
#' Budapest, Akademiai Kiado, pp. 267-281.
#'
#' Schwarz, G. E. (1978). Estimating the dimension of a model. Annals of
#' Statistics 6 (2), 461-464.
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
#' RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2,
#'                    names=c("V1", "V2", "V3", "V4", "V5"))
#'
#' # simulate a sample of size 300 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(300,RVM)
#'
#' # compute AIC and BIC
#' RVineAIC(simdata, RVM)
#' RVineBIC(simdata, RVM)
#'
RVineAIC <- function(data, RVM, par = RVM$par, par2 = RVM$par2) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    ## number of parameters
    npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- RVM$family %in% allfams[onepar] +
        2 * (RVM$family %in% allfams[twopar])

    ## calculate AICs
    like <- RVineLogLik(data, RVM)
    AIC <- -2 * like$loglik + 2 * npar
    pair.AIC <- -2 * like$V$value + 2 * npar_pair

    ## return results
    list(AIC = AIC, pair.AIC = pair.AIC)
}


#'@rdname RVineAIC
RVineBIC <- function(data, RVM, par = RVM$par, par2 = RVM$par2) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- RVM$family %in% allfams[onepar] +
        2 * (RVM$family %in% allfams[twopar])

    like <- RVineLogLik(data, RVM)

    BIC <- -2 * like$loglik + log(n) * npar
    pair.BIC <- -2 * like$V$value + log(n) * npar_pair

    return(list(BIC = BIC, pair.BIC = pair.BIC))
}
