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

    par[is.na(par)] <- 0
    par[upper.tri(par, diag = T)] <- 0
    par2[is.na(par2)] <- 0
    par2[upper.tri(par2, diag = T)] <- 0

    if (any(par != NA) & dim(par)[1] != dim(par)[2])
        stop("Parameter matrix has to be quadratic.")
    if (any(par2 != NA) & dim(par2)[1] != dim(par2)[2])
        stop("Second parameter matrix has to be quadratic.")

    family <- RVM$family

    if (!all(par %in% c(0, NA))) {
        for (i in 2:dim(RVM$Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if ((family[i, j] == 1 || family[i, j] == 2) && abs(par[i, j]) >= 1)
                    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
                if (family[i, j] == 2 && par2[i, j] <= 2)
                    stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
                if ((family[i, j] == 3 || family[i, j] == 13) && par[i, j] <= 0)
                    stop("The parameter of the Clayton copula has to be positive.")
                if ((family[i, j] == 4 || family[i, j] == 14) && par[i, j] < 1)
                    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
                if ((family[i, j] == 6 || family[i, j] == 16) && par[i, j] <= 1)
                    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
                if (family[i, j] == 5 && par[i, j] == 0)
                    stop("The parameter of the Frank copula has to be unequal to 0.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par[i, j] <= 0)
                    stop("The first parameter of the BB1 copula has to be positive.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par2[i,  j] < 1)
                    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par[i,  j] <= 0)
                    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par2[i, j] < 1)
                    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par[i, j] < 1)
                    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par2[i, j] <= 0)
                    stop("The second parameter of the BB7 copula has to be positive.")
                if ((family[i, j] == 10 || family[i, j] == 20) && par[i, j] < 1)
                    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 10 || family[i, j] == 20) && (par2[i, j] <= 0 || par2[i, j] > 1))
                    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
                if ((family[i, j] == 23 || family[i, j] == 33) && par[i, j] >= 0)
                    stop("The parameter of the rotated Clayton copula has to be negative.")
                if ((family[i, j] == 24 || family[i, j] == 34) && par[i, j] > -1)
                    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 26 || family[i, j] == 36) && par[i, j] >= -1)
                    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
                if ((family[i, j] == 27 || family[i, j] == 37) && par[i, j] >= 0)
                    stop("The first parameter of the rotated BB1 copula has to be negative.")
                if ((family[i, j] == 27 || family[i, j] == 37) && par2[i, j] > -1)
                    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par[i, j] >= 0)
                    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par2[i, j] > -1)
                    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par[i, j] > -1)
                    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par2[i, j] >= 0)
                    stop("The second parameter of the rotated BB7 copula has to be negative.")
                if ((family[i, j] == 30 || family[i, j] == 40) && par[i, j] > -1)
                    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 30 || family[i, j] == 40) && (par2[i, j] >= 0 || par2[i, j] < (-1)))
                    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && par[i, j] < 1)
                    stop("Please choose 'par' of the Tawn copula in [1,oo).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && (par2[i, j] < 0 ||  par2[i, j] > 1))
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && par[i, j] > -1)
                    stop("Please choose 'par' of the Tawn copula in (-oo,-1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && (par2[i, j] < 0 ||  par2[i, j] > 1))
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
            }
        }
    }

    npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- RVM$family %in% allfams[onepar] +
        2 * (RVM$family %in% allfams[twopar])

    RVM2 <- RVM
    RVM2$par <- par
    RVM2$par2 <- par2

    like <- RVineLogLik(data, RVM2)

    AIC <- -2 * like$loglik + 2 * npar
    pair.AIC <- -2 * like$V$value + 2 * npar_pair

    return(list(AIC = AIC, pair.AIC = pair.AIC))
}


#'@rdname RVineAIC
RVineBIC <- function(data, RVM, par = RVM$par, par2 = RVM$par2) {

    if (is.vector(data)) {
        data <- t(as.matrix(data))
    } else {
        data <- as.matrix(data)
    }
    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T
    if (n != dim(RVM))
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (is(RVM)[1] != "RVineMatrix")
        stop("'RVM' has to be an RVineMatrix object.")

    par[is.na(par)] <- 0
    par[upper.tri(par, diag = T)] <- 0
    par2[is.na(par2)] <- 0
    par2[upper.tri(par2, diag = T)] <- 0

    if (any(par != NA) & dim(par)[1] != dim(par)[2])
        stop("Parameter matrix has to be quadratic.")
    if (any(par2 != NA) & dim(par2)[1] != dim(par2)[2])
        stop("Second parameter matrix has to be quadratic.")

    family <- RVM$family

    if (!all(par %in% c(0, NA))) {
        for (i in 2:dim(RVM$Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if ((family[i, j] == 1 || family[i, j] == 2) && abs(par[i, j]) >= 1)
                    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
                if (family[i, j] == 2 && par2[i, j] <= 2)
                    stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
                if ((family[i, j] == 3 || family[i, j] == 13) && par[i, j] <= 0)
                    stop("The parameter of the Clayton copula has to be positive.")
                if ((family[i, j] == 4 || family[i, j] == 14) && par[i, j] < 1)
                    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
                if ((family[i, j] == 6 || family[i, j] == 16) && par[i, j] <= 1)
                    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
                if (family[i, j] == 5 && par[i, j] == 0)
                    stop("The parameter of the Frank copula has to be unequal to 0.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par[i, j] <= 0)
                    stop("The first parameter of the BB1 copula has to be positive.")
                if ((family[i, j] == 7 || family[i, j] == 17) && par2[i, j] < 1)
                    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par[i, j] <= 0)
                    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 8 || family[i, j] == 18) && par2[i, j] < 1)
                    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par[i, j] < 1)
                    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 9 || family[i, j] == 19) && par2[i, j] <= 0)
                    stop("The second parameter of the BB7 copula has to be positive.")
                if ((family[i, j] == 10 || family[i, j] == 20) && par[i, j] < 1)
                    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
                if ((family[i, j] == 10 || family[i, j] == 20) && (par2[i, j] <= 0 || par2[i, j] > 1))
                    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
                if ((family[i, j] == 23 || family[i, j] == 33) && par[i, j] >= 0)
                    stop("The parameter of the rotated Clayton copula has to be negative.")
                if ((family[i, j] == 24 || family[i, j] == 34) && par[i, j] > -1)
                    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 26 || family[i, j] == 36) && par[i, j] >= -1)
                    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
                if ((family[i, j] == 27 || family[i, j] == 37) && par[i, j] >= 0)
                    stop("The first parameter of the rotated BB1 copula has to be negative.")
                if ((family[i, j] == 27 || family[i, j] == 37) && par2[i,j] > -1)
                    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par[i, j] >= 0)
                    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 28 || family[i, j] == 38) && par2[i, j] > -1)
                    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par[i, j] > -1)
                    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 29 || family[i, j] == 39) && par2[i, j] >= 0)
                    stop("The second parameter of the rotated BB7 copula has to be negative.")
                if ((family[i, j] == 30 || family[i, j] == 40) && par[i, j] > -1)
                    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
                if ((family[i, j] == 30 || family[i, j] == 40) && (par2[i, j] >= 0 || par2[i, j] < (-1)))
                    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i,  j] == 204 || family[i, j] == 214) && par[i, j] < 1)
                    stop("Please choose 'par' of the Tawn copula in [1,oo).")
                if ((family[i, j] == 104 || family[i, j] == 114 || family[i, j] == 204 || family[i, j] == 214) && (par2[i, j] < 0 || par2[i, j] > 1))
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i, j] == 224 || family[i, j] == 234) && par[i, j] > -1)
                    stop("Please choose 'par' of the Tawn copula in (-oo,-1].")
                if ((family[i, j] == 124 || family[i, j] == 134 || family[i,  j] == 224 || family[i, j] == 234) && (par2[i, j] < 0 || par2[i, j] > 1))
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
            }
        }
    }

    npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
        2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
    npar_pair <- RVM$family %in% allfams[onepar] +
        2 * (RVM$family %in% allfams[twopar])

    RVM2 <- RVM
    RVM2$par <- par
    RVM2$par2 <- par2

    like <- RVineLogLik(data, RVM2)

    BIC <- -2 * like$loglik + log(T) * npar
    pair.BIC <- -2 * like$V$value + log(T) * npar_pair

    return(list(BIC = BIC, pair.BIC = pair.BIC))
}
