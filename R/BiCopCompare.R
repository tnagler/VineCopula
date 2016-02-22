#' Compare Maximum Likelihood Estimates of Bivariate Copula Families
#'
#' This function allows to compare bivariate copula models accross a number of
#' families w.r.t. the fit statistics log-likelihood, AIC, and BIC. For each
#' family, the parameters are estimated by maximum likelihood.
#'
#' First all available copulas are fitted using maximum likelihood estimation.
#' Then the criteria are computed for all available copula families (e.g., if
#' \code{u1} and \code{u2} are negatively
#' dependent, Clayton, Gumbel, Joe, BB1, BB6, BB7 and BB8 and their survival
#' copulas are not considered) and the family with the minimum value is chosen.
#' For observations \eqn{u_{i,j},\ i=1,...,N,\ j=1,2,}{u_{i,j}, i=1,...,N,\
#' j=1,2,} the AIC of a bivariate copula family \eqn{c} with parameter(s)
#' \eqn{\boldsymbol{\theta}} is defined as \deqn{AIC := -2 \sum_{i=1}^N
#' \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + 2k, }{ AIC := -2 \sum_{i=1}^N
#' ln[c(u_{i,1},u_{i,2}|\theta)] + 2k, } where \eqn{k=1} for one parameter
#' copulas and \eqn{k=2} for the two parameter t-, BB1, BB6, BB7 and BB8
#' copulas. Similarly, the BIC is given by \deqn{BIC := -2 \sum_{i=1}^N
#' \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + \ln(N)k. }{ BIC := -2
#' \sum_{i=1}^N ln[c(u_{i,1},u_{i,2}|\theta)] + ln(N)k. } Evidently, if the BIC
#' is chosen, the penalty for two parameter families is stronger than when
#' using the AIC.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param familyset Vector of bivariate copula families to select from.
#' The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' If \code{familyset = NA} (default), selection among all possible families is
#' performed.  Coding of bivariate copula families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param weights Numerical; weights for each observation (optional).
#' @param rotations If \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#' @param ... further arguments passed to \code{\link{BiCopEst}}.
#'
#' @return A list containing
#' \item{models}{a list of \code{\link{BiCop}} objects corresponding to the
#' familyset (only families corresponding to the sign of the empirical
#' Kendall's tau are used),}
#' \item{summary}{a data frame containing the log-likelihoods, AICs, and BICs
#' of all the fitted models.}
#'
#'
#' @author Thomas Nagler
#'
#' @seealso
#' \code{\link{BiCop}},
#' \code{\link{BiCopEst}}
#'
#'
#' @references Akaike, H. (1973). Information theory and an extension of the
#' maximum likelihood principle. In B. N. Petrov and F. Csaki (Eds.),
#' Proceedings of the Second International Symposium on Information Theory
#' Budapest, Akademiai Kiado, pp. 267-281.
#'
#' Schwarz, G. E. (1978). Estimating the dimension of a model. Annals of
#' Statistics 6 (2), 461-464.
#'
#' @examples
#'
#' ## compare models
#' data(daxreturns)
#' comp <- BiCopCompare(daxreturns[, 1], daxreturns[, 4])
#'
#'
BiCopCompare <- function(u1, u2, familyset = NA, weights = NA, rotations = TRUE,
                         ...) {
    if (is.na(familyset[1]))
        familyset <- c(0, allfams)

    ## sanity checks
    if ((is.null(u1) == TRUE) || (is.null(u2) == TRUE))
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2))
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (length(u1) < 2)
        stop("Number of observations has to be at least 2.")
    if (any(u1 > 1) || any(u1 < 0))
        stop("Data has to be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0))
        stop("Data has to be in the interval [0,1].")
    if (!all(familyset %in% c(0, allfams)))
        stop("Copula family not implemented.")

    ## adjust familyset if rotations = TRUE
    if (rotations)
        familyset <- with_rotations(familyset)

    # calculate empirical kendall's tau
    emp_tau <- fasttau(u1, u2, weights)

    ## perform independence test
    p.value.indeptest <- BiCopIndTest(u1, u2)$p.value

    ## sets of families for negative and positive dependence
    negfams <- c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40, 124, 134, 224, 234)
    posfams <- c(1:10, 13, 14, 16:20, 104, 114, 204, 214)

    ## stop if familyset not sufficient
    if (!is.na(familyset[1]) &&
        !(any(familyset %in% negfams) && any(familyset %in% posfams))) {
        txt <- paste0("'familyset' has to include at least one bivariate ",
                      "copula family for positive and one for negative ",
                      "dependence.")
        stop(txt)
    }

    ## find families for which estimation is required
    ## (only families that allow for the empirical kendall's tau)
    if (emp_tau < 0) {
        todo <- c(0, negfams)
    } else {
        todo <- c(0, posfams)
    }
    todo <- todo[which(todo %in% familyset)]

    ## maximum likelihood estimation
    optiout <- list()
    for (i in seq_along(todo)) {
        optiout[[i]] <- BiCopEst.intern(u1, u2,
                                        family = todo[i],
                                        as.BiCop = TRUE,
                                        ...)
    }


    ## store model information in data.frame
    tab <- data.frame(family = todo,
                      logLik = sapply(optiout, function(x) x$logLik),
                      AIC = sapply(optiout, function(x) x$AIC),
                      BIC = sapply(optiout, function(x) x$BIC))

    ## return along with estimated models
    list(models = optiout, summary = round(tab, 2))
}