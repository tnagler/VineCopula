#' Selection and Maximum Likelihood Estimation of Bivariate Copula Families
#'
#' This function selects an appropriate bivariate copula family for given
#' bivariate copula data using one of a range of methods. The corresponding
#' parameter estimates are obtained by maximum likelihood estimation.
#'
#' Copulas can be selected according to the Akaike and Bayesian Information
#' Criteria (AIC and BIC, respectively). First all available copulas are fitted
#' using maximum likelihood estimation. Then the criteria are computed for all
#' available copula families (e.g., if `u1` and `u2` are negatively
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
#' Additionally a test for independence can be performed beforehand.
#'
#' @param u1,u2 Data vectors of equal length with values in \eqn{[0,1]}.
#' @param familyset Vector of bivariate copula families to select from.
#' The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' If `familyset = NA` (default), selection among all possible families is
#' performed. If a vector of negative numbers is provided, selection among all
#' but `abs(familyset)` families is performed. Coding of bivariate copula
#' families: \cr
#' `0` = independence copula \cr
#' `1` = Gaussian copula \cr
#' `2` = Student t copula (t-copula) \cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `7` = BB1 copula \cr
#' `8` = BB6 copula \cr
#' `9` = BB7 copula \cr
#' `10` = BB8 copula \cr
#' `13` = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' `14` = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' `16` = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' `17` = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' `18` = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' `19` = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' `20` = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' `23` = rotated Clayton copula (90 degrees) \cr
#' `24` = rotated Gumbel copula (90 degrees) \cr
#' `26` = rotated Joe copula (90 degrees) \cr
#' `27` = rotated BB1 copula (90 degrees) \cr
#' `28` = rotated BB6 copula (90 degrees) \cr
#' `29` = rotated BB7 copula (90 degrees) \cr
#' `30` = rotated BB8 copula (90 degrees) \cr
#' `33` = rotated Clayton copula (270 degrees) \cr
#' `34` = rotated Gumbel copula (270 degrees) \cr
#' `36` = rotated Joe copula (270 degrees) \cr
#' `37` = rotated BB1 copula (270 degrees) \cr
#' `38` = rotated BB6 copula (270 degrees) \cr
#' `39` = rotated BB7 copula (270 degrees) \cr
#' `40` = rotated BB8 copula (270 degrees) \cr
#' `104` = Tawn type 1 copula \cr
#' `114` = rotated Tawn type 1 copula (180 degrees) \cr
#' `124` = rotated Tawn type 1 copula (90 degrees) \cr
#' `134` = rotated Tawn type 1 copula (270 degrees) \cr
#' `204` = Tawn type 2 copula \cr
#' `214` = rotated Tawn type 2 copula (180 degrees) \cr
#' `224` = rotated Tawn type 2 copula (90 degrees) \cr
#' `234` = rotated Tawn type 2 copula (270 degrees) \cr
#' @param selectioncrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: `selectioncrit = "AIC"` (default),
#' `"BIC"`, or `"logLik"`.
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' `u1` and `u2` is performed before bivariate copula selection
#' (default: `indeptest = FALSE`; see [BiCopIndTest()]).  The
#' independence copula is chosen if the null hypothesis of independence cannot
#' be rejected.
#' @param level Numeric; significance level of the independence test (default:
#' `level = 0.05`).
#' @param weights Numerical; weights for each observation (optional).
#' @param rotations If `TRUE`, all rotations of the families in
#' `familyset` are included (or subtracted).
#' @param se Logical; whether standard error(s) of parameter estimates is/are
#' estimated (default: `se = FALSE`).
#' @param presel Logical; whether to exclude families before fitting based on
#' symmetry properties of the data. Makes the selection about 30% faster
#' (on average), but may yield slightly worse results in few special cases.
#' @param method indicates the estimation method: either maximum
#' likelihood estimation (`method = "mle"`; default) or inversion of
#' Kendall's tau (`method = "itau"`). For `method = "itau"` only
#' one parameter families and the Student t copula can be used (`family =
#' 1,2,3,4,5,6,13,14,16,23,24,26,33,34` or `36`). For the t-copula,
#' `par2` is found by a crude profile likelihood optimization over the
#' interval (2, 10].
#'
#' @return An object of class [BiCop()], augmented with the following
#' entries:
#' \item{se, se2}{standard errors for the parameter estimates (if
#' `se = TRUE`,}
#' \item{nobs}{number of observations,}
#' \item{logLik}{log likelihood}
#' \item{AIC}{Aikaike's Informaton Criterion,}
#' \item{BIC}{Bayesian's Informaton Criterion,}
#' \item{emptau}{empirical value of Kendall's tau,}
#' \item{p.value.indeptest}{p-value of the independence test.}
#'
#' @note For a comprehensive summary of the fitted model, use
#' `summary(object)`; to see all its contents, use `str(object)`.
#'
#' @author Eike Brechmann, Jeffrey Dissmann, Thomas Nagler
#'
#' @note The parameters of the Student t and BB copulas are restricted (see
#' defaults in [BiCopEst()] to avoid being to close to their limiting
#' cases.
#'
#' @seealso
#' [BiCop()],
#' [BiCopEst()],
#' [RVineStructureSelect()],
#' [RVineCopSelect()],
#' [BiCopIndTest()],
#'
#'
#' @references Akaike, H. (1973). Information theory and an extension of the
#' maximum likelihood principle. In B. N. Petrov and F. Csaki (Eds.),
#' Proceedings of the Second International Symposium on Information Theory
#' Budapest, Akademiai Kiado, pp. 267-281.
#'
#' Brechmann, E. C. (2010). Truncated and simplified regular vines and their
#' applications. Diploma thesis, Technische Universitaet Muenchen.\cr
#' <http://mediatum.ub.tum.de/?id=1079285>.
#'
#' Manner, H. (2007). Estimation and model selection of copulas with an
#' application to exchange rates. METEOR research memorandum 07/056, Maastricht
#' University.
#'
#' Schwarz, G. E. (1978). Estimating the dimension of a model. Annals of
#' Statistics 6 (2), 461-464.
#'
#' @examples
#' ## Example 1: Gaussian copula with large dependence parameter
#' par <- 0.7
#' fam <- 1
#' dat1 <- BiCopSim(500, fam, par)
#' # select the bivariate copula family and estimate the parameter(s)
#' cop1 <- BiCopSelect(dat1[, 1], dat1[, 2], familyset = 1:10,
#'                     indeptest = FALSE, level = 0.05)
#' cop1  # short overview
#' summary(cop1)  # comprehensive overview
#' str(cop1)  # see all contents of the object
#'
#' ## Example 2: Gaussian copula with small dependence parameter
#' par <- 0.01
#' fam <- 1
#' dat2 <- BiCopSim(500, fam, par)
#' # select the bivariate copula family and estimate the parameter(s)
#' cop2 <- BiCopSelect(dat2[, 1], dat2[, 2], familyset = 0:10,
#'                     indeptest = TRUE, level = 0.05)
#' summary(cop2)
#'
#' ## Example 3: empirical data
#' data(daxreturns)
#' cop3 <- BiCopSelect(daxreturns[, 1], daxreturns[, 4], familyset = 0:10)
#' summary(cop3)
#'
BiCopSelect <- function(u1, u2, familyset = NA, selectioncrit = "AIC",
                        indeptest = FALSE, level = 0.05, weights = NA,
                        rotations = TRUE, se = FALSE, presel = TRUE,
                        method = "mle") {
    if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
        stop("Selection criterion not implemented.")
    ## preprocessing of arguments
    if (presel)
        todo_fams <- todo_fams_presel
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    prep_familyset,
                    check_twoparams,
                    check_est_pars,
                    check_fam_tau,
                    todo_fams,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    # perform independence test
    p.value.indeptest <- BiCopIndTest(u1, u2)$p.value

    if (indeptest & (p.value.indeptest >= level)) {
        ## select independence copula, if not rejected
        obj <- BiCop(0)
    } else {
        optiout <- list()
        lls  <- rep(Inf, length(familyset))
        AICs <- rep(Inf, length(familyset))
        BICs <- rep(Inf, length(familyset))
        ## maximum likelihood estimation
        for (i in seq_along(familyset)) {
            optiout[[i]] <- BiCopEst.intern(u1, u2,
                                            family = familyset[i],
                                            method = method,
                                            se = se,
                                            weights = weights,
                                            as.BiCop = FALSE)
            if (any(is.na(weights))) {
                lls[i] <- sum(log(BiCopPDF(u1,
                                           u2,
                                           familyset[i],
                                           optiout[[i]]$par,
                                           optiout[[i]]$par2,
                                           check.pars = FALSE)))
            } else {
                lls[i] <- sum(log(BiCopPDF(u1,
                                           u2,
                                           familyset[i],
                                           optiout[[i]]$par,
                                           optiout[[i]]$par2,
                                           check.pars = FALSE)) %*% weights)
            }
            npars <- ifelse(familyset[i] %in% allfams[onepar], 1, 2)
            if (familyset[i] == 0)
                npars <- 0
            AICs[i] <- -2 * lls[i] + 2 * npars
            BICs[i] <- -2 * lls[i] + log(length(u1)) * npars
        }

        ## select the best fitting model
        sel <- switch(selectioncrit,
                      "logLik" = which.max(lls),
                      "AIC"    = which.min(AICs),
                      "BIC"    = which.min(BICs))

        ## store results in BiCop object (dependence measures are calculated)
        obj <- BiCop(optiout[[sel]]$family,
                     optiout[[sel]]$par,
                     optiout[[sel]]$par2,
                     check.pars = FALSE)
    }

    ## add more information about the fit
    if (obj$family == 0) {
        if (se)
            obj$se  <- NA
        obj$nobs   <- length(u1)
        obj$logLik <- 0
        obj$AIC    <- 0
        obj$BIC    <- 0
    } else {
        if (se) {
            obj$se <- optiout[[sel]]$se
            if (obj$family %in% allfams[twopar])
                obj$se2 <- optiout[[sel]]$se2
        }
        obj$nobs   <- length(u1)
        obj$logLik <- lls[sel]
        obj$AIC    <- AICs[sel]
        obj$BIC    <- BICs[sel]
    }
    obj$emptau <- args$emp_tau
    obj$p.value.indeptest <- p.value.indeptest

    ## store the call that created the BiCop object
    obj$call <- match.call()

    ## return final BiCop objectz
    obj
}


##### ----------------------------------------------------------------------
## function for augmenting a familyset with rotations
with_rotations <- function(nums) {
    unique(unlist(lapply(nums, get_rotations)))
}

get_rotations <- function(fam) {
    sgn <- sign(fam)  # indicator for negative selection
    fam <- sgn * fam  # ensure that fam is positive from here on

    if (fam %in% c(0, 1, 2, 5)) {
        # no roations for independence, gaussian, student and frank copulas
        out <- fam
    } else if (fam %in% c(3, 13, 23, 33)) {
        out <- c(3, 13, 23, 33)
    } else if(fam %in% c(4, 14, 24, 34)) {
        out <- c(4, 14, 24, 34)
    } else if(fam %in% c(6, 16, 26, 36)) {
        out <- c(6, 16, 26, 36)
    } else if(fam %in% c(7, 17, 27, 37)) {
        out <- c(7, 17, 27, 37)
    } else if(fam %in% c(8, 18, 28, 38)) {
        out <- c(8, 18, 28, 38)
    } else if(fam %in% c(9, 19, 29, 39)) {
        out <- c(9, 19, 29, 39)
    } else if(fam %in% c(10, 20, 30, 40)) {
        out <- c(10, 20, 30, 40)
    } else if(fam %in% c(104, 114, 124, 134)) {
        out <- c(104, 114, 124, 134)
    } else if(fam %in% c(204, 214, 224, 234)) {
        out <- c(204, 214, 224, 234)
    }

    # adjust for negative selection
    sgn * out
}
