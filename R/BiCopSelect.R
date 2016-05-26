#' Selection and Maximum Likelihood Estimation of Bivariate Copula Families
#'
#' This function selects an appropriate bivariate copula family for given
#' bivariate copula data using one of a range of methods. The corresponding
#' parameter estimates are obtained by maximum likelihood estimation.
#'
#' Copulas can be selected according to the Akaike and Bayesian Information
#' Criteria (AIC and BIC, respectively). First all available copulas are fitted
#' using maximum likelihood estimation. Then the criteria are computed for all
#' available copula families (e.g., if \code{u1} and \code{u2} are negatively
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
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param familyset Vector of bivariate copula families to select from.
#' The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' If \code{familyset = NA} (default), selection among all possible families is
#' performed. If a vector of negative numbers is provided, selection among all
#' but \code{abs(familyset)} families is performed. Coding of bivariate copula
#' families: \cr
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
#' @param selectioncrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{selectioncrit = "AIC"} (default),
#' \code{"BIC"}, or \code{"logLik"}.
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).  The
#' independence copula is chosen if the null hypothesis of independence cannot
#' be rejected.
#' @param level Numeric; significance level of the independence test (default:
#' \code{level = 0.05}).
#' @param weights Numerical; weights for each observation (optional).
#' @param rotations If \code{TRUE}, all rotations of the families in
#' \code{familyset} are included (or substracted).
#' @param se Logical; whether standard error(s) of parameter estimates is/are
#' estimated (default: \code{se = TRUE}).
#'
#' @return An object of class \code{\link{BiCop}}, augmented with the following
#' entries:
#' \item{se, se2}{standard errors for the parameter estimates (if
#' \code{se = TRUE},}
#' \item{nobs}{number of observations,}
#' \item{logLik}{log likelihood}
#' \item{AIC}{Aikaike's Informaton Criterion,}
#' \item{BIC}{Bayesian's Informaton Criterion,}
#' \item{emptau}{empirical value of Kendall's tau,}
#' \item{p.value.indeptest}{p-value of the independence test.}
#'
#' @note For a comprehensive summary of the fitted model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}.
#'
#' @author Thomas Nagler
#'
#' @note The parameters of the Student t and BB copulas are restricted (see
#' defaults in \code{\link{BiCopEst}} to avoid being to close to their limiting
#' cases.
#'
#' @seealso
#' \code{\link{BiCop}},
#' \code{\link{BiCopEst}},
#' \code{\link{RVineStructureSelect}},
#' \code{\link{RVineCopSelect}},
#' \code{\link{BiCopIndTest}},
#'
#'
#' @references Akaike, H. (1973). Information theory and an extension of the
#' maximum likelihood principle. In B. N. Petrov and F. Csaki (Eds.),
#' Proceedings of the Second International Symposium on Information Theory
#' Budapest, Akademiai Kiado, pp. 267-281.
#'
#' Brechmann, E. C. (2010). Truncated and simplified regular vines and their
#' applications. Diploma thesis, Technische Universitaet Muenchen.\cr
#' \url{http://mediatum.ub.tum.de/?id=1079285}.
#'
#' Manner, H. (2007). Estimation and model selection of copulas with an
#' application to exchange rates. METEOR research memorandum 07/056, Maastricht
#' University.
#'
#' Schwarz, G. E. (1978). Estimating the dimension of a model. Annals of
#' Statistics 6 (2), 461-464.
#'
#' @examples
#'
#' ## Example 1: Gaussian copula with large dependence parameter
#' par <- 0.7
#' fam <- 1
#' dat1 <- BiCopSim(500, fam, par)
#'
#' # select the bivariate copula family and estimate the parameter(s)
#' cop1 <- BiCopSelect(dat1[, 1], dat1[, 2], familyset = 1:10,
#'                     indeptest = FALSE, level = 0.05)
#' cop1  # short overview
#' summary(cop1)  # comprehensive overview
#' str(cop1)  # see all contents of the object
#'
#'
#' ## Example 2: Gaussian copula with small dependence parameter
#' par <- 0.01
#' fam <- 1
#' dat2 <- BiCopSim(500, fam, par)
#'
#' # select the bivariate copula family and estimate the parameter(s)
#' cop2 <- BiCopSelect(dat2[, 1], dat2[, 2], familyset = 0:10,
#'                     indeptest = TRUE, level = 0.05)
#' summary(cop2)
#'
#'
#' ## Example 3: empirical data
#' data(daxreturns)
#' cop3 <- BiCopSelect(daxreturns[, 1], daxreturns[, 4], familyset = 0:10)
#' summary(cop3)
#'
BiCopSelect <- function(u1, u2, familyset = NA, selectioncrit = "AIC",
                        indeptest = FALSE, level = 0.05, weights = NA,
                        rotations = TRUE, se = TRUE) {
    if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
        stop("Selection criterion not implemented.")
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    prep_familyset,
                    check_est_pars,
                    check_fam_tau,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    # perform independence test
    p.value.indeptest <- BiCopIndTest(u1, u2)$p.value

    if (indeptest & (p.value.indeptest >= level)) {
        ## select independence copula, if not rejected
        obj <- BiCop(0)
    } else {
        emp_tau <- args$emp_tau
        ## find families for which estimation is required
        ## (only families that allow for the empirical kendall's tau)
        if (emp_tau < 0) {
            todo <- negfams
        } else if (emp_tau < 0) {
            todo <- posfams
        } else {
            todo <- c(negfams, posfams)
        }
        todo <- todo[which(todo %in% familyset)]

        ## maximum likelihood estimation
        optiout <- list()
        for (i in seq_along(todo)) {
            optiout[[i]] <- BiCopEst.intern(u1, u2,
                                            family = todo[i],
                                            se = se,
                                            weights = weights,
                                            as.BiCop = FALSE)
        }

        ## calculate logLik, AIC and BIC
        lls  <- rep(Inf, length(todo))
        AICs <- rep(Inf, length(todo))
        BICs <- rep(Inf, length(todo))
        for (i in seq_along(todo)) {
            if (any(is.na(weights))) {
                lls[i] <- sum(log(BiCopPDF(u1,
                                           u2,
                                           todo[i],
                                           optiout[[i]]$par,
                                           optiout[[i]]$par2,
                                           check.pars = FALSE)))
            } else {
                lls[i] <- sum(log(BiCopPDF(u1,
                                           u2,
                                           todo[i],
                                           optiout[[i]]$par,
                                           optiout[[i]]$par2,
                                           check.pars = FALSE)) %*% weights)
            }
            npars <- ifelse(todo[i] %in% allfams[onepar], 1, 2)
            AICs[i] <- -2 * lls[i] + 2 * npars
            BICs[i] <- -2 * lls[i] + log(length(u1)) * npars

        }

        ## add independence copula
        if (0 %in% familyset) {
            optiout[[length(todo) + 1]] <- list(family = 0, par = 0, par2 = 0)
            lls[length(todo) + 1] <- 0
            AICs[length(todo) + 1] <- 0
            BICs[length(todo) + 1] <- 0
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
    obj$emptau <- emp_tau
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
