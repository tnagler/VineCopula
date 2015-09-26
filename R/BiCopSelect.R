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
#' \eqn{\boldsymbol{\theta}} is defined as \deqn{ }{ AIC := -2 \sum_{i=1}^N
#' ln[c(u_{i,1},u_{i,2}|\theta)] + 2k, }\deqn{AIC := -2 \sum_{i=1}^N
#' \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + 2k, }{ AIC := -2 \sum_{i=1}^N
#' ln[c(u_{i,1},u_{i,2}|\theta)] + 2k, } where \eqn{k=1} for one parameter
#' copulas and \eqn{k=2} for the two parameter t-, BB1, BB6, BB7 and BB8
#' copulas. Similarly, the BIC is given by \deqn{ }{ BIC := -2 \sum_{i=1}^N
#' ln[c(u_{i,1},u_{i,2}|\theta)] + ln(N)k. }\deqn{BIC := -2 \sum_{i=1}^N
#' \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + \ln(N)k. }{ BIC := -2
#' \sum_{i=1}^N ln[c(u_{i,1},u_{i,2}|\theta)] + ln(N)k. } Evidently, if the BIC
#' is chosen, the penalty for two parameter families is stronger than when
#' using the AIC.
#'
#' Additionally a test for independence can be performed beforehand.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param familyset Vector of bivariate copula families to select from (the
#' independence copula MUST NOT be specified in this vector, otherwise it will
#' be selected).  The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' Not listed copula families might be included to better handle limit cases.
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
#' @param selectioncrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{selectioncrit = "AIC"} (default) or
#' \code{"BIC"}.
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).  The
#' independence copula is chosen if the null hypothesis of independence cannot
#' be rejected.
#' @param level Numeric; significance level of the independence test (default:
#' \code{level = 0.05}).
#' @param weights Numerical; weights for each observation (optional).
#' @param rotations If \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#' @param se Logical; whether standard error(s) of parameter estimates is/are
#' estimated (default: \code{se = FALSE}).
#'
#' @return An object of class \code{\link{BiCop}}, i.e., a list containing
#' \item{family}{The selected bivariate copula family.}
#' \item{par, par2}{The estimated bivariate copula parameter(s).}
#' \item{p.value.indeptest}{P-value of the independence test if performed.}
#'
#' @note When the bivariate t-copula is considered and the degrees of freedom
#' are estimated to be larger than 30, then the bivariate Gaussian copula is
#' taken into account instead. Similarly, when BB1 (Clayton-Gumbel), BB6
#' (Joe-Gumbel), BB7 (Joe-Clayton) or BB8 (Joe-Frank) copulas are considered
#' and the parameters are estimated to be very close to one of their boundary
#' cases, the respective one parameter copula is taken into account instead.
#'
#' @author Eike Brechmann, Jeffrey Dissmann
#'
#' @seealso
#' \code{\link{RVineStructureSelect}},
#' \code{\link{RVineCopSelect}},
#' \code{\link{BiCopIndTest}},
#' \code{\link{BiCop}}
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
#' par1 <- 0.7
#' fam1 <- 1
#' dat1 <- BiCopSim(500, fam1, par1)
#'
#' # select the bivariate copula family and estimate the parameter(s)
#' cop1 <- BiCopSelect(dat1[,1], dat1[,2], familyset = c(1:10),
#'                     indeptest = FALSE, level = 0.05)
#' cop1$family
#' cop1$par
#' cop1$par2
#'
#'
#' ## Example 2: Gaussian copula with small dependence parameter
#' par2 <- 0.01
#' fam2 <- 1
#' dat2 <- BiCopSim(500, fam2, par2)
#'
#' # select the bivariate copula family and estimate the parameter(s)
#' cop2 <- BiCopSelect(dat2[,1], dat2[,2], familyset = c(1:10),
#'                     indeptest = TRUE, level = 0.05)
#' cop2$family
#' cop2$par
#' cop2$par2
#'
#'
#' ## Example 3: empirical data
#' data(daxreturns)
#' cop3 <- BiCopSelect(daxreturns[,1], daxreturns[,4],
#'                     familyset = c(1:10, 13, 14, 16,
#'                                   23, 24, 26, 33, 34, 36))
#' cop3$family
#' cop3$par
#' cop3$par2
#'
#' @export BiCopSelect
BiCopSelect <- function(u1, u2, familyset = NA, selectioncrit = "AIC",
                        indeptest = FALSE, level = 0.05, weights = NA,
                        rotations = TRUE, se = TRUE) {
    allfams <- c(1:10,
                 13, 14, 16:20,
                 23, 24, 26:30, 33, 34, 36:40,
                 104, 114, 124, 134, 204, 214, 224, 234)
    if (is.na(familyset[1]))
        familyset <- allfams

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
    if (!all(familyset %in% allfams))
        stop("Copula family not implemented.")
    if ((selectioncrit != "AIC") && (selectioncrit != "BIC"))
        stop("Selection criterion not implemented.")
    if ((level) < 0 & (level > 1))
        stop("Significance level has to be between 0 and 1.")

    ## prepare objects
    out <- list()
    data1 <- u1
    data2 <- u2

    ## adjust familyset if rotations = TRUE
    if (rotations)
        familyset <- with_rotations(familyset)

    if (!is.na(familyset[1]) & any(familyset == 0)) {
        ## select independence if allowed
        out$p.value.indeptest <- NA
        out$family <- 0
        out$par <- out$par2 <- 0
    } else {
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

        # calculate empirical kendall's tau
        emp_tau <- fasttau(data1, data2, weights)

        ## perform independence test (if asked for)
        if (indeptest == TRUE) {
            out$p.value.indeptest <- BiCopIndTest(data1, data2)$p.value
        } else {
            out$p.value.indeptest <- NA
        }

        if ((!is.na(out$p.value.indeptest)) & (out$p.value.indeptest >= level)) {
            ## select independence copula, if not rejected
            out$family <- 0
            out$par <- out$par2 <- 0
        } else {
            ## initial values for parameter optimization
            start <- list()
            start[[1]] <- sin(pi * emp_tau/2)
            start[[2]] <- c(sin(emp_tau * pi/2), 10)
            start[[3]] <- start[[13]] <- 2 * abs(emp_tau)/(1 - abs(emp_tau))
            start[[4]] <- start[[14]] <- 1/(1 - abs(emp_tau))
            if (5 %in% familyset) {
                start[[5]] <- Frank.itau.JJ(emp_tau)
            } else {
                start[[5]] <- 0
            }
            if (any(c(6, 16) %in% familyset)) {
                start[[6]] <- start[[16]] <- Joe.itau.JJ(abs(emp_tau))
            } else {
                start[[6]] <- start[[16]] <- 0
            }
            start[[7]] <- start[[17]] <- c(0.5, 1.5)
            start[[8]] <- start[[18]] <- c(1.5, 1.5)
            start[[9]] <- start[[19]] <- c(1.5, 0.5)
            start[[10]] <- start[[20]] <- c(1.5, 0.5)
            start[[23]] <- start[[33]] <- -2 * abs(emp_tau)/(1 - abs(emp_tau))
            start[[24]] <- start[[34]] <- -1/(1 - abs(emp_tau))
            if (any(c(26, 36) %in% familyset)) {
                start[[26]] <- start[[36]] <- -Joe.itau.JJ(abs(emp_tau))
            } else {
                start[[26]] <- start[[36]] <- 0
            }
            start[[27]] <- start[[37]] <- c(-0.5, -1.5)
            start[[28]] <- start[[38]] <- c(-1.5, -1.5)
            start[[29]] <- start[[39]] <- c(-1.5, -0.5)
            start[[30]] <- start[[40]] <- c(-1.5, -0.5)
            delta <- min(abs(emp_tau) + 0.1, 0.999)
            theta1 <- 1 + 6 * abs(emp_tau)
            start[[104]] <- start[[204]] <- c(theta1, delta)
            start[[114]] <- start[[214]] <- c(theta1, delta)
            start[[124]] <- start[[224]] <- c(-theta1, delta)
            start[[134]] <- start[[234]] <- c(-theta1, delta)

            ## find families for which estimation is required (only families that allow for
            ## the empirical kendall's tau)
            if (emp_tau < 0) {
                todo <- negfams
            } else {
                todo <- posfams
            }
            todo <- todo[which(todo %in% familyset)]


            ## estimate parameters for each of the families (in 'todo')
            optiout <- list()

            # t
            if (any(todo == 2)) {
                optiout[[2]] <- suppressWarnings(BiCopEst(data1,
                                                          data2,
                                                          family = 2,
                                                          max.df = 30,
                                                          weights = weights,
                                                          se = se))
                optiout[[2]]$par <- c(optiout[[2]]$par, optiout[[2]]$par2)
                if (optiout[[2]]$par[2] >= 30) {
                    todo[todo == 2] <- 1
                    todo <- unique(todo)
                    optiout[[2]] <- list()
                }
            }
            # BB1
            if (any(todo == 7)) {
                optiout[[7]] <- MLE_intern(cbind(data1, data2),
                                           start[[7]],
                                           7,
                                           weights = weights,
                                           se = se)
                if (optiout[[7]]$par[1] <= 0.1 | optiout[[7]]$par[2] <= 1.1) {
                    if (optiout[[7]]$par[1] <= 0.1) {
                        todo[todo == 7] <- 4
                        todo <- unique(todo)
                    } else if (optiout[[7]]$par[2] <= 1.1) {
                        todo[todo == 7] <- 3
                        todo <- unique(todo)
                    }
                    optiout[[7]] <- list()
                }
            }
            # BB6
            if (any(todo == 8)) {
                optiout[[8]] <- MLE_intern(cbind(data1, data2),
                                           start[[8]],
                                           8,
                                           weights = weights,
                                           se = se)
                if (optiout[[8]]$par[1] <= 1.1 | optiout[[8]]$par[2] <= 1.1) {
                    if (optiout[[8]]$par[1] <= 1.1) {
                        todo[todo == 8] <- 4
                        todo <- unique(todo)
                    } else if (optiout[[8]]$par[2] <= 1.1) {
                        todo[todo == 8] <- 6
                        todo <- unique(todo)
                    }
                    optiout[[8]] <- list()
                }
            }
            # BB7
            if (any(todo == 9)) {
                optiout[[9]] <- MLE_intern(cbind(data1, data2),
                                           start[[9]],
                                           9,
                                           weights = weights,
                                           se = se)
                if (optiout[[9]]$par[1] <= 1.1 | optiout[[9]]$par[2] <= 0.1) {
                    if (optiout[[9]]$par[1] <= 1.1) {
                        todo[todo == 9] <- 3
                        todo <- unique(todo)
                    } else if (optiout[[9]]$par[2] <= 0.1) {
                        todo[todo == 9] <- 6
                        todo <- unique(todo)
                    }
                    optiout[[9]] <- list()
                }
            }
            # BB8
            if (any(todo == 10)) {
                optiout[[10]] <- MLE_intern(cbind(data1, data2),
                                            start[[10]],
                                            10,
                                            weights = weights,
                                            se = se)
                if (optiout[[10]]$par[2] >= 0.99) {
                    todo[todo == 10] <- 6
                    todo <- unique(todo)
                    optiout[[10]] <- list()
                }
            }
            # SBB1
            if (any(todo == 17)) {
                optiout[[17]] <- MLE_intern(cbind(data1, data2),
                                            start[[17]],
                                            17,
                                            weights = weights,
                                            se = se)
                if (optiout[[17]]$par[1] <= 0.1 | optiout[[17]]$par[2] <= 1.1) {
                    if (optiout[[17]]$par[1] <= 0.1) {
                        todo[todo == 17] <- 14
                        todo <- unique(todo)
                    } else if (optiout[[17]]$par[2] <= 1.1) {
                        todo[todo == 17] <- 13
                        todo <- unique(todo)
                    }
                    optiout[[17]] <- list()
                }
            }
            # SBB6
            if (any(todo == 18)) {
                optiout[[18]] <- MLE_intern(cbind(data1, data2),
                                            start[[18]],
                                            18,
                                            weights = weights,
                                            se = se)
                if (optiout[[18]]$par[1] <= 1.1 | optiout[[18]]$par[2] <= 1.1) {
                    if (optiout[[18]]$par[1] <= 1.1) {
                        todo[todo == 18] <- 14
                        todo <- unique(todo)
                    } else if (optiout[[18]]$par[2] <= 1.1) {
                        todo[todo == 18] <- 16
                        todo <- unique(todo)
                    }
                    optiout[[18]] <- list()
                }
            }
            # SBB7
            if (any(todo == 19)) {
                optiout[[19]] <- MLE_intern(cbind(data1, data2),
                                            start[[19]],
                                            19,
                                            weights = weights,
                                            se = se)
                if (optiout[[19]]$par[1] <= 1.1 | optiout[[19]]$par[2] <= 0.1) {
                    if (optiout[[19]]$par[1] <= 1.1) {
                        todo[todo == 19] <- 13
                        todo <- unique(todo)
                    } else if (optiout[[19]]$par[2] <= 0.1) {
                        todo[todo == 19] <- 16
                        todo <- unique(todo)
                    }
                    optiout[[19]] <- list()
                }
            }
            # SBB8
            if (any(todo == 20)) {
                optiout[[20]] <- MLE_intern(cbind(data1, data2),
                                            start[[20]],
                                            20,
                                            weights = weights,
                                            se = se)
                if (optiout[[20]]$par[2] >= 0.99) {
                    todo[todo == 20] <- 16
                    todo <- unique(todo)
                    optiout[[20]] <- list()
                }
            }
            # BB1_90
            if (any(todo == 27)) {
                optiout[[27]] <- MLE_intern(cbind(data1, data2),
                                            start[[27]],
                                            27,
                                            weights = weights,
                                            se = se)
                if (optiout[[27]]$par[1] >= -0.1 | optiout[[27]]$par[2] >= -1.1) {
                    if (optiout[[27]]$par[1] >= -0.1) {
                        todo[todo == 27] <- 24
                        todo <- unique(todo)
                    } else if (optiout[[27]]$par[2] >= -1.1) {
                        todo[todo == 27] <- 23
                        todo <- unique(todo)
                    }
                    optiout[[27]] <- list()
                }
            }
            # BB6_90
            if (any(todo == 28)) {
                optiout[[28]] <- MLE_intern(cbind(data1, data2),
                                            start[[28]],
                                            28,
                                            weights = weights,
                                            se = se)
                if (optiout[[28]]$par[1] >= -1.1 | optiout[[28]]$par[2] >= -1.1) {
                    if (optiout[[28]]$par[1] >= -1.1) {
                        todo[todo == 28] <- 24
                        todo <- unique(todo)
                    } else if (optiout[[28]]$par[2] >= -1.1) {
                        todo[todo == 28] <- 26
                        todo <- unique(todo)
                    }
                    optiout[[28]] <- list()
                }
            }
            # BB7_90
            if (any(todo == 29)) {
                optiout[[29]] <- MLE_intern(cbind(data1, data2),
                                            start[[29]],
                                            29,
                                            weights = weights,
                                            se = se)
                if (optiout[[29]]$par[1] >= -1.1 | optiout[[29]]$par[2] >= -0.1) {
                    if (optiout[[29]]$par[1] >= -1.1) {
                        todo[todo == 29] <- 23
                        todo <- unique(todo)
                    } else if (optiout[[29]]$par[2] >= -0.1) {
                        todo[todo == 29] <- 26
                        todo <- unique(todo)
                    }
                    optiout[[29]] <- list()
                }
            }

            if (any(todo == 30)) {
                # BB8_90
                optiout[[30]] <- MLE_intern(cbind(data1, data2),
                                            start[[30]],
                                            30,
                                            weights = weights,
                                            se = se)
                if (optiout[[30]]$par[2] <= -0.99) {
                    todo[todo == 30] <- 26
                    todo <- unique(todo)
                    optiout[[30]] <- list()
                }
            }
            # BB1_270
            if (any(todo == 37)) {
                optiout[[37]] <- MLE_intern(cbind(data1, data2),
                                            start[[37]],
                                            37,
                                            weights = weights,
                                            se = se)
                if (optiout[[37]]$par[1] >= -0.1 | optiout[[37]]$par[2] >= -1.1) {
                    if (optiout[[37]]$par[1] >= -0.1) {
                        todo[todo == 37] <- 34
                        todo <- unique(todo)
                    } else if (optiout[[37]]$par[2] >= -1.1) {
                        todo[todo == 37] <- 33
                        todo <- unique(todo)
                    }
                    optiout[[37]] <- list()
                }
            }
            # BB6_270
            if (any(todo == 38)) {
                optiout[[38]] <- MLE_intern(cbind(data1, data2),
                                            start[[38]],
                                            38,
                                            weights = weights,
                                            se = se)
                if (optiout[[38]]$par[1] >= -1.1 | optiout[[38]]$par[2] >= -1.1) {
                    if (optiout[[38]]$par[1] >= -1.1) {
                        todo[todo == 38] <- 34
                        todo <- unique(todo)
                    } else if (optiout[[38]]$par[2] >= -1.1) {
                        todo[todo == 38] <- 36
                        todo <- unique(todo)
                    }
                    optiout[[38]] <- list()
                }
            }
            # BB7_270
            if (any(todo == 39)) {
                optiout[[39]] <- MLE_intern(cbind(data1, data2),
                                            start[[39]],
                                            39,
                                            weights = weights,
                                            se = se)
                if (optiout[[39]]$par[1] >= -1.1 | optiout[[39]]$par[2] >= -0.1) {
                    if (optiout[[39]]$par[1] >= -1.1) {
                        todo[todo == 39] <- 33
                        todo <- unique(todo)
                    } else if (optiout[[39]]$par[2] >= -0.1) {
                        todo[todo == 39] <- 36
                        todo <- unique(todo)
                    }
                    optiout[[39]] <- list()
                }
            }
            # BB8_270
            if (any(todo == 40)) {
                optiout[[40]] <- MLE_intern(cbind(data1, data2),
                                            start[[40]],
                                            40,
                                            weights = weights,
                                            se = se)
                if (optiout[[40]]$par[2] <= -0.99) {
                    todo[todo == 40] <- 36
                    todo <- unique(todo)
                    optiout[[40]] <- list()
                }
            }
            # Tawns
            for (i in todo[(todo %in% allfams[tawns])]) {
                optiout[[i]] <- MLE_intern_Tawn(cbind(data1, data2),
                                                start[[i]],
                                                i,
                                                se = se)
            }
            # one parameter families
            for (i in todo[todo %in% allfams[onepar]]) {
                optiout[[i]] <- MLE_intern(cbind(data1, data2),
                                           start[[i]],
                                           i,
                                           weights = weights,
                                           se = se)
            }

            ## calculate AIC and BIC
            AICs <- rep(Inf, max(todo))
            BICs <- rep(Inf, max(todo))
            lls  <- rep(Inf, max(todo))
            for (i in todo) {
                if (i %in% allfams[twopar]) {
                    if (any(is.na(weights))) {
                        lls[i] <- sum(log(BiCopPDF(data1,
                                                   data2,
                                                   i,
                                                   optiout[[i]]$par[1],
                                                   optiout[[i]]$par[2],
                                                   check.pars = FALSE)))
                    } else {
                        lls[i] <- sum(log(BiCopPDF(data1,
                                                   data2,
                                                   i,
                                                   optiout[[i]]$par[1],
                                                   optiout[[i]]$par[2],
                                                   check.pars = FALSE)) %*% weights)
                    }
                    AICs[i] <- -2 * lls[i] + 4
                    BICs[i] <- -2 * lls[i] + 2 * log(length(data1))
                } else {
                    if (any(is.na(weights))) {
                        lls[i] <- sum(log(BiCopPDF(data1,
                                                   data2,
                                                   i,
                                                   optiout[[i]]$par,
                                                   check.pars = FALSE)))
                    } else {
                        lls[i] <- sum(log(BiCopPDF(data1,
                                                   data2,
                                                   i,
                                                   optiout[[i]]$par,
                                                   check.pars = FALSE)) %*% weights)
                    }
                    AICs[i] <- -2 * lls[i] + 2
                    BICs[i] <- -2 * lls[i] + 2 * log(length(data1))
                }
            }

            ## select the best fitting model
            if (selectioncrit == "AIC") {
                out$family <- todo[which.min(AICs[todo])][1]
            } else {
                out$family <- todo[which.min(BICs[todo])][1]
            }

            ## for one-parameter families, set par2 = 0 (default)
            out$par <- optiout[[out$family]]$par[1]
            if (out$family %in% allfams[onepar]) {
                out$par2 <- 0
            } else {
                out$par2 <- optiout[[out$family]]$par[2]
            }
        }
    }

    ## store results in BiCop object (dependence measures are calculated)
    out <- BiCop(out$family, out$par, out$par2, check.pars = FALSE)

    ## add more information about the fit
    if (out$family == 0) {
        if (se)
            out$se <- NA
        out$nobs   <- length(u1)
        out$logLik <- 0
        out$AIC    <- 0
        out$BIC    <- 0
    } else {
        if (se) {
            out$se <- optiout[[out$family]]$se[1]
            if (out$family %in% allfams[twopar])
                out$se2 <- optiout[[out$family]]$se[2]
        }
        out$nobs   <- length(u1)
        out$logLik <- lls[out$family]
        out$AIC    <- AICs[out$family]
        out$BIC    <- BICs[out$family]
    }

    ## store the call that created the BiCop object
    out$call <- match.call()

    ## return final BiCop objectz
    out
}


##### ----------------------------------------------------------------------
## function for augmenting a familyset with rotations
with_rotations <- function(nums) {
    unique(unlist(lapply(nums, get_rotations)))
}

get_rotations <- function(i) {
    # no roations for independence, gaussian, student and frank copulas
    out <- i

    ## rotations for other families
    if(i %in% c(3, 13, 23, 33)) out <- c(3, 13, 23, 33)
    if(i %in% c(4, 14, 24, 34)) out <- c(4, 14, 24, 34)
    if(i %in% c(6, 16, 26, 36)) out <- c(6, 16, 26, 36)
    if(i %in% c(7, 17, 27, 37)) out <- c(7, 17, 27, 37)
    if(i %in% c(8, 18, 28, 38)) out <- c(8, 18, 28, 38)
    if(i %in% c(9, 19, 29, 39)) out <- c(9, 19, 29, 39)
    if(i %in% c(10, 20, 30, 40)) out <- c(10, 20, 30, 40)
    if(i %in% c(104, 114, 124, 134)) out <- c(104, 114, 124, 134)
    if(i %in% c(204, 214, 224, 234)) out <- c(204, 214, 224, 234)

    out
}
