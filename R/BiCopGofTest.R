#' Goodness-of-Fit Test for Bivariate Copulas
#' 
#' This function performs a goodness-of-fit test for bivariate copulas, either
#' based on White's information matrix equality (White 1982) as introduced by
#' Huang and Prokhorov (2011) or based on Kendall's process. It computes the
#' test statistics and p-values.
#' 
#' \code{method = "white"}:\cr This goodness-of fit test uses the information
#' matrix equality of White (1982) and was investigated by Huang and Prokhorov
#' (2011). The main contribution is that under correct model specification the
#' Fisher Information can be equivalently calculated as minus the expected
#' Hessian matrix or as the expected outer product of the score function. The
#' null hypothesis is \deqn{ H_0: \boldsymbol{H}(\theta) +
#' \boldsymbol{C}(\theta) = 0 } against the alternative \deqn{ H_0:
#' \boldsymbol{H}(\theta) + \boldsymbol{C}(\theta) \neq 0 , } where
#' \eqn{\boldsymbol{H}(\theta)} is the expected Hessian matrix and
#' \eqn{\boldsymbol{C}(\theta)} is the expected outer product of the score
#' function. For the calculation of the test statistic we use the consistent
#' maximum likelihood estimator \eqn{\hat{\theta}} and the sample counter parts
#' of \eqn{\boldsymbol{H}(\theta)} and \eqn{\boldsymbol{C}(\theta)}. The
#' correction of the covariance-matrix in the test statistic for the
#' uncertainty in the margins is skipped. The implemented tests assumes that
#' where is no uncertainty in the margins. The correction can be found in Huang
#' and Prokhorov (2011). It involves two-dimensional integrals.\cr WARNING: For
#' the t-copula the test may be instable. The results for the t-copula
#' therefore have to be treated carefully.\cr \cr \code{method = "kendall"}:\cr
#' This copula goodness-of-fit test is based on Kendall's process as
#' investigated by Genest and Rivest (1993) and Wang and Wells (2000). For
#' rotated copulas the input arguments are transformed and the goodness-of-fit
#' procedure for the corresponding non-rotated copula is used.
#' 
#' @param u1,u2 Numeric vectors of equal length with values in [0,1].
#' @param family An integer defining the bivariate copula family: \cr \code{0}
#' = independence copula \cr \code{1} = Gaussian copula \cr \code{2} = Student
#' t copula (t-copula) (only for \code{method = "white"}; see details)\cr
#' \code{3} = Clayton copula \cr \code{4} = Gumbel copula \cr \code{5} = Frank
#' copula \cr \code{6} = Joe copula (only for \code{method = "kendall"}) \cr
#' \code{7} = BB1 copula (only for \code{method = "kendall"})\cr \code{8} = BB6
#' copula (only for \code{method = "kendall"})\cr \code{9} = BB7 copula (only
#' for \code{method = "kendall"})\cr \code{10} = BB8 copula (only for
#' \code{method ="kendall"})\cr \code{13} = rotated Clayton copula (180
#' degrees; ``survival Clayton'') \cr \code{14} = rotated Gumbel copula (180
#' degrees; ``survival Gumbel'') \cr \code{16} = rotated Joe copula (180
#' degrees; ``survival Joe'') \cr \code{17} = rotated BB1 copula (180 degrees;
#' ``survival BB1''; only for \code{method = "kendall"})\cr \code{18} = rotated
#' BB6 copula (180 degrees; ``survival BB6''; only for \code{method =
#' "kendall"})\cr \code{19} = rotated BB7 copula (180 degrees; ``survival
#' BB7''; only for \code{method = "kendall"})\cr \code{20} = rotated BB8 copula
#' (180 degrees; ``survival BB8''; only for \code{method = "kendall"})\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr \code{24} = rotated
#' Gumbel copula (90 degrees) \cr \code{26} = rotated Joe copula (90 degrees)
#' \cr \code{27} = rotated BB1 copula (90 degrees; only for \code{method =
#' "kendall"}) \cr \code{28} = rotated BB6 copula (90 degrees; only for
#' \code{method = "kendall"}) \cr \code{29} = rotated BB7 copula (90 degrees;
#' only for \code{method = "kendall"}) \cr \code{30} = rotated BB8 copula (90
#' degrees; only for \code{method = "kendall"}) \cr \code{33} = rotated Clayton
#' copula (270 degrees) \cr \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr \code{37} = rotated BB1
#' copula (270 degrees; only for \code{method = "kendall"}) \cr \code{38} =
#' rotated BB6 copula (270 degrees; only for \code{method = "kendall"}) \cr
#' \code{39} = rotated BB7 copula (270 degrees; only for \code{method =
#' "kendall"}) \cr \code{40} = rotated BB8 copula (270 degrees; only for
#' \code{method = "kendall"})
#' @param par Copula parameter (optional).
#' @param par2 Second parameter for bivariate t-copula (optional); default:
#' \code{par2 = 0}.
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: \code{max.df = 30}).
#' @param method A string indicating the goodness-of-fit method:\cr
#' \code{"white"} = goodness-of-fit test based on White's information matrix
#' equality (default) \cr \code{"kendall"} = goodness-of-fit test based on
#' Kendall's process
#' @param B Integer; number of bootstrap samples (default: \code{B = 100}).
#' For \code{B = 0} only the the test statistics are returned.\cr WARNING: If
#' \code{B} is chosen too large, computations will take very long.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @return For \code{method = "white"}: \item{p.value}{Asymptotic p-value.}
#' \item{statistic}{The observed test statistic.}\cr For \code{method =
#' "kendall"} \item{p.value.CvM}{Bootstrapped p-value of the goodness-of-fit
#' test using the Cramer-von Mises statistic (if \code{B > 0}).}
#' \item{p.value.KS}{Bootstrapped p-value of the goodness-of-fit test using the
#' Kolmogorov-Smirnov statistic (if \code{B > 0}).} \item{statistic.CvM}{The
#' observed Cramer-von Mises test statistic.} \item{statistic.KS}{The observed
#' Kolmogorov-Smirnov test statistic.}
#' @author Ulf Schepsmeier, Wanling Huang, Jiying Luo, Eike Brechmann
#' @seealso \code{\link{BiCopDeriv2}}, \code{\link{BiCopDeriv}},
#' \code{\link{BiCopIndTest}}, \code{\link{BiCopVuongClarke}}
#' @references Genest, C. and L.-P. Rivest (1993). Statistical inference
#' procedures for bivariate Archimedean copulas. Journal of the American
#' Statistical Association, 88 (423), 1034-1043.
#' 
#' Huang, w. and A. Prokhorov (2011). A goodness-of-fit test for copulas. to
#' appear in Econometric Reviews
#' 
#' Luo J. (2011). Stepwise estimation of D-vines with arbitrary specified
#' copula pairs and EDA tools. Diploma thesis, Technische Universitaet
#' Muenchen.\cr \url{http://mediatum.ub.tum.de/?id=1079291}.
#' 
#' Wang, W. and M. T. Wells (2000). Model selection and semiparametric
#' inference for bivariate failure-time data. Journal of the American
#' Statistical Association, 95 (449), 62-72.
#' 
#' White, H. (1982) Maximum likelihood estimation of misspecified models,
#' Econometrica, 50, 1-26.
#' @examples
#' 
#' # simulate from a bivariate Clayton copula
#' set.seed(123)
#' simdata <- BiCopSim(300, 3, 2)
#' u1 <- simdata[,1]
#' u2 <- simdata[,2]
#' 
#' # perform White's goodness-of-fit test for the true copula
#' BiCopGofTest(u1, u2, family = 3)
#' 
#' # perform White's goodness-of-fit test for the Frank copula
#' BiCopGofTest(u1, u2, family = 5)
#' 
#' # perform Kendall's goodness-of-fit test for the true copula
#' gof <- BiCopGofTest(u1, u2, family = 3, method = "kendall", B=50)
#' gof$p.value.CvM
#' gof$p.value.KS
#' 
#' # perform Kendall's goodness-of-fit test for the Frank copula
#' gof <- BiCopGofTest(u1, u2, family = 5, method = "kendall", B=50)
#' gof$p.value.CvM
#' gof$p.value.KS
#' 
#' @export BiCopGofTest
BiCopGofTest <- function(u1, u2, family, par = 0, par2 = 0, method = "white", max.df = 30, 
                         B = 100, obj = NULL) {
    if (method == "White") 
        method <- "white"
    if (method == "Kendall") 
        method <- "kendall"
    
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
        stop("Data has be in the interval [0,1].")
    
    ## extract family and parameters if BiCop object is provided
    if (missing(family))
        family <- NA
    # for short hand usage extract obj from family
    if (class(family) == "BiCop")
        obj <- family
    if (!is.null(obj)) {
        stopifnot(class(obj) == "BiCop")
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }
    
    ## sanity checks for family and parameters
    if (is.na(family)) 
        stop("Provide either 'family' and 'par' or 'obj'")
    if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 
                        20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 43, 44))) 
        stop("Copula family not implemented.")
    if (c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) %in% family && par2 == 0) 
        stop("For t-, BB1, BB6, BB7 and BB8 copulas, 'par2' must be set.")
    if (c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36) %in% family && length(par) < 1) 
        stop("'par' not set.")
    if (par != 0) 
        BiCopCheck(family, par, par2)
    if (family == 2 && method == "kendall") 
        stop("The goodness-of-fit test based on Kendall's process is not implemented for the t-copula.")
    if (family %in% c(7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) && 
            method == "white") 
        stop("The goodness-of-fit test based on White's information matrix equality is not implemented for the BB copulas.")
    # if((level < 0 || level > 1) && method=='kendall') stop('Significance level has
    # to be between 0 and 1.')
    
    T <- length(u1)
    
    if (method == "white") {
        # Step 1: maximum likelihood estimation
        if (family == 2) {
            if (par == 0) {
                pars <- BiCopEst(u1, u2, family = family, method = "mle", max.df = max.df)
                theta <- pars$par
                nu <- pars$par2
                print(theta)
                print(nu)
            } else {
                theta <- par
                nu <- par2
            }
        } else {
            nu <- 0
            theta <- BiCopEst(u1, u2, family = family, method = "mle")$par
        }
        
        # Step 2: Calculation of Hesse and gradient and D
        
        if (family == 2) {
            Dprime <- matrix(0, 3, T)
            Vt <- array(0, dim = c(3, 3, T))
            grad <- c(0, 0)
            for (t in 1:T) {
                rho_teil <- f_rho(u1[t], u2[t], theta, nu)
                nu_teil <- f_nu(u1[t], u2[t], theta, nu)
                rho_nu_teil <- f_rho_nu(u1[t], u2[t], theta, nu)
                H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, nu_teil), 2, 2)  # Hesse matrix
                Hprime <- as.vector(H[lower.tri(H, diag = TRUE)])
                grad[1] <- BiCopDeriv(u1[t],
                                      u2[t], 
                                      family = family,
                                      par = theta, 
                                      par2 = nu,
                                      deriv = "par",
                                      log = TRUE)
                grad[2] <- BiCopDeriv(u1[t], 
                                      u2[t],
                                      family = family,
                                      par = theta, 
                                      par2 = nu,
                                      deriv = "par2",
                                      log = TRUE)
                C <- grad %*% t(grad)
                Cprime <- as.vector(C[lower.tri(C, diag = TRUE)])
                Dprime[, t] <- Hprime + Cprime
                
                Vt[, , t] <- Dprime %*% t(Dprime)
            }
            D <- apply(Dprime, 1, mean)
            V0 <- apply(Vt, c(1, 2), mean)
        } else {
            d <- rep(0, T)
            for (t in 1:T) {
                b <- BiCopPDF(u1[t],
                              u2[t], 
                              family, 
                              theta,
                              nu)
                d[t] <- BiCopDeriv2(u1[t],
                                    u2[t],
                                    family = family,
                                    par = theta, 
                                    par2 = nu, 
                                    deriv = "par")/b
            }
            D <- mean(d)
            Vt <- d^2
            V0 <- mean(Vt)
        }
        
        
        # Teststatistik und p-Wert
        if (family == 2) {
            handle <- try(solve(V0), TRUE)
            if (is.null(dim(handle))) 
                handle <- ginv(V0)
            test <- T * (t(D) %*% handle %*% D)
            pvalue <- 1 - pchisq(test, df = length(D))
        } else {
            test <- T * D * solve(V0) * D
            pvalue <- 1 - pchisq(test, df = 1)
        }
        out <- list(p.value = pvalue, statistic = test)
    } else if (method == "IR") {
        # Information ratio GOF Step 1: maximum likelihood estimation
        if (family == 2) {
            if (par == 0) {
                pars <- BiCopEst(u1, 
                                 u2, 
                                 family = family,
                                 method = "mle",
                                 max.df = max.df)
                theta <- pars$par
                nu <- pars$par2
                print(theta)
                print(nu)
            } else {
                theta <- par
                nu <- par2
            }
        } else {
            nu <- 0
            theta <- BiCopEst(u1, 
                              u2, 
                              family = family,
                              method = "mle")$par
        }
        
        # Step 2: Calculation of Hesse and gradient
        if (family == 2) {
            grad <- c(0, 0)
            rho_teil <- f_rho(u1, u2, theta, nu)
            nu_teil <- f_nu(u1, u2, theta, nu)
            rho_nu_teil <- f_rho_nu(u1, u2, theta, nu)
            H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, nu_teil), 2, 2)  # Hesse matrix
            grad[1] <- BiCopDeriv(u1,
                                  u2, 
                                  family = family, 
                                  par = theta, 
                                  par2 = nu, 
                                  deriv = "par", 
                                  log = TRUE)
            grad[2] <- BiCopDeriv(u1, 
                                  u2, 
                                  family = family,
                                  par = theta, 
                                  par2 = nu, 
                                  deriv = "par2",
                                  log = TRUE)
            C <- grad %*% t(grad)
        } else {
            d <- rep(0, T)
            for (t in 1:T) {
                b <- BiCopPDF(u1[t], 
                              u2[t],
                              family,
                              theta,
                              nu)
                d[t] <- BiCopDeriv2(u1[t],
                                    u2[t],
                                    family = family, 
                                    par = theta,
                                    par2 = nu, 
                                    deriv = "par")/b
            }
            H <- mean(d)
            C <- BiCopDeriv(u1, 
                            u2, 
                            family = family,
                            par = theta,
                            par2 = nu,
                            deriv = "par", 
                            log = TRUE)
        }
        Phi <- -solve(H) %*% C
        IR <- trace(Phi)/dim(H)[1]  #Zwischenergebnis
        
        # Bootstrap procedure
        if (B == 0) {
            out <- list(IR = IR, p.value = NULL)
        } else {
            IR_boot <- boot.IR(family, theta, nu, B, length(u1))
            sigma2 <- var(IR_boot)
            IR_new <- ((IR - 1)/sqrt(sigma2))^2
            IR_boot <- ((IR_boot - 1)/sqrt(sigma2))^2
            p.value <- mean(IR_boot >= IR_new)
            
            out <- list(IR = IR, p.value = p.value)
        }
        
    } else if (method == "kendall") {
        if (family %in% c(13, 14, 16, 17, 18, 19, 20)) {
            u1 <- 1 - u1
            u2 <- 1 - u2
            family <- family - 10
        } else if (family %in% c(23, 24, 26, 27, 28, 29, 30)) {
            u1 <- 1 - u1
            family <- family - 20
        } else if (family %in% c(33, 34, 36, 37, 38, 39, 40)) {
            u2 <- 1 - u2
            family <- family - 30
        }
        
        ostat <- obs.stat(u1, u2, family)
        
        if (B == 0) {
            # no bootstrap
            
            sn.obs <- ostat$Sn
            tn.obs <- ostat$Tn
            out <- list(Sn = sn.obs, Tn = tn.obs)
            
        } else {
            bstat <- list()
            for (i in 1:B) bstat[[i]] <- boot.stat(u1, u2, family)
            
            sn.boot <- rep(0, B)
            tn.boot <- rep(0, B)
            for (i in 1:B) {
                sn.boot[i] <- bstat[[i]]$sn
                tn.boot[i] <- bstat[[i]]$tn
            }
            
            sn.obs <- ostat$Sn
            tn.obs <- ostat$Tn
            
            # k<-as.integer((1-level)*B) sn.critical <- sn.boot[k] # critical value of test
            # at level 0.05 tn.critical <- tn.boot[k] # critical value of test at level 0.05
            
            
            pv.sn <- sapply(sn.obs,
                            function(x) (1/B) * length(which(sn.boot[1:B] >= x)))  # P-value of Sn
            pv.tn <- sapply(tn.obs,
                            function(x) (1/B) * length(which(tn.boot[1:B] >= x)))  # P-value of Tn
            
            out <- list(p.value.CvM = pv.sn, 
                        p.value.KS = pv.tn,
                        statistic.CvM = sn.obs, 
                        statistic.KS = tn.obs)
            
        }
    } else {
        stop("Method not implemented")
    }
    
    
    
    return(out)
}


f_rho <- function(u1, u2, par, par2) {
    a <- .C("diff2lPDF_rho_tCopula",
            as.double(u1),
            as.double(u2),
            as.integer(length(u1)), 
            as.double(c(par, par2)),
            as.integer(2),
            as.double(rep(0, length(u1))),
            PACKAGE = "VineCopula")[[6]]
    
    return(sum(a))
}

f_nu <- function(u1, u2, par, par2) {
    a <- .C("diff2lPDF_nu_tCopula_new",
            as.double(u1),
            as.double(u2),
            as.integer(length(u1)), 
            as.double(c(par, par2)),
            as.integer(2),
            as.double(rep(0, length(u1))),
            PACKAGE = "VineCopula")[[6]]
    
    return(sum(a))
}

f_rho_nu <- function(u1, u2, par, par2) {
    a <- .C("diff2lPDF_rho_nu_tCopula_new",
            as.double(u1),
            as.double(u2),
            as.integer(length(u1)), 
            as.double(c(par, par2)), 
            as.integer(2), 
            as.double(rep(0, length(u1))),
            PACKAGE = "VineCopula")[[6]]
    
    return(sum(a))
}

boot.stat <- function(u, v, fam) {
    
    n <- length(u)
    t <- seq(1, n)/(n + 1e-04)
    kt <- rep(0, n)
    
    # estimate paramemter for different copula family from (u,v)
    param <- suppressWarnings({
        BiCopEst(u, v, family = fam)
    })
    # calulate k(t) and kn(t) of bootstrap sample data
    if (fam == 1) {
        # normal
        sam <- BiCopSim(n, 1, param$par, param$par2)  # generate data for the simulation of K(t)
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # parameter estimation of sample data
        sim <- BiCopSim(10000, 1, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
        cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
        dcop <- rep(0, 10000)
        for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]),
                                                        qnorm(sim[i, 2])),
                                              corr = cormat)
        kt <- sapply(t,
                     function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
    }
    if (fam == 2) {
        # t
        sam <- BiCopSim(n, fam, param$par, param$par2)  # generate data for the simulation of K(t)
        sam.par <- suppressWarnings({
            BiCopEst(sam[, 1], sam[, 2], family = fam)
        })  # parameter estimation of sample data
        sim <- BiCopSim(10000, fam, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
        
        # par2 muss auf einen Integer gesetzt werden f?r mvtnorm
        param$par2 <- round(param$par2)
        
        cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
        dcop <- rep(0, 10000)
        for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2), 
                                                     qt(sim[i, 2], df = param$par2)),
                                           corr = cormat,
                                           df = param$par2)
        kt <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
    } else if (fam == 3) {
        # Clayton
        sam <- BiCopSim(n, 3, param$par)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
        kt <- t + t * (1 - t^sam.par)/sam.par
    } else if (fam == 4) {
        # gumbel
        sam <- BiCopSim(n, 4, param$par)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
        kt <- t - t * log(t)/(sam.par)
    } else if (fam == 5) {
        # frank
        sam <- BiCopSim(n, 5, param$par)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
        kt <- t + log((1 - exp(-sam.par))/(1 - exp(-sam.par * t))) * (1 - exp(-sam.par * t))/(sam.par * exp(-sam.par * t))
    } else if (fam == 6) {
        sam <- BiCopSim(n, 6, param$par)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data        
        kt <- t - (log(1 - (1 - t)^sam.par) * (1 - (1 - t))^sam.par)/(sam.par * (1 - t)^(sam.par - 1))
    } else if (fam == 7) {
        # BB1
        sam <- BiCopSim(n, 7, param$par, param$par2)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
        theta <- sam.par$par
        delta <- sam.par$par2
        kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
    } else if (fam == 8) {
        # BB6
        sam <- BiCopSim(n, 8, param$par, param$par2)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
        theta <- sam.par$par
        delta <- sam.par$par2
        kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
    } else if (fam == 9) {
        # BB7
        sam <- BiCopSim(n, 9, param$par, param$par2)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
        theta <- sam.par$par
        delta <- sam.par$par2
        kt <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
    } else if (fam == 10) {
        # BB8
        sam <- BiCopSim(n, 10, param$par, param$par2)  # generate sample data
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
        theta <- sam.par$par
        delta <- sam.par$par2
        kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
    }
    
    # calculate emp. Kn
    w <- rep(0, n)
    w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > sam[, 1] & y > sam[, 2])),
                     sam[, 1],
                     sam[, 2])
    w <- sort(w)
    kn <- rep(0, n)
    kn <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
    
    # calculate test statistic Sn
    Sn1 <- 0
    Sn2 <- 0
    for (j in 1:(n - 1)) {
        Sn1 <- Sn1 + ((kn[j])^2 * (kt[j + 1] - kt[j]))
        Sn2 <- Sn2 + (kn[j]) * ((kt[j + 1])^2 - (kt[j])^2)
    }
    sn <- n/3 + n * Sn1 - n * Sn2
    # calculation of test statistics Tn
    tm <- matrix(0, n - 1, 2)
    # mit i=0
    for (j in 1:(n - 1)) {
        tm[j, 1] <- abs(kn[j] - kt[j])
    }
    # mit i=1
    for (j in 1:(n - 1)) {
        tm[j, 2] <- abs(kn[j] - kt[j + 1])
    }
    tn <- max(tm) * sqrt(n)
    
    sn <- sort(sn)  # vector of ordered statistic Sn
    tn <- sort(tn)  # vector of ordered statistic Tn
    out <- list(sn = sn, tn = tn)
}



obs.stat <- function(u, v, fam) {
    
    n <- length(u)
    t <- seq(1, n)/(n + 1e-04)
    kt <- rep(0, n)
    
    # estimate paramemter for different copula family from (u,v)
    param <- suppressWarnings({
        BiCopEst(u, v, family = fam)
    })
    
    # calculate observed K(t) of (u,v)
    kt.obs <- rep(0, n)
    if (fam == 1) {
        sim <- BiCopSim(10000, 1, param$par)  # generate data for the simulation of K(t)
        cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
        dcop <- rep(0, 10000)
        for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]), 
                                                        qnorm(sim[i, 2])),
                                              corr = cormat)
        kt.obs <- sapply(t, 
                         function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
    } else if (fam == 2) {
        sim <- BiCopSim(10000, 2, param$par, param$par2)  # generate data for the simulation of K(t)
        cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
        dcop <- rep(0, 10000)
        for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2), 
                                                     qt(sim[i, 2], df = param$par2)), 
                                           corr = cormat,
                                           df = param$par2)
        kt.obs <- sapply(t, 
                         function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
    } else if (fam == 3) {
        kt.obs <- t + t * (1 - t^param$par)/param$par
    } else if (fam == 4) {
        kt.obs <- t - t * log(t)/(param$par)
    } else if (fam == 5) {
        kt.obs <- t + log((1 - exp(-param$par))/(1 - exp(-param$par * t))) * (1 - exp(-param$par * t))/(param$par * exp(-param$par * t))
    } else if (fam == 6) {
        kt.obs <- t - (log(1 - (1 - t)^param$par) * (1 - (1 - t))^param$par)/(param$par * (1 - t)^(param$par - 1))
    } else if (fam == 7) {
        theta <- param$par
        delta <- param$par2
        kt.obs <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
    } else if (fam == 8) {
        theta <- param$par
        delta <- param$par2
        kt.obs <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
    } else if (fam == 9) {
        theta <- param$par
        delta <- param$par2
        kt.obs <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
    } else if (fam == 10) {
        theta <- param$par
        delta <- param$par2
        kt.obs <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
    }
    # calculation of observed Kn
    w <- rep(0, n)
    w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > u & y > v)), u, v)
    
    w <- sort(w)
    kn.obs <- rep(0, n)
    kn.obs <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
    
    # calculation of observed value Sn
    Sn1 <- 0
    Sn2 <- 0
    for (j in 1:(n - 1)) {
        Sn1 <- Sn1 + ((kn.obs[j])^2 * (kt.obs[j + 1] - kt.obs[j]))
        Sn2 <- Sn2 + (kn.obs[j]) * ((kt.obs[j + 1])^2 - (kt.obs[j])^2)
    }
    Sn <- n/3 + n * Sn1 - n * Sn2  # observed Sn
    
    # calculation of observed Tn
    tn.obs <- matrix(0, n - 1, 2)
    # mit i=0
    for (j in 1:(n - 1)) {
        tn.obs[j, 1] <- abs(kn.obs[j] - kt.obs[j])
    }
    # mit i=1
    for (j in 1:(n - 1)) {
        tn.obs[j, 2] <- abs(kn.obs[j] - kt.obs[j + 1])
    }
    Tn <- max(tn.obs) * sqrt(n)
    out <- list(Sn = Sn, Tn = Tn)
    return(out)
}


############################ 

# bootstrap for IR

boot.IR <- function(family, theta, nu, B, n) {
    # theta und nu sind die geschaetzten Parameter
    IR <- rep(0, B)
    for (i in 1:B) {
        sam <- BiCopSim(n, family, theta, nu)
        sam.par <- BiCopEst(sam[, 1], sam[, 2], family = family)  # parameter estimation of sample data
        if (family == 2) {
            theta2 <- sam.par[1]
            nu2 <- sam.par[2]
            grad <- c(0, 0)
            rho_teil <- f_rho(sam[, 1], sam[, 2], theta2, nu2)
            nu_teil <- f_nu(sam[, 1], sam[, 2], theta2, nu2)
            rho_nu_teil <- f_rho_nu(sam[, 1], sam[, 2], theta2, nu2)
            H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, nu_teil), 2, 2)  # Hesse matrix
            grad[1] <- BiCopDeriv(sam[, 1], 
                                  sam[, 2],
                                  family = family,
                                  par = theta2, 
                                  par2 = nu2, 
                                  deriv = "par",
                                  log = TRUE)
            grad[2] <- BiCopDeriv(sam[, 1],
                                  sam[, 2],
                                  family = family,
                                  par = theta2, 
                                  par2 = nu2,
                                  deriv = "par2", 
                                  log = TRUE)
            C <- grad %*% t(grad)
        } else {
            theta2 <- sam.par
            nu2 <- 0
            d <- rep(0, T)
            for (t in 1:T) {
                b <- BiCopPDF(sam[t, 1], sam[t, 2], family, theta, nu)
                d[t] <- BiCopDeriv2(sam[t, 1],
                                    sam[t, 2],
                                    family = family, 
                                    par = theta, 
                                    par2 = nu,
                                    deriv = "par")/b
            }
            H <- mean(d)
            C <- BiCopDeriv(sam[, 1],
                            sam[, 2], 
                            family = family,
                            par = theta2, 
                            par2 = nu2, 
                            deriv = "par", 
                            log = TRUE)
        }
        Phi <- -solve(H) %*% C
        IR[i] <- trace(Phi)/dim(H)[1]
    }
    
    return(IR)
}
