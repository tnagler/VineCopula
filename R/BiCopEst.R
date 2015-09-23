BiCopEst <- function(u1, u2, family, method = "mle", se = FALSE, max.df = 30, max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), 
                     weights = NA) {
    # Function that estimates the parameter(s) of the bivatiate copula
    #---------------------------------------------------------
    # INPUT:
    #   u1,u2      Data for which to estimate parameter
    #   family            The array definig the copulas in the pcc copula construction
    # OUTPUT:
    #   theta      Estimated Parameters
    #----------------------------------------------------------
    # Author: Carlos Almeida  <calmeida at ma.tum.de>
    # Update: Ulf Schepsmeier <schepsmeier at ma.tum.de>
    # Date: 2008-12-08
    # Update date: 2011-05-27
    # Version: 1.1
    #---------------------------------------------------------------
    
    # sanity checks
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (length(u1) < 2) 
        stop("Number of observations has to be at least 2.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20,
                        23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40,
                        41, 51, 61,  71, 104, 114, 124, 134, 204, 214, 224, 234))) 
        stop("Copula family not implemented.")
    
    if (max.df <= 2) 
        stop("The upper bound for the degrees of freedom parameter has to be larger than 2.")
    if (!is.list(max.BB)) 
        stop("'max.BB' has to be a list.")
    if (max.BB$BB1[1] < 0.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB1[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[2] < 0.001) 
        stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
        stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
    
    if (method != "mle" && method != "itau") 
        stop("Estimation method has to be either 'mle' or 'itau'.")
    if (method == "itau" && family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234)) {
        message("For two parameter copulas the estimation method 'itau' cannot be used. The method is automatically set to 'mle'.")
        method <- "mle"
    }
    if (is.logical(se) == FALSE) 
        stop("'se' has to be a logical variable (TRUE or FALSE).")
    
    
    ## calculate empirical kendall's tau
    if (family != 0) {
        # tau <- cor(u1,u2,method='kendall')
        tau <- fasttau(u1, u2)
    }
    
    ## inversion of kendall's tau
    theta <- 0
    if (family == 0) {
        # independent
        theta <- 0
    } else if (family == 1) {
        ## Gaussian
        theta <- sin(tau * pi/2)
    } else if (family == 3 || family == 13) {
        ## Clayton
        if (tau <= 0) {
            warning("Clayton copula cannot be used for negatively dependent data.")
            tau <- 0.05
        }
        theta <- max(0, 2 * tau/(1 - tau))
    } else if (family == 4 || family == 14) {
        ## Gumbel
        if (tau < 0) {
            warning("Gumbel copula cannot be used for negatively dependent data.")
            tau <- 0.05
        }
        theta <- max(1, 1/(1 - tau))
    } else if (family == 5) {
        ## Frank
        theta <- Frank.itau.JJ(tau)
    } else if (family == 6 || family == 16) {
        ## Joe
        if (tau <= 0) {
            warning("Joe copula cannot be used for negatively dependent data.")
            tau <- 0.05
        }
        theta <- Joe.itau.JJ(tau)
    } else if (family == 23 || family == 33) {
        if (tau >= 0) {
            warning("Rotated Clayton copula cannot be used for positively dependent data.")
            tau <- -0.05
        }
        theta <- (2 * tau/(1 + tau))
    } else if (family == 24 || family == 34) {
        if (tau > 0) {
            warning("Rotated Gumbel copula cannot be used for positively dependent data.")
            tau <- -0.05
        }
        theta <- -(1/(1 + tau))
    } else if (family == 26 || family == 36) {
        if (tau >= 0) {
            warning("Rotated Joe copula cannot be used for positively dependent data.")
            tau <- -0.05
        }
        theta <- -Joe.itau.JJ(-tau)
    } else if (family %in% c(41, 51)) {
        theta <- ipsA.tau2cpar(tau)
    } else if (family %in% c(61, 71)) {
        theta <- -ipsA.tau2cpar(-tau)
    }
    
    ## standard errors for method itau
    se1 <- 0
    if (method == "itau" && se == TRUE) {
        p <- 2
        n <- length(u1)
        ec <- numeric(n)
        u <- cbind(u1, u2)
        v <- matrix(0, n, p * (p - 1)/2)
        
        if (family == 1) 
            tauder <- function(x) 2/(pi * sqrt(1 - x^2)) else if (family %in% c(3, 13, 23, 33)) 
                tauder <- function(x) 2 * (2 + x)^(-2) else if (family %in% c(4, 14, 24, 34)) 
                    tauder <- function(x) x^(-2) else if (family == 5) {
                        tauder <- function(x) {
                            f <- function(x) x/(exp(x) - 1)
                            4/x^2 - 8/x^3 * integrate(f, lower = 0 + .Machine$double.eps^0.5, upper = x)$value + 4/(x * (exp(x) - 1))
                        }
                    } else if (family %in% c(6, 16, 26, 36)) {
                        tauder <- function(x) {
                            euler <- 0.577215664901533
                            -((-2 + 2 * euler + 2 * log(2) + digamma(1/x) + digamma(1/2 * (2 + x)/x) + x)/(-2 + x)^2) + ((-trigamma(1/x)/x^2 + trigamma(1/2 * (2 + 
                                                                                                                                                                   x)/x) * (1/(2 + x) - (2 + x)/(2 * x^2)) + 1)/(-2 + x))
                        }
                    } else if (family %in% c(41, 51, 61, 71)) {
                        tauder <- function(x) {
                            2 * sqrt(pi) * gamma(0.5 + x) * (digamma(1 + x) - digamma(0.5 + x))/gamma(1 + x)
                        }
                    }
        
        l <- 1
        for (j in 1:(p - 1)) {
            for (i in (j + 1):p) {
                for (k in 1:n) ec[k] <- sum(u[, i] <= u[k, i] & u[, j] <= u[k, j])/n
                v[, l] <- 2 * ec - u[, i] - u[, j]
                l <- l + 1
            }
        }
        
        if (family == 0) 
            D <- 0 else if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 41, 51)) 
                D <- 1/tauder(theta) else if (family %in% c(23, 33, 24, 34, 26, 36, 61, 71)) 
                    D <- 1/tauder(-theta)
        
        
        se1 <- as.numeric(sqrt(16/n * var(v %*% D)))
    }
    
    ## set starting parameters for maximum likelihood estimation
    if (method == "mle") {
        theta1 <- 0
        delta <- 0
        
        if (!(family %in% c(2, 6, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234))) {
            theta1 <- theta
        }
        if (family == 2) {
            ## t
            theta1 <- sin(tau * pi/2)
            delta1 <- min(10, (max.df + 2)/2)  # Take the middle between 2 and max.df
            delta <- MLE_intern(cbind(u1, u2),
                                c(theta1, delta1),
                                family = family,
                                se = FALSE,
                                max.df, 
                                max.BB, 
                                cor.fixed = TRUE,
                                weights)$par[2]
        } else if (family == 7 || family == 17) {
            ## BB1
            if (tau < 0) {
                print("The BB1 or survival BB1 copula cannot be used for negatively dependent data.")
                delta <- 1.001
                theta1 <- 0.001
            } else {
                delta <- min(1.5, max((max.BB$BB1[2] + 1.001)/2, 1.001))
                theta1 <- min(0.5, max((max.BB$BB1[1] + 0.001)/2, 0.001))
            }
        } else if (family == 27 || family == 37) {
            ## BB1
            if (tau > 0) {
                print("The rotated BB1 copulas cannot be used for positively dependent data.")
                delta <- -1.001
                theta1 <- -0.001
            } else {
                delta <- max(-1.5, -max((max.BB$BB1[2] + 1.001)/2, 1.001))
                theta1 <- max(-0.5, -max((max.BB$BB1[1] + 0.001)/2, 0.001))
            }
        } else if (family == 8 || family == 18) {
            ## BB6
            if (tau < 0) {
                print("The BB6 or survival BB6 copula cannot be used for negatively dependent data.")
                delta <- 1.001
                theta1 <- 1.001
            } else {
                delta <- min(1.5, max((max.BB$BB6[2] + 1.001)/2, 1.001))
                theta1 <- min(1.5, max((max.BB$BB6[1] + 1.001)/2, 1.001))
            }
        } else if (family == 28 || family == 38) {
            ## BB6
            if (tau > 0) {
                print("The rotated BB6 copulas cannot be used for positively dependent data.")
                delta <- -1.001
                theta1 <- -1.001
            } else {
                delta <- max(-1.5, -max((max.BB$BB6[2] + 1.001)/2, 1.001))
                theta1 <- max(-1.5, -max((max.BB$BB6[1] + 1.001)/2, 1.001))
            }
        } else if (family == 9 || family == 19) {
            ## BB7
            if (tau < 0) {
                print("The BB7 or survival BB7 copula cannot be used for negatively dependent data.")
                delta <- 0.001
                theta <- 1.001
            } else {
                delta <- min(0.5, max((max.BB$BB7[2] + 0.001)/2, 0.001))
                theta1 <- min(1.5, max((max.BB$BB7[1] + 1.001)/2, 1.001))
            }
        } else if (family == 29 || family == 39) {
            ## BB7
            if (tau > 0) {
                print("The rotated BB7 copulas cannot be used for positively dependent data.")
                delta <- -0.001
                theta1 <- -1.001
            } else {
                delta <- max(-0.5, -max((max.BB$BB7[2] + 0.001)/2, 0.001))
                theta1 <- max(-1.5, -max((max.BB$BB7[1] + 1.001)/2, 1.001))
            }
        } else if (family == 10 || family == 20) {
            ## BB8
            if (tau < 0) {
                print("The BB8 or survival BB8 copula cannot be used for negatively dependent data.")
                delta <- 0.001
                theta <- 1.001
            } else {
                delta <- min(0.5, max((max.BB$BB8[2] + 0.001)/2, 0.001))
                theta1 <- min(1.5, max((max.BB$BB8[1] + 1.001)/2, 1.001))
            }
        } else if (family == 30 || family == 40) {
            ## BB8
            if (tau > 0) {
                print("The rotated BB8 copulas cannot be used for positively dependent data.")
                delta <- -0.001
                theta1 <- -1.001
            } else {
                delta <- max(-0.5, -max((max.BB$BB8[2] + 0.001)/2, 0.001))
                theta1 <- max(-1.5, -max((max.BB$BB8[1] + 1.001)/2, 1.001))
            }
        } else if (family %in% c(104, 114, 124, 134, 204, 214, 224, 234)) {
            ## Tawn
            
            # the folllowing gives a theoretical kendall's tau close to the empirical one
            delta <- min(abs(tau) + 0.1, 0.999)
            theta1 <- 1 + 6 * abs(tau)
            
            # check if data can be modeled by selected family
            if (family %in% c(104, 114)) {
                if (tau < 0) {
                    print("The Tawn or survival Tawn copula cannot be used for negatively dependent data.")
                    delta <- 1
                    theta1 <- 1.001
                }
            } else if (family %in% c(124, 134)) {
                if (tau > 0) {
                    print("The rotated Tawn copula cannot be used for positively dependent data.")
                    delta <- 1
                    theta1 <- -1.001
                } else theta1 <- -theta1
                
            } else if (family %in% c(204, 214)) {
                if (tau < 0) {
                    print("The Tawn2 or survival Tawn2 copula cannot be used for negatively dependent data.")
                    delta <- 1
                    theta1 <- 1.001
                }
            } else if (family %in% c(224, 234)) {
                if (tau > 0) {
                    print("The rotated Tawn2 copula cannot be used for positively dependent data.")
                    delta <- 1
                    theta1 <- -1.001
                } else theta1 <- -theta1
            }
        }
        
        ## likelihood optimization
        if (family != 0 && family < 100) {
            out <- MLE_intern(cbind(u1, u2),
                              c(theta1, delta),
                              family = family,
                              se,
                              max.df, 
                              max.BB
                              , weights)
            theta <- out$par
            if (se == TRUE) 
                se1 <- out$se
        } else if (family != 0 && family > 100) {
            # New
            out <- MLE_intern_Tawn(cbind(u1, u2),
                                   c(theta1, delta),
                                   family = family,
                                   se)
            theta <- out$par
            if (se == TRUE) 
                se1 <- out$se
        }
    }
    
    ## store estimated parameters
    out2 <- list(family = family)
    if (length(theta) == 2) {
        out2$par <- theta[1]
        out2$par2 <- theta[2]
    } else {
        out2$par <- theta
        out2$par2 <- 0
    }
    
    ## store standard errors (if asked for)
    if (se == TRUE) {
        if (length(se1) == 2) {
            out2$se <- se1[1]
            out2$se2 <- se1[2]
        } else {
            out2$se <- se1
            out2$se2 <- 0
        }
    }
    
    ## return results
    class(out2) <- "BiCop"
    out2
}



Frank.itau.JJ <- function(tau) {
    a <- 1
    if (tau < 0) {
        a <- -1
        tau <- -tau
    }
    v <- uniroot(function(x) tau - (1 - 4/x + 4/x * debye1(x)), 
                 lower = 0 + .Machine$double.eps^0.5, upper = 5e5, 
                 tol = .Machine$double.eps^0.5)$root
    return(a*v)
}


Joe.itau.JJ <- function(tau) {
    if (tau < 0) {
        return(1.000001)
    } else {
        tauF <- function(par) {
            param1 <- 2/par + 1
            tem <- digamma(2) - digamma(param1)
            tau <- 1 + tem * 2/(2 - par)
            tau[par == 2] <- 1 - trigamma(2)
            tau
        }
        
        v <- uniroot(function(x) tau - tauF(x),
                     lower = 1, 
                     upper = 5e5,
                     tol = .Machine$double.eps^0.5)$root
        return(v)
    }
}

ipsA.tau2cpar <- function(tau, mxiter = 20, eps = 1e-06, dstart = 0, iprint = FALSE) {
    con <- log((1 - tau) * sqrt(pi)/2)
    de <- dstart
    if (dstart <= 0) 
        de <- tau + 1
    iter <- 0
    diff <- 1
    while (iter < mxiter & max(abs(diff)) > eps) {
        g <- con + lgamma(1 + de) - lgamma(de + 0.5)
        gp <- digamma(1 + de) - digamma(de + 0.5)
        iter <- iter + 1
        diff <- g/gp
        de <- de - diff
        while (min(de) <= 0) {
            diff <- diff/2
            de <- de + diff
        }
        if (iprint) 
            cat(iter, " ", de, " ", diff, "\n")
    }
    if (iter >= mxiter) 
        cat("did not converge\n")
    de
}


#############################################################
# bivariate MLE function
#
#------------------------------------------------------------
# INPUT:
#   data    Data for which to estimate parameter
#   start.parm	Start parameter for the MLE
#   Maxiter	max number of iterations
#   se		TRUE or FALSE
# OUTPUT:
#   out     Estimated Parameters and standard error (if se==TRUE)
#--------------------------------------------------------------
# Author: Ulf Schepsmeier
# Date: 2011-02-04
# Version: 1.1
#---------------------------------------------------------------

MLE_intern <- function(data, start.parm, family, se = FALSE, max.df = 30,
                       max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), 
                       weights = NULL, cor.fixed = FALSE) {
    
    n <- dim(data)[1]
    if (any(is.na(weights))) 
        weights <- NULL
    
    if (family %in% c(7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        t_LL <- function(param) {
            
            if (is.null(weights)) {
                ll <- .C("LL_mod2", 
                         as.integer(family),
                         as.integer(n),
                         as.double(data[, 1]),
                         as.double(data[, 2]), 
                         as.double(param[1]),
                         as.double(param[2]), 
                         as.double(0), 
                         PACKAGE = "VineCopula")[[7]]
            } else {
                ll <- .C("LL_mod_seperate", 
                         as.integer(family), 
                         as.integer(n),
                         as.double(data[, 1]), 
                         as.double(data[, 2]),
                         as.double(param[1]),
                         as.double(param[2]), 
                         as.double(rep(0, n)),
                         PACKAGE = "VineCopula")[[7]] %*% weights
            }
            
            if (is.infinite(ll) || is.na(ll) || ll < -10^250) 
                ll <- -10^250
            
            return(ll)
        }
        
        if (family == 7 || family == 17) {
            low <- c(0.001, 1.001)
            up <- max.BB$BB1
        } else if (family == 8 || family == 18) {
            low <- c(1.001, 1.001)
            up <- max.BB$BB6
        } else if (family == 9 | family == 19) {
            low <- c(1.001, 0.001)
            up <- max.BB$BB7
        } else if (family == 10 | family == 20) {
            low <- c(1.001, 0.001)
            up <- max.BB$BB8
        } else if (family == 27 | family == 37) {
            up <- c(-0.001, -1.001)
            low <- -max.BB$BB1
        } else if (family == 28 | family == 38) {
            up <- c(-1.001, -1.001)
            low <- -max.BB$BB6
        } else if (family == 29 | family == 39) {
            up <- c(-1.001, -0.001)
            low <- -max.BB$BB7
        } else if (family == 30 | family == 40) {
            up <- c(-1.001, -0.001)
            low <- -max.BB$BB8
        }
        
        if (se == TRUE) {
            optimout <- optim(par = start.parm, 
                              fn = t_LL,
                              method = "L-BFGS-B",
                              lower = low, 
                              upper = up,
                              control = list(fnscale = -1, maxit = 500), 
                              hessian = TRUE)
        } else {
            optimout <- optim(par = start.parm, 
                              fn = t_LL,
                              method = "L-BFGS-B", 
                              lower = low, 
                              upper = up, 
                              control = list(fnscale = -1, maxit = 500))
        }
        
    } else if (family == 2) {
        
        if (cor.fixed == FALSE) {
            
            t_LL <- function(param) {
                if (param[1] < -0.9999 | param[1] > 0.9999 | param[2] < 2.0001 | param[2] > max.df) {
                    ll <- -10^10
                } else {
                    if (is.null(weights)) {
                        ll <- .C("LL_mod2",
                                 as.integer(family),
                                 as.integer(n),
                                 as.double(data[, 1]), 
                                 as.double(data[, 2]), 
                                 as.double(param[1]), 
                                 as.double(param[2]), 
                                 as.double(0),
                                 PACKAGE = "VineCopula")[[7]]
                    } else {
                        ll <- .C("LL_mod_seperate",
                                 as.integer(family),
                                 as.integer(n), 
                                 as.double(data[, 1]), 
                                 as.double(data[, 2]), 
                                 as.double(param[1]),
                                 as.double(param[2]), 
                                 as.double(rep(0, n)),
                                 PACKAGE = "VineCopula")[[7]] %*% weights
                    }
                    
                    if (is.infinite(ll) || is.na(ll) || ll < -10^10) 
                        ll <- -10^10
                }
                return(ll)
            }
            
            gr_LL <- function(param) {
                gr <- rep(0, 2)
                gr[1] <- sum(BiCopDeriv(data[, 1], 
                                        data[, 2],
                                        family = 2,
                                        par = param[1],
                                        par2 = param[2], 
                                        deriv = "par",
                                        log = TRUE,
                                        check.pars = FALSE))
                gr[2] <- sum(BiCopDeriv(data[, 1], 
                                        data[, 2], 
                                        family = 2, 
                                        par = param[1],
                                        par2 = param[2], 
                                        deriv = "par2",
                                        log = TRUE,
                                        check.pars = FALSE))
                return(gr)
            }
            
            if (is.null(weights)) {
                if (se == TRUE) {
                    optimout <- optim(par = start.parm, 
                                      fn = t_LL, 
                                      gr = gr_LL, 
                                      method = "L-BFGS-B", 
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE, 
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                } else {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL,
                                      gr = gr_LL, 
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                }
            } else {
                if (se == TRUE) {
                    optimout <- optim(par = start.parm,
                                      fn = t_LL, 
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500), 
                                      hessian = TRUE, 
                                      lower = c(-0.9999, 2.0001),
                                      upper = c(0.9999, max.df))
                } else {
                    optimout <- optim(par = start.parm, 
                                      fn = t_LL, 
                                      method = "L-BFGS-B", 
                                      control = list(fnscale = -1, maxit = 500), 
                                      lower = c(-0.9999, 2.0001), 
                                      upper = c(0.9999, max.df))
                }
            }
            
            if (optimout$par[2] >= (max.df - 1e-04)) 
                warning(paste("Degrees of freedom of the t-copula estimated to be larger than ", 
                              max.df, ". Consider using the Gaussian copula instead.", 
                              sep = ""))
            
        } else {
            t_LL <- function(param) {
                if (is.null(weights)) {
                    ll <- .C("LL_mod2",
                             as.integer(family),
                             as.integer(n), 
                             as.double(data[, 1]),
                             as.double(data[, 2]), 
                             as.double(start.parm[1]),
                             as.double(param[1]), 
                             as.double(0), 
                             PACKAGE = "VineCopula")[[7]]
                } else {
                    ll <- .C("LL_mod_seperate",
                             as.integer(family),
                             as.integer(n), 
                             as.double(data[, 1]), 
                             as.double(data[, 2]), 
                             as.double(start.parm[1]), 
                             as.double(param[1]),
                             as.double(rep(0, n)), 
                             PACKAGE = "VineCopula")[[7]] %*% weights
                }

                if (is.infinite(ll) || is.na(ll) || ll < -10^250) 
                    ll <- -10^250
                
                return(ll)
            }
            
            gr_LL <- function(param) {
                gr <- sum(BiCopDeriv(data[, 1],
                                     data[, 2],
                                     family = 2,
                                     par = start.parm[1],
                                     par2 = param[1], 
                                     deriv = "par2", 
                                     log = TRUE,
                                     check.pars = FALSE))
                return(gr)
            }
            
            if (se == TRUE) {
                if (is.null(weights)) {
                    optimout <- optim(par = start.parm[2], 
                                      fn = t_LL, 
                                      gr = gr_LL, 
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE, 
                                      lower = 2.0001, 
                                      upper = max.df)
                } else {
                    optimout <- optim(par = start.parm[2], 
                                      fn = t_LL, 
                                      method = "L-BFGS-B", 
                                      control = list(fnscale = -1, maxit = 500),
                                      hessian = TRUE,
                                      lower = 2.0001, 
                                      upper = max.df)
                }
            } else {
                optimout <- optimize(f = t_LL, 
                                     maximum = TRUE, 
                                     interval = c(2.0001, max.df))
                optimout$par <- optimout$maximum
            }
            optimout$par <- c(0, optimout$par)
            
        }
        
    } else {
        
        t_LL <- function(param) {
            if (is.null(weights)) {
                ll <- .C("LL_mod2", as.integer(family), 
                         as.integer(n), 
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param),
                         as.double(0), as.double(0), 
                         PACKAGE = "VineCopula")[[7]]
            } else {
                ll <- .C("LL_mod_seperate", 
                         as.integer(family), 
                         as.integer(n), 
                         as.double(data[, 1]),
                         as.double(data[, 2]),
                         as.double(param[1]), 
                         as.double(0), 
                         as.double(rep(0, n)), 
                         PACKAGE = "VineCopula")[[7]] %*% weights
            }
            if (is.infinite(ll) || is.na(ll) || ll < -10^250) 
                ll <- -10^250
            
            return(ll)
        }
        
        gr_LL <- function(param) {
            gr <- sum(BiCopDeriv(data[, 1],
                                 data[, 2],
                                 family = family, 
                                 par = param,
                                 deriv = "par", 
                                 log = TRUE,
                                 check.pars = FALSE))
            return(gr)
        }
        
        low <- -Inf
        up <- Inf
        
        if (family == 1) {
            low <- -0.9999
            up <- 0.9999
        } else if (family %in% c(3, 13)) {
            low <- 1e-04
            up <- BiCopTau2Par(family, 0.99)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.95)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.9)
        } else if (family %in% c(4, 14)) {
            low <- 1.0001
            up <- BiCopTau2Par(family, 0.99)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.95)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.9)
        } else if (family %in% c(5)) {
            low <- BiCopTau2Par(family, -0.99)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.9)
            up <- BiCopTau2Par(family, 0.99)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.95)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.9)
        } else if (family %in% c(6, 16)) {
            low <- 1.0001
            up <- BiCopTau2Par(family, 0.99)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.95)
            if (t_LL(up) == -10^250) 
                up <- BiCopTau2Par(family, 0.9)
        } else if (family %in% c(23, 33)) {
            up <- -1e-04
            low <- BiCopTau2Par(family, -0.99)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.9)
        } else if (family %in% c(24, 34)) {
            up <- -1.0001
            low <- BiCopTau2Par(family, -0.99)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.9)
        } else if (family %in% c(26, 36)) {
            up <- -1.0001
            low <- BiCopTau2Par(family, -0.99)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.9)
        } else if (family %in% c(41, 51)) {
            low <- 1e-04
            up <- BiCopTau2Par(family, 0.85)
            # if(t_LL(up)==-10^250) up=BiCopTau2Par(family,0.95) if(t_LL(up)==-10^250) up=BiCopTau2Par(family,0.9)
        } else if (family %in% c(61, 71)) {
            up <- -1e-04
            low <- BiCopTau2Par(family, -0.85)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.95)
            if (t_LL(low) == -10^250) 
                low <- BiCopTau2Par(family, -0.9)
        }
        
        pscale <- ifelse(family == 1, 0.001, 1)
        
        if (se == TRUE) {
            if (is.null(weights)) {
                optimout <- try(expr = optim(par = start.parm[1],
                                             fn = t_LL,
                                             gr = gr_LL,
                                             method = "L-BFGS-B", 
                                             control = list(fnscale = -1,
                                                            maxit = 500,
                                                            parscale = pscale),
                                             lower = low,
                                             upper = up, 
                                             hessian = TRUE), 
                                silent = TRUE)
                if (class(optimout) == "try-error") {
                    optimout <- optim(par = start.parm[1],
                                      fn = t_LL, 
                                      gr = NULL, 
                                      method = "L-BFGS-B",
                                      control = list(fnscale = -1, 
                                                     maxit = 500,
                                                     parscale = pscale), 
                                      lower = low,
                                      upper = up,
                                      hessian = TRUE)
                }
            } else {
                optimout <- optim(par = start.parm[1],
                                  fn = t_LL, 
                                  gr = NULL, 
                                  method = "L-BFGS-B",
                                  control = list(fnscale = -1,
                                                 maxit = 500,
                                                 parscale = pscale), 
                                  lower = low,
                                  upper = up,
                                  hessian = TRUE)
            }
        } else {
            optimout <- optimize(f = t_LL, interval = c(low, up), maximum = TRUE)
            optimout$par <- optimout$maximum
        }
        optimout$par <- c(optimout$par, 0)
        
    }
    
    out <- list()
    
    if (se == TRUE) {
        if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
            out$par <- optimout$par
            
            if (det(optimout$hessian) == 0) {
                var <- diag(1, dim(optimout$hessian)[1])
            } else {
                var <- (-solve(optimout$hessian))
            }
            
            out$se <- sqrt(diag(var))
            
            if (family == 2 && out$par[2] >= (max.df - 1e-04)) 
                out$se[2] <- NA
            
        } else {
            out$par <- optimout$par[1]
            
            if (optimout$hessian == 0) {
                var <- 1
            } else {
                if (optimout$hessian < 0) 
                    optimout$hessian <- -optimout$hessian
                var <- 1/optimout$hessian
            }
            
            out$se <- as.numeric(sqrt(var))
        }
        
    } else {
        if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
            
            out$par <- optimout$par
            
        } else {
            
            out$par[1] <- optimout$par[1]
        }
        
    }
    out$value <- optimout$value
    return(out)
}


# New for Tawn

MLE_intern_Tawn <- function(data, start.parm, family, se = FALSE) {
    
    n <- dim(data)[1]
    
    ## set bounds for optimization
    tau <- fasttau(data[, 1], data[, 2])
    if (family == 104 || family == 114 || family == 204 || family == 214) {
        parlower <- c(1.001, max(tau - 0.1, 1e-04))
        parupper <- c(20, min(tau + 0.2, 0.99))
    } else if (family == 124 || family == 134 || family == 224 || family == 234) {
        parlower <- c(-20, max(-tau - 0.1, 1e-04))
        parupper <- c(-1.001, min(-tau + 0.2, 0.99))
    }
    
    ## log-liklihood function
    loglikfunc <- function(param) {
        ll <- .C("LL_mod2",
                 as.integer(family), 
                 as.integer(n), 
                 as.double(data[, 1]),
                 as.double(data[, 2]), 
                 as.double(param[1]), 
                 as.double(param[2]), 
                 as.double(0), 
                 PACKAGE = "VineCopula")[[7]]
        if (is.infinite(ll) || is.na(ll) || ll < -10^250) 
            ll <- -10^250
        # print(param) print(ll)
        return(ll)
    }
    
    ## optimize log-likelihood
    out <- list()
    # print(start.parm)
    if (se == TRUE) {
        optimout <- optim(par = start.parm,
                          fn = loglikfunc,
                          method = c("L-BFGS-B"), 
                          lower = parlower,
                          upper = parupper, 
                          control = list(fnscale = -1, maxit = 500), 
                          hessian = TRUE)
        if (det(optimout$hessian) == 0) {
            var <- diag(1, dim(optimout$hessian)[1])
        } else {
            var <- (-solve(optimout$hessian))
        }
        
        out$se <- sqrt(diag(var))
    } else {
        optimout <- optim(par = start.parm,
                          fn = loglikfunc,
                          method = c("L-BFGS-B"),
                          lower = parlower,
                          upper = parupper, 
                          control = list(fnscale = -1, maxit = 500))
    }
    
    ## return results
    out$par <- optimout$par
    out$value <- optimout$value
    return(out)
}



fasttau <- function(x, y, weights = NA) {
    if (any(is.na(weights))) {
        m <- length(x)
        n <- length(y)
        if (m == 0 || n == 0) 
            stop("both 'x' and 'y' must be non-empty")
        if (m != n) 
            stop("'x' and 'y' must have the same length")
        out <- .C("ktau",
                  x = as.double(x),
                  y = as.double(y),
                  N = as.integer(n),
                  tau = as.double(0),
                  S = as.double(0),
                  D = as.double(0),
                  T = as.integer(0), 
                  U = as.integer(0), 
                  V = as.integer(0), 
                  PACKAGE = "VineCopula")
        ktau <- out$tau
    } else {
        ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
    }
    return(ktau)
}
