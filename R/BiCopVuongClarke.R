#########################################################################
# Functions of Natalia Djunushalieva for Vuong and Clarke test    	#
# (copula selection - goodness-of-fit)					#
# vuong.clarke.sim							#
# vuong.test								#
# clarke.test								#
# score									#
# gofVC			goodness-of-fit on the basis of Vuong and Clarke#
#########################################################################


#################################################################
# vuong.clarke.sim						#
#								#
# Input:							#
# xy		2-dim matrix of data				#
# cop.model1	character, copula on the left side		#
# cop.model2	vector of character, copula on the righ side	#
# Output:							#
# sim.result	Simulation of the models			#
#################################################################

vuong.clarke.sim <- function(xy, cop.model1, cop.model2, correction, level) {
    n <- dim(xy)[1]
    
    loglik.model1 <- my.loglik(xy, cop.model1)
    loglik.model2 <- data.frame(sapply(cop.model2, function(x) {
        my.loglik(fam = x, data = xy)
    }, USE.NAMES = TRUE))
    vuong.result.temp <- list()
    clarke.result.temp <- list()
    
    for (i in 1:length(cop.model2)) {
        p1 <- 1
        p2 <- 1
        
        
        if (cop.model1 == cop.model2[i]) {
            test.result <- data.frame(NA, NA, NA, NA)
            names(test.result) <- c("model", "nu", "p.value", "kurtosis")
            vuong.result.temp[[i]] <- test.result
            names(vuong.result.temp)[i] <- cop.model2[i]
            clarke.result.temp[[i]] <- test.result
            names(clarke.result.temp)[i] <- cop.model2[i]
        } else if (cop.model1 != cop.model2[i]) {
            if (cop.model1 %in% c(2, 7:10, 17:20, 27:30, 37:40,
                                  104, 114, 124, 134, 204, 214, 224, 234)) {
                p1 <- 2
            }
            if (cop.model2[i] %in% c(2, 7:10, 17:20, 27:30, 37:40,
                                     104, 114, 124, 134, 204, 214, 224, 234)) {
                p2 <- 2
            }
            
            vuong.result.temp[[i]] <- vuong.test(loglik.model1,
                                                 loglik.model2[, i], 
                                                 alpha = level,
                                                 p1 = p1, 
                                                 p2 = p2,
                                                 print.result = FALSE, 
                                                 correction)
            names(vuong.result.temp)[i] <- cop.model2[i]
            clarke.result.temp[[i]] <- clarke.test(loglik.model1,
                                                   loglik.model2[, i],
                                                   alpha = level, 
                                                   p1 = p1, 
                                                   p2 = p2, 
                                                   print.result = FALSE,
                                                   correction)
            names(clarke.result.temp)[i] <- cop.model2[i]
        }
    }
    sim.result <- list(vuong.result.temp, clarke.result.temp)
    names(sim.result) <- c("vuong", "clarke")
    return(sim.result)
}


#########################################
# my.loglik				#
#					#
# Input:				#
# data		2-dim matrix		#
# fam		copula family		#
# Output:				#
# loglik.temp	Log likelihood		#
#########################################

my.loglik <- function(data, fam) {
    T <- dim(data)[1]
    
    # Parameter sch?tzen
    par <- BiCopEst(data[, 1], data[, 2], fam)
    if (fam %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234)) {
        theta <- par$par
        nu <- par$par2
    } else {
        theta <- par$par
        nu <- 0
    }
    
    loglik.temp <- rep(0, T)
    for (i in 1:T) {
        out <- .C("LL_mod2", 
                  as.integer(fam), 
                  as.integer(1), 
                  as.double(data[i, 1]), 
                  as.double(data[i, 2]), 
                  as.double(theta), 
                  as.double(nu), 
                  as.double(0), 
                  PACKAGE = "VineCopula")
        loglik.temp[i] <- out[[7]]
        if (loglik.temp[i] %in% c(NA, NaN, Inf))  loglik.temp[i] <- 1e+10
    }
    
    return(loglik.temp)
}




#################################################################
# vuong.test							#
#								#
# Input:							#
# loglik.model1	vector of log likelihood of the first model	#
# loglik.model1	vector of log likelihood of the second model	#
# alpha		alpha-level for p-value				#
# p1		number of parameters (model 1)			#
# p2		number of parameters (model 2)			#
# Output:							#
# Vuong test							#
#################################################################

vuong.test <- function(loglik.model1, loglik.model2, alpha = 0.05, p1 = 0, p2 = 0, print.result = TRUE, correction) {
    # cat('H0: model (1) is eqivalent to model (2)', '\n')
    
    # Calculate test statistic
    n <- length(loglik.model1)
    if (correction == FALSE) {
        m.i <- loglik.model1 - loglik.model2 
    } else if (correction == "Akaike") { 
        m.i <- loglik.model1 - loglik.model2 - (p1 - p2)/n 
    } else if (correction == "Schwarz") {
        m.i <- loglik.model1 - loglik.model2 - (p1 - p2) * log(n)/(2 * n)
    }
    kurt.ratios <- kurtosis(loglik.model1 - loglik.model2)
    nu <- (sqrt(n) * mean(m.i))/(sqrt((n - 1)/n * var(m.i)))
    
    if (abs(nu) < qnorm(1 - alpha/2)) {
        decision <- "Decision: none of the models is favored."
        result <- 0
    }
    if (nu >= qnorm(1 - alpha/2)) {
        decision <- "Decision: favor model 1."
        result <- 1
    }
    if (nu <= -qnorm(1 - alpha/2)) {
        decision <- "Decision: favor model 2."
        result <- 2
    }
    # cat(decision,'\n')
    
    pvalue <- 2 * pnorm(-abs(nu))
    result <- data.frame(result,
                         round(nu, digits = 3),
                         round(pvalue, digits = 3), 
                         round(kurt.ratios, digits = 3))
    names(result) <- c("model", "nu", "p.value", "kurtosis")
    
    if (print.result) print(result)
    
    return(result)
    rm(n, m.i, kurt.ratios, nu, pvalue, result, decision)
}

#################################################################
# clarke.test							#
#								#
# Input:							#
# loglik.model1	vector of log likelihood of the first model	#
# loglik.model1	vector of log likelihood of the second model	#
# alpha		alpha-level for p-value				#
# p1		number of parameters (model 1)			#
# p2		number of parameters (model 2)			#
# Output:							#
# Clarke test							#
#################################################################


clarke.test <- function(loglik.model1, loglik.model2, alpha = 0.05, p1 = 0, p2 = 0, print.result = TRUE, correction) {
    # cat('H0: model (1) is eqivalent to model (2)', '\n')
    
    # Calculate test statistic
    n <- length(loglik.model1)
    if (correction == FALSE) {
        m.i <- loglik.model1 - loglik.model2 
    } else if (correction == "Akaike") { 
        m.i <- loglik.model1 - loglik.model2 - (p1 - p2)/n 
    } else if (correction == "Schwarz") {
        m.i <- loglik.model1 - loglik.model2 - (p1 - p2) * log(n)/(2 * n)
    }
    kurt.ratios <- kurtosis(loglik.model1 - loglik.model2)
    B <- sum(m.i > 0)
    
    # Calculate critical value
    decision <- "Decision: non of the models is favoured"
    result <- 0
    if (B >= n/2) {
        pvalue <- 2 * (1 - pbinom(B - 1, n, 0.5))
        if (pvalue <= alpha) {
            decision <- "Decision: favour model 1"
            result <- 1
        }
    }
    if (B < n/2) {
        pvalue <- 2 * (pbinom(B, n, 0.5))
        if (pvalue <= alpha) {
            decision <- "Decision: favour model 2"
            result <- 2
        }
    }
    # cat(decision,'\n')
    
    result <- data.frame(result, 
                         round(B, digits = 3), 
                         round(pvalue, digits = 3), 
                         round(kurt.ratios, digits = 3))
    names(result) <- c("model", "B", "p.value", "kurtosis")
    
    if (print.result) print(result)
    
    return(result)
    rm(n, m.i, kurt.ratios, B, result, decision, pvalue)
}


#################################
# Score Funktion		#
#				#
# Input:			#
# dat		data matrix	#
# Output:			#
# score				#
#################################

score <- function(dat) {
    if (sum(is.na(dat)) != length(dat)) {
        yes <- sum(dat == 1, na.rm = TRUE)  # zaehle wieviel mal kommt yes vor
        no <- sum(dat == 2, na.rm = TRUE) * (-1)  # zaehle wieviel mal kommt no vor
        score.out <- yes + no
    }
    if (sum(is.na(dat)) == length(dat)) {
        score.out <- NA
    }
    return(score.out)
}



#################################################
# Goodness-of-fit test with Vuong and Clarke	#
#						#
# Input:					#
# xy		data matrix (2-dim)		#
# family	vector of copula families	#
# Output:					#
# score.gof	Goodness-of-fit score		#
#################################################

BiCopVuongClarke <- function(u1, u2, familyset = NA, correction = FALSE, level = 0.05) {
    
    if (is.na(familyset[1])) 
        familyset <- c(1:10, 13, 14, 16:20, 23, 24, 26:30,
                       33, 34, 36:40, 41, 51,  61, 71,
                       104, 114, 124, 134, 204, 214, 224, 234)
    # Sicherheitsabfragen
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
    for (i in 1:length(familyset)) {
        if (!(familyset[i] %in% c(1:10, 13, 14, 16:20, 23, 24, 26:30,
                                  33, 34, 36:40, 41, 51, 61, 71,
                                  104, 114, 124, 134, 204, 214, 224, 234))) 
            stop("Copula family not implemented.")
    }
    
    xy <- cbind(u1, u2)
    score.gof <- matrix(rep(NA, 2 * length(familyset)), nrow = 2)
    dimnames(score.gof) <- list(c("Vuong", "Clarke"), familyset)
    
    for (i in 1:length(familyset)) {
        temp <- vuong.clarke.sim(xy, familyset[i], familyset, correction, level)
        temp.vuong <- c()
        temp.clarke <- c()
        for (j in 1:length(temp$vuong)) {
            temp.vuong <- cbind(temp.vuong, t(temp$vuong[[j]]))
            dimnames(temp.vuong)[[2]][j] <- names(temp$vuong)[j]
            temp.clarke <- cbind(temp.clarke, t(temp$clarke[[j]]))
            dimnames(temp.clarke)[[2]][j] <- names(temp$clarke)[j]
        }
        score.temp <- rbind(score(temp.vuong[1, ]), score(temp.clarke[1, ]))
        score.gof[, i] <- score.temp
    }
    
    return(score.gof)
}


#########################
# kurtosis		#
#########################

kurtosis <- function(x) {
    x <- x[!is.na(x)]
    kurtosis <- sum((x - mean(x))^4/as.numeric(var(x))^2)/length(x)
    return(kurtosis)
}
