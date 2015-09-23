# ===============================================================================
# -------------------- CHI-PLOT FOR BIVARIATE DATA
# -----------------------------
# ===============================================================================
# Author: Natalia Djunushalieva, TU Muenchen, April 2010 Update: Ulf
# Schepsmeier, TU Muenchen, June 2010 For more detail see 'Everything you
# always wanted to now about copula modeling but were afraid to ask',
# Christian Genest, Anne-Catherine Favre NOTE: It is also possible to
# calculate chi-plot for righ upper and left lower quadrant of the data for
# determining upper or lower tail dependence. For more details see 'A simple
# graphical method to explore tail-dependence in stock-return pairs', Klaus
# Abberger, University of Konstanz, Germany
# -------------------------------------------------------------------------------
BiCopChiPlot <- function(u1, u2, PLOT = TRUE, mode = "NULL", ...) {
    # -----------------------------------------------------------------------------
    # PARAMETER DESCRIPTION: u1,u2 -> numeric, two vectors u1 and u2 of the same
    # length or u1 is a data frame of dimension d=2 PLOT -> logical, should the
    # results be plotted? If FALSE, the values W.in and Hi and contron bounds
    # will be returned mode -> character or NULL, generall chi plot or
    # lower/upper chi plot. Possible values are NULL, 'upper','lower'
    # -----------------------------------------------------------------------------
    
    # validation of input parameter
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (length(u1) < 2) 
        stop("Number of observations has to be at least 2.")
    
    if (PLOT != TRUE && PLOT != FALSE) 
        stop("The parameter 'PLOT' has to be set to 'TRUE' or 'FALSE'.")
    
    # Computing of results
    n <- length(u1)
    Hi <- H(u1, u2)
    Fi <- F(u1)
    Gi <- F(u2)
    
    lambda <- 4 * sign((Fi - 0.5) * (Gi - 0.5)) * apply(data.frame((Fi - 0.5)^2, 
                                                                   (Gi - 0.5)^2), 1, max)
    control.bounds <- c(1.54/sqrt(n), -1.54/sqrt(n))
    
    if (mode == "upper") {
        to.keep <- intersect(which(u1 > mean(u1)), which(u2 > mean(u2)))
        # to.keep<-intersect(which(u1>0),which(u2>0))
        to.keep <- intersect(to.keep, which(lambda > 0))
        Hi <- Hi[to.keep]
        Fi <- Fi[to.keep]
        Gi <- Gi[to.keep]
        lambda <- lambda[to.keep]
    }
    if (mode == "lower") {
        to.keep <- intersect(which(u1 < mean(u1)), which(u2 < mean(u2)))
        # to.keep<-intersect(which(u1<0),which(u2<0))
        to.keep <- intersect(to.keep, which(lambda > 0))
        Hi <- Hi[to.keep]
        Fi <- Fi[to.keep]
        Gi <- Gi[to.keep]
        lambda <- lambda[to.keep]
    }
    # remote entries with null values
    to.remove <- c(which(Fi == 0), which(Gi == 0), which(Fi == 1), which(Gi == 
                                                                             1))
    if (length(to.remove) != 0) {
        Hi <- Hi[-to.remove]
        Fi <- Fi[-to.remove]
        Gi <- Gi[-to.remove]
        lambda <- lambda[-to.remove]
    }
    
    chi <- (Hi - Fi * Gi)/sqrt(Fi * (1 - Fi) * Gi * (1 - Gi))
    control.bounds <- c(1.54/sqrt(n), -1.54/sqrt(n))
    
    if (PLOT) {
        # plotting of results
        plot(lambda, chi, xlab = expression(lambda), ylab = expression(chi), 
             ...)
        abline(h = control.bounds[1], col = "gray", lty = "dashed")
        abline(h = control.bounds[2], col = "gray", lty = "dashed")
        abline(h = 0, col = "gray")
        abline(v = 0, col = "gray")
    } else {
        # create output data
        chi.plot.output <- list(lambda, chi, control.bounds)
        names(chi.plot.output) <- c("lambda", "chi", "control.bounds")
        return(chi.plot.output)
    }
}  # end of cho.plot-fnction


# ===============================================================================
# ----------------- KENDALL-PLOT FOR BIVARIATE DATA
# ----------------------------
# ===============================================================================
# Author: Natalia Djunushalieva, TU Muenchen, April 2010 Update: Ulf
# Schepsmeier, TU Muenchen, June 2010 For more detail see 'Everything you
# always wanted to now about copula modeling but were afraid to ask',
# Christian Genest, Anne-Catherine Favre
# -------------------------------------------------------------------------------
BiCopKPlot <- function(u1, u2, PLOT = TRUE, ...) {
    # -----------------------------------------------------------------------------
    # INPUT: PARAMETER DESCRIPTION: u1,u2 - numeric, two vectors u1 and u2 of
    # the same length or u1 is a data frame of dimension d=2 PLOT - logical,
    # should the results be plotted? If FALSE, the values W.in and Hi will be
    # return.
    # -----------------------------------------------------------------------------
    
    # validation of input data
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (length(u1) < 2) 
        stop("Number of observations has to be at least 2.")
    
    if (PLOT != TRUE && PLOT != FALSE) 
        stop("The parameter 'PLOT' has to be set to 'TRUE' or 'FALSE'.")
    
    # Computing of results
    Wi <- W(u1, u2)
    Hi <- H(u1, u2)
    Hi.sort <- sort(Hi)
    
    n <- length(u1)
    
    W.in <- rep(NA, n)
    for (i in 1:n) {
        f <- function(w) {
            w * (-log(w)) * (w - w * log(w))^(i - 1) * (1 - w + w * log(w))^(n - 
                                                                                 i)
        }  # function to be integrated
        W.in[i] <- n * choose(n - 1, i - 1) * (integrate(f, lower = 0, upper = 1)$value)  # W_{i:n} for i=1:n
    }
    g <- function(w) {
        w - w * log(w)
    }  # K_{0}(w)=P(UV<=w)
    
    if (PLOT) {
        # should the results be plotted?
        plot(g, xlim = c(0, 1), ylim = c(0, 1), pch = "x", xlab = expression(W[1:n]), 
             ylab = "H", ...)  #Kurve K_{0}(w)
        points(W.in, Hi.sort, pch = "x", cex = 0.4, ...)
        abline(a = 0, b = 1)  # angle bisector
    } else {
        # create output data
        kendall.plot.output <- list(W.in, Hi.sort)
        names(kendall.plot.output) <- c("W.in", "Hi.sort")
        return(kendall.plot.output)
    }
}  # end of kendall.plot-function



# ===============================================================================
# --------------- HELP FUNCTIONS FOR CHI- AND KENDALL-PLOTS
# --------------------
# ===============================================================================
F <- function(t) {
    # help function for chi-plot
    n <- length(t)
    result <- (rank(t) - 1)/(n - 1)
    return(result)
}
# -------------------------------------------------------------------------------
H <- function(t, s) {
    # help function for chi-plot
    n <- length(t)
    H.result <- c()
    
    for (i in 1:n) {
        rank.smaller.t <- setdiff(which(rank(t) <= rank(t)[i]), i)
        rank.smaller.s <- setdiff(which(rank(s) <= rank(s)[i]), i)
        H.result[i] <- length(intersect(rank.smaller.t, rank.smaller.s))
    }
    
    H.result <- H.result/(n - 1)
    return(H.result)
}
# -------------------------------------------------------------------------------
W <- function(t, s) {
    # Help function for kendall.plot
    n <- length(t)
    result <- ((n - 1) * H(t, s) + 1)/n
    return(result)
}
# -------------------------------------------------------------------------------
