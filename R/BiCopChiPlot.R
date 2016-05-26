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


#' Chi-plot for Bivariate Copula Data
#'
#' This function creates a chi-plot of given bivariate copula data.
#'
#' For observations \eqn{u_{i,j},\ i=1,...,N,\ j=1,2,}{u_{i,j}, i=1,...,N,
#' j=1,2,} the chi-plot is based on the following two quantities: the
#' chi-statistics
#' \deqn{\chi_i = \frac{\hat{F}_{U_1U_2}(u_{i,1},u_{i,2})
#' - \hat{F}_{U_1}(u_{i,1})\hat{F}_{U_2}(u_{i,2})}{
#' \sqrt{\hat{F}_{U_1}(u_{i,1})(1-\hat{F}_{U_1}(u_{i,1}))
#' \hat{F}_{U_2}(u_{i,2})(1-\hat{F}_{U_2}(u_{i,2}))}}, }{
#' \chi_i = F_{U_1,U_2}(u_{i,1},u_{i,2}) - F_{U_1}(u_{i,1})F_{U_2}(u_{i,2})
#' / (F_{U_1}(u_{i,1}) (1-F_{U_1}(u_{i,1}))F_{U_2}(u_{i,2})
#' (1-F_{U_2}(u_{i,2}))^0.5, } and the lambda-statistics
#' \deqn{\lambda_i = 4 sgn\left( \tilde{F}_{U_1}(u_{i,1}),\tilde{F}_{U_2}(u_{i,2}) \right)
#' \cdot \max\left( \tilde{F}_{U_1}(u_{i,1})^2,\tilde{F}_{U_2}(u_{i,2})^2 \right), }{
#' \lambda_i = 4 sgn( tildeF_{U_1}(u_{i,1}),tildeF_{U_2}(u_{i,2}) )
#' * max( tildeF_{U_1}(u_{i,1})^2,tildeF_{U_2}(u_{i,2})^2 ), }
#' where \eqn{\hat{F}_{U_1}}{F_{U_1}}, \eqn{\hat{F}_{U_2}}{F_{U_2}} and
#' \eqn{\hat{F}_{U_1U_2}}{F_{U_1U_2}} are the empirical distribution functions
#' of the uniform random variables \eqn{U_1} and \eqn{U_2} and of
#' \eqn{(U_1,U_2)}, respectively. Further,
#' \eqn{\tilde{F}_{U_1}=\hat{F}_{U_1}-0.5}{tildeF_{U_1}=F_{U_1}-0.5} and
#' \eqn{\tilde{F}_{U_2}=\hat{F}_{U_2}-0.5}{tildeF_{U_2}=F_{U_2}-0.5}.
#'
#' These quantities only depend on the ranks of the data and are scaled to the
#' interval \eqn{[0,1]}. \eqn{\lambda_i} measures a distance of a data point
#' \eqn{\left(u_{i,1},u_{i,2}\right)}{(u_{i,1},u_{i,2})} to the center of the
#' bivariate data set, while \eqn{\chi_i} corresponds to a correlation
#' coefficient between dichotomized values of \eqn{U_1} and \eqn{U_2}. Under
#' independence it holds that \eqn{\chi_i \sim
#' \mathcal{N}(0,\frac{1}{N})}{\chi_i~N(0,1/N)} and \eqn{\lambda_i \sim
#' \mathcal{U}[-1,1]}{\lambda_i~U[0,1]} asymptotically, i.e., values of
#' \eqn{\chi_i} close to zero indicate independence---corresponding to
#' \eqn{F_{U_1U_2}=F_{U_1}F_{U_2}}.
#'
#' When plotting these quantities, the pairs of \eqn{\left(\lambda_i, \chi_i
#' \right)}{(\lambda_i,\chi_i)} will tend to be located above zero for
#' positively dependent margins and vice versa for negatively dependent
#' margins. Control bounds around zero indicate whether there is significant
#' dependence present.
#'
#' If \code{mode = "lower"} or \code{"upper"}, the above quantities are
#' calculated only for those \eqn{u_{i,1}}'s and \eqn{u_{i,2}}'s which are
#' smaller/larger than the respective means of
#' \code{u1}\eqn{=(u_{1,1},...,u_{N,1})} and
#' \code{u2}\eqn{=(u_{1,2},...,u_{N,2})}.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param PLOT Logical; whether the results are plotted. If \code{PLOT =
#' FALSE}, the values \code{lambda}, \code{chi} and \code{control.bounds} are
#' returned (see below; default: \code{PLOT = TRUE}).
#' @param mode Character; whether a general, lower or upper chi-plot is
#' calculated.  Possible values are \code{mode = "NULL"}, \code{"upper"} and
#' \code{"lower"}. \cr \code{"NULL"} = general chi-plot (default)\cr
#' \code{"upper"} = upper chi-plot\cr \code{"lower"} = lower chi-plot
#' @param ... Additional plot arguments.
#' @return \item{lambda}{Lambda-statistics (x-axis).} \item{chi}{Chi-statistics
#' (y-axis).} \item{control.bounds}{A 2-dimensional vector of bounds
#' \eqn{((1.54/\sqrt{n},-1.54/\sqrt{n})}, where \eqn{n} is the length of
#' \code{u1} and where the chosen values correspond to an approximate
#' significance level of 10\%.}
#' @author Natalia Belgorodski, Ulf Schepsmeier
#' @seealso \code{\link{BiCopMetaContour}}, \code{\link{BiCopKPlot}},
#' \code{\link{BiCopLambda}}
#' @references Abberger, K. (2004). A simple graphical method to explore
#' tail-dependence in stock-return pairs. Discussion Paper, University of
#' Konstanz, Germany.
#'
#' Genest, C. and A. C. Favre (2007). Everything you always wanted to know
#' about copula modeling but were afraid to ask. Journal of Hydrologic
#' Engineering, 12 (4), 347-368.
#' @examples
#'
#' ## chi-plots for bivariate Gaussian copula data
#'
#' # simulate copula data
#' fam <- 1
#' tau <- 0.5
#' par <- BiCopTau2Par(fam, tau)
#' cop <- BiCop(fam, par)
#' set.seed(123)
#' dat <- BiCopSim(500, cop)
#'
#' # create chi-plots
#' op <- par(mfrow = c(1, 3))
#' BiCopChiPlot(dat[,1], dat[,2], xlim = c(-1,1), ylim = c(-1,1),
#'              main="General chi-plot")
#' BiCopChiPlot(dat[,1], dat[,2], mode = "lower", xlim = c(-1,1),
#'              ylim = c(-1,1), main = "Lower chi-plot")
#' BiCopChiPlot(dat[,1], dat[,2], mode = "upper", xlim = c(-1,1),
#'              ylim = c(-1,1), main = "Upper chi-plot")
#' par(op)
#'
#' @export BiCopChiPlot
BiCopChiPlot <- function(u1, u2, PLOT = TRUE, mode = "NULL", ...) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    ## sanity checks
    if (length(u1) < 2)
        stop("Number of observations has to be at least 2.")
    stopifnot(is.logical(PLOT))

    ## Computations
    n <- length(u1)
    Hi <- Hn(u1, u2)
    Fi <- Fn(u1)
    Gi <- Fn(u2)
    lambda <- 4 * sign((Fi - 0.5) * (Gi - 0.5)) * apply(data.frame((Fi - 0.5)^2,
                                                                   (Gi - 0.5)^2),
                                                        1,
                                                        max)
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


#' Kendall's Plot for Bivariate Copula Data
#'
#' This function creates a Kendall's plot (K-plot) of given bivariate copula
#' data.
#'
#' For observations \eqn{u_{i,j},\ i=1,...,N,\ j=1,2,}{u_{i,j}, i=1,...,N,
#' j=1,2,} the K-plot considers two quantities: First, the ordered values of
#' the empirical bivariate distribution function
#' \eqn{H_i:=\hat{F}_{U_1U_2}(u_{i,1},u_{i,2})} and, second, \eqn{W_{i:N}},
#' which are the expected values of the order statistics from a random sample
#' of size \eqn{N} of the random variable \eqn{W=C(U_1,U_2)} under the null
#' hypothesis of independence between \eqn{U_1} and \eqn{U_2}. \eqn{W_{i:N}}
#' can be calculated as follows \deqn{ W_{i:n}= N {N-1 \choose i-1}
#' \int\limits_{0}^1 \omega k_0(\omega) ( K_0(\omega) )^{i-1} ( 1-K_0(\omega)
#' )^{N-i} d\omega, } where
#' \deqn{K_0(\omega) = \omega - \omega \log(\omega), }{
#' K_=(\omega) = \omega - \omega log(\omega) }
#' and \eqn{k_0(\cdot)}{k_0()} is the corresponding density.
#'
#' K-plots can be seen as the bivariate copula equivalent to QQ-plots. If the
#' points of a K-plot lie approximately on the diagonal \eqn{y=x}, then
#' \eqn{U_1} and \eqn{U_2} are approximately independent. Any deviation from
#' the diagonal line points towards dependence. In case of positive dependence,
#' the points of the K-plot should be located above the diagonal line, and vice
#' versa for negative dependence. The larger the deviation from the diagonal,
#' the stronger is the degree of dependency. There is a perfect positive
#' dependence if points \eqn{\left(W_{i:N},H_i\right)} lie on the curve
#' \eqn{K_0(\omega)} located above the main diagonal. If points
#' \eqn{\left(W_{i:N},H_i\right)}{(W_{i:N},H_i)} however lie on the x-axis,
#' this indicates a perfect negative dependence between \eqn{U_1} and
#' \eqn{U_2}.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @param PLOT Logical; whether the results are plotted. If \code{PLOT =
#' FALSE}, the values \code{W.in} and \code{Hi.sort} are returned (see below;
#' default: \code{PLOT = TRUE}).
#' @param ... Additional plot arguments.
#' @return \item{W.in}{W-statistics (x-axis).} \item{Hi.sort}{H-statistics
#' (y-axis).}
#' @author Natalia Belgorodski, Ulf Schepsmeier
#' @seealso \code{\link{BiCopMetaContour}}, \code{\link{BiCopChiPlot}},
#' \code{\link{BiCopLambda}}, \code{\link{BiCopGofTest}}
#' @references Genest, C. and A. C. Favre (2007). Everything you always wanted
#' to know about copula modeling but were afraid to ask. Journal of Hydrologic
#' Engineering, 12 (4), 347-368.
#' @examples
#'
#' ## Gaussian and Clayton copulas
#' n <- 500
#' tau <- 0.5
#'
#' # simulate from Gaussian copula
#' fam <- 1
#' par <- BiCopTau2Par(fam, tau)
#' cop1 <- BiCop(fam, par)
#' set.seed(123)
#' dat1 <- BiCopSim(n, cop1)
#'
#' # simulate from Clayton copula
#' fam <- 3
#' par <- BiCopTau2Par(fam, tau)
#' cop2 <- BiCop(fam, par)
#' set.seed(123)
#' dat2 <- BiCopSim(n, cop2)
#'
#' # create K-plots
#' op <- par(mfrow = c(1, 2))
#' BiCopKPlot(dat1[,1], dat1[,2], main = "Gaussian copula")
#' BiCopKPlot(dat2[,1], dat2[,2], main = "Clayton copula")
#' par(op)
#'
#' @export BiCopKPlot
BiCopKPlot <- function(u1, u2, PLOT = TRUE, ...) {
    stopifnot(is.logical(PLOT))

    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    remove_nas,
                    check_nobs,
                    check_if_01,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    # Computations
    Wi <- Wn(u1, u2)
    Hi <- Hn(u1, u2)
    Hi.sort <- sort(Hi)
    n <- length(u1)
    W.in <- rep(NA, n)
    for (i in 1:n) {
        # function to be integrated
        f <- function(w) {
            w * (-log(w)) * (w - w * log(w))^(i - 1) * (1 - w + w * log(w))^(n - i)
        }
        W.in[i] <- n * choose(n - 1, i - 1) * (integrate(f, lower = 0, upper = 1)$value)
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
Fn <- function(t) {
    # help function for chi-plot
    n <- length(t)
    result <- (rank(t) - 1)/(n - 1)
    return(result)
}
# -------------------------------------------------------------------------------
Hn <- function(t, s) {
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
Wn <- function(t, s) {
    # Help function for kendall.plot
    n <- length(t)
    result <- ((n - 1) * Hn(t, s) + 1)/n
    return(result)
}
# -------------------------------------------------------------------------------
