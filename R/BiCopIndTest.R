#' Independence Test for Bivariate Copula Data
#'
#' This function returns the p-value of a bivariate asymptotic independence
#' test based on Kendall's tau.
#'
#' The test exploits the asymptotic normality of the test statistic
#' \deqn{\texttt{statistic} := T =
#' \sqrt{\frac{9N(N - 1)}{2(2N + 5)}} \times |\hat{\tau}|, }{
#' statistic := T = ( (9N(N-1)) / (2(2N+5)) )^0.5 * |\tau|, }
#' where \eqn{N} is the number of observations (length of \code{u1}) and
#' \eqn{\hat{\tau}} the empirical Kendall's tau of the data vectors \code{u1}
#' and \code{u2}. The p-value of the null hypothesis of bivariate independence
#' hence is asymptotically
#' \deqn{\texttt{p.value} = 2 \times \left(1 - \Phi\left(T\right)\right), }{
#' p.value = 2*(1-\Phi(T)), }
#' where \eqn{\Phi} is the standard normal distribution function.
#'
#' @param u1,u2 Data vectors of equal length with values in [0,1].
#' @return \item{statistic}{Test statistic of the independence test.}
#' \item{p.value}{P-value of the independence test.}
#' @author Jeffrey Dissmann
#' @seealso \code{\link{BiCopGofTest}}, \code{\link{BiCopPar2Tau}},
#' \code{\link{BiCopTau2Par}}, \code{\link{BiCopSelect}},\cr
#' \code{\link{RVineCopSelect}}, \code{\link{RVineStructureSelect}}
#' @references Genest, C. and A. C. Favre (2007). Everything you always wanted
#' to know about copula modeling but were afraid to ask. Journal of Hydrologic
#' Engineering, 12 (4), 347-368.
#' @examples
#'
#' ## Example 1: Gaussian copula with large dependence parameter
#' par <- 0.7
#' fam <- 1
#' cop <- BiCop(fam, par)
#' set.seed(123)
#' dat <- BiCopSim(500, cop)
#'
#' # perform the asymptotic independence test
#' BiCopIndTest(dat[,1], dat[,2])
#'
#'
#' ## Example 2: Gaussian copula with small dependence parameter
#' par <- 0.01
#' fam <- 1
#' cop <- BiCop(fam, par)
#' set.seed(123)
#' dat <- BiCopSim(500, cop)
#'
#' # perform the asymptotic independence test
#' BiCopIndTest(dat[,1], dat[,2])
#'
#' @export BiCopIndTest
BiCopIndTest <- function(u1, u2) {
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

    # tau = cor(u1,u2,method='kendall')
    tau <- fasttau(u1, u2)

    N <- length(u1)
    f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)

    return(list(statistic = f, p.value = 2 * (1 - pnorm(f))))
}
