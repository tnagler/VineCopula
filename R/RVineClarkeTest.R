#' Clarke Test Comparing Two R-Vine Copula Models
#' 
#' This function performs a Clarke test between two d-dimensional R-vine copula
#' models as specified by their \code{\link{RVineMatrix}} objects.
#' 
#' The test proposed by Clarke (2007) allows to compare non-nested models. For
#' this let \eqn{c_1} and \eqn{c_2} be two competing vine copulas in terms of
#' their densities and with estimated parameter sets
#' \eqn{\hat{\boldsymbol{\theta}}_1}{\theta_1} and
#' \eqn{\hat{\boldsymbol{\theta}}_2}{\theta_2}. The null hypothesis of
#' statistical indistinguishability of the two models is \deqn{ }{ H_0: P(m_i >
#' 0) = 0.5 forall i=1,..,N, }\deqn{H_0: P(m_i > 0) = 0.5\ \forall i=1,..,N, }{
#' H_0: P(m_i > 0) = 0.5 forall i=1,..,N, } where
#' \eqn{m_i:=\log\left[\frac{c_1(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_1)}{c_2(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_2)}\right]}{m_i:=log[
#' c_1(u_i|\theta_1) / c_2(u_i|\theta_2) ]} for observations
#' \eqn{\boldsymbol{u}_i,\ i=1,...,N}{u_i, i=1,...,N}.
#' 
#' Since under statistical equivalence of the two models the log likelihood
#' ratios of the single observations are uniformly distributed around zero and
#' in expectation \eqn{50\%} of the log likelihood ratios greater than zero,
#' the tets statistic \deqn{ }{ statistic := B = \sum_{i=1}^N
#' 1_{(0,\infty)}(m_i), }\deqn{\texttt{statistic} := B = \sum_{i=1}^N
#' \mathbf{1}_{(0,\infty)}(m_i), }{ statistic := B = \sum_{i=1}^N
#' 1_{(0,\infty)}(m_i), } where \eqn{\mathbf{1}}{1} is the indicator function,
#' is distributed Binomial with parameters \eqn{N} and \eqn{p=0.5}, and
#' critical values can easily be obtained. Model 1 is interpreted as
#' statistically equivalent to model 2 if \eqn{B} is not significantly
#' different from the expected value \eqn{Np = \frac{N}{2}}{np=N/2}.
#' 
#' Like AIC and BIC, the Clarke test statistic may be corrected for the number
#' of parameters used in the models. There are two possible corrections; the
#' Akaike and the Schwarz corrections, which correspond to the penalty terms in
#' the AIC and the BIC, respectively.
#' 
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM1,RVM2 \code{\link{RVineMatrix}} objects of models 1 and 2.
#' @return \item{statistic, statistic.Akaike, statistic.Schwarz}{Test
#' statistics without correction, with Akaike correction and with Schwarz
#' correction.} \item{p.value, p.value.Akaike, p.value.Schwarz}{P-values of
#' tests without correction, with Akaike correction and with Schwarz
#' correction.}
#' @author Jeffrey Dissmann, Eike Brechmann
#' @seealso \code{\link{RVineVuongTest}}, \code{\link{RVineAIC}},
#' \code{\link{RVineBIC}}
#' @references Clarke, K. A. (2007). A Simple Distribution-Free Test for
#' Nonnested Model Selection. Political Analysis, 15, 347-363.
#' @examples
#' 
#' \dontrun{
#' # vine structure selection time-consuming (~ 20 sec)
#' 
#' # load data set
#' data(daxreturns)
#' 
#' # select the R-vine structure, families and parameters
#' RVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6))
#' RVM$Matrix
#' RVM$par
#' RVM$par2
#' 
#' # select the C-vine structure, families and parameters
#' CVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6), type = "CVine")
#' CVM$Matrix
#' CVM$par
#' CVM$par2
#' 
#' # compare the two models based on the data
#' clarke <- RVineClarkeTest(daxreturns[,1:5], RVM, CVM)
#' clarke$statistic
#' clarke$statistic.Schwarz
#' clarke$p.value
#' clarke$p.value.Schwarz
#' }
#' 
#' @export RVineClarkeTest
RVineClarkeTest <- function(data, RVM1, RVM2) {
    
    N <- dim(data)[1]
    
    if (dim(data)[2] < 2) 
        stop("Dimension has to be at least 2.")
    if (N < 2) 
        stop("Number of observations has to be at least 2.")
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    if (is(RVM1)[1] != "RVineMatrix") 
        stop("'RVM1' has to be an RVineMatrix object.")
    if (is(RVM2)[1] != "RVineMatrix") 
        stop("'RVM2' has to be an RVineMatrix object.")
    
    Model1.ll <- RVineLogLik(data, RVM1, separate = TRUE)$loglik
    Model2.ll <- RVineLogLik(data, RVM2, separate = TRUE)$loglik
    
    anz.1 <- sum(RVM1$family >= 1, na.rm = TRUE) + sum(RVM1$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    anz.2 <- sum(RVM2$family >= 1, na.rm = TRUE) + sum(RVM2$family %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    
    B <- sum(Model1.ll - Model2.ll > 0)
    B.Schwarz <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2) * log(N)/(2 * N) > 0)
    B.Akaike <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2)/N > 0)
    
    if (B == 0 | B == N/2) {
        p <- 1 
    } else { 
        p <- 2 * min(pbinom(B, N, 0.5), 1 - pbinom(B - 1, N, 0.5))
    }
    
    if (B.Schwarz == 0 | B.Schwarz == N/2) {
        p.Schwarz <- 1
    } else {
        p.Schwarz <- 2 * min(pbinom(B.Schwarz, N, 0.5), 1 - pbinom(B.Schwarz - 1, N, 0.5))
    }
    
    if (B.Akaike == 0 | B.Akaike == N/2) {
        p.Akaike <- 1 
    }   else {
        p.Akaike <- 2 * min(pbinom(B.Akaike, N, 0.5), 1 - pbinom(B.Akaike - 1, N, 0.5))
    }
    
    return(list(statistic = B, 
                statistic.Akaike = B.Akaike,
                statistic.Schwarz = B.Schwarz,
                p.value = p, 
                p.value.Akaike = p.Akaike,
                p.value.Schwarz = p.Schwarz))
}
