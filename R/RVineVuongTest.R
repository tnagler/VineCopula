#' Vuong Test Comparing Two R-Vine Copula Models
#'
#' This function performs a Vuong test between two d-dimensional R-vine copula
#' models as specified by their \code{\link{RVineMatrix}} objects.
#'
#' The likelihood-ratio based test proposed by Vuong (1989) can be used for
#' comparing non-nested models. For this let \eqn{c_1} and \eqn{c_2} be two
#' competing vine copulas in terms of their densities and with estimated
#' parameter sets \eqn{\hat{\boldsymbol{\theta}}_1}{\theta_1} and
#' \eqn{\hat{\boldsymbol{\theta}}_2}{\theta_2}. We then compute the
#' standardized sum, \eqn{\nu}, of the log differences of their pointwise
#' likelihoods
#' \eqn{m_i:=\log\left[\frac{c_1(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_1)}{c_2(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_2)}\right]}{m_i:=log[c_1(u_i|\theta_1)
#' / c_2(u_i|\theta_2) ]} for observations \eqn{\boldsymbol{u}_i\in[0,1],\
#' i=1,...,N}{u_i in [0,1],i=1,...,N} , i.e.,
#' \deqn{\texttt{statistic} := \nu = \frac{\frac{1}{n}\sum_{i=1}^N
#' m_i}{\sqrt{\sum_{i=1}^N\left(m_i - \bar{m} \right)^2}}. }{ statistic := \nu
#' = (1/n\sum_{i=1}^N m_i) / ((\sum_{i=1}^N (m_i - \bar{m} )^2)^0.5). } Vuong
#' (1989) shows that \eqn{\nu} is asymptotically standard normal. According to
#' the null-hypothesis \deqn{H_0:
#' E[m_i] = 0\ \forall i=1,...,N, }{ H_0: E[m_i] = 0 forall i=1,...,N, } we
#' hence prefer vine model 1 to vine model 2 at level \eqn{\alpha} if
#' \deqn{\nu>\Phi^{-1}\left(1-\frac{\alpha}{2}\right), }{ \nu >
#' \Phi^{-1}(1-\alpha/2), } where \eqn{\Phi^{-1}} denotes the inverse of the
#' standard normal distribution function. If
#' \eqn{\nu<-\Phi^{-1}\left(1-\frac{\alpha}{2}\right)}{\nu<-\Phi^{-1}(1-\alpha/2)}
#' we choose model 2.  If, however,
#' \eqn{|\nu|\leq\Phi^{-1}\left(1-\frac{\alpha}{2}\right)}{|\nu| <=
#' \Phi^{-1}(1-\alpha/2)}, no decision among the models is possible.
#'
#' Like AIC and BIC, the Vuong test statistic may be corrected for the number
#' of parameters used in the models. There are two possible corrections; the
#' Akaike and the Schwarz corrections, which correspond to the penalty terms in
#' the AIC and the BIC, respectively.
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM1,RVM2 \code{\link{RVineMatrix}} objects of models 1 and 2.
#'
#' @return \item{statistic, statistic.Akaike, statistic.Schwarz}{Test
#' statistics without correction, with Akaike correction and with Schwarz
#' correction.} \item{p.value, p.value.Akaike, p.value.Schwarz}{P-values of
#' tests without correction, with Akaike correction and with Schwarz
#' correction.}
#'
#' @author Jeffrey Dissmann, Eike Brechmann
#'
#' @seealso \code{\link{RVineClarkeTest}}, \code{\link{RVineAIC}},
#' \code{\link{RVineBIC}}
#'
#' @references Vuong, Q. H. (1989). Ratio tests for model selection and
#' non-nested hypotheses. Econometrica 57 (2), 307-333.
#'
#' @examples
#' \donttest{
#' # vine structure selection time-consuming (~ 20 sec)
#'
#' # load data set
#' data(daxreturns)
#'
#' # select the R-vine structure, families and parameters
#' RVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6))
#'
#' # select the C-vine structure, families and parameters
#' CVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6), type = "CVine")
#'
#' # compare the two models based on the data
#' vuong <- RVineVuongTest(daxreturns[,1:5], RVM, CVM)
#' vuong$statistic
#' vuong$statistic.Schwarz
#' vuong$p.value
#' vuong$p.value.Schwarz
#' }
#'
#' @export RVineVuongTest
RVineVuongTest <- function(data, RVM1, RVM2) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_if_01,
                    check_nobs,
                    check_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())
    N <- args$n

    Model1.ll <- RVineLogLik(data, RVM1, separate = TRUE)$loglik
    Model2.ll <- RVineLogLik(data, RVM2, separate = TRUE)$loglik

    anz.1 <- sum(RVM1$family >= 1, na.rm = TRUE) + sum(RVM1$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)
    anz.2 <- sum(RVM2$family >= 1, na.rm = TRUE) + sum(RVM2$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234), na.rm = TRUE)

    if (all(Model1.ll - Model2.ll == 0)) {
        # models are the same
        V <- 0
        V.Schwarz <- 0
        V.Akaike <- 0

        p <- 1
        p.Schwarz <- 1
        p.Akaike <- 1
    } else {
        w <- 1/N * sum((Model1.ll - Model2.ll)^2) + (1/N * sum(Model1.ll - Model2.ll))^2
        w <- sqrt(w)

        LR <- sum(Model1.ll) - sum(Model2.ll)
        LR.Schwarz <- LR - ((anz.1/2 * log(N) - anz.2/2 * log(N)))
        LR.Akaike <- LR - (anz.1 - anz.2)

        V <- LR/(sqrt(N) * w)
        V.Schwarz <- LR.Schwarz/(sqrt(N) * w)
        V.Akaike <- LR.Akaike/(sqrt(N) * w)

        p <- 2 * min(pnorm(V), 1 - pnorm(V))
        p.Schwarz <- 2 * min(pnorm(V.Schwarz), 1 - pnorm(V.Schwarz))
        p.Akaike <- 2 * min(pnorm(V.Akaike), 1 - pnorm(V.Akaike))
    }

    return(list(statistic = V,
                statistic.Akaike = V.Akaike,
                statistic.Schwarz = V.Schwarz,
                p.value = p,
                p.value.Akaike = p.Akaike,
                p.value.Schwarz = p.Schwarz))
}
