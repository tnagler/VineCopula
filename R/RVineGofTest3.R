#' Goodness-of-Fit Tests for R-Vine Copula Models
#'
#' This function performs a goodness-of-fit test for R-vine copula models.
#' There are 15 different goodness-of-fit tests implemented, described in
#' Schepsmeier (2013).
#'
#' \code{method = "White"}: \cr
#' This goodness-of fit test uses the information
#' matrix equality of White (1982) and was original investigated by Huang and
#' Prokhorov (2011) for copulas. \cr
#' Schepsmeier (2012) enhanced their approach
#' to the vine copula case. \cr
#' The main contribution is that under correct
#' model specification the Fisher Information can be equivalently calculated as
#' minus the expected Hessian matrix or as the expected outer product of the
#' score function.
#' The null hypothesis is
#' \deqn{ H_0: \boldsymbol{H}(\theta) + \boldsymbol{C}(\theta) = 0 }{H_0:
#' H(\theta) + C(\theta) = 0 }
#' against the alternative
#' \deqn{ H_1: \boldsymbol{H}(\theta) +
#' \boldsymbol{C}(\theta) \neq 0 , }{H_1: H(\theta) + C(\theta) != 0, }
#' where
#' \eqn{\boldsymbol{H}(\theta)}{H(\theta)} is the expected Hessian matrix and
#' \eqn{\boldsymbol{C}(\theta)}{C(\theta)} is the expected outer product of the
#' score function. \cr
#' For the calculation of the test statistic we use the
#' consistent maximum likelihood estimator \eqn{\hat{\theta}} and the sample
#' counter parts of \eqn{\boldsymbol{H}(\theta)}{H(\theta)} and
#' \eqn{\boldsymbol{C}(\theta)}{C(\theta)}. \cr
#' The correction of the
#' Covariance-Matrix in the test statistic for the uncertainty in the margins
#' is skipped. The implemented test assumes that there is no uncertainty in the
#' margins. The correction can be found in Huang and Prokhorov (2011) for
#' bivariate copulas and in Schepsmeier (2013) for vine copulas. It involves
#' multi-dimensional integrals. \cr
#'
#' \code{method = "IR"}: \cr
#' As the White test the information matrix ratio
#' test is based on the expected Hessian matrix
#' \eqn{\boldsymbol{H}(\theta)}{H(\theta)}
#' and the expected outer product of
#' the score function \eqn{\boldsymbol{C}(\theta)}{C(\theta)}. \cr
#' \deqn{ H_0:-\boldsymbol{H}(\theta)^{-1}\boldsymbol{C}(\theta) = I_{p} }{H_0:
#' -H(\theta)^{-1}C(\theta) = I_p}
#' against the alternative
#' \deqn{ H_1:
#' -\boldsymbol{H}(\theta)^{-1}\boldsymbol{C}(\theta) \neq I_{p} . }{H_1:
#' -H(\theta)^{-1}C(\theta) != I_p.}
#' The test statistic can then be calculated as
#' \deqn{ IR_n:=tr(\Phi(\theta))/p } with
#' \eqn{\Phi(\theta)=-\boldsymbol{H}(\theta)^{-1}\boldsymbol{C}(\theta)}{\Phi(\theta)=-H(\theta)^{-1}C(\theta)},
#' \eqn{p} is the number of parameters, i.e. the length of \eqn{\theta}, and
#' \eqn{tr(A)} is the trace of the matrix \eqn{A} \cr
#' For details see Schepsmeier (2013) \cr
#'
#' \code{method = "Breymann"}, \code{method = "Berg"} and \code{method =
#' "Berg2"}: \cr
#' These tests are based on the multivariate probability integral
#' transform (PIT) applied in \code{\link{RVinePIT}}. The multivariate data
#' \eqn{y_{i}} returned form the PIT are aggregated to univariate data by
#' different aggregation functions \eqn{\Gamma(\cdot)} in the sum \deqn{
#' s_t=\sum_{i=1}^d \Gamma(y_{it}), t=1,...,n }.
#' In Breymann et al. (2003) the weight function is suggested as
#' \eqn{\Gamma(\cdot)=\Phi^{-1}(\cdot)^2}{\Gamma(y)=\Phi^{-1}y^2}, while in
#' Berg and Bakken (2007) the weight function is either
#' \eqn{\Gamma(\cdot)=|\cdot-0.5|}{\Gamma(y)=|y-0.5|} (\code{method="Berg"}) or
#' \eqn{\Gamma(\cdot)=(\cdot-0.5)^{\alpha},\alpha=2,4,6,...}{\Gamma(y)=(y-0.5)^{\alpha},\alpha=2,4,6,...}
#' (\code{method="Berg2"}). \cr Furthermore, the \code{"Berg"} and
#' \code{"Berg2"} test are based on the order statistics of the PIT returns.
#' \cr See Berg and Bakken (2007) or Schepsmeier (2013) for details. \cr
#'
#' \code{method = "ECP"} and \code{method = "ECP2"}: \cr
#' Both tests are test
#' for \eqn{H_0: C \in C_0} against \eqn{H_1: C \notin C_0} where C denotes the
#' (vine) copula distribution function and \eqn{C_0} is a class of parametric
#' (vine) copulas with \eqn{\Theta\subseteq R^p} being the parameter space of
#' dimension p. They are based on the empirical copula process (ECP)
#' \deqn{
#' \hat{C}_n(u)-C_{\hat{\theta}_n}(u), }
#' with
#' \eqn{u=(u_1,\ldots,u_d)\in[0,1]^d} and
#' \eqn{\hat{C}_n(u) =
#' \frac{1}{n+1}\sum_{t=1}^n \boldsymbol{1}_{\{U_{t1}\leq u_1,\ldots,U_{td}\leq
#' u_d \}} }{\hat{C}_n(u) = 1/(n+1)\sum_{t=1}^n 1_{U_{t1}\leq
#' u_1,\ldots,U_{td}\leq u_d} }.
#' The ECP is utilized in a multivariate
#' Cramer-von Mises (CvM) or multivariate Kolmogorov-Smirnov (KS) based test
#' statistic. An extension of the ECP-test is the combination of the
#' multivariate PIT approach with the ECP. The general idea is that the
#' transformed data of a multivariate PIT should be "close" to the independence
#' copula Genest et al. (2009). Thus a distance of CvM or KS type between them
#' is considered. This approach is called ECP2. Again we refer to Schepsmeier
#' (2013) for details.
#'
#' @param data An N x d data matrix (with uniform margins).
#' @param RVM \code{\link{RVineMatrix}} objects of the R-vine model under the
#' null hypothesis. \cr
#' Only the following copula families are allowed in
#' \code{RVM$family} due to restrictions in \code{\link{RVineGrad}} and
#' \code{\link{RVineHessian}} \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula)\cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' @param method A string indicating the goodness-of-fit method:\cr
#' \code{"White"} = goodness-of-fit test based on White's information matrix
#' equality (default) \cr
#' \code{"IR"} = goodness-of-fit test based on the
#' information ratio \cr
#' \code{"Breymann"} = goodness-of-fit test based on the
#' probability integral transform (PIT) and the aggregation to univariate data
#' by Breymann et al. (2003). \cr
#' \code{"Berg"} = goodness-of-fit test based on
#' the probability integral transform (PIT) and the aggregation to univariate
#' data by Berg and Bakken (2007). \cr
#' \code{"Berg2"} = second goodness-of-fit
#' test based on the probability integral transform (PIT) and the aggregation
#' to univariate data by Berg and Bakken (2007). \cr
#' \code{"ECP"} =
#' goodness-of-fit test based on the empirical copula process (ECP) \cr
#' \code{"ECP2"} = goodness-of-fit test based on the combination of probability
#' integral transform (PIT) and empirical copula process (ECP) (Genest et al.
#' 2009) \cr
#' @param statistic A string indicating the goodness-of-fit test statistic
#' type:\cr \code{"CvM"} = Cramer-von Mises test statistic (univariate for
#' \code{"Breymann"}, \code{"Berg"} and \code{"Berg2"}, multivariate for
#' \code{"ECP"} and \code{"ECP2"}) \cr \code{"KS"} = Kolmogorov-Smirnov test
#' statistic (univariate for \code{"Breymann"}, \code{"Berg"} and
#' \code{"Berg2"}, multivariate for \code{"ECP"} and \code{"ECP2"}) \cr
#' \code{"AD"} = Anderson-Darling test statistic (only univariate for
#' \code{"Breymann"}, \code{"Berg"} and \code{"Berg2"})
#' @param B an integer for the number of bootstrap steps (default \code{B =
#' 200})\cr
#' For \code{B = 0} the asymptotic p-value is returned if available,
#' otherwise only the test statistic is returned.\cr
#' WARNING: If \code{B} is chosen too large, computations will take very long.
#' @param alpha an integer of the set \code{2,4,6,...} for the \code{"Berg2"}
#' goodness-of-fit test (default \code{alpha = 2})
#'
#' @return For \code{method = "White"}:
#' \item{White}{test statistic}
#' \item{p.value}{p-value, either asymptotic for \code{B = 0} or bootstrapped
#' for \code{B > 0}}
#' For \code{method = "IR"}:
#' \item{IR}{test statistic (raw version as stated above)}
#' \item{p.value}{So far no p-value is returned nigher a asymptotic nor a
#' bootstrapped one. How to calculated a bootstrapped p-value is explained in
#' Schepsmeier (2013). Be aware, that the test statistics than have to be adjusted
#' with the empirical variance.}
#' For \code{method = "Breymann"}, \code{method = "Berg"}
#' and \code{method = "Berg2"}:
#' \item{CvM, KS, AD}{test statistic according to
#' the choice of \code{statistic}}
#' \item{p.value}{p-value, either asymptotic
#' for \code{B = 0} or bootstrapped for \code{B > 0}.  A asymptotic p-value is
#' only available for the Anderson-Darling test statistic if the R-package
#' \code{ADGofTest} is loaded. \cr
#' Furthermore, a asymptotic p-value can be
#' calculated for the Kolmogorov-Smirnov test statistic. For the Cramer-von
#' Mises no asymptotic p-value is available so far.}
#' For \code{method = "ECP"} and \code{method = "ECP2"}:
#' \item{CvM, KS}{test statistic according to the
#' choice of \code{statistic}}
#' \item{p.value}{bootstrapped p-value} \cr
#' Warning: The code for all the p-values are not yet approved since some of them are
#' moved from R-code to C-code. If you need p-values the best way is to write your own
#' algorithm as suggested in Schepsmeier (2013) to get bootstrapped p-values.
#'
#' @author Ulf Schepsmeier
#'
#' @seealso \code{\link{BiCopGofTest}}, \code{\link{RVinePIT}}
#'
#' @references Berg, D. and H. Bakken (2007) A copula goodness-of-fit apprach
#' based on the conditional probability integral transformation.
#' \url{http://www.danielberg.no/publications/Btest.pdf}
#'
#' Breymann, W., A. Dias and P. Embrechts (2003) Dependence structures for
#' multivariate high-frequence data in finance. Quantitative Finance 3, 1-14
#'
#' Genest, C., B. Remillard, and D. Beaudoin (2009) Goodness-of-fit tests for
#' copulas: a review and power study.  Insur. Math. Econ. 44, 199-213.
#'
#' Huang, w. and A. Prokhorov (2011). A goodness-of-fit test for copulas. to
#' appear in Econometric Reviews
#'
#' Schepsmeier, U. (2013) A goodness-of-fit test for regular vine copula
#' models.  Preprint \url{http://arxiv.org/abs/1306.0818}
#'
#' Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for
#' vine copula models with fixed margins. Journal of Multivariate Analysis 138,
#' 34-52.
#'
#' White, H. (1982) Maximum likelihood estimation of misspecified models,
#' Econometrica, 50, 1-26.
#'
#' @examples
#' \donttest{## time-consuming example
#'
#' # load data set
#' data(daxreturns)
#'
#' # select the R-vine structure, families and parameters
#' RVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6))
#'
#' # White test with asymptotic p-value
#' RVineGofTest(daxreturns[,1:5], RVM, B = 0)
#'
#' # ECP2 test with Cramer-von-Mises test statistic and a bootstrap
#' # with 200 replications for the calculation of the p-value
#' RVineGofTest(daxreturns[,1:5], RVM, method = "ECP2",
#'              statistic = "CvM", B = 200)
#' }
#'
RVineGofTest <- function(data, RVM, method = "White", statistic = "CvM", B = 200, alpha = 2) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    remove_nas,
                    check_if_01,
                    check_nobs,
                    check_RVMs,
                    prep_RVMs,
                    na.txt = " Only complete observations are used.")
    list2env(args, environment())

    if (any(!(RVM$family %in% c(0, 1:6, 13, 14, 16, 23, 24, 26, 33, 34, 36))))
        stop("Copula family not implemented.")
    if (statistic == "Cramer-von Mises") {
        statistic <- "CvM"
    } else if (statistic == "Kolmogorov-Smirnov") {
        statistic <- "KS"
    } else if (statistic == "Anderson-Darling") {
        statistic <- "AD"
    }
    T <- dim(data)[1]
    d <- dim(data)[2]

    if (method == "White") {
        out <- gof_White(data, RVM, B)
    } else if (method == "Breymann" || method == "Berg" || method == "Berg2") {
        out <- gof_PIT(data, RVM, method, B, statistic, alpha)
    } else if (method == "ECP" || method == "ECP2") {
        if (statistic == "AD")
            stop("The Anderson-Darling statistic is not available for the empirical copula process based goodness-of-fit tests.")
        out <- gof_ECP(data, RVM, B, method, statistic)
    } else if (method == "IR") {
        out2 <- RVineHessian(data, RVM)
        C <- out2$der
        H <- out2$hessian
        p <- dim(C)[1]
        Z <- solve(-C, H)
        IR <- sum(diag(Z))/p

        # TODO: p-value via bootstrap (?)
        out <- list(IR = IR, pvalue = NULL)
    }

    return(out)
}
