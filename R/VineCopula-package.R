
#' Statistical Inference of Vine Copulas
#'
#' Vine copulas are a flexible class of dependence models consisting of
#' bivariate building blocks (see e.g., Aas et al., 2009). This package is
#' primarily made for the statistical analysis of vine copula
#' models. The package includes tools for parameter estimation, model selection,
#' simulation, goodness-of-fit tests, and visualization. Tools for estimation,
#' selection and exploratory data analysis of bivariate copula models are
#' also provided.
#'
#' \tabular{ll}{
#' Package: \tab VineCopula\cr
#' Type: \tab Package\cr
#' Version: \tab 2.0.1 \cr
#' Date: \tab 2016-06-09\cr
#' License: \tab GPL (>=2)\cr
#' Depends: \tab R (\eqn{\geq 2.11.0}{>= 2.11.0})\cr
#' Imports: \tab graphics, grDevices,
#' stats, utils, MASS, mvtnorm, network, methods, copula (>= 0.999-10),
#' kdecopula (>= 0.6.0), ADGofTest, lattice, doParallel, parallel, foreach \cr
#' Suggests: \tab CDVine, TSP, shiny\cr
#' LazyLoad: \tab yes }
#'
#' @name VineCopula-package
#' @aliases VineCopula-package VineCopula
#' @docType package
#'
#' @section Remark: The package \code{VineCopula} is a continuation of the
#' package \code{CDVine} by U. Schepsmeier and E. C. Brechmann (see Brechmann
#' and Schepsmeier (2013)). It includes all functions implemented in CDVine for
#' the bivariate case (BiCop-functions).
#'
#' @author Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt
#' Graeler, Thomas Nagler, Tobias Erhardt
#'
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
#' Pair-copula constructions of multiple dependence. Insurance: Mathematics and
#' Economics 44 (2), 182-198.
#'
#' Bedford, T. and R. M. Cooke (2001). Probability density decomposition for
#' conditionally dependent random variables modeled by vines. Annals of
#' Mathematics and Artificial intelligence 32, 245-268.
#'
#' Bedford, T. and R. M. Cooke (2002). Vines - a new graphical model for
#' dependent random variables. Annals of Statistics 30, 1031-1068.
#'
#' Brechmann, E. C., C. Czado, and K. Aas (2012). Truncated regular vines in
#' high dimensions with applications to financial data. Canadian Journal of
#' Statistics 40 (1), 68-85.
#'
#' Brechmann, E. C. and C. Czado (2011). Risk management with high-dimensional
#' vine copulas: An analysis of the Euro Stoxx 50. Statistics & Risk Modeling,
#' 30 (4), 307-342.
#'
#' Brechmann, E. C. and U. Schepsmeier (2013). Modeling Dependence with C- and
#' D-Vine Copulas: The R Package CDVine. Journal of Statistical Software, 52
#' (3), 1-27. \url{http://www.jstatsoft.org/v52/i03/}.
#'
#' Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation
#' of mixed C-vines with application to exchange rates. Statistical Modelling,
#' 12(3), 229-255.
#'
#' Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
#' Selecting and estimating regular vine copulae and application to financial
#' returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
#' Eschenburg, P. (2013). Properties of extreme-value copulas Diploma thesis,
#' Technische Universitaet Muenchen
#' \url{http://mediatum.ub.tum.de/node?id=1145695}
#'
#' Joe, H. (1996). Families of m-variate distributions with given margins and
#' m(m-1)/2 bivariate dependence parameters. In L. Rueschendorf, B. Schweizer,
#' and M. D. Taylor (Eds.), Distributions with fixed marginals and related
#' topics, pp. 120-141. Hayward: Institute of Mathematical Statistics.
#'
#' Joe, H. (1997). Multivariate Models and Dependence Concepts. London: Chapman
#' and Hall.
#'
#' Knight, W. R. (1966). A computer method for calculating Kendall's tau with
#' ungrouped data. Journal of the American Statistical Association 61 (314),
#' 436-439.
#'
#' Kurowicka, D. and R. M. Cooke (2006). Uncertainty Analysis with High
#' Dimensional Dependence Modelling. Chichester: John Wiley.
#'
#' Kurowicka, D. and H. Joe (Eds.) (2011). Dependence Modeling: Vine Copula
#' Handbook. Singapore: World Scientific Publishing Co.
#'
#' Nelsen, R. (2006).  An introduction to copulas.  Springer
#'
#' Schepsmeier, U. and J. Stoeber (2014). Derivatives and Fisher information of
#' bivariate copulas. Statistical Papers, 55 (2), 525-542. \cr
#' \url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.
#'
#' Schepsmeier, U. (2013) A goodness-of-fit test for regular vine copula
#' models.  Preprint \url{http://arxiv.org/abs/1306.0818}
#'
#' Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for
#' vine copula models with fixed margins. Journal of Multivariate Analysis 138,
#' 34-52.
#'
#' Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular
#' vine copula models. Computational Statistics, 28 (6), 2679-2707 \cr
#' \url{http://link.springer.com/article/10.1007/s00180-013-0423-8#}.
#'
#' White, H. (1982) Maximum likelihood estimation of misspecified models,
#' Econometrica, 50, 1-26.
NULL


#' Major German Stocks
#'
#' This data set contains transformed standardized residuals of daily log
#' returns of 15 major German stocks represented in the index DAX observed from
#' January 2005 to August 2009. Each time series is filtered using a GARCH(1,1)
#' model with Student t innovations.
#'
#'
#' @name daxreturns
#' @docType data
#' @format A data frame with 1158 observations on 15 variables. Column names
#' correspond to ticker symbols of the stocks.
#' @seealso \code{\link{RVineStructureSelect}}
#' @source Yahoo! Finance
#' @examples
#'
#' # load the data set
#' data(daxreturns)
#'
#' # compute the empirical Kendall's tau matrix
#' TauMatrix(daxreturns)
#'
NULL






