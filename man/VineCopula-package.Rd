\name{VineCopula-package}
\alias{VineCopula-package}        
\alias{VineCopula}
\docType{package}

\title{Statistical Inference of Vine Copulas}

\description{
Tools for bivariate exploratory data analysis, bivariate copula selection and (vine) tree construction are provided.
Vine copula models can be estimated either sequentially or by joint maximum likelihood estimation.
Sampling algorithms and plotting methods are included.
Data is assumed to lie in the unit hypercube (so-called copula data).
For C- and D-vines links to the package 'CDVine' are provided.
}

\section{Remark}{
The package \code{VineCopula} is a continuation of the package \code{CDVine} by U. Schepsmeier and E. C. Brechmann (see Brechmann and Schepsmeier (2013)).
It includes all functions implemented in CDVine for the bivariate case (BiCop-functions).
}

\section{Bivariate copula families}{
In this package several bivariate copula families are included for bivariate analysis as well as for multivariate analysis using vine copulas.
It provides functionality of elliptical (Gaussian and Student-t) as well as Archimedean (Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8) copulas
to cover a large bandwidth of possible dependence structures.
For the Archimedean copula families rotated versions are included to cover negative dependence too.
The two parameter BB1, BB6, BB7 and BB8 copulas are however numerically instable for large parameters,
in particular, if BB6, BB7 and BB8 copulas are close to the Joe copula which is a boundary case of these three copula families.
In general, the user should be careful with extreme parameter choices.

As an asymmetric extension of the Gumbel copula, the Tawn copula with three parameters is also included in the package.
Both the Gumbel and the Tawn copula are extreme-value copulas, which can be defined in terms of their corresponding Pickands dependence functions.
For simplicity, we implemented two versions of the Tawn copula with two parameters each.
Each type has one of the asymmetry parameters fixed to 1, so that the corresponding Pickands dependence is either left- or right-skewed. In the manual we will call these two new copulas "Tawn type 1" and "Tawn type 2".

The following table shows the parameter ranges of bivariate copula families with parameters \code{par} and \code{par2}:
\tabular{lll}{
Copula family \tab \code{par} \tab \code{par2} \cr
Gaussian \tab \eqn{(-1,1)} \tab - \cr
Student t \tab \eqn{(-1,1)} \tab \eqn{(2,\infty)} \cr
(Survival) Clayton \tab \eqn{(0,\infty)} \tab - \cr
(Survival) Gumbel \tab \eqn{[1,\infty)} \tab - \cr
Frank \tab \eqn{R\backslash\{0\}} \tab - \cr
(Survival) Joe \tab \eqn{(1,\infty)} \tab - \cr
Rotated Clayton (90 and 270 degrees) \tab \eqn{(-\infty,0)} \tab - \cr
Rotated Gumbel (90 and 270 degrees)\tab \eqn{(-\infty,-1]} \tab - \cr
Rotated Joe (90 and 270 degrees) \tab \eqn{(-\infty,-1)} \tab - \cr
(Survival) Clayton-Gumbel (BB1) \tab \eqn{(0,\infty)} \tab \eqn{[1,\infty)} \cr
(Survival) Joe-Gumbel (BB6) \tab \eqn{[1,\infty)} \tab \eqn{[1,\infty)} \cr
(Survival) Joe-Clayton (BB7) \tab \eqn{[1,\infty)} \tab \eqn{(0,\infty)} \cr
(Survival) Joe-Frank (BB8) \tab \eqn{[1,\infty)} \tab \eqn{(0,1]} \cr
Rotated Clayton-Gumbel (90 and 270 degrees) \tab \eqn{(-\infty,0)} \tab \eqn{(-\infty,-1]} \cr
Rotated Joe-Gumbel (90 and 270 degrees) \tab \eqn{(-\infty,-1]} \tab \eqn{(-\infty,-1]} \cr
Rotated Joe-Clayton (90 and 270 degrees) \tab \eqn{(-\infty,-1]} \tab \eqn{(-\infty,0)} \cr
Rotated Joe-Frank (90 and 270 degrees) \tab \eqn{(-\infty,-1]} \tab \eqn{[-1,0)} \cr
(Survival) Tawn type 1 and type 2 \tab \eqn{[1,\infty)} \tab \eqn{[0,1]} \cr
Rotated Tawn type 1 and type 2 (90 and 270 degrees) \tab \eqn{(-\infty,-1]} \tab \eqn{[0,1]} \cr
}  
}

\section{R-vine copula models}{
The specification of an R-vine is done in matrix notation, introduced by Dissmann et al. (2013). One matrix contains the R-vine tree structure,
one the copula families utilized and two matrices corresponding parameter values. 
These four matrices are stored in an \code{RVineMatrix} object created by the function \code{RVineMatrix}. Each matrix is a d x d lower triangular matrix.
Since C- and D-vines are special cases, boundary cases, of R-vines one can write each C- or D-vine in R-vine notation. The transformation
of notation to an R-vine can be done via \code{C2RVine} and \code{D2RVine}, which provide an interface to the package CDVine.
For more details see the documentation of the functions.
}


\section{Acknowledgment}{
We acknowledge substantial contributions by our working group at Technische Universitaet Muenchen,
in particular by Carlos Almeida and Aleksey Min.
In addition, we like to thank Shing (Eric) Fu, Feng Zhu, Guang (Jack) Yang, and Harry Joe for providing their implementation
of the method by Knight (1966) for efficiently computing the empirical Kendall's tau.
We are especially grateful to Harry Joe for his contributions to the implementation of the bivariate Archimedean copulas.
}

\details{
\tabular{ll}{
Package: \tab VineCopula\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2015-07-30\cr
License: \tab GPL (>=2)\cr
Depends: \tab R (\eqn{\geq 2.11.0}{>= 2.11.0})\cr
Imports: \tab graphics, grDevices, stats, utils, MASS, mvtnorm, network, methods, copula, ADGofTest, lattice\cr
Suggests: \tab CDVine, TSP\cr
LazyLoad: \tab yes
}
}

\author{
Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler, Thomas Nagler, Tobias Erhardt
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.

Bedford, T. and R. M. Cooke (2001).
Probability density decomposition for conditionally dependent random variables modeled by vines.
Annals of Mathematics and Artificial intelligence 32, 245-268.

Bedford, T. and R. M. Cooke (2002).
Vines - a new graphical model for dependent random variables.
Annals of Statistics 30, 1031-1068.

Brechmann, E. C., C. Czado, and K. Aas (2012).
Truncated regular vines in high dimensions with applications to financial data.
Canadian Journal of Statistics 40 (1), 68-85.

Brechmann, E. C. and C. Czado (2011).
Risk management with high-dimensional vine copulas: An analysis of the Euro Stoxx 50.
Statistics & Risk Modeling, 30 (4), 307-342.

Brechmann, E. C. and U. Schepsmeier (2013).
Modeling Dependence with C- and D-Vine Copulas: The R Package CDVine.
Journal of Statistical Software, 52 (3), 1-27.
\url{http://www.jstatsoft.org/v52/i03/}.

Czado, C., U. Schepsmeier, and A. Min (2012).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Statistical Modelling, 12(3), 229-255.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.

Eschenburg, P. (2013).
Properties of extreme-value copulas
Diploma thesis, Technische Universitaet Muenchen
\url{http://mediatum.ub.tum.de/node?id=1145695}

Joe, H. (1996).
Families of m-variate distributions with given margins and m(m-1)/2 bivariate dependence parameters.
In L. Rueschendorf, B. Schweizer, and M. D. Taylor (Eds.), Distributions with fixed marginals and related topics, pp. 120-141.
Hayward: Institute of Mathematical Statistics.

Joe, H. (1997).
Multivariate Models and Dependence Concepts.
London: Chapman and Hall.

Knight, W. R. (1966).
A computer method for calculating Kendall's tau with ungrouped data.
Journal of the American Statistical Association 61 (314), 436-439.

Kurowicka, D. and R. M. Cooke (2006).
Uncertainty Analysis with High Dimensional Dependence Modelling.
Chichester: John Wiley.

Kurowicka, D. and H. Joe (Eds.) (2011).
Dependence Modeling: Vine Copula Handbook.
Singapore: World Scientific Publishing Co.

Nelsen, R. (2006). 
An introduction to copulas. 
Springer

Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Statistical Papers, 55 (2), 525-542. \cr
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.

Schepsmeier, U. (2013)
A goodness-of-fit test for regular vine copula models. 
Preprint
\url{http://arxiv.org/abs/1306.0818}

Schepsmeier, U. (2015)
Efficient information based goodness-of-fit tests for vine copula models with fixed margins.
Journal of Multivariate Analysis 138, 34-52.

Stoeber, J. and U. Schepsmeier (2013).
Estimating standard errors in regular vine copula models. Computational Statistics, 28 (6), 2679-2707 \cr
\url{http://link.springer.com/article/10.1007/s00180-013-0423-8#}.

White, H. (1982)
Maximum likelihood estimation of misspecified models,
Econometrica, 50, 1-26.
}
