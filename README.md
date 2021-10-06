VineCopula
==========

[![R build status](https://github.com/tnagler/VineCopula/workflows/R-CMD-check/badge.svg)](https://github.com/tnagler/VineCopula)
[![CRAN version](https://www.r-pkg.org/badges/version/VineCopula)](https://cran.r-project.org/package=VineCopula) 
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/VineCopula)](https://cran.r-project.org/package=VineCopula)

Vine copulas are a flexible class of dependence models consisting of bivariate 
building blocks (see e.g., Aas et al., 2009). You can find a comprehensive 
list of publications and other materials on [vine-copula.org](https://www.statistics.ma.tum.de/en/research/vine-copula-models/).

This package is primarily made for the statistical analysis of **vine copula
models**. The package includes tools for parameter estimation, model selection,
simulation, goodness-of-fit tests, and visualization. Tools for estimation,
selection and exploratory data analysis of **bivariate copula** models are also
provided. Please see the [API documentation](https://tnagler.github.io/VineCopula/)
for a detailed description of all functions.


Status
-----------------

The library is no longer actively developed, but will continued to be maintained.
Check out the [rvinecopulib](https://github.com/vinecopulib/rvinecopulib) package
for an alternative with several benefits:  

  * a sleaker and more modern API,  
  
  * shorter runtimes, especially in high dimensions,  
  
  * nonparametric and multi-parameter families,  
  
  * ability to model discrete variables.  


Table of contents
-----------------

- [How to install](#how-to-install)
- [Package overview](#package-overview)
	- [Bivariate copula modeling: the BiCop-family](#bivariate-copula-modeling-the-bicop-family)
	- [Vine copula modeling: the RVine-family](#vine-copula-modeling-the-rvine-family)
	- [Additional features](#additional-features)
	- [Bivariate copula families](#bivariate-copula-families)
- [Associated shiny apps](#associated-shiny-apps)
- [References](#references)


------------------------------------------------------------------------


How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("VineCopula")
    ```

-   the latest development version:

    ``` r
    # install.packages("remotes")
    remotes::install_github("tnagler/VineCopula")
    ```

------------------------------------------------------------------------

Package overview
----------------

Below, we list most functions and features you should know about. As usual in 
copula models, data are assumed to be serially independent and lie in the unit
hypercube. 

### Bivariate copula modeling: the BiCop-family

  * `BiCop`: Creates a bivariate copula by specifying the family and parameters
    (or Kendall's tau). Returns an object of class `BiCop`. The class has the
    following methods:
     
     * `print`, `summary`: a brief or comprehensive overview of the bivariate
        copula, respectively. 
            
     * `plot`, `contour`: surface/perspective and contour plots of the copula
        density. Possibly coupled with standard normal margins (default for
        `contour`). 
        
  * `BiCopSim`: Simulates from a bivariate copula.

  * `BiCopEst`: Estimates parameters of a bivariate copula with a prespecified
    family. Estimation can be done by maximum likelihood (`method = "mle"`) or
    inversion of the empirical Kendall's tau (`method = "itau"`, only available
    for one-parameter families). Returns an object of class `BiCop`.
     
  * `BiCopSelect`: Estimates the parameters of a bivariate copula for a set 
    of families and selects the best fitting model (using either AIC or BIC). 
    Returns an object of class `BiCop`.
    
  * `BiCopGofTest`: Goodness-of-Fit tests for bivariate copulas.
    
  * `BiCopVuongClarke`: Vuong and Clarke tests for model comparison within a 
    prespecified set of copula families.
     
  * `BiCopPar2Tau`, `BiCopTau2Par`, `BiCopPar2Beta`, `BiCopPar2TailDep`: 
    Conversion between dependence measures and parameters (for a given family).
    Functions are vectorized in all arguments.
     
  * Evaluate functions related to a bivariate copula: `BiCopPDF`, `BiCopCDF`, 
    `BiCopDeriv`, `BiCopDeriv2`, `BiCopHfunc`, `BiCopHfuncDeriv`, 
    `BiCopHfuncDeriv2`, `BiCopHinv`. Functions are vectorized in the `family`,
    `par`, and `par2` arguments. 
    
  * `BiCopKDE`: Kernel density plots for copula data.
    
  * `BiCopLambda`, `BiCopKPlot`, `BiCopChiPlot`: Further plot types for the 
    analysis of bivariate copulas.
    
For most functions, you can provide an object of class `BiCop` instead of
specifying `family`, `par` and `par2` manually.


### Vine copula modeling: the RVine-family

  * `RVineMatrix`: Creates a vine copula model by specifying structure, family
    and parameter matrices. Such matrices have been introduced by Dissman et al. 
    (2013). Returns an object of class `RVineMatrix`. The class has
    the following methods:
    
    * `plot`: Plots the trees of the the R-vine tree structure. Optionally,
      you can annotate the edges with pair-copula families and parameters.
      
    * `contour`: Creates a matrix of contour plots associated with the
      pair-copulas.

  * `RVineSim`: Simulates from a vine copula model.
      
  * `RVineSeqEst`: Estimates the parameters of a vine copula model with 
    prespecified structure and families.
    
  * `RVineCopSelect`: Estimates the parameters and selects the best family for a
    vine copula model with prespecified structure matrix.
    
  * `RVineStructureSelect`: Fits a vine copula model assuming no prior knowledge.
    It selects the R-vine structure using Dissmann et al. (2013)'s 
    method, estimates parameters for various families, and selects the best 
    family for each pair.

  * `RVineGoFTest`: Goodness-of-Fit tests for a vine copula model (c.f., 
    Schepsmeier, 2013, 2015). Related functions are `RVineGrad`, 
    `RVineHessian`, `RVineStdError`, and `RVinePIT`.
    
  * `RVineVoungTest`, `RVineClarkeTest`: Vuong and Clarke tests for comparing
    two vine copula models.

  * `RVinePar2Tau`, `RVinePar2Beta`: Calculate dependence measures 
    corresponding to a vine copula model.
    
  * `RVinePDF`, `RVineLogLik`, `RVineAIC`, `RVineBIC`: Calculate the density, 
    log-likelihood, AIC, and BIC of a vine copula.


### Additional features
  
The functions `C2RVine` and `D2RVine` create `RVineMatrix` objects for C- and 
D-vine copulas. This is particularly useful for former users of the CDVine
package. 

Furthermore, bivariate and vine copula models from this packages can be used
with the [copula](https://cran.r-project.org/package=copula) package 
(Hofert et al., 2015). For example, `vineCopula` transforms an `RVineMatrix` 
object into an object of class `vineCopula` which provides methods for
`dCopula`, `pCopula`, and `rCopula`. For more details, we refer to the package 
manual.


### Bivariate copula families

In this package several bivariate copula families are included for bivariate 
and multivariate analysis using vine copulas. It provides 
functionality of elliptical (Gaussian and Student-t) as well as Archimedean 
(Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8) copulas to cover a large
range of dependence patterns. For Archimedean copula families,
rotated versions are included to cover negative dependence as well.

The Tawn copula is a non-exchangable extension of the Gumbel copula with three 
parameters. For simplicity, we implemented two versions of the Tawn copula with
two parameters each. Each type has one of the asymmetry parameters fixed to 1, 
so that the corresponding copula density is either left- or right-skewed 
(relative to the main diagonal). In the manual we will call these two new 
copulas "Tawn type 1" and "Tawn type 2".

The following table shows the parameter ranges of bivariate copula families with 
parameters `par` and `par2` and internal coding `family`:

| Copula family                               | `family`     | `par`        | `par2`       |
|:--------------------------------------------|:-------------|:-------------|:-------------|
| Gaussian                                    | `1`          | `(-1, 1)`    | -            |
| Student t                                   | `2`          | `(-1, 1)`    | `(2,Inf)`    |
| (Survival) Clayton                          | `3`, `13`    | `(0, Inf)`   | -            |
| Rotated Clayton (90 and 270 degrees)        | `23`, `33`   | `(-Inf, 0)`  | -            |
| (Survival) Gumbel                           | `4`, `14`    | `[1, Inf)`   | -            |
| Rotated Gumbel (90 and 270 degrees)         | `24`, `34`   | `(-Inf, -1]` | -            |
| Frank                                       | `5`          | `R \ {0}`    | -            |
| (Survival) Joe                              | `6`, `16`    | `(1, Inf)`   | -            |
| Rotated Joe (90 and 270 degrees)            | `26`, `36`   | `(-Inf, -1)` | -            |
| (Survival) Clayton-Gumbel (BB1)             | `7`, `17`    | `(0, Inf)`   | `[1, Inf)`   |
| Rotated Clayton-Gumbel (90 and 270 degrees) | `27`, `37`   | `(-Inf, 0)`  | `(-Inf, -1]` |
| (Survival) Joe-Gumbel (BB6)                 | `8`, `18`    | `[1 ,Inf)`   | `[1, Inf)`   |
| Rotated Joe-Gumbel (90 and 270 degrees)     | `28`, `38`   | `(-Inf, -1]` | `(-Inf, -1]` |
| (Survival) Joe-Clayton (BB7)                | `9`, `19`    | `[1, Inf)`   | `(0, Inf)`   |
| Rotated Joe-Clayton (90 and 270 degrees)    | `29`, `39`   | `(-Inf, -1]` | `(-Inf, 0)`  |
| (Survival) Joe-Frank (BB8)                  | `10`, `20`   | `[1, Inf)`   | `(0, 1]`     |
| Rotated Joe-Frank (90 and 270 degrees)      | `30`, `40`   | `(-Inf, -1]` | `[-1, 0)`    |
| (Survival) Tawn type 1                      | `104`, `114` | `[1, Inf)`   | `[0, 1]`     |
| Rotated Tawn type 1(90 and 270 degrees)     | `124`, `134` | `(-Inf, -1]` | `[0, 1]`     |
| (Survival) Tawn type 2                      | `204`, `214` | `[1, Inf)`   | `[0, 1]`     |
| Rotated Tawn type 2 (90 and 270 degrees)    | `224`, `234` | `(-Inf, -1]` | `[0, 1]`     |


------------------------------------------------------------------------

Associated shiny apps
---------------------

### rvinegraph
This small shiny app enables the user to draw nice tree plots of an R-Vine
copula model using the package 
[d3Network](https://cran.r-project.org/package=d3Network). Models have to be
set up locally in an `RVineMatrix` object and uploaded as .RData file. The page
is still under construction.  
Author: Ulf Schepsmeier

https://rvinegraph.shinyapps.io/rvinegraph

------------------------------------------------------------------------

References
----------

Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.

Bedford, T. and R. M. Cooke (2001). Probability density decomposition for conditionally dependent random variables modeled by vines. Annals of Mathematics and Artificial intelligence 32, 245-268.

Bedford, T. and R. M. Cooke (2002). Vines - a new graphical model for dependent random variables. Annals of Statistics 30, 1031-1068.

Brechmann, E. C., C. Czado, and K. Aas (2012). Truncated regular vines in high dimensions with applications to financial data. Canadian Journal of Statistics 40 (1), 68-85.

Brechmann, E. C. and C. Czado (2011). Risk management with high-dimensional vine copulas: An analysis of the Euro Stoxx 50. Statistics & Risk Modeling, 30 (4), 307-342.

Brechmann, E. C. and U. Schepsmeier (2013). Modeling Dependence with C- and D-Vine Copulas: The R Package CDVine. Journal of Statistical Software, 52 (3), 1-27. <https://doi.org/10.18637/jss.v052.i03>.

Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation of mixed C-vines with application to exchange rates. Statistical Modelling, 12(3), 229-255.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013). Selecting and estimating regular vine copulae and application to financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.

Eschenburg, P. (2013). Properties of extreme-value copulas Diploma thesis, Technische Universitaet Muenchen <https://mediatum.ub.tum.de/node?id=1145695>.

Hofert, M., I. Kojadinovic, M. Maechler, and J. Yan (2015). copula: Multivariate
Dependence with Copulas. R package version 0.999-13 
<https://cran.r-project.org/package=VineCopula>

Joe, H. (1996). Families of m-variate distributions with given margins and m(m-1)/2 bivariate dependence parameters. In L. Rueschendorf, B. Schweizer, and M. D. Taylor (Eds.), Distributions with fixed marginals and related topics, pp. 120-141. Hayward: Institute of Mathematical Statistics.

Joe, H. (1997). Multivariate Models and Dependence Concepts. London: Chapman and Hall.

Knight, W. R. (1966). A computer method for calculating Kendall's tau with ungrouped data. Journal of the American Statistical Association 61 (314), 436-439.

Kurowicka, D. and R. M. Cooke (2006). Uncertainty Analysis with High Dimensional Dependence Modelling. Chichester: John Wiley.

Kurowicka, D. and H. Joe (Eds.) (2011). Dependence Modeling: Vine Copula Handbook. Singapore: World Scientific Publishing Co.

Nelsen, R. (2006). An introduction to copulas. Springer

Nagler, T. (2015). kdecopula: Kernel Smoothing for Bivariate Copula Densities. R package
version 0.6.0. <https://cran.r-project.org/package=kdecopula>

Schepsmeier, U. and J. Stoeber (2012). Derivatives and Fisher information of bivariate copulas. Statistical Papers, 55 (2), 525-542. <https://link.springer.com/article/10.1007/s00362-013-0498-x>.

Schepsmeier, U. (2013) A goodness-of-fit test for regular vine copula models. Preprint. <https://arxiv.org/abs/1306.0818>.

Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for vine copula models with fixed margins. Journal of Multivariate Analysis 138, 34-52.

Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular vine copula models. Computational Statistics, 28 (6), 2679-2707 <https://link.springer.com/article/10.1007/s00180-013-0423-8>.

White, H. (1982) Maximum likelihood estimation of misspecified models, Econometrica, 50, 1-26.
