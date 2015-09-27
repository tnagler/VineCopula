VineCopula
==========

> Statistical inference of vine copulas

[![Build status Linux](https://travis-ci.org/tnagler/VineCopula.svg?branch=master)](https://travis-ci.org/tnagler/VineCopula)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/tnagler/VineCopula?svg=true)](https://ci.appveyor.com/project/tnagler/VineCopula)
[![CRAN version](http://www.r-pkg.org/badges/version/VineCopula)](https://cran.r-project.org/web/packages/VineCopula/index.html) 
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/VineCopula)](https://cran.r-project.org/web/packages/VineCopula/index.html)

This package is made primarily for the statistical analysis with **vine copula** 
models. Vine copulas are a flexible class of dependence models consisting of
bivariate building blocks (see e.g., Aas et al., 2009).

The package includes model selection, model verification via goodness-of-fit
tests, parameter estimation, simulation and visualization tools. Vine copula 
models can be estimated either sequentially or by joint maximum likelihood
estimation. Sampling algorithms and plotting methods are included.

Tools for bivariate exploratory data analysis, **bivariate copula** selection 
and estimation are also provided. 

For all functions data is assumed to lie in the unit hypercube (so-called copula
data).

You can install:

-   the stable release on CRAN:

    ``` r
    install.pacakges("VineCopula")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("tnagler/VineCopula")
    ```

------------------------------------------------------------------------

Package overview
----------------

### Functionality

The package provides functions for dependence modeling with bivariate and R-vine
copulas. Below, we list most functions you should know.

#### Bivariate copula modeling: the BiCop-family

  * `BiCop`: Creates a bivariate copula by specifying the family and parameters
    (or Kendall's tau). Returns an object of class `BiCop`. The class has the
    following methods:
     
     * `print`, `summary`: a brief or comprehensive overview of the bivariate
        copula, respectively. 
            
     * `plot`, `contour`: surface/perspective and contour plots of the copula
        density. Possibly coupled with standard normal margins (default for
        `contour`.) 
     
  * `BiCopEst`: Estimates a bivariate copula to given data and a prespecified
    family. Estimation can be done by maximum likelihood (`method = "mle"`) or
    inversion of the empirical Kendall's tau (`method = "itau"`, only available
    for one-parameter families). Returns an object of class `BiCop`.
     
  * `BiCopSelect`: Estimates a bivariate copula for a selection of families and
    selects the best fitting model (using either AIC or BIC). 
    Returns an object of class `BiCop`.
    
  * `BiCopGofTest`: Goodness-of-Fit test for bivariate copulas.
    
  * `BiCopVuongClarke`: Vuong and Clarke tests for model comparison within a 
    prespecified set of copula families.
     
  * `BiCopSim`: Simulates from a bivariate copula.
   
  * `BiCopPar2Tau`, `BiCopTau2Par`, `BiCopPar2Beta`, `BiCopPar2TailDep`: 
    Conversion between dependence measures and parameters (for a given fmaily).
    Functions are vectorized in all arguments.
     
  * Evaluate functions related to a bivariate copula: `BiCopPDF`, `BiCopCDF`, 
    `BiCopDeriv`, `BiCopDeriv2`, `BiCopHfunc`, `BiCopHfuncDeriv`, 
    `BiCopHfuncDeriv2`, `BiCopHinv`. Functions are vectorized in the `family`,
    `par`, `par2` arguments. 
    
  * `BiCopMetaContour`: Contour plots for a bivariate copula. If data is
    provided, you can create a plot of a kernel estimate of the copula density.
    In the latter case, we recommed to use the 
    [kdecopula](https://github.com/tnagler/kdecopula) package (Nagler, 2015)
    which implements state-of-the art kernel smoothing techniques for copula
    densities.
    
  * `BiCopLambda`, `BiCopKPlot`, `BiCopChiPlot`: Further plot types for the 
    analysis of bivariate copulas.
    
For most functions, you can provide an object of class `BiCop` instead of
specifying `family`, `par` and `par2` manually.


#### R-vine copula modeling: the RVine-family

  * `RVineMatrix`: Creates an R-vine copula by specifying structure, family and
    parameter matrices. Such matrices have been introduced by Dissman et al. 
    (2013). Returns an object of class `RVineMatrix`. The class has
    the following methods:
    
    * `plot`: Plots the trees of the the R-vine tree structure. Optionally,
      you can annotate the edges with pair-copula families and parameters.
      
    * `contour`: Creates a matrix of contour plots associated with the
      pair-copulas in the model.
      
  * `RVineSeqEst`: Estimates the parameters of a vine copula model with 
    prespecified structure and family matrices.
    
  * `RVineCopSelect`: Estimates the parameters and selects the best fitting
    families in a vine copula model with prespecified strucutre matrix.
    
  * `RVineStructureSelect`: Fits an vine model assuming no prior knowledge.
    It selects the R-vine structure using Dissmann et al. (2013)'s 
    method, estimates parameters for selected families, and selects the best 
    family for each pair-copula.

  * `RVineGoFTest`: Goodness-of-Fit tests for a vine copula moel (c.f., 
    Schepsmeier, 2013, 2015). Related functions are `RVineGrad`, 
    `RVineHessian`, `RVineStdError`, and `RVinePIT`.
    
  * `RVineVoungTest`, `RVineClarkeTest`: Vuong and Clarke tests for comparison
    between two vine copula models.

  * `RVineSim`: Simulates from a vine copula model.
  
  * `RVinePar2Tau`, `RVinePar2Beta`: Calculates dependence measures 
    corresponding to a vine copula model.
    
  * `RVinePDF`, `RVineAIC`, `RVineBIC`: Calculates the PDF, AIC, and BIC of a 
    vine copula model, respectively.


#### Additional features
  
The functions `C2RVine` and `D2RVine` create `RVineMatrix` objects for C- and 
D-vine copula models. This is particularly useful for former users of the CDVine
package. 

Furthermore, all bivariate and vine copula models from this packages can be used
within the [copula package](https://r-forge.r-project.org/R/?group_id=600) 
(Hofert et al., 2015). For more details, we refer to the package manual.
    
    

  
### Bivariate copula families

In this package several bivariate copula families are included for bivariate 
analysis as well as for multivariate analysis using vine copulas. It provides 
functionality of elliptical (Gaussian and Student-t) as well as Archimedean 
(Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8) copulas to cover a large
bandwidth of possible dependence structures. For the Archimedean copula families
rotated versions are included to cover negative dependence too. The two 
parameter BB1, BB6, BB7 and BB8 copulas are however numerically instable for 
large parameters, in particular, if BB6, BB7 and BB8 copulas are close to the
Joe copula which is a boundary case of these three copula families. In general,
the user should be careful with extreme parameter choices.

As an asymmetric extension of the Gumbel copula, the Tawn copula with three 
parameters is also included in the package. Both the Gumbel and the Tawn copula 
are extreme-value copulas, which can be defined in terms of their corresponding
Pickands dependence functions. For simplicity, we implemented two versions of 
the Tawn copula with two parameters each. Each type has one of the asymmetry 
parameters fixed to 1, so that the corresponding Pickands dependence is either
left- or right-skewed. In the manual we will call these two new copulas 
"Tawn type 1" and "Tawn type 2".

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

### R-vine copula models

The specification of an R-vine is done in matrix notation, introduced by 
Dissmann et al. (2013). One matrix contains the R-vine tree structure, one the
copula families utilized and two matrices corresponding parameter values. These 
four matrices are stored in an `RVineMatrix` object created by the function 
`RVineMatrix`. Each matrix is a d x d lower triangular matrix. Since C- and 
D-vines are special cases, boundary cases, of R-vines one can write each C- or 
D-vine in R-vine notation. The transformation of notation to an R-vine can be 
done via `C2RVine` and `D2RVine`, which provide an interface to the package 
**CDVine** (https://cran.r-project.org/package=CDVine). For more details see
the documentation of the functions.

------------------------------------------------------------------------

Associated shiny apps
---------------------

### Copulatheque
This small shiny app illustrates a couple of bivariate copula families 
implemented in the [copula](http://cran.r-project.org/web/packages/copula/index.html), [VineCopula](http://cran.r-project.org/web/packages/VineCopula/index.html) and [spcopula](http://r-forge.r-project.org/projects/spcopula/) R packages. Density
and pairs plots are drawn as well as Kendall's tau and tail dependence 
coefficients. This side and script is provided by Benedikt Gräler (Universtity 
of Münster) 

http://ifgi.uni-muenster.de/~b_grae02/indexCopulatheque.html


### rvinegraph
This small shiny app enables the user to draw nice tree plots of an R-Vine
copula model using the package 
[d3Network](https://cran.r-project.org/package=d3Network). Models have to be
set up in an `RVineMatrix` object. The page is still under construction.
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

Brechmann, E. C. and U. Schepsmeier (2013). Modeling Dependence with C- and D-Vine Copulas: The R Package CDVine. Journal of Statistical Software, 52 (3), 1-27. <http://www.jstatsoft.org/v52/i03/>.

Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation of mixed C-vines with application to exchange rates. Statistical Modelling, 12(3), 229-255.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013). Selecting and estimating regular vine copulae and application to financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.

Eschenburg, P. (2013). Properties of extreme-value copulas Diploma thesis, Technische Universitaet Muenchen <http://mediatum.ub.tum.de/node?id=1145695>.

Hofert, M., I. Kojadinovic, M. Maechler, and J. Yan (2015). copula: Multivariate
Dependence with Copulas. R package version 0.999-13 
<http://CRAN.R-project.org/package=copula>

Joe, H. (1996). Families of m-variate distributions with given margins and m(m-1)/2 bivariate dependence parameters. In L. Rueschendorf, B. Schweizer, and M. D. Taylor (Eds.), Distributions with fixed marginals and related topics, pp. 120-141. Hayward: Institute of Mathematical Statistics.

Joe, H. (1997). Multivariate Models and Dependence Concepts. London: Chapman and Hall.

Knight, W. R. (1966). A computer method for calculating Kendall's tau with ungrouped data. Journal of the American Statistical Association 61 (314), 436-439.

Kurowicka, D. and R. M. Cooke (2006). Uncertainty Analysis with High Dimensional Dependence Modelling. Chichester: John Wiley.

Kurowicka, D. and H. Joe (Eds.) (2011). Dependence Modeling: Vine Copula Handbook. Singapore: World Scientific Publishing Co.

Nelsen, R. (2006). An introduction to copulas. Springer

Nagler, T. (2015). kdecopula: Kernel Smoothing for Bivariate Copula Densities. R package
version 0.2.0. https://github.com/tnagler/kdecopula

Schepsmeier, U. and J. Stoeber (2012). Derivatives and Fisher information of bivariate copulas. Statistical Papers, 55 (2), 525-542. <http://link.springer.com/article/10.1007/s00362-013-0498-x>.

Schepsmeier, U. (2013) A goodness-of-fit test for regular vine copula models. Preprint. <http://arxiv.org/abs/1306.0818>.

Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for vine copula models with fixed margins. Journal of Multivariate Analysis 138, 34-52.

Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular vine copula models. Computational Statistics, 28 (6), 2679-2707 <http://link.springer.com/article/10.1007/s00180-013-0423-8>.

White, H. (1982) Maximum likelihood estimation of misspecified models, Econometrica, 50, 1-26.
