<!-- README.md is generated from README.Rmd. Please edit that file -->
VineCopula
==========

> Statistical inference of vine copulas

[![CRAN version](http://www.r-pkg.org/badges/version/VineCopula)](https://cran.r-project.org/web/packages/VineCopula/index.html) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/VineCopula)](https://cran.r-project.org/web/packages/VineCopula/index.html)

This package is made primarily for the statistical analysis with **vine copula** models.
Vine copulas are a flexible class of dependence models consisting of bivariate building blocks.

The package includes model selection, model verification via goodness-of-fit tests, parameter estimation, simulation and visualization tools.
Vine copula models can be estimated either sequentially or by joint maximum likelihood estimation. Sampling algorithms and plotting methods are included.

Tools for bivariate exploratory data analysis, **bivariate copula** selection and estimation are also provided. 

For all functions data is assumed to lie in the unit hypercube (so-called copula data).

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

### Bivariate copula families

In this package several bivariate copula families are included for bivariate analysis as well as for multivariate analysis using vine copulas. It provides functionality of elliptical (Gaussian and Student-t) as well as Archimedean (Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8) copulas to cover a large bandwidth of possible dependence structures. For the Archimedean copula families rotated versions are included to cover negative dependence too. The two parameter BB1, BB6, BB7 and BB8 copulas are however numerically instable for large parameters, in particular, if BB6, BB7 and BB8 copulas are close to the Joe copula which is a boundary case of these three copula families. In general, the user should be careful with extreme parameter choices.

As an asymmetric extension of the Gumbel copula, the Tawn copula with three parameters is also included in the package. Both the Gumbel and the Tawn copula are extreme-value copulas, which can be defined in terms of their corresponding Pickands dependence functions. For simplicity, we implemented two versions of the Tawn copula with two parameters each. Each type has one of the asymmetry parameters fixed to 1, so that the corresponding Pickands dependence is either left- or right-skewed. In the manual we will call these two new copulas "Tawn type 1" and "Tawn type 2".

The following table shows the parameter ranges of bivariate copula families with parameters `par` and `par2` and internal coding `family`:

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

The specification of an R-vine is done in matrix notation, introduced by Dissmann et al. (2013). One matrix contains the R-vine tree structure, one the copula families utilized and two matrices corresponding parameter values. These four matrices are stored in an `RVineMatrix` object created by the function `RVineMatrix`. Each matrix is a d x d lower triangular matrix. Since C- and D-vines are special cases, boundary cases, of R-vines one can write each C- or D-vine in R-vine notation. The transformation of notation to an R-vine can be done via `C2RVine` and `D2RVine`, which provide an interface to the package **CDVine** (https://cran.r-project.org/package=CDVine). For more details see the documentation of the functions.

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

Joe, H. (1996). Families of m-variate distributions with given margins and m(m-1)/2 bivariate dependence parameters. In L. Rueschendorf, B. Schweizer, and M. D. Taylor (Eds.), Distributions with fixed marginals and related topics, pp. 120-141. Hayward: Institute of Mathematical Statistics.

Joe, H. (1997). Multivariate Models and Dependence Concepts. London: Chapman and Hall.

Knight, W. R. (1966). A computer method for calculating Kendall's tau with ungrouped data. Journal of the American Statistical Association 61 (314), 436-439.

Kurowicka, D. and R. M. Cooke (2006). Uncertainty Analysis with High Dimensional Dependence Modelling. Chichester: John Wiley.

Kurowicka, D. and H. Joe (Eds.) (2011). Dependence Modeling: Vine Copula Handbook. Singapore: World Scientific Publishing Co.

Nelsen, R. (2006). An introduction to copulas. Springer

Schepsmeier, U. and J. Stoeber (2012). Derivatives and Fisher information of bivariate copulas. Statistical Papers, 55 (2), 525-542. <http://link.springer.com/article/10.1007/s00362-013-0498-x>.

Schepsmeier, U. (2013) A goodness-of-fit test for regular vine copula models. Preprint. <http://arxiv.org/abs/1306.0818>.

Schepsmeier, U. (2015) Efficient information based goodness-of-fit tests for vine copula models with fixed margins. Journal of Multivariate Analysis 138, 34-52.

Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular vine copula models. Computational Statistics, 28 (6), 2679-2707 <http://link.springer.com/article/10.1007/s00180-013-0423-8>.

White, H. (1982) Maximum likelihood estimation of misspecified models, Econometrica, 50, 1-26.
