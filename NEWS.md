VineCopula 2.1.3 (August 15, 2017)
----------------------------------------------------------------

NEW FEATURES

  * `summary.RVineMatrix()` invisibly returns a `data.frame` containg most of
    what is printed as output.
    
  * Less restrictive conditions on what is considered an appropriate `treecrit` 
    function in `RVineStructureSelect()` (thanks to Thibault Vatter, #40).
    
  * New option `verbose` in `RVinePDF()` (thanks to @cag51, #36).

BUG FIXES

  * Fix parameter bounds in `RVineMLE()` corresponding to updated requirements 
    in `BiCop`-functions.
   
  * Adapt to re-naming of `tailIndex` to `lambda` in the copula package (thanks
    to Benedikt Graeler, #41, #42).
    
  * Fix bug in detecting C-vine copulas for `summary/print.RVineMatrix()` (#38).
    

VineCopula 2.1.2 (April 24, 2017)
----------------------------------------------------------------

NEW FEATURES

  * Online API documentation on https://tnagler.github.io/VineCopula/.
  
  * Faster `BiCopCDF()`.

BUG FIXES

  * Fixed bug in preprocessing of `weights` argument.

  * More informative error message when family is unknown.

  * Safer `BiCopSelect()` with `presel = TRUE` and insufficient data.
  
  * Fixed `RVineMatrix()` (output dimension was `d-1` and `naturalOrder = TRUE`
    wasn't working, thanks @tvatter).
    
  * Safer `BiCopEst()` for `method = "itau"`.


VineCopula 2.1.1 (January 11, 2017)
----------------------------------------------------------------

IMPORTS

  * Package now requires `copula (>= 0.999-16)`. The new version of copula
    requires VineCopula to be re-installed, because the old `fitCopula()` 
    method doesn't work any longer, but was re-exported by VineCopula.
    
NEW FEATURES

  * Package developers can use the VineCopula C-functionality by linking against
    VineCopula through the `LinkingTo` field.
  


VineCopula 2.1.0 (December 23, 2016)
----------------------------------------------------------------

# DEPENDS

  * Now depends explicitly on `R (>= 3.1.0)`. So far, this dependence was
    implicit trhough our dependence on the copula package.


NEW FEATURES

  * All estimation functions now have a `method` argument. The default is 
    `method = "mle"` and corresponds to the old behavior. The other option, 
    `method = "itau"`, estimates the parameters by inversion of Kendall's. It is
    much faster than `method = "mle"`, but is only available for one-parameter 
    families and the t-copula. Big thanks to Thibault Vatter who did most of the
    work (PR #25).
    
  * New function `RVineMatrixSample` that randomly generates valid R-vine 
    matrices using the algorithm of Joe et al. (2011). Contributed by 
    Thibault Vatter (PR #27).

  * Faster versions of `RVineMatrix` and `RVineCopSelect`, and 
   `RVineStructureSelect` by avoiding unnecessary computations (thanks to 
    Thibault Vatter, PRs #29 and #31).
    
  * All `-Select` functions now have a `presel` argument. If `TRUE` (default) 
    the familyset is reduced to families that have asymmetry charactistics that
    conform with the observed data.
    
  * `RVineSim`, `RVineLogLik`, and `RVinePDF` have been implemented in a way
    that demands less memory. For `RVineLogLik`, the option `calculate.V` has
    to be set to `TRUE`.
    
    
BUG FIXES

  * Fixed bug in upper tail dependence coefficient for survival BB1 (reported
    by Viviana Fernandez, thanks!).
    
  * `RVineStructureSelect` now works in dimension two as well.
  
  * `RVineMatrixCheck` returns the new error code `-4` when the matrix is 
    neither lower nor upper triangular (threw an actual error before).
    
  * `contour.RVineMatrix` now arranges the contour matrix conforming with the
    family and parameter matrices.
    


VineCopula 2.0.5 (September 25, 2016)
----------------------------------------------------------------

DEPENDS

  * Require higher version of the copula package (>= 0.999-15).


BUG FIXES

  * Use `se = FALSE` as default throughout the package for faster estimation.

  * Fixed errors in `BiCopGofTest` for the Student t copula and `par2` close to
    2 and for most other families for `tau` close to 0 (reported by Thong Huy 
    Nguyen, thanks!).
  
  * Fix sign in starting parameters for Tawn MLE.
  
  * Fixed storage positions of `par2` in output matrix of `RVineMLE` (reported
    by Robin Evans, thanks!).



VineCopula 2.0.4 (August 8, 2016)
----------------------------------------------------------------

NEW FEATURES

  * Option for kernel contours in `contour.RVineMAtrix`.


BUG FIXES

  * Return scalar instead of 1-dim array in tau <-> par conversion.
  
  * Fix vectorized call of `BiCopCDF`.
  
  * Negative selection of families is now working properly.
  
  * Correct logLik calculation in output of `RVineMLE`.
  
  * Dependence measures are updated in output of `RVineCor2pcor`.
  
  * Fix bug in annotation of edges with Kendall's when using `plot.RVineMatrix`.


VineCopula 2.0.1 (June 9, 2016)
----------------------------------------------------------------

BUG FIXES

  * fixed small memory leak (reported by Prof. Ripley, thanks!).


VineCopula 2.0.0 (June 8, 2016)
----------------------------------------------------------------

MAINTAINER

 * changed from Tobias Erhardt to Thomas Nagler (thomas.nagler@tum.de).


DEPENDS

  * igraph has been removed from `Imports`.
  
  * network has been added to `Imports`.
  
  
NEW FEATURES

  * All functions in the package can now handle NAs in the data. A warning
    message indicates presence of NAs and explains how they are treated.

  * Extend `BiCop` and `RVineMatrix` objects (now include: associated dependence
    measures; fit statistics and p-values, if available).

  * New methods `print` and `summary` for objects of class `BiCop` and 
    `RVineMatrix`.

  * New plotting generics:
  
    * `plot.RVineMatrix` for plotting vine trees,
  
    * `contour.RVineMatrix` for a matrix of contour plots,
  
    * `contour.BiCop` as short hand for `plot.BiCop(..., type = "contour")`.

  * Vectorize `BiCopXyz`-functions w.r.t. `family`, `par`, `par2`:
  
    * in C: `BiCopPDF`, `BiCopHfunc`, `BiCopHinv`, `BiCopDeriv`, `BiCopDeriv2`,
      `BiCopHfuncDeriv`, `BiCopHfuncDeriv2`.
      
    * in R: `BiCopCDF`, `BiCopPar2Tau`, `BiCopPar2Beta`, `BiCopPar2TailDep`.

  * Etimation of vine copulas 
    (`RVineStructureSelect`, `RVineCopSelect`, `RVineSeqEst`) can now be done in 
    parallel using the `cores` argument (based on `foreach`) and is more memory
    efficient (pseudo-observations are discarded as soon they are useless).
    
  * `RVineStructureSelect` now takes an argument `treecrit` allowing for 
    several preimplemented (and custom) choices of the edge weight used in 
    Dissmann's algorithm.
  
  * Treatment of `familyset` in the -Select functions:
  
    * independence copula is handled as a regular family,
    
    * negative integers can be used to select from all but a subset of 
      families.
     
  * New function `BiCopCompare`: A shiny app where the user can visually assess 
    how well several families fit the data. 

  * New function `BiCopKDE` for kernel density plots (based on kdecopula 
    package). 
    
  * New function `BiCopCondSim` for conditional simulation from a bivariate 
    copula.
    
  * New function `BiCopHinv` for computation of inverse h-functions.
  
  * New functions `BiCopHfunc1`, `BiCopHfunc2`, `BiCopHinv1`, and `BiCopHinv2`
    that only compute one of the two h-functions (or inverse h-functions).

  * New function `BiCopCheck`  for checking of family/parameter
    consistency.
    
  * Add `check.pars`/`check.taus` argument to the above functions for the option
    to omit family/parameter consistency checks (for internal usage). When
    `FALSE`, the Clayton and Frank copulas can be used with `par = 0`.
    
  * Faster implementations of `BiCopPar2Tau`/`BiCopTau2Par` for Frank copula and
    `BiCopTau2Par` conversion for Joe copula. 

    
BUG FIXES

  * Correct call for non-t families in `BiCopHfuncDeriv2(..., deriv = "par1u2")`.

  * The S4-class objets of the Tawn copulas pointed to Archimedean CDFs, now 
    corrected to true CDFs based on C-code.

  * `TauMatrix`: restriction for input data to be in [0,1] removed.

  * `RVineCopSelect`: no printing of family matrix.

  * Added methods for Pickand's dependence function "A" for `tawnT1Copula`, 
    `surTawnT1Copula`, `tawnT2Copula` and `surTawnT2Copula`.

  * Use C-code instead of R-code and remove redundant C-code of
    Tawn copulas.
    
  * Small bug fix in log-likelihood for families 3, 4, 7, 17.
  
  * Increased upper limit for uniroot in `Joe.itau.JJ`.
  
  * Fixed rotations of Tawns (they were actually reflection w.r.t. the axes 
    u = 0.5 and u2 = 0.5)
	
  * Correct calculation of the goodness-of-fit test based on Whites Information 
    matrix test for bivariate copulas `BiCopGofTest(..., method = "white")`.
	The variance matrix needed for the test statistic had a poor approximation.
	Thereby the asymptotic p-values are corrected.
	
  * Bound parameter ranges for Archimedean copulas to avoid numerical 
    instabilities in -PDF and -Sim functions.
   

VineCopula 1.6-1 (November 9, 2015)
----------------------------------------------------------------

DEPENDS:
  
  * Removed CDVine from Suggests.

BUG FIXES:

  * Fix code/documentation mismatch in function 'pobs' following a 
    change in the copula package.
  

VineCopula 1.6 (July 16, 2015)
----------------------------------------------------------------

NEW FEATURES

  * `RVineTreePlot`: option for a legend (and numbered nodes and edges).

BUG FIXES

  * Definition of "C" in BiCopCDF for tawn copulas used constants `u1` and `u2` 
    instead of arguments `u` and `v`.

  * `RVineStructureSelect`: Adjust to new version of igraph. Tree structure was 
    not selected correctly. igraph function names changed to the names used in
    the new version. Some small modifications to avoid some for loops and make 
    the code easier to read.

IMPORTS

  * Extend Imports to avoid undefined globals (CRAN E-mail 02.07.2015).

  * New version requires `igraph (>= 1.0.0)`.


VineCopula 1.5 (June 2, 2015)
----------------------------------------------------------------

NEW FEATURES

  * `as.copuladata`: coerce to class `copuladata`.

  * `pairs.copuladata`: `pairs` plots for objects of class `copuladata`.

  * `RVinePDF`: PDF of an R-Vine Copula Model.

  * `BiCopSelect`, `RVineCopSelect`, `RVineStructureSelect`: add option 
    `rotations = TRUE` which augments the familyset with all rotations to a 
    given family.

  * `RVineMatrix`, `RVineStructureSelect`: allow upper triangular matrices as 
    input (output remains lower triangular).

  * `BiCop` objects for bivariate copulas:

  * add constructor `BiCop` and plotting generic `plot.BiCop`.

  * define results of `BiCopEst`/`BiCopSelect` as `BiCop` objects.

  * add compatibility with other `BiCopXyz` functions (`BiCopPDF`, 
    `BiCopPar2Tau`, etc.).

BUG FIXES

  * `BiCopEst`: extend search interval for Tawn MLE to avoid `optim`-errors.

  * `BiCopEst`: fix for optim error ('non-finite value supplied').

  * `RVineSim`: reorder `U` so that it corresponds to the order of `RVM`.

  * `RVineCor2pcor`: include normalization step for a more intuitive behavior, 
    bug fix for $d = 2, 3$ and $d > 9$.

  * `RVinePcor2cor`: bug fixes for $d = 2$ and $d > 9$.

  * `RVineCopSelect`: `RVM` object now uses variable names as provided by data.


VineCopula 1.4 (January 26, 2015)
----------------------------------------------------------------

NEW FEATURES

  * `BiCopTau2Par` and `BiCopPar2Tau`: fully vectorized (parameter/tau input), 
    and sanity checks extended. Before vector input was not prohibited. 
    However, both functions were not intended to be used for vectorized input.


VineCopula 1.3-2 (January 19, 2015)
----------------------------------------------------------------

NEW AUTHOR

  * Thomas Nagler

NEW FEATURES

  * Import/Export of function 'pobs' from 'copula' package.

BUG FIXES

  * `RVineStructureSelect`: Bug concerning the dimensions of input data/security 
    queries fixed (Reported by Sarka Cerna, Radek Solnicky and Ludovic Theate. 
    Thanks a lot!)

  * `RVineStructureSelect`: Correct handling of rotated BBs and Tawns.

  * `BiCopSelect`, `BiCopEst`: Improved starting values for Tawn MLE.

  * `hfunc.c`: 

    * Correct `Hfunc1` for Tawns.

    * Bound all results to lie in [0,1] (`Hfunc1` and `Hfunc2`)

    * Extension of `Hinv1` and `Hinv2` in analogy to `Hfunc1` and `Hfunc2`.

  * `incompleteBeta.c`: Misuse of the C function abs (as reported by CRAN)
    corrected to `fabs`.

  * `gof_PIT.R`: Use of `require()` replaced by `requireNamespace` according to
    'Writing R Extensions'.

  * Package `ADGofTest` removed from `Suggests` (see 'Writing R Extensions' 
    for usage of Suggests).

  * Import of function `ad.test` from `ADGofTest` for `gof_PIT.R`.


VineCopula 1.3-1 (September 10, 2014)
----------------------------------------------------------------

BUG FIXES

  * Bootstrap procedure for the White test (`RVineGofTest`, `gof_White`) was 
    incorrect. 
    (Reported by Piotr Zuraniewski and Daniel Worm. Thanks!)

  * Bootstrap procedure for the PIT based and the ECP based test were incorrect. 
    First, C starts to count at 0 not 1. This could result in zero entries in
    the bootstrapped data matrix. Second, forget to permute vdirect and vindirect 
    according to the permutation of data.

  * `BiCopSelect`: For the rotated BB7 and BB8 (`family = 37, 38`) the limiting
    cases were incorrect for very small parameters (copy&paste error).
    (Reported by Radek Solnicky. Thanks!)


VineCopula 1.3 (March 26, 2014)
----------------------------------------------------------------

MAINTAINER

 * changed from Ulf Schepsmeier to Tobias Erhardt (tobias.erhardt@tum.de).


VineCopula 1.2-1 (March 21, 2014)
----------------------------------------------------------------

NEW FEATURES

  * Added tests generated from example code.

IMPORTS

  * Moved copula from `Depends` to the more appropriate `Import` field.


VineCopula 1.2-1 (March 4, 2014)
----------------------------------------------------------------

NEW FEATURES

  * `RVineSim` allows to commit a `(N x d)`-matrix of $U[0,1]$ random variates
    to be transformed to the copula sample. For example if you want to use quasi
    random variables instead of the pseudo random variables implemented in R. 
    (Thanks to Marius Hofert)

  * The package now contains class wrappers that are compatible with the `copula`
    class from the `copula` R-package. These include all bivariate families 
    currently implemented: The class representation for different rotated 
    families of e.g the BB6 family are represented as `BB6Copula`, `r90B6Copula`,
    `surBB6Copula` and `r270BB6Copula`. These bivariate classes are fully 
    compatible with the standard copula methods such as `dCopula`, `pCopula`, 
    `rCopula` or `fitCopula` including `persp` and `contour`. A vine copula can
    as well be coerced into a class representation of `vineCopula`. However,
    the support of the standard methods is limited. See the corresponding help 
    pages for details. Earlier introduced R-wrapper of C-functions have been 
    removed, as they are no longer needed by the `spcopula` R-package.

  * Added parameter `verbose` to `RVineLogLik` to allow to suppress some debug
    output.

BUG FIXES

  * `RVineMLE`: the `optim` argument `parscale` was not correctly defined for 
    all cases.
    
  * `RVineAIC`/`RVineBIC`: Instead of the function arguments `par` and `par2` 
    the calculation was based on `RVM$par` and `RVM$par2`. This is corrected now.
    (reported by Marcel Duellmann; thanks)
    
  * `RVineStructureSelect`: The new igraph version returned a different 
    variable type causing an error in the second and higher order tree selection.


VineCopula 1.2 (October 09, 2013)
----------------------------------------------------------------

NEW FEATURES

  * `RVinePIT`: calculation of the probability integral transform (PIT) for 
    R-vines.

  * `RVineGofTest`: 15 different goodness-of-fit tests for R-vine copulas 
    (Schepsmeier 2013).

  * `print.RVM`: A more detailed summary is printed if `print(RVM, detail = TRUE)`
    is set.

  * `BetaMatrix`: Matrix of empirical Blomqvist's beta values.

  * `BiCopPar2Beta`: Blomqvist's beta value of a bivariate copula.

  * `RVinePar2Beta`: Blomqvist's beta values of an R-vine copula model.

  * `RVineCor2pcor`: correlations to partial correlations for R-vines.

  * `RVinePcor2cor`: partial correlations to correlations for R-vines.

  * New copula families for most of the BiCop as well as for the
    `RVine`-functions: 

  * As an asymmetric extension of the Gumbel copula, the Tawn copula with
    three parameters is now also included in the package. Both the Gumbel and
    the Tawn copula are extreme-value copulas, which can be  defined in terms
    of their corresponding Pickands dependence functions.

  * For simplicity, we implemented two versions of the Tawn copula with two
    parameters each. Each type has one of the asymmetry parameters fixed to 1,
    so that the corresponding Pickands dependence is either left- or 
    right-skewed. In the manual we will call these two new copulas
    "Tawn type 1" and "Tawn type 2". 

  * The families `104`, `114`, `124`, `134` denote the Tawn copula and their 
    rotated versions in the case of left skewness (Tawn type 1). 

  * The families `204`, `214`, `224`, `234` denote the Tawn copula and their 
    rotated versions in the case of right skewness (Tawn type 2).

BUG FIXES

  * `BiCopPar2Tau`: corrected calculation of Kendall's tau of rotated BB7.
(Reported by Giampiero Marra. Thanks!)

  * `RVineStructureSelect`: Corrected code for the igraph package.

  * `RVineTreePlot`: Now a 3-dimensional R-vine can be plotted too.

  * Corrected upper tail dependence coefficient for the survival BB1 copula 
    (`BiCopPar2TailDep`).

  * Minor improvement in `BiCopSelect` regarding the starting values for 
    parameter estimation.

DOCUMENTAION

  * Updated manual files.


VineCopula 1.1-2 (July 09, 2013)
----------------------------------------------------------------

NEW AUTHOR

  * Benedikt Graeler  

NEW FEATURES

  * Changed dependency from igraph0 to igraph since the support for igraph0
    will be quit soon.

  * Additional validity check of the R-vine matrix in `RVineMatrix` (Code
    provided by Harry Joe). Also available as separate function 
    `RVineMatrixCheck`.

  * New bivariate copula: Reflection asymmetric Archimedean copula. In our 
    functions it is `family = 41, 51, 61, 71` ( with rotated versions).
    So far only implemented in some bivariate functions (not documented so 
    far; experimental).


BUG FIXES

  * New (correct) examples for the Clarke and Vuong test.
  
  * Fixed memory problem in the C-function `ktau` (`TauMatrix`).
  

VineCopula 1.1-1 (February 7, 2013)
----------------------------------------------------------------

BUG FIXES

  * Fixed issue with the inverse h-function of the Gumbel copula.


VineCopula 1.1 (February 4, 2013)
----------------------------------------------------------------

NEW FEATURES

  * `BiCopGofTest`: Goodness-of-fit test for bivariate copulas based on White's 
    information matrix equality as introduced by Wanling and Prokhorov (2011).
    The formally included function `BiCopGofKendall` is now part of 
    `BiCopGofTest` (`method = "kendall"`).

  * Additional edge label `pair` in `RVineTreePlot` to display the indices of
    the (conditioned) pairs of variables identified by the edges.

  * In `RVineStructureSelect` and `RVineCopSelect` a truncation level can be set.

  * Improved inverse h-functions for the Gumbel and Joe copulas (thanks to 
    Harry Joe).

  * C to R wrapping functions for the h-functions (`Hfunc1`, `Hfunc2`), the 
    bivariate log-likelihood function (`LL_mod_seperate`), the bivariate
    Archimedean copula CDF (`archCDF`) and the simulation function for C- and
    D-vines (`pcc`) (request of Benedikt Graeler for the R-package `spcopula`).

  * The functions R2CVine and R2Dine were removed, since they were only correct
    in special cases.
    
BUG FIXES

  * Work around for a problem with `optim` and the analytical gradient in 
    `BiCopEst`.
    
  * Improvement of the bivariate maximum likelihood estimation (`BiCopEst`).
  
  * In the functions `BiCopCDF` and `BiCopGoFTest(..., method = "Kendall")` the
    t-copula is not implemented any more. The calculation of the CDF was 
    incorrect for non-integer values of the degrees-of-freedom parameter. The 
    implemented algorithm in the `mvtnorm` package only works for integer values 
    of the degrees-of-freedom parameter.
    
  * Improvement in the calculation of the cdf of the Frank copula (`BiCopCDF`).
