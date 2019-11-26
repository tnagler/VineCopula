#' Blomqvist's Beta Value of a Bivariate Copula
#'
#' This function computes the theoretical Blomqvist's beta value of a bivariate
#' copula for given parameter values.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version \cr
#' \preformatted{BiCopPar2Beta(obj)} can be used.
#'
#' @param family integer; single number or vector of size `n`; defines the
#' bivariate copula family: \cr
#' `0` = independence copula \cr
#' `2` = Student t copula (t-copula) \cr
#' `1` = Gaussian copula \cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `7` = BB1 copula \cr
#' `8` = BB6 copula \cr
#' `9` = BB7 copula \cr
#' `10` = BB8 copula \cr
#' `13` = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' `14` = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' `16` = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' `17` = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' `18` = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' `19` = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' `20` = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' `23` = rotated Clayton copula (90 degrees) \cr
#' `24` = rotated Gumbel copula (90 degrees) \cr
#' `26` = rotated Joe copula (90 degrees) \cr
#' `27` = rotated BB1 copula (90 degrees) \cr
#' `28` = rotated BB6 copula (90 degrees) \cr
#' `29` = rotated BB7 copula (90 degrees) \cr
#' `30` = rotated BB8 copula (90 degrees) \cr
#' `33` = rotated Clayton copula (270 degrees) \cr
#' `34` = rotated Gumbel copula (270 degrees) \cr
#' `36` = rotated Joe copula (270 degrees) \cr
#' `37` = rotated BB1 copula (270 degrees) \cr
#' `38` = rotated BB6 copula (270 degrees) \cr
#' `39` = rotated BB7 copula (270 degrees) \cr
#' `40` = rotated BB8 copula (270 degrees) \cr
#' `104` = Tawn type 1 copula \cr
#' `114` = rotated Tawn type 1 copula (180 degrees) \cr
#' `124` = rotated Tawn type 1 copula (90 degrees) \cr
#' `134` = rotated Tawn type 1 copula (270 degrees) \cr
#' `204` = Tawn type 2 copula \cr
#' `214` = rotated Tawn type 2 copula (180 degrees) \cr
#' `224` = rotated Tawn type 2 copula (90 degrees) \cr
#' `234` =  rotated Tawn type 2 copula (270 degrees) \cr
#' Note that the Student's t-copula is not allowed since the CDF of the t-copula
#' is not implemented (see [BiCopCDF()]).
#' @param par numeric; single number or vector of size `n`; copula
#' parameter.
#' @param par2 numeric; single number or vector of size `n`; second
#' parameter for the two parameter BB1, BB6, BB7, BB8, Tawn type 1 and type 2
#' copulas (default: `par2 = 0`).
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#'
#' @return Theoretical value of Blomqvist's beta corresponding to the bivariate
#' copula `family` and parameter(s) `par`, `par2`.
#'
#' @note The number `n` can be chosen arbitrarily, but must agree across
#' arguments.
#'
#' @author Ulf Schepsmeier
#'
#' @references Blomqvist, N. (1950).  On a measure of dependence between two
#' random variables. The Annals of Mathematical Statistics, 21(4), 593-600.
#'
#' Nelsen, R. (2006). An introduction to copulas.  Springer
#'
#' @examples
#' ## Example 1: Gaussian copula
#' BiCopPar2Beta(family = 1, par = 0.7)
#' BiCop(1, 0.7)$beta  # alternative
#'
#' ## Example 2: Clayton copula
#' BiCopPar2Beta(family = 3, par = 2)
#'
#' ## Example 3: different copula families
#' BiCopPar2Beta(family = c(3,4,6), par = 2:4)
#'
BiCopPar2Beta <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate beta
    Cuv <- BiCopCDF(rep(0.5, length(par)),
                    rep(0.5, length(par)),
                    family,
                    par,
                    par2,
                    check.pars = check.pars)
    4 * Cuv - 1
}
