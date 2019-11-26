#' Conditional simulation from a Bivariate Copula
#'
#' This function simulates from a parametric bivariate copula, where on of
#' the variables is fixed. I.e., we simulate either from
#' \eqn{C_{2|1}(u_2|u_1;\theta)} or \eqn{C_{1|2}(u_1|u_2;\theta)}, which are both
#' conditional distribution functions of one variable given another.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version
#' \preformatted{BiCopCondSim(N, cond.val, cond.var, obj)}
#' can be used.
#'
#'
#' @param N Number of observations simulated.
#' @param cond.val numeric vector of length `N` containing the values to
#' condition on.
#' @param cond.var either `1` or `2`; the variable to condition on.
#' @param family integer; single number or vector of size `N`; defines the
#' bivariate copula family: \cr
#' `0` = independence copula \cr
#' `1` = Gaussian copula \cr
#' `2` = Student t copula (t-copula) \cr
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
#' `234` = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par numeric; single number or vector of size `N`; copula
#' parameter.
#' @param par2 numeric; single number or vector of size `N`; second
#' parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8,
#' Tawn type 1 and type 2; default: `par2 = 0`). `par2` should be a
#' positive integer for the Students's t copula `family = 2`.
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#'
#' @return A length `N` vector of simulated from conditional distributions
#' related to bivariate copula with `family` and parameter(s) `par`,
#' `par2`.
#'
#' @author Thomas Nagler
#'
#' @seealso [BiCopCDF()], [BiCopPDF()],
#' [RVineSim()]
#'
#' @examples
#' # create bivariate t-copula
#' obj <- BiCop(family = 2, par = -0.7, par2 = 4)
#'
#' # simulate 500 observations of (U1, U2)
#' sim <- BiCopSim(500, obj)
#' hist(sim[, 1])  # data have uniform distribution
#' hist(sim[, 2])  # data have uniform distribution
#'
#' # simulate 500 observations of (U2 | U1 = 0.7)
#' sim1 <- BiCopCondSim(500, cond.val = 0.7, cond.var = 1, obj)
#' hist(sim1)  # not uniform!
#'
#' # simulate 500 observations of (U1 | U2 = 0.1)
#' sim2 <- BiCopCondSim(500, cond.val = 0.1, cond.var = 2, obj)
#' hist(sim2)  # not uniform!
#'
BiCopCondSim <- function(N, cond.val, cond.var, family, par, par2 = 0,
                         obj = NULL, check.pars = TRUE) {
    if (length(cond.val) == 1)
        cond.val <- rep(cond.val, N)
    if (length(cond.val) != N)
        stop("cond.val must be a numeric vector of length 1 or N")
    if (any(cond.val <= 0)  || any(cond.val >= 1))
        stop("cond.val must be in the interval (0, 1)")
    if (!all(cond.var %in% 1:2))
        stop("cond.var has to be either 1 or 2")
    if (missing(family))
        family <- NA
    if (missing(par))
        par <- NA

    ## start with independent uniforms (byrow for backwards compatibility)
    w <- runif(N)

    ## simulate from copula conditional on U_cond.var = cond.val
    switch(cond.var,
           "1" = BiCopHinv1(cond.val, w, family, par, par2, obj, check.pars),
           "2" = BiCopHinv2(w, cond.val, family, par, par2, obj, check.pars))
}
