#' Density of a Bivariate Copula
#'
#' This function evaluates the probability density function (PDF) of a given
#' parametric bivariate copula.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version \cr
#' \preformatted{BiCopPDF(u1, u2, obj)}
#' can be used.
#'
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param family integer; single number or vector of size `length(u1)`;
#' defines the bivariate copula family: \cr
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
#' @param par numeric; single number or vector of size `length(u1)`;
#' copula parameter.
#' @param par2 numeric; single number or vector of size `length(u1)`;
#' second parameter for bivariate copulas with two parameters (t, BB1, BB6,
#' BB7, BB8, Tawn type 1 and type 2; default: `par2 = 0`). `par2`
#' should be an positive integer for the Students's t copula `family = 2`.
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#'
#' @return A numeric vector of the bivariate copula density
#' \itemize{
#' \item of the copula `family`
#' \item with parameter(s) `par`, `par2`
#' \item evaluated at `u1` and `u2`.
#' }
#'
#' @author Eike Brechmann
#'
#' @seealso [BiCopCDF()], [BiCopHfunc()],
#' [BiCopSim()], [BiCop()]
#'
#' @examples
#' \dontshow{set.seed(123)}
#' ## simulate from a bivariate Student-t copula
#' cop <- BiCop(family = 2, par = -0.7, par2 = 4)
#' simdata <- BiCopSim(100, cop)
#'
#' ## evaluate the density of the bivariate t-copula
#' u1 <- simdata[,1]
#' u2 <- simdata[,2]
#' BiCopPDF(u1, u2, cop)
#'
#' ## select a bivariate copula for the simulated data
#' fit <- BiCopSelect(u1, u2)
#' summary(fit)
#' ## and evaluate its PDF
#' round(BiCopPDF(u1, u2, fit), 3)
#'
BiCopPDF <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## evaluate log-density
    n <- args$n
    if (length(par) == 1) {
        if (family == 1004) {
            family <- ifelse(par >= 0, 4, 24)
            par <- sign(par) * (1 + abs(par))
        }

        # unvectorized call
        coplik <- .C("PDF_seperate",
                     as.integer(family),
                     as.integer(n),
                     as.double(u1),
                     as.double(u2),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        coplik <- .C("PDF_seperate_vec",
                     as.integer(family),
                     as.integer(n),
                     as.double(u1),
                     as.double(u2),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
    }

    # reset NAs
    out <- reset_nas(coplik, args)
    # return result
    out
}
