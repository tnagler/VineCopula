#' Conditional Distribution Function of a Bivariate Copula
#'
#' Evaluate the conditional distribution function (h-function)
#' of a given parametric bivariate copula.
#'
#' The h-function is defined as the conditional distribution function of a
#' bivariate copula, i.e.,
#' \deqn{h_1(u_2|u_1;\boldsymbol{\theta}) := P(U_2 \le u_2 | U_1 = u_1)
#' = \frac{\partial C(u_1, u_2; \boldsymbol{\theta})}{\partial u_1}, }{
#' h_1(u_2|u_1,\theta) :=  P(U_2 \le u_2 | U_1 = u_1)
#' = \partial C(u_1,u_2) / \partial u_1, }
#' \deqn{h_2(u_1|u_2;\boldsymbol{\theta}) := P(U_1 \le u_1 | U_2 = u_2)
#'  = \frac{\partial C(u_1, u_2; \boldsymbol{\theta})}{\partial u_2}, }{
#' h_2(u_1|u_2,\theta) := P(U_1 \le u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2, }
#' where \eqn{(U_1, U_2) \sim C}, and \eqn{C} is a bivariate copula distribution
#' function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}.
#' For more details see Aas et al. (2009). \cr \cr
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative versions
#' \preformatted{BiCopHfunc(u1, u2, obj)
#' BiCopHfunc1(u1, u2, obj)
#' BiCopHfunc2(u1, u2, obj)}
#' can be used.
#'
#' @aliases BiCopHfunc1 BiCopHfunc2
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
#' @return `BiCopHfunc` returns a list with
#' \item{hfunc1}{Numeric vector of the conditional distribution
#' function (h-function) of the copula `family` with parameter(s)
#' `par`, `par2` evaluated at `u2` given `u1`, i.e.,
#' \eqn{h_1(u_2|u_1;\boldsymbol{\theta})}{h_1(u_2|u_1;\theta)}.}
#' \item{hfunc2}{Numeric vector of the conditional distribution function
#' (h-function) of the copula `family` with parameter(s) `par`,
#' `par2` evaluated at `u1` given `u2`, i.e.,
#' \eqn{h_2(u_1|u_2;\boldsymbol{\theta})}{h_2(u_1|u_2; \theta)}.}
#' `BiCopHfunc1` is a faster version that only calculates `hfunc1`;
#' `BiCopHfunc2` only calculates `hfunc2`.
#'
#' @author Ulf Schepsmeier
#'
#' @seealso [BiCopHinv()], [BiCopPDF()], [BiCopCDF()],
#' [RVineLogLik()], [RVineSeqEst()], [BiCop()]
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
#' Pair-copula constructions of multiple dependence. Insurance: Mathematics and
#' Economics 44 (2), 182-198.
#'
#' @examples
#' data(daxreturns)
#'
#' # h-functions of the Gaussian copula
#' cop <- BiCop(family = 1, par = 0.5)
#' h <- BiCopHfunc(daxreturns[, 2], daxreturns[, 1], cop)
#' \dontshow{h}
#' # or using the fast versions
#' h1 <- BiCopHfunc1(daxreturns[, 2], daxreturns[, 1], cop)
#' h2 <- BiCopHfunc2(daxreturns[, 2], daxreturns[, 1], cop)
#' all.equal(h$hfunc1, h1)
#' all.equal(h$hfunc2, h2)
#'
BiCopHfunc <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate h-functions
    n <- args$n
    if (length(par) == 1) {
        # call for single parameters
        hfunc1 <- .C("Hfunc1",                      # h(u2|u1)
                     as.integer(family),
                     as.integer(n),
                     as.double(u2),
                     as.double(u1),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
        hfunc2 <- .C("Hfunc2",                      # h(u1|u2)
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
        hfunc1 <- .C("Hfunc1_vec",                      # h(u2|u1)
                     as.integer(family),
                     as.integer(n),
                     as.double(u2),
                     as.double(u1),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
        hfunc2 <- .C("Hfunc2_vec",                      # h(u1|u2)
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
    hfunc1 <- reset_nas(hfunc1, args)
    suppressWarnings(hfunc2 <- reset_nas(hfunc2, args))  # suppress second warning

    # return result
    list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}


#' @rdname BiCopHfunc
BiCopHfunc1 <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate h-function
    n <- args$n
    if (length(par) == 1) {
        if (family == 1004) {
            family <- ifelse(par >= 0, 4, 24)
            par <- sign(par) * (1 + abs(par))
        }
        # call for single parameters
        hfunc1 <- .C("Hfunc1",                      # h(u2|u1)
                     as.integer(family),
                     as.integer(n),
                     as.double(u2),
                     as.double(u1),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        hfunc1 <- .C("Hfunc1_vec",                      # h(u2|u1)
                     as.integer(family),
                     as.integer(n),
                     as.double(u2),
                     as.double(u1),
                     as.double(par),
                     as.double(par2),
                     as.double(rep(0, n)),
                     PACKAGE = "VineCopula")[[7]]
    }

    # reset NAs
    hfunc1 <- reset_nas(hfunc1, args)
    # return result
    hfunc1
}

#' @rdname BiCopHfunc
BiCopHfunc2 <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate h-function
    n <- args$n
    if (length(par) == 1) {
        if (family == 1004) {
            family <- ifelse(par >= 0, 4, 24)
            par <- sign(par) * (1 + abs(par))
        }
        # call for single parameters
        hfunc2 <- .C("Hfunc2",                      # h(u1|u2)
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
        hfunc2 <- .C("Hfunc2_vec",                      # h(u1|u2)
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
    hfunc2 <- reset_nas(hfunc2, args)
    # return results
    hfunc2
}

