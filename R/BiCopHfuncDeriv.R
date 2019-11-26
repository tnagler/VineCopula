#' Derivatives of the h-Function of a Bivariate Copula
#'
#' This function evaluates the derivative of a given conditional parametric
#' bivariate copula (h-function) with respect to its parameter(s) or one of its
#' arguments.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version \cr
#' \preformatted{BiCopHfuncDeriv(u1, u2, obj, deriv = "par")}
#' can be used.
#'
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param family integer; single number or vector of size `length(u1)`;
#' defines the bivariate copula family: \ \cr
#' `0` = independence copula \cr
#' `1` = Gaussian copula \cr
#' `2` = Student t copula (t-copula) \cr
#' `3` = Clayton copula \cr
#' `4` = Gumbel copula \cr
#' `5` = Frank copula \cr
#' `6` = Joe copula \cr
#' `13` = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' `14` = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' `16` = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' `23` = rotated Clayton copula (90 degrees) \cr
#' `24` = rotated Gumbel copula (90 degrees) \cr
#' `26` = rotated Joe copula (90 degrees) \cr
#' `33` = rotated Clayton copula (270 degrees) \cr
#' `34` = rotated Gumbel copula (270 degrees) \cr
#' `36` = rotated Joe copula (270 degrees) \cr
#' @param par numeric; single number or vector of size `length(u1)`;
#' copula parameter.
#' @param par2 integer; single number or vector of size `length(u1)`;
#' second parameter for the t-Copula; default is `par2 = 0`, should be an
#' positive integer for the Students's t copula `family = 2`.
#' @param deriv Derivative argument \cr
#' `"par"` = derivative with respect to the first parameter (default)\cr
#' `"par2"` = derivative with respect to the second parameter
#' (only available for the t-copula) \cr
#' `"u2"` = derivative with respect to the second argument `u2` \cr
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#' @return A numeric vector of the conditional bivariate copula derivative
#' \itemize{
#' \item of the copula `family`,
#' \item with parameter(s) `par`, `par2`,
#' \item with respect to `deriv`,
#' \item evaluated at `u1` and `u2`.
#' }
#' @author Ulf Schepsmeier
#' @seealso [RVineGrad()], [RVineHessian()],
#' [BiCopDeriv2()], [BiCopDeriv2()],
#' [BiCopHfuncDeriv()], [BiCop()]
#' @references Schepsmeier, U. and J. Stoeber (2014). Derivatives and Fisher
#' information of bivariate copulas. Statistical Papers, 55 (2), 525-542. \cr
#' <http://link.springer.com/article/10.1007/s00362-013-0498-x>.
#' @examples
#'
#' ## simulate from a bivariate Student-t copula
#' set.seed(123)
#' cop <- BiCop(family = 2, par = -0.7, par2 = 4)
#' simdata <- BiCopSim(100, cop)
#'
#' ## derivative of the conditional Student-t copula
#' ## with respect to the first parameter
#' u1 <- simdata[,1]
#' u2 <- simdata[,2]
#' BiCopHfuncDeriv(u1, u2, cop, deriv = "par")
#'
#' ## estimate a Student-t copula for the simulated data
#' cop <- BiCopEst(u1, u2, family = 2)
#' ## and evaluate the derivative of the conditional copula
#' ## w.r.t. the second argument u2
#' BiCopHfuncDeriv(u1, u2, cop, deriv = "u2")
#'
BiCopHfuncDeriv <- function(u1, u2, family, par, par2 = 0, deriv = "par", obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## check if specification is admissible for this function
    if (!all(family %in% c(0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)))
        stop("Copula family not implemented.")
    if ((deriv == "par2") && any(family != 2))
        stop("The derivative with respect to the second parameter can only be derived for the t-copula.")

    ## call C routines for specified 'deriv' case
    n <- args$n
    if (length(par) == 1) {
        ## call for single parameters
        if (deriv == "par") {
            if (family == 2) {
                out <- .C("diffhfunc_rho_tCopula",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(2),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diffhfunc_mod",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(family),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2") {
            out <- .C("diffhfunc_nu_tCopula_new",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(2),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "u2") {
            out <- .C("diffhfunc_v_mod",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else {
            stop("This kind of derivative is not implemented")
        }
    } else {
        # vectorized call
        if (deriv == "par") {
                out <- .C("diffhfunc_mod_vec",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(par),
                          as.double(par2),
                          as.integer(family),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2") {
            out <- .C("diffhfunc_nu_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "u2") {
            out <- .C("diffhfunc_v_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else {
            stop("This kind of derivative is not implemented")
        }
    }

    # reset NAs
    out <- reset_nas(out, args)
    # return results
    out
}
