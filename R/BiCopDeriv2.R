#' Second Derivatives of a Bivariate Copula Density
#'
#' This function evaluates the second derivative of a given parametric
#' bivariate copula density with respect to its parameter(s) and/or its
#' arguments.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version \cr
#' \preformatted{BiCopDeriv2(u1, u2, obj, deriv = "par")}
#' can be used.
#'
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param family integer; single number or vector of size `length(u1)`;
#' defines the bivariate copula family:  \cr
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
#' @param par Copula parameter.
#' @param par2 integer; single number or vector of size `length(u1)`;
#' second parameter for the t-Copula; default is `par2 = 0`, should be an
#' positive integer for the Students's t copula `family = 2`.
#' @param deriv Derivative argument \cr
#' `"par"` = second derivative with respect to
#' the first parameter (default)\cr
#' `"par2"` = second derivative with respect to
#' the second parameter (only available for the t-copula) \cr
#' `"u1"` = second derivative with respect to
#' the first argument `u1` \cr
#' `"u2"` = second derivative with respect to
#' the second argument `u2` \cr
#' `"par1par2"` = second derivative with respect to
#' the first and second parameter (only available for the t-copula)
#' \cr `"par1u1"` = second derivative with respect to
#' the first parameter and the first argument \cr
#' `"par2u1"` = second derivative with respect to the
#' second parameter and the first argument (only available for the t-copula) \cr
#' `"par1u2"` = second derivative with respect to
#' the first parameter and the second argument \cr
#' `"par2u2"` = second derivative with respect to
#' the second parameter and the second argument
#' (only available for the t-copula) \cr
#' @param obj `BiCop` object containing the family and parameter
#' specification.
#' @param check.pars logical; default is `TRUE`; if `FALSE`, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#' @return A numeric vector of the second-order bivariate copula derivative
#' \itemize{
#' \item of the copula `family`
#' \item with parameter(s) `par`, `par2`
#' \item with respect to `deriv`
#' \item evaluated at `u1` and `u2`.
#' }
#' @author Ulf Schepsmeier, Jakob Stoeber
#' @seealso [RVineGrad()], [RVineHessian()],
#' [BiCopDeriv()], [BiCopHfuncDeriv()], [BiCop()]
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
#' ## second derivative of the Student-t copula w.r.t. the first parameter
#' u1 <- simdata[,1]
#' u2 <- simdata[,2]
#' BiCopDeriv2(u1, u2, cop, deriv = "par")
#'
#' ## estimate a Student-t copula for the simulated data
#' cop <- BiCopEst(u1, u2, family = 2)
#' ## and evaluate its second derivative w.r.t. the second argument u2
#' BiCopDeriv2(u1, u2, cop, deriv = "u2")
#'
BiCopDeriv2 <- function(u1, u2, family, par, par2 = 0, deriv = "par", obj = NULL, check.pars = TRUE) {
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
    if ((deriv %in% c("par2", "par1par2", "par2u1", "par2u2")) && any(family != 2))
        stop("The derivative with respect to the second parameter can only be derived for the t-copula.")

    ## call C routines for specified 'deriv' case
    n <- args$n
    if (length(par) == 1) {
        ## call for single parameters
        if (deriv == "par") {
            if (family == 2) {
                out <- .C("diff2PDF_rho_tCopula",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(2),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diff2PDF_mod",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(family),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2") {
            out <- .C("diff2PDF_nu_tCopula_new",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(2),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "u1") {
            out <- .C("diff2PDF_u_mod",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "u2") {
            out <- .C("diff2PDF_v_mod",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "par1par2") {
            out <- .C("diff2PDF_rho_nu_tCopula_new",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(2),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "par1u1") {
            if (family == 2) {
                out <- .C("diff2PDF_rho_u_tCopula_new",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(2),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diff2PDF_par_u_mod",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(family),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2u1") {
            out <- .C("diff2PDF_nu_u_tCopula_new",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(2),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else if (deriv == "par1u2") {
            if (family == 2) {
                out <- .C("diff2PDF_rho_v_tCopula_new",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(2),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            } else {
                out <- .C("diff2PDF_par_v_mod",
                          as.double(u1),
                          as.double(u2),
                          as.integer(n),
                          as.double(c(par, par2)),
                          as.integer(family),
                          as.double(rep(0, n)),
                          PACKAGE = "VineCopula")[[6]]
            }
        } else if (deriv == "par2u2") {
            out <- .C("diff2PDF_nu_v_tCopula_new",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(c(par, par2)),
                      as.integer(2),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[6]]
        } else {
            stop("This kind of derivative is not implemented")
        }
    } else {
        ## vectorized calls
        if (deriv == "par") {
            out <- .C("diff2PDF_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2") {
            out <- .C("diff2PDF_nu_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "u1") {
            out <- .C("diff2PDF_u_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "u2") {
            out <- .C("diff2PDF_v_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par1par2") {
            out <- .C("diff2PDF_rho_nu_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par1u1") {
            out <- .C("diff2PDF_par_u_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2u1") {
            out <- .C("diff2PDF_nu_u_tCopula_new_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par1u2") {
            out <- .C("diff2PDF_par_v_mod_vec",
                      as.double(u1),
                      as.double(u2),
                      as.integer(n),
                      as.double(par),
                      as.double(par2),
                      as.integer(family),
                      as.double(rep(0, n)),
                      PACKAGE = "VineCopula")[[7]]
        } else if (deriv == "par2u2") {
            out <- .C("diff2PDF_nu_v_tCopula_new_vec",
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
    # return result
    out
}
