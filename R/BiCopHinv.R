#' Inverse Conditional Distribution Function of a Bivariate Copula
#'
#' Evaluate the inverse conditional distribution function
#' (inverse h-function) of a given parametric bivariate copula.
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
#' If the family and parameter specification is stored in a \code{\link{BiCop}}
#' object \code{obj}, the alternative version
#' \preformatted{BiCopHinv(u1, u2, obj),
#' BiCopHinv1(u1, u2, obj),
#' BiCopHinv2(u1, u2, obj)}
#' can be used.
#'
#' @aliases BiCopHinv1 BiCopHinv2
#'
#' @param u1,u2 numeric vectors of equal length with values in [0,1].
#' @param family integer; single number or vector of size \code{length(u1)};
#' defines the bivariate copula family: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par numeric; single number or vector of size \code{length(u1)};
#' copula parameter.
#' @param par2 numeric; single number or vector of size \code{length(u1)};
#' second parameter for bivariate copulas with two parameters (t, BB1, BB6,
#' BB7, BB8, Tawn type 1 and type 2; default: \code{par2 = 0}). \code{par2}
#' should be an positive integer for the Students's t copula \code{family = 2}.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#'
#' @return \code{BiCopHinv} returns a list with
#' \item{hinv1}{Numeric vector of the inverse conditional distribution function
#' (inverse h-function) of the copula \code{family} with parameter(s)
#' \code{par}, \code{par2} evaluated at \code{u2} given \code{u1}, i.e.,
#' \eqn{h_1^{-1}(u_2|u_1;\boldsymbol{\theta})}{h_1^{-1}(u_2|u_1;\theta)}.}
#' \item{hinv2}{Numeric vector of the inverse conditional distribution function
#' (inverse h-function) of the copula \code{family} with parameter(s) \code{par},
#' \code{par2} evaluated at \code{u1} given \code{u2}, i.e.,
#' \eqn{h_2^{-1}(u_1|u_2;\boldsymbol{\theta})}{h_2^{-1}(u_1|u_2; \theta)}.}
#' \code{BiCopHinv1} is a faster version that only calculates \code{hinv1};
#' \code{BiCopHinv2} only calculates \code{hinv2}.
#'
#' @author Ulf Schepsmeier, Thomas Nagler
#'
#' @seealso \code{\link{BiCopHfunc}}, \code{\link{BiCopPDF}}, \code{\link{BiCopCDF}},
#' \code{\link{RVineLogLik}}, \code{\link{RVineSeqEst}}, \code{\link{BiCop}}
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
#' Pair-copula constructions of multiple dependence. Insurance: Mathematics and
#' Economics 44 (2), 182-198.
#'
#' @examples
#' # inverse h-functions of the Gaussian copula
#' cop <- BiCop(1, 0.5)
#' hi <- BiCopHinv(0.1, 0.2, cop)
#' \dontshow{hi}
#' # or using the fast versions
#' hi1 <- BiCopHinv1(0.1, 0.2, cop)
#' hi2 <- BiCopHinv2(0.1, 0.2, cop)
#' all.equal(hi$hinv1, hi1)
#' all.equal(hi$hinv2, hi2)
#'
#' # check if it is actually the inverse
#' cop <- BiCop(3, 3)
#' all.equal(0.2, BiCopHfunc1(0.1, BiCopHinv1(0.1, 0.2, cop), cop))
#' all.equal(0.1, BiCopHfunc2(BiCopHinv2(0.1, 0.2, cop), 0.2, cop))
#'
#' @export BiCopHinv
BiCopHinv <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate inverse h-functions
    n <- args$n
    if (length(par) == 1) {
        # call for single parameters
        hinv1 <- .C("Hinv1",                      # h(u2|u1)
                    as.integer(family),
                    as.integer(n),
                    as.double(u2),
                    as.double(u1),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2",                      # h(u1|u2)
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
        hinv1 <- .C("Hinv1_vec",                      # h(u2|u1)
                    as.integer(family),
                    as.integer(n),
                    as.double(u2),
                    as.double(u1),
                    as.double(par),
                    as.double(par2),
                    as.double(rep(0, n)),
                    PACKAGE = "VineCopula")[[7]]
        hinv2 <- .C("Hinv2_vec",                      # h(u1|u2)
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
    hinv1 <- reset_nas(hinv1, args)
    suppressWarnings(hinv2 <- reset_nas(hinv2, args))  # suppress second warning

    # return result
    list(hinv1 = hinv1, hinv2 = hinv2)
}

#' @rdname BiCopHinv
BiCopHinv1 <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate inverse h-function
    n <- args$n
    if (length(par) == 1) {
        # call for single parameters
        hinv1 <- .C("Hinv1",                      # h(u2|u1)
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
        hinv1 <- .C("Hinv1_vec",                      # h(u2|u1)
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
    hinv1 <- reset_nas(hinv1, args)
    # return result
    hinv1
}

#' @rdname BiCopHinv
BiCopHinv2 <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_u,
                    fix_nas,
                    check_if_01,
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate inverse h-functions
    n <- args$n
    if (length(par) == 1) {
        # call for single parameters
        hinv2 <- .C("Hinv2",                      # h(u1|u2)
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
        hinv2 <- .C("Hinv2_vec",                      # h(u1|u2)
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
    hinv2 <- reset_nas(hinv2, args)
    # return result
    hinv2
}

