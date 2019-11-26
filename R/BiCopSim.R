#' Simulation from a Bivariate Copula
#'
#' This function simulates from a given parametric bivariate copula.
#'
#' If the family and parameter specification is stored in a [BiCop()]
#' object `obj`, the alternative version
#' \preformatted{BiCopSim(N, obj)}
#' can be used.
#'
#' @param N Number of bivariate observations simulated.
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
#' @return An `N` x 2 matrix of data simulated from the bivariate copula
#' with `family` and parameter(s) `par`, `par2`.
#'
#' @author Ulf Schepsmeier
#'
#' @seealso
#' [BiCop()],
#' [RVineSim()]
#' @examples
#' # simulate from a bivariate t-copula
#' simdata <- BiCopSim(100, 2, -0.7, par2 = 4)
#'
#' # or alternatively
#' obj <- BiCop(family = 2, par = -0.7, par2 = 4)
#' simdata2 <- BiCopSim(100, obj)
#'
BiCopSim <- function(N, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## stop if lengths are not 1 or N
    if (!(length(family) %in% c(1, N)))
        stop("'family' has to be a single number or a size N vector")
    if (!(length(par) %in% c(1, N)))
        stop("'par' has to be a single number or a size N vector")
    if (!(length(par2) %in% c(1, N)))
        stop("'par2' has to be a single number or a size N vector")

    ## start with independent uniforms (byrow for backwards compatibility)
    w <- matrix(runif(2*N), ncol = 2, byrow = TRUE)

    ## simulate from copula by inverse rosenblatt transform
    if (length(par) == 1) {
        # call for single parameters
        tmp <- .C("Hinv1",
                  as.integer(family),
                  as.integer(N),
                  as.double(w[, 2]),
                  as.double(w[, 1]),
                  as.double(par),
                  as.double(par2),
                  as.double(rep(0, N)),
                  PACKAGE = "VineCopula")[[7]]
    } else {
        # vectorized call
        tmp <- .C("Hinv1_vec",
                  as.integer(family),
                  as.integer(N),
                  as.double(w[, 2]),
                  as.double(w[, 1]),
                  as.double(par),
                  as.double(par2),
                  as.double(rep(0, N)),
                  PACKAGE = "VineCopula")[[7]]
    }

    ## return results
    U <- matrix(c(w[, 1], tmp), ncol = 2)
    U
}
