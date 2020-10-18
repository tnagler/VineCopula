#' Transform C-Vine to R-Vine Structure
#'
#' This function transforms a C-vine structure from the package CDVine to the
#' corresponding R-vine structure.
#'
#'
#' @param order A d-dimensional vector specifying the order of the root nodes
#' in the C-vine.
#' @param family A d*(d-1)/2 vector of pair-copula families with values\cr
#' `0` = independence copula \cr `1` = Gaussian copula \cr `2` =
#' Student t copula (t-copula) \cr `3` = Clayton copula \cr `4` =
#' Gumbel copula \cr `5` = Frank copula \cr `6` = Joe copula \cr
#' `7` = BB1 copula \cr `8` = BB6 copula \cr `9` = BB7 copula
#' \cr `10` = BB8 copula \cr `13` = rotated Clayton copula (180
#' degrees; ``survival Clayton'') \cr `14` = rotated Gumbel copula (180
#' degrees; ``survival Gumbel'') \cr `16` = rotated Joe copula (180
#' degrees; ``survival Joe'') \cr `17` = rotated BB1 copula (180 degrees;
#' ``survival BB1'')\cr `18` = rotated BB6 copula (180 degrees; ``survival
#' BB6'')\cr `19` = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' `20` = rotated BB8 copula (180 degrees; ``survival BB8'')\cr `23`
#' = rotated Clayton copula (90 degrees) \cr `24` = rotated Gumbel copula
#' (90 degrees) \cr `26` = rotated Joe copula (90 degrees) \cr `27` =
#' rotated BB1 copula (90 degrees) \cr `28` = rotated BB6 copula (90
#' degrees) \cr `29` = rotated BB7 copula (90 degrees) \cr `30` =
#' rotated BB8 copula (90 degrees) \cr `33` = rotated Clayton copula (270
#' degrees) \cr `34` = rotated Gumbel copula (270 degrees) \cr `36` =
#' rotated Joe copula (270 degrees) \cr `37` = rotated BB1 copula (270
#' degrees) \cr `38` = rotated BB6 copula (270 degrees) \cr `39` =
#' rotated BB7 copula (270 degrees) \cr `40` = rotated BB8 copula (270
#' degrees) \cr `104` = Tawn type 1 copula \cr `114` = rotated Tawn
#' type 1 copula (180 degrees) \cr `124` = rotated Tawn type 1 copula (90
#' degrees) \cr `134` = rotated Tawn type 1 copula (270 degrees) \cr
#' `204` = Tawn type 2 copula \cr `214` = rotated Tawn type 2 copula
#' (180 degrees) \cr `224` = rotated Tawn type 2 copula (90 degrees) \cr
#' `234` = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par A d*(d-1)/2 vector of pair-copula parameters.
#' @param par2 A d*(d-1)/2 vector of second pair-copula parameters (optional;
#' default:\cr `par2 = rep(0,length(family))`), necessary for the t-, BB1,
#' BB6, BB7, BB8, Tawn type 1 and type 2 copulas.
#'
#' @return An [RVineMatrix()] object.
#'
#' @author Ulf Schepsmeier, Eike Brechmann
#'
#' @seealso [RVineMatrix()], [D2RVine()]
#'
#' @examples
#' # set up C-vine copula model with mixed pair-copulas
#' d <- 4
#' dd <- d*(d-1)/2
#' order <- 1:d
#' family <- c(1, 2, 3, 4, 7, 3)
#' par <- c(0.5, 0.4, 2, 1.5, 1.2, 1.5)
#' par2 <- c(0, 5, 0, 0, 2, 0)
#'
#' # transform to R-vine matrix notation
#' RVM <- C2RVine(order, family, par, par2)
C2RVine <- function(order, family, par, par2 = rep(0, length(family))) {
    # check input length
    dd <- length(family)
    d <- (1 + sqrt(1 + 8 * dd))/2
    if (dd < 1)
        stop("Length of 'family' has to be positive.")
    if (length(par) != length(par2))
        stop("Lengths of 'par' and 'par2' do not match.")
    if (length(par) != dd)
        stop("Lengths of 'family' and 'par' do not match.")
    if (length(order) != d)
        stop("Length of 'order' and dimension of the D-vine do not match.")

    # check parameters
    BiCopCheck(family, par, par2)

    Matrix <- matrix(rep(0, d * d), d, d)
    Copula.Params <- matrix(rep(0, d * d), d, d)
    Copula.Params2 <- matrix(rep(0, d * d), d, d)
    Copula.Types <- matrix(rep(0, d * d), d, d)

    # structure
    for (i in 1:d) {
        for (j in 1:(d - i + 1)) {
            Matrix[(d - i + 1), j] <- order[i]
        }
    }

    # copula properties
    k <- 1
    for (i in 1:(d - 1)) {
        for (j in 1:(d - i)) {
            Copula.Types[(d - i + 1), (d - j - i + 1)] <- family[k]
            Copula.Params[(d - i + 1), (d - j - i + 1)] <- par[k]
            Copula.Params2[(d - i + 1), (d - j - i + 1)] <- par2[k]
            k <- k + 1
        }
    }

    # return as RVineMatrix object
    RVineMatrix(Matrix = Matrix,
                family = Copula.Types,
                par = Copula.Params,
                par2 = Copula.Params2)
}
