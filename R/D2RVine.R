#' Transform D-Vine to R-Vine Structure
#'
#' This function transforms a D-vine structure from the package CDVine to the
#' corresponding R-vine structure.
#'
#'
#' @param order A d-dimensional vector specifying the order of the nodes in the
#' D-vine.
#' @param family A d*(d-1)/2 vector of pair-copula families with values\cr
#' \code{0} = independence copula \cr \code{1} = Gaussian copula \cr \code{2} =
#' Student t copula (t-copula) \cr \code{3} = Clayton copula \cr \code{4} =
#' Gumbel copula \cr \code{5} = Frank copula \cr \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr \code{8} = BB6 copula \cr \code{9} = BB7 copula
#' \cr \code{10} = BB8 copula \cr \code{13} = rotated Clayton copula (180
#' degrees; ``survival Clayton'') \cr \code{14} = rotated Gumbel copula (180
#' degrees; ``survival Gumbel'') \cr \code{16} = rotated Joe copula (180
#' degrees; ``survival Joe'') \cr \code{17} = rotated BB1 copula (180 degrees;
#' ``survival BB1'')\cr \code{18} = rotated BB6 copula (180 degrees; ``survival
#' BB6'')\cr \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr \code{23}
#' = rotated Clayton copula (90 degrees) \cr \code{24} = rotated Gumbel copula
#' (90 degrees) \cr \code{26} = rotated Joe copula (90 degrees) \cr \code{27} =
#' rotated BB1 copula (90 degrees) \cr \code{28} = rotated BB6 copula (90
#' degrees) \cr \code{29} = rotated BB7 copula (90 degrees) \cr \code{30} =
#' rotated BB8 copula (90 degrees) \cr \code{33} = rotated Clayton copula (270
#' degrees) \cr \code{34} = rotated Gumbel copula (270 degrees) \cr \code{36} =
#' rotated Joe copula (270 degrees) \cr \code{37} = rotated BB1 copula (270
#' degrees) \cr \code{38} = rotated BB6 copula (270 degrees) \cr \code{39} =
#' rotated BB7 copula (270 degrees) \cr \code{40} = rotated BB8 copula (270
#' degrees) \cr \code{104} = Tawn type 1 copula \cr \code{114} = rotated Tawn
#' type 1 copula (180 degrees) \cr \code{124} = rotated Tawn type 1 copula (90
#' degrees) \cr \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr \code{214} = rotated Tawn type 2 copula
#' (180 degrees) \cr \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par A d*(d-1)/2 vector of pair-copula parameters.
#' @param par2 A d*(d-1)/2 vector of second pair-copula parameters (optional;
#' default:\cr \code{par2 = rep(0,length(family))}), necessary for the t-, BB1,
#' BB6, BB7, BB8, Tawn type 1 and type 2 copulas.
#' @return An \code{\link{RVineMatrix}} object.
#' @author Ulf Schepsmeier
#' @seealso \code{\link{RVineMatrix}}, \code{\link{C2RVine}}
#' @examples
#'
#' # set up D-vine copula model with mixed pair-copulas
#' d <- 4
#' dd <- d*(d-1)/2
#' order <- 1:d
#' family <- c(1, 2, 3, 4, 7, 3)
#' par <- c(0.5, 0.4, 2, 1.5, 1.2, 1.5)
#' par2 <- c(0, 5, 0, 0, 2, 0)
#'
#' # transform to R-vine matrix notation
#' RVM <- D2RVine(order, family, par, par2)
#'
#' \dontrun{
#' # load package CDVine for comparison
#' library(CDVine)
#'
#' # simulate a sample of size 500 from a 4-dimensional D-vine
#' type <- 2  # D-vine
#' simdata <- CDVineSim(500, family, par, par2, type)
#'
#' # determine log-likelihood
#' out <- CDVineLogLik(simdata, family, par, par2, type)
#' out$loglik
#'
#' # check that log-likelihood stays the same
#' out2 <- RVineLogLik(simdata, RVM)
#' out2$loglik
#' }
#'
#' @export D2RVine
D2RVine <- function(order, family, par, par2 = rep(0, length(family))) {
    dd <- length(family)
    if (dd < 1)
        stop("Length of 'family' has to be positive.")
    if (length(par) != length(par2))
        stop("Lengths of 'par' and 'par2' do not match.")
    if (length(par) != dd)
        stop("Lengths of 'family' and 'par' do not match.")

    for (i in 1:dd) {
        if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40)))
            stop("Copula family not implemented.")
        # Parameterbereiche abfragen
        if ((family[i] == 1 || family[i] == 2) && abs(par[i]) >= 1)
            stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
        if (family[i] == 2 && par2[i] <= 2)
            stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
        if ((family[i] == 3 || family[i] == 13) && par[i] <= 0)
            stop("The parameter of the Clayton copula has to be positive.")
        if ((family[i] == 4 || family[i] == 14) && par[i] < 1)
            stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
        if ((family[i] == 6 || family[i] == 16) && par[i] <= 1)
            stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
        if (family[i] == 5 && par[i] == 0)
            stop("The parameter of the Frank copula has to be unequal to 0.")
        if ((family[i] == 7 || family[i] == 17) && par[i] <= 0)
            stop("The first parameter of the BB1 copula has to be positive.")
        if ((family[i] == 7 || family[i] == 17) && par2[i] < 1)
            stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
        if ((family[i] == 8 || family[i] == 18) && par[i] <= 0)
            stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
        if ((family[i] == 8 || family[i] == 18) && par2[i] < 1)
            stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
        if ((family[i] == 9 || family[i] == 19) && par[i] < 1)
            stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
        if ((family[i] == 9 || family[i] == 19) && par2[i] <= 0)
            stop("The second parameter of the BB7 copula has to be positive.")
        if ((family[i] == 10 || family[i] == 20) && par[i] < 1)
            stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
        if ((family[i] == 10 || family[i] == 20) && (par2[i] <= 0 || par2[i] > 1))
            stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
        if ((family[i] == 23 || family[i] == 33) && par[i] >= 0)
            stop("The parameter of the rotated Clayton copula has to be negative.")
        if ((family[i] == 24 || family[i] == 34) && par[i] > -1)
            stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
        if ((family[i] == 26 || family[i] == 36) && par[i] >= -1)
            stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
        if ((family[i] == 27 || family[i] == 37) && par[i] >= 0)
            stop("The first parameter of the rotated BB1 copula has to be negative.")
        if ((family[i] == 27 || family[i] == 37) && par2[i] > -1)
            stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 28 || family[i] == 38) && par[i] >= 0)
            stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 28 || family[i] == 38) && par2[i] > -1)
            stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 29 || family[i] == 39) && par[i] > -1)
            stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 29 || family[i] == 39) && par2[i] >= 0)
            stop("The second parameter of the rotated BB7 copula has to be negative.")
        if ((family[i] == 30 || family[i] == 40) && par[i] > -1)
            stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
        if ((family[i] == 30 || family[i] == 40) && (par2[i] >= 0 || par2[i] < (-1)))
            stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
    }

    d <- (1 + sqrt(1 + 8 * dd))/2

    if (length(order) != d)
        stop("Length of 'order' and dimension of the D-vine do not match.")

    Matrix <- matrix(rep(0, d * d), d, d)
    Copula.Params <- matrix(rep(0, d * d), d, d)
    Copula.Params2 <- matrix(rep(0, d * d), d, d)
    Copula.Types <- matrix(rep(0, d * d), d, d)

    # structure diagonale
    for (i in 1:d) {
        Matrix[(d - i + 1), (d - i + 1)] <- order[i]
    }

    # below the diagonal + copula properties
    k <- 1
    for (i in 1:(d - 1)) {
        for (j in 1:(d - i)) {
            Matrix[(d - i + 1), (d - j - i + 1)] <- order[j]
            Copula.Types[(d - i + 1), (d - j - i + 1)] <- family[k]
            Copula.Params[(d - i + 1), (d - j - i + 1)] <- par[k]
            Copula.Params2[(d - i + 1), (d - j - i + 1)] <- par2[k]
            k <- k + 1
        }
    }

    RVM <- RVineMatrix(Matrix = Matrix,
                       family = Copula.Types,
                       par = Copula.Params,
                       par2 = Copula.Params2)

    return(RVM)
}
