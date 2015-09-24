#' Kendall's Tau Values of an R-Vine Copula Model
#' 
#' This function computes the values of Kendall's tau corresponding to the
#' parameters of an R-vine copula model.
#' 
#' 
#' @param RVM An \code{\link{RVineMatrix}} object.
#' @return Matrix with the same structure as the family and parameter matrices
#' of the \code{\link{RVineMatrix}} object \code{RVM} where the entries are
#' values of Kendall's tau corresponding to the families and parameters of the
#' R-vine copula model given by \code{RVM}.
#' @author Jeffrey Dissmann
#' @seealso \code{\link{RVineMatrix}}, \code{\link{BiCopPar2Tau}}
#' @examples
#' 
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#' 
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#' 
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#' 
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#' 
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#' 
#' # compute the Kendall's tau values
#' tau <- RVinePar2Tau(RVM)
#' 
#' @export RVinePar2Tau
RVinePar2Tau <- function(RVM) {
    
    if (!is(RVM, "RVineMatrix")) 
        stop("'RVM' has to be an RVineMatrix object.")
    
    taus <- RVM$par
    n <- dim(RVM)
    
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            taus[i, j] <- BiCopPar2Tau(RVM$family[i, j],
                                       RVM$par[i, j],
                                       RVM$par2[i, j])
        }
    }
    
    return(taus)
}
