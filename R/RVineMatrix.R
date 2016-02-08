#' R-Vine Copula Model in Matrix Notation
#'
#' This function creates an \code{\link{RVineMatrix}} object which encodes an
#' R-vine copula model. It contains the matrix identifying the R-vine tree
#' structure, the matrix identifying the copula families utilized and two
#' matrices for corresponding parameter values.
#'
#'
#' @param Matrix Lower (or upper) triangular d x d matrix that defines the
#' R-vine tree structure.
#' @param family Lower (or upper) triangular d x d matrix with zero diagonal
#' entries that assigns the pair-copula families to each (conditional) pair
#' defined by \code{Matrix} (default: \code{family =
#' array(0,dim=dim(Matrix))}).  The bivariate copula families are defined as
#' follows:\cr
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
#' @param par Lower (or upper) triangular d x d matrix with zero diagonal
#' entries that assigns the (first) pair-copula parameter to each (conditional)
#' pair defined by \code{Matrix} \cr (default: \code{par = array(NA, dim =
#' dim(Matrix))}).
#' @param par2 Lower (or upper) triangular d x d matrix with zero diagonal
#' entries that assigns the second parameter for pair-copula families with two
#' parameters to each (conditional) pair defined by \code{Matrix} (default:
#' \code{par2 = array(NA, dim = dim(Matrix))}).
#' @param names A vector of names for the d variables; default: \code{names =
#' NULL}.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#'
#' @return An \code{\link{RVineMatrix}} object with the following matrix
#' components: \item{Matrix}{R-vine tree structure matrix.}
#' \item{family}{Pair-copula family matrix with values as above.}
#' \item{par}{Pair-copula parameter matrix.} \item{par2}{Second pair-copula
#' parameter matrix with parameters necessary for pair-copula families with two
#' parameters.}
#'
#' @note The \code{print} function writes the R-vine matrix defined by
#' \code{Matrix}. A detailed output is given by \code{print(RVM, detail=TRUE)},
#' where \code{RVM} is the \code{\link{RVineMatrix}} object. \cr The
#' \code{\link{RVineMatrix}} function automatically checks if the given matrix
#' is a valid R-vine matrix (see \code{\link{RVineMatrixCheck}}). \cr Although
#' the function allows upper triangular matrices as its input, it will always
#' store them as lower triangular matrices.
#'
#' @author Jeffrey Dissmann, Thomas Nagler
#'
#' @seealso
#' \code{\link{RVineMatrixCheck}},
#' \code{\link{RVineMLE}},
#' \code{\link{RVineSim}},
#' \code{\link{C2RVine}},
#' \code{\link{D2RVine}}
#'
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#'
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
#' # Print detailed information
#' print(RVM, detail = TRUE)
#'
RVineMatrix <- function(Matrix,
                        family = array(0, dim = dim(Matrix)),
                        par = array(NA, dim = dim(Matrix)),
                        par2 = array(NA, dim = dim(Matrix)),
                        names = NULL, check.pars = TRUE) {

    ## set NAs to zero
    Matrix[is.na(Matrix)] <- 0
    family[is.na(family)] <- 0
    par[is.na(par)] <- 0
    par2[is.na(par2)] <- 0

    ## convert to lower triangular matrix if necessary
    Matrix <- ToLowerTri(Matrix)
    family <- ToLowerTri(family)
    par    <- ToLowerTri(par)
    par2   <- ToLowerTri(par2)

    ## set upper triangle to zero
    family[upper.tri(family, diag = T)] <- 0
    par[upper.tri(par, diag = T)] <- 0
    par2[upper.tri(par2, diag = T)] <- 0

    ## check matrices
    if (dim(Matrix)[1] != dim(Matrix)[2])
        stop("Structure matrix has to be quadratic.")
    if (any(par != NA) & dim(par)[1] != dim(par)[2])
        stop("Parameter matrix has to be quadratic.")
    if (any(par2 != NA) & dim(par2)[1] != dim(par2)[2])
        stop("Second parameter matrix has to be quadratic.")
    if (any(family != 0) & dim(family)[1] != dim(family)[2])
        stop("Copula family matrix has to be quadratic.")
    if (max(Matrix) > dim(Matrix)[1])
        stop("Error in the structure matrix.")
    if (length(names) > 0 & length(names) != dim(Matrix)[1])
        stop("Length of the vector 'names' is not correct.")
    if (RVineMatrixCheck(Matrix) != 1)
        stop("'Matrix' is not a valid R-vine matrix")

    ## check for family/parameter consistency
    sel <- lower.tri(family)  # selector for lower triangular matrix
    if (check.pars & (any(family != 0) | any(!is.na(par)))) {
        BiCopCheck(family[sel], par[sel], par2[sel])
    }

    ## create help matrices
    MaxMat <- createMaxMat(Matrix)
    CondDistr <- neededCondDistr(Matrix)

    ## create RVineMatrix object
    RVM <- list(Matrix = Matrix,
                family = family,
                par = par,
                par2 = par2,
                names = names,
                MaxMat = MaxMat,
                CondDistr = CondDistr)
    class(RVM) <- "RVineMatrix"

    ## add dependence measures
    # create list of BiCop ojbects
    objlst <- apply(cbind(family[sel], par[sel], par2[sel]),
                    1,
                    function(x) BiCop(x[1], x[2], x[3], check.pars = FALSE))
    # construct dependence measure matrices
    taus <- utds <- ltds <- bets <- matrix(0, nrow(Matrix), ncol(Matrix))
    taus[sel] <- sapply(objlst, function(x) x$tau)
    utds[sel] <- sapply(objlst, function(x) x$taildep$upper)
    ltds[sel] <- sapply(objlst, function(x) x$taildep$lower)
    bets[sel] <- sapply(objlst, function(x) x$beta)
    RVM$tau <- taus
    RVM$taildep$upper <- utds
    RVM$taildep$lower <- ltds
    RVM$beta <- bets

    ## return results
    RVM
}

normalizeRVineMatrix <- function(RVM) {

    oldOrder <- diag(RVM$Matrix)
    Matrix <- reorderRVineMatrix(RVM$Matrix)

    return(RVineMatrix(Matrix,
                       RVM$family,
                       RVM$par,
                       RVM$par2,
                       names = rev(RVM$names[oldOrder]),
                       check.pars = FALSE))
}

reorderRVineMatrix <- function(Matrix) {
    oldOrder <- diag(Matrix)

    O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)

    for (i in 1:nrow(Matrix)) {
        Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
    }

    return(Matrix)
}

# exported version of normalizeRVineMatrix


#' Normalization of R-Vine Matrix
#'
#' An \code{\link{RVineMatrix}} is permuted to achieve a natural ordering (i.e.
#' \code{diag(RVM$Matrix) == d:1})
#'
#'
#' @param RVM \code{\link{RVineMatrix}} defining the R-vine structure
#' @return \item{RVM}{An \code{\link{RVineMatrix}} in natural ordering with
#' entries in \code{RVM$names} keeping track of the reordering.}
#' @keywords vine
#' @examples
#'
#' Matrix <- matrix(c(5, 2, 3, 1, 4,
#'                    0, 2, 3, 4, 1,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 4, 1,
#'                    0, 0, 0, 0, 1), 5, 5)
#' family <- matrix(1,5,5)
#'
#' par <- matrix(c(0, 0.2, 0.9, 0.5, 0.8,
#'                 0,   0, 0.1, 0.6, 0.9,
#'                 0,   0,   0, 0.7, 0.5,
#'                 0,   0,   0,   0, 0.8,
#'                 0,   0,   0,   0,   0), 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix, family, par)
#'
#' # normalise the RVine
#' RVineMatrixNormalize(RVM)
#'
#' @export RVineMatrixNormalize
RVineMatrixNormalize <- function(RVM) {
    stopifnot(is(RVM, "RVineMatrix"))

    if (is.null(RVM$names))
        RVM$names <- paste("V", 1:nrow(RVM$Matrix), sep = "")
    oldOrder <- diag(RVM$Matrix)

    return(normalizeRVineMatrix(RVM))
}

dim.RVineMatrix <- function(x) {
    RVine <- x
    return(dim(RVine$Matrix)[1])
    NextMethod("dim")
}

print.RVineMatrix <- function(x, ...) {
    RVine <- x
    cat("R-vine copula with the following pair-copulas:\n")
    d <- dim(RVine)
    cat("")
    cat("Tree 1:\n")
    for (i in 1:(d - 1)) {
        a <- paste(RVine$Matrix[i, i],
                   ",",
                   RVine$Matrix[d, i],
                   sep = "")
        a <- paste(a,
                   "   ",
                   BiCopName(RVine$family[d, i], short = FALSE),
                   sep = "")
        if (RVine$family[d, i] != 0) {
            a <- paste(a, " (par = ", round(RVine$par[d, i], 2), sep = "")
            if (RVine$family[d, i] %in% c(2, 7, 8, 9, 10,
                                          17, 18, 19, 20,
                                          27, 28, 29, 30,
                                          37, 38, 39, 40,
                                          104, 114, 124, 134,
                                          204, 214, 224, 234)) {
                a <- paste(a, ", par2 = ", round(RVine$par2[d, i], 2), sep = "")
            }
            a <- paste(a,
                       ", tau = ",
                       round(BiCopPar2Tau(RVine$family[d, i],
                                          RVine$par[d, i],
                                          RVine$par2[d, i]), 2),
                       ")\n",
                       sep = "")
        }
        cat(a)
    }
    cat("\n")
    for (j in 2:(d - 1)) {
        cat("")
        a <- paste("Tree ", j, ":\n", sep = "")
        cat(a)
        for (i in 1:(d - j)) {
            a <- paste(RVine$Matrix[i, i],
                       ",",
                       RVine$Matrix[d - j + 1, i],
                       sep = "")
            a <- paste(a, ";", sep = "")
            conditioningSet <- (d - j + 2):d
            for (k in 1:length(conditioningSet)) {
                if (k > 1) {
                    a <- paste(a, ",", sep = "")
                }
                a <- paste(a,
                           RVine$Matrix[conditioningSet[k], i],
                           sep = "")
            }
            a <- paste(a,
                       "   ",
                       BiCopName(RVine$family[d - j + 1, i], short = FALSE),
                       sep = "")
            if (RVine$family[d - j + 1, i] != 0) {
                a <- paste(a,
                           " (par = ",
                           round(RVine$par[d - j + 1, i], 2),
                           sep = "")
                if (RVine$family[d - j + 1, i] %in% c(2, 7, 8, 9, 10,
                                                      17, 18, 19, 20,
                                                      27, 28, 29, 30,
                                                      37, 38, 39, 40,
                                                      104, 114, 124, 134,
                                                      204, 214, 224, 234)) {
                    a <- paste(a,
                               ", par2 = ",
                               round(RVine$par2[d - j + 1, i], 2),
                               sep = "")
                }
                a <- paste(a,
                           ", tau = ",
                           round(RVine$tau[d - j + 1, i], 2),
                           ")\n",
                           sep = "")
            }
            cat(a)
        }
        if (j < d - 1) cat("\n")
    }
    # show names if provided
    if (!is.null(RVine$names)) {
        cat("\n")
        cat("Where  ")
        for (i in 1:length(RVine$names)) {
            cat(i, "<->", RVine$names[[i]])
            if (i < length(RVine$names))
                cat(",   ")
        }
    }

}


createMaxMat <- function(Matrix) {

    if (dim(Matrix)[1] != dim(Matrix)[2])
        stop("Structure matrix has to be quadratic.")

    MaxMat <- reorderRVineMatrix(Matrix)

    n <- nrow(MaxMat)

    for (j in 1:(n - 1)) {
        for (i in (n - 1):j) {
            MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
        }
    }

    tMaxMat <- MaxMat
    tMaxMat[is.na(tMaxMat)] <- 0

    oldSort <- diag(Matrix)
    oldSort <- oldSort[n:1]

    for (i in 1:n) {
        MaxMat[tMaxMat == i] <- oldSort[i]
    }

    return(MaxMat)
}

neededCondDistr <- function(Vine) {

    if (dim(Vine)[1] != dim(Vine)[2])
        stop("Structure matrix has to be quadratic.")

    Vine <- reorderRVineMatrix(Vine)

    MaxMat <- createMaxMat(Vine)

    d <- nrow(Vine)

    M <- list()
    M$direct <- matrix(FALSE, d, d)
    M$indirect <- matrix(FALSE, d, d)

    M$direct[2:d, 1] <- TRUE

    for (i in 2:(d - 1)) {
        v <- d - i + 1

        bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v

        direct <- Vine[i:d, 1:(i - 1)] == v

        M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)

        M$direct[i:d, i] <- TRUE

        M$direct[i, i] <- any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
    }

    return(M)
}

as.RVineMatrix <- function(RVine) as.RVM2(RVine)


###########################################################################
# Code from Harry Joe (Thanks for that)


# varray2NO:  vine array to natural order
# irev=F means A1[d,d]=A[d,d]
# irev=T means A1[d,d]=A[d-1,d] (this option used to check if A is in
#                   equivalence class of size 1 or 2).
# A is a vine array; 1:d on diagonal is not necessary

varray2NO <- function(A, irev = FALSE, iprint = FALSE) {
    d <- nrow(A)
    d2 <- d - 2
    d1 <- d - 1
    A1 <- matrix(0, d, d)
    T <- vpartner(A)
    if (irev) {
        A1[d, d] <- A[d1, d]
    } else {
        A1[d, d] <- A[d, d]
    }
    for (k in d:2) {
        x <- A1[k, k]
        for (ell in 1:(k - 1)) A1[ell, k] <- which(T[x, ] == ell)
        T[x, ] <- 0
        T[, x] <- 0
        A1[k - 1, k - 1] <- A1[k - 1, k]
    }
    # A1 satisfies A[i,i]=A[i,i+1]
    if (iprint)
        print(A1)
    # now apply permutation
    iorder <- order(diag(A1))
    A2 <- A1
    for (i in 1:d) {
        for (j in i:d) A2[i, j] <- iorder[A1[i, j]]
    }
    if (iprint)
        print(A2)
    list(NOa = A1, NO = A2, perm = iorder, diag = diag(A1))
}

# This function is not in NAMESPACE (not for direct use).
# Function with T(x,y)=k if x-y are conditioned partners in tree k
vpartner <- function(A) {
    d <- nrow(A)
    tree <- matrix(0, d, d)
    for (j in 2:d) {
        x <- A[j, j]
        for (k in 1:(j - 1)) {
            y <- A[k, j]
            tree[x, y] <- k
            tree[y, x] <- k
        }
    }
    tree
}

# Check whether A in natural order is a valid vine array if it has
#  permutation of 1:d on diagonal and permutation of 1;j in A[1:j,j]
# Natural order also means A[j-1,j]=j-1 and A[j,j]=j
# Function with A vine array (assumed dxd) and calls to vinvstepb
# if(b==0) on input,
#    columns 4 to d, binary elements of b are randomly generated
# output is b matrix with NA in lower triangle if A is OK
#   otherwise returns -1
#varraycheck=function(A)

varray2bin <- function(A) {
    d <- nrow(A)
    b <- matrix(NA, d, d)
    b[1, ] <- 1
    diag(b) <- 1
    for (i in 3:d) b[i - 1, i] <- 1
    for (i in 4:d) {
        b0 <- vinvstepb(A, i)
        # print(b0)
        if (min(b0) == -1)
            return(-1)
        b[1:i, i] <- b0
    }
    b
}


# inverse for column i:
# input A has dimension at least ixi
# output b has length i
# This function is not in NAMESPACE (not for direct use);
# it is used by varraycheck

vinvstepb <- function(A, i, ichk0 = 0) {
    # do these basic checks first
    if (ichk0 > 0) {
        diagA <- diag(A[1:i, 1:i])
        if (max(diagA - (1:i)) != 0)
            return(-1)
        for (k in 2:i) {
            if (A[k - 1, k] != k - 1)
                return(-1)
        }
        if (A[1][3] != 1)
            return(-1)
    }

    b <- rep(1, i)
    itaken <- rep(0, i)
    itaken[i] <- 1
    itaken[i - 1] <- 1
    ac <- i - 2  # active column
    for (k in (i - 2):1) {
        if (A[k, i] == A[ac, ac]) {
            b[k] <- 1
            tem <- A[ac, ac]
            itaken[tem] <- 1
            if (k > 1) {
                ac <- max((1:i)[itaken == 0])
            }
        } else if (A[k, i] == A[k - 1, ac]) {
            b[k] <- 0
            tem <- A[k - 1, ac]
            itaken[tem] <- 1
        } else return(-1)  # not valid A in NO(i)
    }
    b
}


# various checks for validity of a vine array A with dimension d>=4
# return 1 for OK
# return -3 for diagonal not 1:d
# return -2 for not permutation of 1:j in column j
# return -1 if cannot find proper binary array from array in natural order



#' R-Vine Matrix Check
#'
#' The given matrix is tested to be a valid R-vine matrix.
#'
#'
#' @param M A dxd vine matrix: only lower triangle is used; For the check, M is
#' assumed to be in natural order, i.e. d:1 on diagonal. Further M[j+1,j]=d-j
#' and M[j,j]=d-j
#' @return \item{code}{ \code{1} for OK; \cr \code{-3} diagonal can not be put
#' in order d:1; \cr \code{-2} for not permutation of j:d in column d-j; \cr
#' \code{-1} if cannot find proper binary array from array in natural order.  }
#' @note The matrix M do not have to be given in natural order or the diagonal
#' in order d:1. The test checks if it can be done in order to be a valid
#' R-vine matrix. \cr If a function in this package needs the natural order the
#' \code{RVineMatrix} object is automatically "normalized". \cr The function
#' \code{\link{RVineMatrix}} automatically checks if the given R-vine matrix is
#' valid.
#' @author Harry Joe
#' @seealso \code{\link{RVineMatrix}}
#' @references Joe H, Cooke RM and Kurowicka D (2011). Regular vines:
#' generation algorithm and number of equivalence classes. In Dependence
#' Modeling: Vine Copula Handbook, pp 219--231. World Scientific, Singapore.
#' @examples
#'
#' A1 <- matrix(c(6, 0, 0, 0, 0, 0,
#' 			         5, 5, 0, 0, 0, 0,
#' 			         3, 4, 4, 0, 0, 0,
#' 			         4, 3, 3, 3, 0, 0,
#' 			         1, 1, 2, 2, 2, 0,
#' 			         2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
#' b1 <- RVineMatrixCheck(A1)
#' print(b1)
#' # improper vine matrix, code=-1
#' A2 <- matrix(c(6, 0, 0, 0, 0, 0,
#' 			         5, 5, 0, 0, 0, 0,
#' 			         4, 4, 4, 0, 0, 0,
#' 			         1, 3, 3, 3, 0, 0,
#' 			         3, 1, 2, 2, 2, 0,
#' 			         2, 2, 1, 1, 1,1 ), 6, 6, byrow = TRUE)
#' b2 <- RVineMatrixCheck(A2)
#' print(b2)
#' # improper vine matrix, code=-2
#' A3 <- matrix(c(6, 0, 0, 0, 0, 0,
#' 			         3, 5, 0, 0, 0, 0,
#' 			         3, 4, 4, 0, 0, 0,
#' 			         4, 3, 3, 3, 0, 0,
#' 			         1, 1, 2, 2, 2, 0,
#' 			         2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
#' b3 <- RVineMatrixCheck(A3)
#' print(b3)
#'
#' @export RVineMatrixCheck
RVineMatrixCheck <- function(M) {
    A <- M
    d <- nrow(A)
    if (d != ncol(A))
        return(-1)

    A <- A[d:1, d:1]  # unsere Notation <-> Harrys Notation

    if (sum(abs(sort(diag(A)) - (1:d))) != 0)
        return(-3)
    # convert to 1:d on diagonal
    iorder <- order(diag(A))
    A2 <- A
    for (i in 1:d) {
        for (j in i:d) A2[i, j] <- iorder[A[i, j]]
    }
    # print(A2)
    for (j in 2:d) {
        if (sum(abs(sort(A2[1:(j - 1), j]) - (1:(j - 1)))) != 0)
            return(-2)
    }
    # next convert to natural order for more checks
    if (d <= 3)
        return(1)
    ANOobj <- varray2NO(A2)
    # print(ANOobj)
    b <- varray2bin(ANOobj$NO)  # if OK, a binary matrix is returned here
    if (is.matrix(b))
        return(1) else return(-1)
}

#### -------------------------------------------------------------
## function that converts upper triagonal matrix to lower triagonal
ToLowerTri <- function(x) {
    ## only change matrix if not already lower triagonal
    if(all(x[lower.tri(x)] == 0)) {
        x[nrow(x):1, ncol(x):1]
    } else {
        x
    }
}

ToUpperTri <- function(x) {
    ## only change matrix if not already upper triagonal
    if(all(x[upper.tri(x)] == 0)) {
        x[nrow(x):1, ncol(x):1]
    } else {
        x
    }
}
