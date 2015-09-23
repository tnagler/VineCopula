RVineMatrix <- function(Matrix, family = array(0, dim = dim(Matrix)), par = array(NA, dim = dim(Matrix)), par2 = array(NA, dim = dim(Matrix)), names = NULL, check.pars = TRUE) {
    
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
    if (any(!(family %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40, 43, 44, 104, 114, 124, 134, 204, 214, 224, 234)))) 
        stop("Copula family not implemented.")
    if (length(names) > 0 & length(names) != dim(Matrix)[1]) 
        stop("Length of the vector 'names' is not correct.")
    if (RVineMatrixCheck(Matrix) != 1) 
        stop("'Matrix' is not a valid R-vine matrix")
    
    ## check for family/parameter consistency
    if (check.pars) {
        if (!all(par %in% c(0, NA))) {
            for (i in 2:dim(Matrix)[1]) {
                for (j in 1:(i - 1)) {
                    BiCopCheck(family[i, j], par[i, j], par2[i, j])
                }
            }
        }
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
                       names = rev(RVM$names[oldOrder])))
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

print.RVineMatrix <- function(x, detail = FALSE, ...) {
    RVine <- x
    message("R-vine matrix:")
    print(RVine$Matrix, ...)
    
    # Falls namen diese auch ausgeben
    if (!is.null(RVine$names)) {
        message("")
        message("Where")
        for (i in 1:length(RVine$names)) {
            message(i, " <-> ", RVine$names[[i]])
        }
    }
    # NextMethod('print')
    
    d <- dim(RVine)
    if (detail == TRUE || detail == T) {
        message("")
        message("Tree 1:")
        for (i in 1:(d - 1)) {
            a <- paste(RVine$names[[RVine$Matrix[i, i]]], 
                       ",",
                       RVine$names[[RVine$Matrix[d, i]]],
                       sep = "")
            a <- paste(a,
                       ": ",
                       BiCopName(RVine$family[d, i], short = FALSE),
                       sep = "")
            if (RVine$family[d, i] != 0) {
                a <- paste(a, " with par=", round(RVine$par[d, i], 2), sep = "")
                if (RVine$family[d, i] %in% c(2, 7, 8, 9, 10,
                                              17, 18, 19, 20,
                                              27, 28, 29, 30, 
                                              37, 38, 39, 40,
                                              104, 114, 124, 134, 
                                              204, 214, 224, 234)) {
                    a <- paste(a, " and par2=", round(RVine$par2[d, i], 2), sep = "")
                }
                a <- paste(a,
                           " (tau=", 
                           round(BiCopPar2Tau(RVine$family[d, i],
                                              RVine$par[d, i], 
                                              RVine$par2[d, i]), 2),
                           ")", 
                           sep = "")
            }
            message(a)
        }
        for (j in 2:(d - 1)) {
            message("")
            a <- paste("Tree ", j, ":", sep = "")
            message(a)
            for (i in 1:(d - j)) {
                a <- paste(RVine$names[[RVine$Matrix[i, i]]], 
                           ",", 
                           RVine$names[[RVine$Matrix[d - j + 1, i]]],
                           sep = "")
                a <- paste(a, "|", sep = "")
                conditioningSet <- (d - j + 2):d
                for (k in 1:length(conditioningSet)) {
                    if (k > 1) {
                        a <- paste(a, ",", sep = "")
                    }
                    a <- paste(a,
                               RVine$names[[RVine$Matrix[conditioningSet[k], i]]], 
                               sep = "")
                }
                a <- paste(a,
                           ": ", 
                           BiCopName(RVine$family[d - j + 1, i], short = FALSE), 
                           sep = "")
                if (RVine$family[d - j + 1, i] != 0) {
                    a <- paste(a, 
                               " with par=",
                               round(RVine$par[d - j + 1, i], 2),
                               sep = "")
                    if (RVine$family[d - j + 1, i] %in% c(2, 7, 8, 9, 10,
                                                          17, 18, 19, 20,
                                                          27, 28, 29, 30, 
                                                          37, 38, 39, 40,
                                                          104, 114, 124, 134, 
                                                          204, 214, 224, 234)) {
                        a <- paste(a, 
                                   " and par2=",
                                   round(RVine$par2[d - j + 1, i], 2), 
                                   sep = "")
                    }
                    a <- paste(a,
                               " (tau=",
                               round(BiCopPar2Tau(RVine$family[d - j + 1, i], 
                                                  RVine$par[d - j + 1, i], 
                                                  RVine$par2[d - j + 1, i]), 2),
                               ")",
                               sep = "")
                }
                message(a)
            }
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
