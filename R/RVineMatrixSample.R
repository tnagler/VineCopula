#' R-Vine Matrix Sample
#'
#' Sample R-Vine matrices.
#'
#' @param d Dimension of the R-Vine matrices.
#' @param size Number of matrices to sample.
#' @param naturalOrder Should the matrices be in the natural order
#' (default: \code{naturalOrder = FALSE}).
#' @return A list of length \code{size} with each element containing one
#' R-Vine matrix.
#' @author Thibault Vatter
#' @seealso \code{\link{RVineMatrix}}, \code{\link{RVineMatrixCheck}}
#' @references Joe H, Cooke RM and Kurowicka D (2011). Regular vines:
#' generation algorithm and number of equivalence classes. In Dependence
#' Modeling: Vine Copula Handbook, pp 219--231. World Scientific, Singapore.
#' @examples
#'
#' # Matrix and sample sizes
#' d <- 10
#' size <- 5
#'
#' # Sample in the natural order
#' RVM <- RVineMatrixSample(d, size)
#' sapply(RVM, RVineMatrixCheck)
#'
#' # Sample without natural order
#' RVM <- RVineMatrixSample(d, size, naturalOrder = TRUE)
#' sapply(RVM, RVineMatrixCheck)
#'
#' @export RVineMatrixSample
RVineMatrixSample <- function(d, size = 1, naturalOrder = FALSE) {
    stopifnot(d > 1)

    ## Sample the required binary vectors
    if (d > 3) {
        sampleBvect <- lapply(4:d,
                              function(j) sampleBinaryVector(j, size, TRUE))
        sampleRVM <- vector("list", size)
    }

    ## Initialize RVM
    initRVM <- diag(1:d)
    delta <- col(initRVM) - row(initRVM)
    initRVM[delta == 1] <- 1:(d - 1)
    if (d > 2) {
        initRVM[1, 3] <- 1
    }

    ## Part of the RVM that needs to be sampled (for d > 3)
    selUpper <- delta > 1

    ## The output
    sampleRVM <- vector("list", size)
    for (k in 1:size) {

        ## Use the initRVM
        RVM <- initRVM

        if (d > 3) {
            ## Get the required binary sample
            b <- lapply(4:d, function(j) sampleBvect[[j - 3]][k, ])

            ## Call the C code
            RVM[selUpper] <- getRVineMatrix(b)
        }

        ## Permute nodes if required
        if (naturalOrder == FALSE) {
            RVM <- reorderRVineMatrix(RVM, sample.int(d, d))
        }

        sampleRVM[[k]] <- ToLowerTri(RVM)
    }

    sampleRVM
}

getRVineMatrix <- function(b) {

    if (is.list(b)) {
        d <- length(b) + 3

        ## Add the first tree (trivial elements)
        b <- append(lapply(1:3, function(j) rep(1, j)), b)

        ## Transform to binary matrix
        b <- sapply(b, function(x) c(rep(0, d - length(x)), x))[d:1, ]
    } else {
        d <- nrow(b)
    }

    RVM <- rep(0, d * (d-1) / 2 - (d - 1))
    RVM <- .C("getRVM",
              as.integer(b),
              as.integer(d),
              as.integer(RVM),
              PACKAGE = "VineCopula")[[3]]

    RVM
}

sampleBinaryVector <- function(d, size = d, constraint = FALSE) {
    if (constraint == FALSE) {
        return(matrix(rbinom(d * size, 1, 0.5), size, d))
    } else {
        if (d < 4) {
            return(matrix(1, size, d))
        } else {
            z <- sampleBinaryVector(d-3, size, constraint = FALSE)
            return(cbind(1, z, 1, 1))
        }
    }
}
