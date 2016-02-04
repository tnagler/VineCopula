#' Sequential Pair-Copula Selection and Estimation for R-Vine Copula Models
#'
#' This function fits a R-vine copula model to a d-dimensional copula data set.
#' Pair-copula families are selected using \code{\link{BiCopSelect}} and
#' estimated sequentially.
#'
#' R-vine copula models with unknown structure can be specified using
#' \code{\link{RVineStructureSelect}}.
#'
#' @param data N x d data matrix (with uniform margins).
#' @param familyset integer vector of pair-copula families to select from
#' (the independence copula MUST NOT be specified in this vector unless one
#' wants to fit an independence vine!).  The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases.  If \code{familyset = NA} (default), selection among all
#' possible families is performed.  The coding of pair-copula families is shown
#' below.
#' @param Matrix lower or upper triangular d x d matrix that defines the R-vine
#' tree structure.
#' @param selectioncrit Character indicating the criterion for pair-copula
#' selection. Possible choices: \code{selectioncrit = "AIC"} (default) or
#' \code{"BIC"} (see \code{\link{BiCopSelect}}).
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).  The
#' independence copula is chosen for a (conditional) pair if the null
#' hypothesis of independence cannot be rejected.
#' @param level numeric; significance level of the independence test (default:
#' \code{level = 0.05}).
#' @param trunclevel integer; level of truncation.
#' @param rotations logical; if \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#' @param cores integer; if \code{cores > 1}, estimation will be parallized
#' within each tree (using \code{\link[foreach]{foreach}}).
#'
#' @return An \code{\link{RVineMatrix}} object with the following matrix
#' components \item{Matrix}{R-vine tree structure matrix as given by the
#' argument \code{Matrix}.}
#' \item{family}{Selected pair-copula family matrix
#' with values corresponding to\cr
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
#' \code{234} = rotated Tawn type 2 copula (270 degrees)
#' }
#'
#' @author Eike Brechmann, Thomas Nagler
#'
#' @seealso
#' \code{\link{RVineStructureSelect}},
#' \code{\link{BiCopSelect}},
#' \code{\link{RVineSeqEst}},
#' \code{\link[foreach]{foreach}}
#'
#' @references Brechmann, E. C., C. Czado, and K. Aas (2012). Truncated regular
#' vines in high dimensions with applications to financial data. Canadian
#' Journal of Statistics 40 (1), 68-85.
#'
#' Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
#' Selecting and estimating regular vine copulae and application to financial
#' returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
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
#' # simulate a sample of size 1000 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(1000, RVM)
#'
#' # determine the pair-copula families and parameters
#' RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4, 5 ,6), Matrix)
#' RVM1$family
#' round(RVM1$par, 2)
#'
RVineCopSelect <- function(data, familyset = NA, Matrix, selectioncrit = "AIC", indeptest = FALSE,
                           level = 0.05, trunclevel = NA, rotations = TRUE, cores = 1) {
    d <- n <- ncol(data)
    N <- nrow(data)

    ## sanity checks
    if (nrow(Matrix) != ncol(Matrix))
        stop("Structure matrix has to be quadratic.")
    if (max(Matrix) > nrow(Matrix))
        stop("Error in the structure matrix.")
    if (N < 2)
        stop("Number of observations has to be at least 2.")
    if (d < 2)
        stop("Dimension has to be at least 2.")
    if (any(data > 1) || any(data < 0))
        stop("Data has be in the interval [0,1].")
    if (!is.na(familyset[1])) {
        if (any(!(familyset %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40,
                                   104, 114, 124, 134, 204, 214, 224, 234))))
            stop("Copula family not implemented.")
    }
    if (selectioncrit != "AIC" && selectioncrit != "BIC")
        stop("Selection criterion not implemented.")
    if (level < 0 & level > 1)
        stop("Significance level has to be between 0 and 1.")

    ## set variable names and trunclevel if not provided
    if (is.null(colnames(data)))
        colnames(data) <- paste("V", 1:d, sep = "")
    varnames <- colnames(data)
    if (is.na(trunclevel))
        trunclevel <- d

    ## adjust familyset
    types <- familyset
    if (trunclevel == 0)
        types <- 0

    ## reorder matrix to natural order
    Matrix <- ToLowerTri(Matrix)
    M <- Matrix
    Mold <- M
    o <- diag(M)
    M <- reorderRVineMatrix(M)
    data <- data[, o[length(o):1]]

    ## create matrices required for selection of h-functions
    MaxMat <- createMaxMat(M)
    CondDistr <- neededCondDistr(M)

    ## create objects for results
    Types   <- matrix(0, d, d)
    Params  <- matrix(0, d, d)
    Params2 <- matrix(0, d, d)
    V <- list()
    V$direct <- array(NA, dim = c(d, N))
    V$indirect <- array(NA, dim = c(d, N))
    V$direct <- t(data[, d:1])

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    ## loop over all trees and pair-copulas
    for (k in d:2) {
        doEst <- function(i) {
            if (k > i) {
                ## get pseudo-observaions
                m <- MaxMat[k, i]
                zr1 <- V$direct[i, ]

                zr2 <- if (m == M[k, i]) {
                    V$direct[(d - m + 1), ]
                } else {
                    V$indirect[(d - m + 1), ]
                }

                ## select pair-copula
                if (trunclevel <= (d-k)) {
                    cfit <- BiCop(0, 0, check.pars = FALSE)
                } else {
                    cfit <- BiCopSelect(zr2,
                                        zr1,
                                        familyset,
                                        selectioncrit,
                                        indeptest,
                                        level,
                                        weights = NA,
                                        rotations)
                }

                ## transform data to pseudo-oberstavions in next tree
                direct <- indirect <- NULL
                if (CondDistr$direct[k - 1, i])
                    direct <- BiCopHfunc(zr1,
                                         zr2,
                                         cfit,
                                         check.pars = FALSE)$hfunc2
                if (CondDistr$indirect[k - 1, i])
                    indirect <- BiCopHfunc(zr1,
                                           zr2,
                                           cfit,
                                           check.pars = FALSE)$hfunc1

                ## return results
                list(direct = direct, indirect = indirect, cfit = cfit)
            } else {
                list(cfit = BiCop(0, 0))
            }
        }

        ## run pair-copula selection for tree k
        res.k <- if (cores > 1) {
            foreach(i = 1:(k-1),
                    .packages = c("VineCopula"),
                    .export = ) %dopar% doEst(i)
        } else {
            lapply(1:(k-1), doEst)
        }

        for (i in 1:(k-1)) {
            ## store results in matrices
            Types[k, i] <- res.k[[i]]$cfit$family
            Params[k, i] <- res.k[[i]]$cfit$par
            Params2[k, i] <- res.k[[i]]$cfit$par2
            if (!is.null(res.k[[i]]$direct))
                V$direct[i, ] <- res.k[[i]]$direct
            if (!is.null(res.k[[i]]$indirect))
                V$indirect[i, ] <- res.k[[i]]$indirect
        } # end i = 1:(d-1)
    } # end k = d:2

    ## free memory and return results
    .RVM <- RVineMatrix(Mold,
                        family = Types,
                        par = Params,
                        par2 = Params2,
                        names = varnames)
    rm(list = ls())
    .RVM
}
