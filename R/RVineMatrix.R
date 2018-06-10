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
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#'
#' @return An object of class \code{\link{RVineMatrix}}, i.e., a list with the
#' following components:
#' \item{Matrix}{R-vine tree structure matrix.}
#' \item{family}{pair-copula family matrix with values as above.}
#' \item{par}{pair-copula parameter matrix.}
#' \item{par2}{second pair-copula parameter matrix with parameters necessary for
#'  pair-copula families with two parameters.}
#' \item{names}{variable names (defaults to \code{V1, V2, ...}).}
#' \item{MaxMat, CondDistr}{additional matrices required internally for
#' evaluating the density etc.,}
#' \item{type}{the type of the vine copula structure; possible types are:
#' \itemize{
#' \item{\code{"C-vine": }}{all trees consist of a star,}
#' \item{\code{"D-vine": }}{all trees consist of a path,}
#' \item{\code{"R-vine": }}{all structures that are neither a C- nor D-vine,}
#' }}
#' \item{tau}{Kendall's tau matrix,}
#' \item{taildep}{matrices of lower and upper tail dependence coefficients,}
#' \item{beta}{Blomqvist's beta matrix.}
#' Objects of this class are also returned by the \code{\link{RVineSeqEst}},
#' \code{\link{RVineCopSelect}}, and \code{\link{RVineStructureSelect}}
#' functions. In this case, further information about the fit is added.
#'
#'
#' @note For a comprehensive summary of the vine copula model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}.\cr
#' The \code{\link{RVineMatrix}} function automatically checks if the given
#' matrix is a valid R-vine matrix (see \code{\link{RVineMatrixCheck}}). \cr
#' Although the function allows upper triangular matrices as its input, it will
#' always store them as lower triangular matrices.
#'
#' @author Jeffrey Dissmann, Thomas Nagler
#'
#' @seealso
#' \code{\link{RVineMatrixCheck}},
#' \code{\link{RVineSeqEst}},
#' \code{\link{RVineCopSelect}},
#' \code{\link{RVineStructureSelect}},
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
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' ## define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' ## see the object's content or a summary
#' str(RVM)
#' summary(RVM)
#'
#' ## inspect the model using plots
#' \dontrun{plot(RVM)  # tree structure}
#' contour(RVM)  # contour plots of all pair-copulas
#'
#' ## simulate from the vine copula model
#' plot(RVineSim(500, RVM))
#'
RVineMatrix <- function(Matrix,
                        family = array(0, dim = dim(Matrix)),
                        par = array(NA, dim = dim(Matrix)),
                        par2 = array(NA, dim = dim(Matrix)),
                        names = NULL, check.pars = TRUE) {
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_matrix,
                    check_fammat,
                    check_parmat,
                    check_par2mat)
    list2env(args, environment())

    ## check matrices
    if (length(names) > 0 & length(names) != dim(Matrix)[1])
        stop("Length of the vector 'names' is not correct.")

    ## check for family/parameter consistency
    sel <- lower.tri(family)  # selector for lower triangular matrix
    if (check.pars & (any(family != 0) | any(!is.na(par)))) {
        BiCopCheck(family[sel], par[sel], par2[sel], call = match.call())
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

    ## add vine type
    if (is.CVine(RVM)) {
        RVM$type <- "C-vine"
    } else if (is.DVine(RVM)) {
        RVM$type <- "D-vine"
    } else {
        RVM$type <- "R-vine"
    }

    ## add dependence measures

    # create list of BiCop ojbects
    objlst <- apply(cbind(family[sel], par[sel], par2[sel]), 1, function(x)
        if (x[1] == 0) NA else BiCop(x[1], x[2], x[3], check.pars = FALSE))

    # construct dependence measure matrices
    taus <- utds <- ltds <- bets <- matrix(0, nrow(Matrix), ncol(Matrix))
    taus[sel] <- vapply(objlst, function(x)
        ifelse(inherits(x, "BiCop"), x$tau, 0), numeric(1))
    utds[sel] <- vapply(objlst, function(x)
        ifelse(inherits(x, "BiCop"), x$taildep$upper, 0), numeric(1))
    ltds[sel] <- vapply(objlst, function(x)
        ifelse(inherits(x, "BiCop"), x$taildep$lower, 0), numeric(1))
    bets[sel] <- vapply(objlst, function(x)
        ifelse(inherits(x, "BiCop"), x$beta, 0), numeric(1))

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

reorderRVineMatrix <- function(Matrix, oldOrder = NULL) {

    if (length(oldOrder) == 0) {
        oldOrder <- diag(Matrix)
    }
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
    cat(x$type, "copula with the following pair-copulas:\n")
    d <- dim(RVine)
    for (j in 1:(d - 1)) {
        cat("")
        a <- paste("Tree ", j, ":\n", sep = "")
        cat(a)

        pc.nums <- sapply(1:(d - j), get_num, tree = j, RVM = RVine)
        pc.nums <- sapply(pc.nums, function(x) gsub(" ", "", x))
        pc.nums.len <- nchar(pc.nums)
        pc.nums.space <- max(pc.nums.len) - pc.nums.len
        maxa <- 0
        for (i in 1:(d - j)) {
            a <- draw_blanks(pc.nums.space[i])
            a <- paste0(a, pc.nums[i])
            a <- paste(a,
                       "  ",
                       BiCopName(RVine$family[d - j + 1, i], short = FALSE),
                       sep = "")
            if (RVine$family[d - j + 1, i] != 0) {
                a <- paste(a,
                           " (par = ",
                           round(RVine$par[d - j + 1, i], 2),
                           sep = "")
                if (RVine$family[d - j + 1, i] %in% allfams[twopar]) {
                    a <- paste(a,
                               ", par2 = ",
                               round(RVine$par2[d - j + 1, i], 2),
                               sep = "")
                }
                a <- paste(a,
                           ", tau = ",
                           round(RVine$tau[d - j + 1, i], 2),
                           ")",
                           sep = "")
            }
            a <- paste(a, "\n")
            maxa <- max(maxa, nchar(a))
            cat(a)
        }
        if (j < d - 1) cat("\n")
    }
    # show names if provided
    if (!is.null(RVine$names)) {
        linelen <- maxa
        cat("\n")
        cat("---\n")
        txt <- paste0(1, " <-> ", RVine$names[[1]])
        for (i in 2:(d - 1)) {
            if (nchar(txt) > linelen) {
                cat(txt, ",\n", sep = "")
                txt <- paste0(i, " <-> ", RVine$names[[i]])
            } else {
                txt <- paste0(txt, ",   ", i, " <-> ", RVine$names[[i]])
            }
        }
        if (nchar(txt) > linelen) {
            cat(txt, ",\n", sep = "")
            txt <- paste0(d, " <-> ", RVine$names[[d]])
        } else {
            txt <- paste0(txt, ",   ", d, " <-> ", RVine$names[[d]])
        }
        cat(txt)
    }
}

summary.RVineMatrix <- function(object, with.se = TRUE, ...) {

    ## create character matrices with pair-copula info
    #     cat("Pair-copulas:\n")
    d <- nrow(object$Matrix)
    fammat  <- matrix("", d, d)
    parmat  <- formatC(object$par, 2, format = "f")
    par2mat <- formatC(object$par2, 2, format = "f")
    taumat  <- formatC(object$tau, 2, format = "f")
    utdmat  <- formatC(object$taildep$upper, 2, format = "f")
    ltdmat  <- formatC(object$taildep$lower, 2, format = "f")
    nammat  <- matrix("", d, d)
    nummat  <- matrix(0, d, d)
    with.se <- with.se & !is.null(object$se)
    if (with.se) {
        semat  <- formatC(object$se, 2, format = "f")
        se2mat <- formatC(object$se2, 2, format = "f")
    }

    ## get names and clean matrices
    for (i in 2:d) {
        for (j in 1:(i - 1)) {
            fammat[i, j] <- BiCopName(object$family[i, j])
            nummat[i, j] <- object$family[i, j]
            nammat[i, j] <- gsub(" ", "", get_num(j, d - i + 1, object))
            if (fammat[i, j] == "I") {
                parmat[i, j] <- "-"
                par2mat[i, j] <- "-"
            } else {
                if (with.se) {
                    parmat[i, j] <- paste0(parmat[i, j],
                                           " (",
                                           semat[i, j],
                                           ")")
                    if (object$family[i, j] %in% allfams[twopar]) {
                        par2mat[i, j] <- paste0(par2mat[i, j],
                                                " (",
                                                se2mat[i, j],
                                                ")")
                    } else {
                        par2mat[i, j] <- "-"
                    }
                }
            }
            if (object$taildep$upper[i, j] == 0)
                utdmat[i, j] <- "-"
            if (object$taildep$lower[i, j] == 0)
                ltdmat[i, j] <- "-"
        }
    }

    ## maximal number of characters for each category
    ltree <- nchar("tree")
    lfam  <- nchar("family")
    lfname <- max(nchar("cop"), max(sapply(fammat, nchar)))
    lpar  <- max(nchar("par"), max(sapply(parmat, nchar)))
    lpar2 <- max(nchar("par2"), max(sapply(par2mat, nchar)))
    ltau  <- max(nchar("tau"), max(sapply(taumat, nchar)))
    lutd  <- max(nchar("UTD"), max(sapply(utdmat, nchar)))
    lltd  <- max(nchar("LTD"), max(sapply(ltdmat, nchar)))
    lnam  <- max(nchar("edge"), max(sapply(nammat, nchar)))


    ## line with headings
    txt <- "tree "
    # substract nchar(edge) - 1 (for space) = 3
    txt <- paste0(txt, draw_blanks(max(1, lnam - 3)), "edge ")
    txt <- paste0(txt, "| family ")
    txt <- paste0(txt, draw_blanks(max(1, lfname - 2)), "cop ")
    txt <- paste0(txt, draw_blanks(max(1, lpar - 2)), "par ")
    txt <- paste0(txt, draw_blanks(max(1, lpar2 - 3)), "par2 |")
    txt <- paste0(txt, draw_blanks(max(1, ltau - 2)), "tau ")
    txt <- paste0(txt, draw_blanks(max(1, lutd - 2)), "utd ")
    txt <- paste0(txt, draw_blanks(max(1, lltd - 2)), "ltd")
    cat(txt, "\n")
    linelen <- nchar(txt)
    cat(draw_lines(linelen), "\n")

    for (tree in 1:(d - 1)) {
        for (edge in 1:(d - tree)) {
            ## print tree number
            if (edge == 1) {
                cat(draw_blanks(max(0, ltree - nchar(tree))))
                cat(tree, "")
            } else {
                cat("     ")
            }

            ## print edge label
            tmpch <- nammat[d + 1 - tree, edge]
            cat(draw_blanks(max(0, lnam - nchar(tmpch))), tmpch)

            ## print copula family
            cat(" |")
            cat(formatC(nummat[d + 1 - tree, edge], lfam))
            tmpch <- fammat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, lfname - nchar(tmpch))) + 1), tmpch)

            ## print parameters
            tmpch <- parmat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, lpar - nchar(tmpch)) + 1)), tmpch)
            tmpch <- par2mat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, lpar2 - nchar(tmpch)) + 1)), tmpch)

            ## print dependence measures
            cat(" |")
            tmpch <- taumat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, ltau - nchar(tmpch)))), tmpch)
            tmpch <- utdmat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, lutd - nchar(tmpch)) + 1)), tmpch)
            tmpch <- ltdmat[d + 1 - tree, edge]
            cat(draw_blanks(min(max(0, lltd - nchar(tmpch)) + 1)), tmpch)


            cat("\n")

        }
    }

    ## print general info
    cat("---\n")
    cat("type:", object$type, "   ")
    if (!is.null(object$logLik)) {
        cat("logLik:", round(object$logLik, 2), "   ")
        cat("AIC:", round(object$AIC, 2), "   ")
        cat("BIC:", round(object$BIC, 2), "   ")
    }
    # show names if provided
    if (!is.null(object$names)) {
        linelen <- min(linelen, 90)
        cat("\n")
        cat("---\n")
        txt <- paste0(1, " <-> ", object$names[[1]])
        for (i in 2:(d - 1)) {
            if (nchar(txt) > linelen) {
                cat(txt, ",\n", sep = "")
                txt <- paste0(i, " <-> ", object$names[[i]])
            } else {
                txt <- paste0(txt, ",   ", i, " <-> ", object$names[[i]])
            }
        }
        if (nchar(txt) > linelen) {
            cat(txt, ",\n", sep = "")
            txt <- paste0(d, " <-> ", object$names[[d]])
        } else {
            txt <- paste0(txt, ",   ", d, " <-> ", object$names[[d]])
        }
        cat(txt)
    }

    sel <- upper.tri(nammat)
    tab <- data.frame(
        tree  = do.call(c, lapply(1:(d - 1), function(i) rep(i, d - i))),
        edge  = rev(t(nammat)[sel]),
        family = rev(t(object$family)[sel]),
        cop = rev(t(fammat)[sel]),
        par  = rev(t(object$par)[sel]),
        par2 = rev(t(object$par2)[sel]),
        tau = rev(t(object$tau)[sel]),
        utd = rev(t(object$taildep$upper)[sel]),
        ltd = rev(t(object$taildep$lower)[sel])
    )
    invisible(tab)
}

draw_blanks <- function(len) {
    do.call(paste0, as.list(rep(" ", len)))
}

draw_lines <- function(len) {
    do.call(paste0, as.list(rep("-", len)))
}

## A D-vine has a path in the first tree (and thus in all trees)
is.DVine <- function(Matrix) {
    if (inherits(Matrix, "RVineMatrix"))
        Matrix <- Matrix$Matrix
    Matrix <- reorderRVineMatrix(Matrix)
    d <- nrow(Matrix)
    nums <- c(diag(Matrix)[seq_len(d - 1)], Matrix[d, seq_len(d - 1)])
    repcount <- table(table(nums))
    # no node in more than two edges, only 2 nodes in one edge
    (names(repcount) == 1:2) && (repcount[1] == 2)
}

## A C-vine has a star in each tree
is.CVine <- function(Matrix) {
    if (inherits(Matrix, "RVineMatrix"))
        Matrix <- Matrix$Matrix
    Matrix <- reorderRVineMatrix(Matrix)
    d <- nrow(Matrix)

    # a vine in less then 4 dimensions is always a C-vine
    if (d < 4)
        return(TRUE)

    # check conditioning sets of each tree (same number has to enter at all
    # edges)
    all.trees.star <- TRUE
    for (tree in 2:(d - 2)) {
        all.trees.star <- all.trees.star &
            (length(unique(Matrix[d - tree + 2, 1:(d - tree)])) == 1)
    }

    all.trees.star
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
    TT <- vpartner(A)
    if (irev) {
        A1[d, d] <- A[d1, d]
    } else {
        A1[d, d] <- A[d, d]
    }
    for (k in d:2) {
        x <- A1[k, k]
        for (ell in 1:(k - 1)) A1[ell, k] <- which(TT[x, ] == ell)
        TT[x, ] <- 0
        TT[, x] <- 0
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


#' R-Vine Matrix Check
#'
#' The given matrix is tested to be a valid R-vine matrix.
#'
#'
#' @param M A \eqn{dxd} vine matrix.
#' @return \item{code}{ \code{1} for OK; \cr
#' \code{-4} matrix is neither lower nor upper triangular;\cr
#' \code{-3} diagonal can not be put in order d:1;\cr
#' \code{-2} for not permutation of j:d in column d-j; \cr
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
#'                5, 5, 0, 0, 0, 0,
#'                3, 4, 4, 0, 0, 0,
#'                4, 3, 3, 3, 0, 0,
#'                1, 1, 2, 2, 2, 0,
#'                2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
#' b1 <- RVineMatrixCheck(A1)
#' print(b1)
#' # improper vine matrix, code=-1
#' A2 <- matrix(c(6, 0, 0, 0, 0, 0,
#'                5, 5, 0, 0, 0, 0,
#'                4, 4, 4, 0, 0, 0,
#'                1, 3, 3, 3, 0, 0,
#'                3, 1, 2, 2, 2, 0,
#'                2, 2, 1, 1, 1,1 ), 6, 6, byrow = TRUE)
#' b2 <- RVineMatrixCheck(A2)
#' print(b2)
#' # improper vine matrix, code=-2
#' A3 <- matrix(c(6, 0, 0, 0, 0, 0,
#'                3, 5, 0, 0, 0, 0,
#'                3, 4, 4, 0, 0, 0,
#'                4, 3, 3, 3, 0, 0,
#'                1, 1, 2, 2, 2, 0,
#'                2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
#' b3 <- RVineMatrixCheck(A3)
#' print(b3)
#'
RVineMatrixCheck <- function(M) {
    lmat <- M[lower.tri(M)]
    umat <- M[upper.tri(M)]
    if (!(all(lmat == 0) | all(umat == 0)))
        return(-4)
    A <- ToLowerTri(M)
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
    ANOobj <- tryCatch(varray2NO(A2), error = function(e) NULL)
    if (is.null(ANOobj))
        return(-1)
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
