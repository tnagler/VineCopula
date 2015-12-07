#' Distribution Function of a Bivariate Copula
#'
#' This function evaluates the cumulative distribution function (CDF) of a
#' given parametric bivariate copula.
#'
#' If the family and parameter specification is stored in a \code{\link{BiCop}}
#' object \code{obj}, the alternative version \cr \preformatted{BiCopCDF(u1,
#' u2, obj)} can be used.
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
#' second parameter for bivariate copulas with two parameters (BB1, BB6, BB7,
#' BB8, Tawn type 1 and type 2; default: \code{par2 = 0}).
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#'
#' @return A numeric vector of the bivariate copula distribution function
#' \itemize{
#' \item of the copula \code{family}
#' \item with parameter(s) \code{par}, \code{par2}
#' \item evaluated at \code{u1} and \code{u2}.
#' }
#'
#' @note The calculation of the cumulative distribution function (CDF) of the
#' Student's t copula (\code{family = 2}) is not implemented any more since the
#' calculation was wrong for non-integer degrees-of-freedom.
#'
#' @author Eike Brechmann
#'
#' @seealso
#' \code{\link{BiCopPDF}},
#' \code{\link{BiCopHfunc}},
#' \code{\link{BiCopSim}},
#' \code{\link{BiCop}}
#'
#' @examples
#'
#' ## simulate from a bivariate Clayton
#' simdata <- BiCopSim(300, 3, 3.4)
#'
#' ## evaluate the distribution function of the bivariate Clayton copula
#' u1 <- simdata[,1]
#' u2 <- simdata[,2]
#' BiCopCDF(u1, u2, 3, 3.4)
#'
#' ## estimate a bivariate copula from the data and evaluate its CDF
#' cop <- BiCopSelect(u1, u2)
#' BiCopCDF(u1, u2, cop)
#'
#' @export BiCopCDF
BiCopCDF <- function(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE)
        stop("u1 and/or u2 are not set or have length zero.")
    if (any(u1 > 1) || any(u1 < 0))
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0))
        stop("Data has be in the interval [0,1].")
    if (length(u1) != length(u2))
        stop("Lengths of 'u1' and 'u2' do not match.")
    n <- length(u1)

    ## extract family and parameters if BiCop object is provided
    if (missing(family))
        family <- NA
    if (missing(par))
        par <- NA
    # for short hand usage extract obj from family
    if (class(family) == "BiCop")
        obj <- family
    if (!is.null(obj)) {
        stopifnot(class(obj) == "BiCop")
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }

    ## check for reasonable input
    if (any(is.na(family)) | any(is.na(par)))
        stop("Provide either 'family' and 'par' or 'obj'")
    if (any(family == 2))
        stop("The CDF of the t-copula is not implemented.")

    ## adjust length for parameter vectors; stop if not matching
    if (any(c(length(family), length(par), length(par2)) == n)) {
        if (length(family) == 1)
            family <- rep(family, n)
        if (length(par) == 1)
            par <- rep(par, n)
        if (length(par2) == 1)
            par2 <- rep(par2, n)
    }
    if (!(length(family) %in% c(1, n)))
        stop("'family' has to be a single number or a size n vector")
    if (!(length(par) %in% c(1, n)))
        stop("'par' has to be a single number or a size n vector")
    if (!(length(par2) %in% c(1, n)))
        stop("'par2' has to be a single number or a size n vector")

    ## sanity checks for family and parameters
    if (check.pars) {
        BiCopCheck(family, par, par2)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        family[(family %in% c(3, 13, 23, 33)) & (par == 0)] <- 0
        family[(family == 5) & (par == 0)] <- 0
    }

    ## calculate CDF
    if (length(par) == 1) {
        # call for single parameters
        out <- calcCDF(u1, u2, family, par, par2)
    } else {
        # vectorized call
        out <- vapply(1:length(par),
                      function(i) calcCDF(u1[i],
                                          u2[i],
                                          family[i],
                                          par[i],
                                          par2[i]),
                      numeric(1))
    }

    ## return result
    out
}

calcCDF <- function(u1, u2, family, par, par2) {
    if (family == 0) {
        res <- u1 * u2
    } else if (family == 1) {
        cdf <- function(u, v) pmvnorm(upper = c(qnorm(u), qnorm(v)),
                                      corr = matrix(c(1,   par, par, 1), 2, 2))
        res <- mapply(cdf, u1, u2, SIMPLIFY = TRUE)
        # }else if(family == 2){ par2=round(par2) cdf = function(u,v)
        # pmvt(upper=c(qt(u,df=par2),qt(v,df=par2)), corr=matrix(c(1,par,par,1),2,2),
        # df=par2) res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)
    } else if (family %in% c(3:10, 41)) {
        res <- .C("archCDF",
                  as.double(u1),
                  as.double(u2),
                  as.integer(length(u1)),
                  as.double(c(par, par2)),
                  as.integer(family),
                  as.double(rep(0, length(u1))),
                  PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(13, 14, 16:20, 51)) {
        res <- u1 + u2 - 1 + .C("archCDF",
                                as.double(1 - u1),
                                as.double(1 - u2),
                                as.integer(length(u1)),
                                as.double(c(par, par2)),
                                as.integer(family - 10),
                                as.double(rep(0, length(u1))),
                                PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(23, 24, 26:30, 61)) {
        res <- u2 - .C("archCDF",
                       as.double(1 - u1),
                       as.double(u2),
                       as.integer(length(u1)),
                       as.double(c(-par, -par2)),
                       as.integer(family - 20),
                       as.double(rep(0, length(u1))),
                       PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(33, 34, 36:40, 71)) {
        res <- u1 - .C("archCDF",
                       as.double(u1),
                       as.double(1 - u2),
                       as.integer(length(u1)),
                       as.double(c(-par, -par2)),
                       as.integer(family - 30),
                       as.double(rep(0, length(u1))),
                       PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(104, 114, 124, 134, 204, 214, 224, 234)) {

        if (family == 104) {
            par3 <- 1
            res <- .C("TawnC",
                      as.double(u1),
                      as.double(u2),
                      as.integer(length(u1)),
                      as.double(par),
                      as.double(par2),
                      as.double(par3),
                      as.double(rep(0, length(u1))),
                      PACKAGE = "VineCopula")[[7]]
        }
        if (family == 114) {
            par3 <- 1
            res <- u1 + u2 - 1 + .C("TawnC",
                                    as.double(1-u1),
                                    as.double(1-u2),
                                    as.integer(length(u1)),
                                    as.double(par),
                                    as.double(par2),
                                    as.double(par3),
                                    as.double(rep(0, length(u1))),
                                    PACKAGE = "VineCopula")[[7]]
        }
        if (family == 124) {
            par3 <- par2
            par2 <- 1
            res <- u2 - .C("TawnC",
                           as.double(1-u1),
                           as.double(u2),
                           as.integer(length(u1)),
                           as.double(-par),
                           as.double(par2),
                           as.double(par3),
                           as.double(rep(0, length(u1))),
                           PACKAGE = "VineCopula")[[7]]
        }
        if (family == 134) {
            par3 <- par2
            par2 <- 1
            res <- u1 - .C("TawnC",
                           as.double(u1),
                           as.double(1-u2),
                           as.integer(length(u1)),
                           as.double(-par),
                           as.double(par2),
                           as.double(par3),
                           as.double(rep(0, length(u1))),
                           PACKAGE = "VineCopula")[[7]]
        }
        if (family == 204) {
            par3 <- par2
            par2 <- 1
            res <- .C("TawnC",
                      as.double(u1),
                      as.double(u2),
                      as.integer(length(u1)),
                      as.double(par),
                      as.double(par2),
                      as.double(par3),
                      as.double(rep(0, length(u1))),
                      PACKAGE = "VineCopula")[[7]]
        }
        if (family == 214) {
            par3 <- par2
            par2 <- 1
            res <- u1 + u2 - 1 + .C("TawnC",
                                    as.double(1-u1),
                                    as.double(1-u2),
                                    as.integer(length(u1)),
                                    as.double(par),
                                    as.double(par2),
                                    as.double(par3),
                                    as.double(rep(0, length(u1))),
                                    PACKAGE = "VineCopula")[[7]]
        }
        if (family == 224) {
            par3 <- 1
            res <- u2 - .C("TawnC",
                           as.double(1-u1),
                           as.double(u2),
                           as.integer(length(u1)),
                           as.double(-par),
                           as.double(par2),
                           as.double(par3),
                           as.double(rep(0, length(u1))),
                           PACKAGE = "VineCopula")[[7]]
        }
        if (family == 234) {
            par3 <- 1
            res <- u1 - .C("TawnC",
                           as.double(u1),
                           as.double(1-u2),
                           as.integer(length(u1)),
                           as.double(-par),
                           as.double(par2),
                           as.double(par3),
                           as.double(rep(0, length(u1))),
                           PACKAGE = "VineCopula")[[7]]
        }
    } else {
        res <- rep(NA, length(u1))
    }

    ## return results
    res
}
