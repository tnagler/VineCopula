#' Check for family/parameter consistency in bivariate copula models
#'
#' The function checks if a certain combination of copula family and parameters
#' can be used within other functions of this package.
#'
#' @param family An integer defining the bivariate copula family: \cr
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
#' @param par Copula parameter.
#' @param par2 Second parameter for bivariate copulas with two parameters (t,
#' BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default is \code{par2 = 0}).
#' @param \dots used internally.
#'
#' @return A logical indicating whether the family can be used with the parameter
#' specification.
#'
#' @author Thomas Nagler
#'
#' @examples
#' ## check parameter of Clayton copula
#' BiCopCheck(3, 1)  # works
#'
#' \dontrun{BiCopCheck(3, -1)  # does not work (only positive parameter is allowed)}
#'
BiCopCheck <- function(family, par, par2 = 0, ...) {
    # see if call has been passed via ...
    cl <- ifelse(!is.null(list(...)$call), list(...)$call[1], "BiCopCheck")

    ## check if all required parameters are set
    if (!(all(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19,
                            20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39,
                            40, 41, 42, 51, 52,  61, 62, 71, 72,
                            104, 114, 124, 134, 204, 214, 224, 234))))
        stop("\n In ", cl, ": ",
             "Copula family not implemented.",
             call. = FALSE)
    if (any((family %in% allfams[twopar]) & (par2 == 0)))
        stop("\n In ", cl, ": ",
             "For t-, BB1, BB6, BB7, BB8 and Tawn copulas, 'par2' must be set.",
             call. = FALSE)
    if (length(par) < 1)
        stop("\n In ", cl, ": ",
             "'par' not set.",
             call. = FALSE)
    stopifnot(length(par) == length(par2))

    apply(cbind(family, par, par2), 1, checkPars, cl = cl)

    ## return TRUE if all checks pass
    TRUE
}

## check for family/parameter consistency
checkPars <- function(x, cl) {
    family <- x[1]
    par <- x[2]
    par2 <- x[3]
    if ((family == 1 || family == 2)) {
        if (abs(par) >= 1)
            stop("\n In ", cl, ": ",
                 "The parameter of the Gaussian and t-copula has to be in the interval (-1,1).",
                 call. = FALSE)
        if (any((family == 2) & (par2 == 0)))
            stop("For t-copulas, 'par2' must be set.")
        if (family == 2 && par2 <= 2)
            stop("\n In ", cl, ": ",
                 "The degrees of freedom parameter of the t-copula has to be larger than 2.",
                 call. = FALSE)
    } else if ((family == 3 || family == 13) && (par <= 0 || par > 100)) {
        stop("\n In ", cl, ": ",
             "The parameter of the Clayton copula has to be in the interval (0,100]",
             call. = FALSE)
    } else if ((family == 4 || family == 14) && (par < 1 || par > 100)) {
        stop("\n In ", cl, ": ",
             "The parameter of the Gumbel copula has to be in the interval [1,100].",
             call. = FALSE)
    } else if (family == 5) {
        if (par == 0)
            stop("\n In ", cl, ": ",
                 "The parameter of the Frank copula has to be unequal to 0.",
                 call. = FALSE)
        if (abs(par) > 100)
            stop("\n In ", cl, ": ",
                 "The parameter of the Frank copula has to be in the interval [-100, 100].",
                 call. = FALSE)
    } else if ((family == 6 || family == 16) && (par <= 1 || par > 50)) {
        stop("\n In ", cl, ": ",
             "The parameter of the Joe copula has to be in the interval (1,50].",
             call. = FALSE)
    } else if ((family == 7 || family == 17)) {
        if (par <= 0 || par > 7)
            stop("\n In ", cl, ": ",
                 "The first parameter of the BB1 copula has to be in the interval [0,7]",
                 call. = FALSE)
        if (par2 < 1 || par2 > 7)
            stop("\n In ", cl, ": ",
                 "The second parameter of the BB1 copula has to be in the interval [1,7].",
                 call. = FALSE)
    } else if ((family == 8 || family == 18)) {
        if (par <= 0 || par > 6)
            stop("\n In ", cl, ": ",
                 "The first parameter of the BB6 copula has to be in the interval [1,6].",
                 call. = FALSE)
        if (par2 < 1 || par2 > 8)
            stop("\n In ", cl, ": ",
                 "The second parameter of the BB6 copula has to be in the interval [1,8].",
                 call. = FALSE)
    } else if ((family == 9 || family == 19)) {
        if (par < 1 || par > 6)
            stop("\n In ", cl, ": ",
                 "The first parameter of the BB7 copula has to be in the interval [1,6].",
                 call. = FALSE)
        if (par2 <= 0 || par2 > 75)
            stop("\n In ", cl, ": ",
                 "The second parameter of the BB7 copula has to be in the interval [0,75]",
                 call. = FALSE)
    } else if ((family == 10 || family == 20)) {
        if (par < 1 || par > 8)
            stop("\n In ", cl, ": ",
                 "The first parameter of the BB8 copula has to be in the interval [1,8].",
                 call. = FALSE)
        if (par2 < 1e-4 || par2 > 1)
            stop("\n In ", cl, ": ",
                 "The second parameter of the BB8 copula has to be in the interval [1e-4,1].",
                 call. = FALSE)
    } else if ((family == 23 || family == 33) && (par >= 0 || par < -100)) {
        stop("\n In ", cl, ": ",
             "The parameter of the rotated Clayton copula has to be be in the interval [-100,0)",
             call. = FALSE)
    } else if ((family == 24 || family == 34) && (par > -1 || par < -100)) {
        stop("\n In ", cl, ": ",
             "The parameter of the rotated Gumbel copula has to be in the interval [-100,-1].",
             call. = FALSE)
    } else if ((family == 26 || family == 36) && (par >= -1 || par < -50)) {
        stop("\n In ", cl, ": ",
             "The parameter of the rotated Joe copula has to be in the interval [-50,-1).",
             call. = FALSE)
    } else if ((family == 27 || family == 37)) {
        if (par >= 0 || par < -7)
            stop("\n In ", cl, ": ",
                 "The first parameter of the rotated BB1 copula has to be in the interval [-7,0]",
                 call. = FALSE)
        if (par2 > -1 || par < -7)
            stop("\n In ", cl, ": ",
                 "The second parameter of the rotated BB1 copula has to be in the interval [-7,-1].",
                 call. = FALSE)
    } else if ((family == 28 || family == 38)) {
        if (par >= 0 || par < -6)
            stop("\n In ", cl, ": ",
                 "The first parameter of the rotated BB6 copula has to be in the interval [-6,-1].",
                 call. = FALSE)
        if (par2 > -1 || par2 < -8)
            stop("\n In ", cl, ": ",
                 "The second parameter of the rotated BB6 copula has to be in the interval [-8,-1].",
                 call. = FALSE)
    } else if ((family == 29 || family == 39)) {
        if (par > -1 || par < -6)
            stop("\n In ", cl, ": ",
                 "The first parameter of the rotated BB7 copula has to be in the interval [-6,-1].",
                 call. = FALSE)
        if (par2 >= 0 || par2 < -75)
            stop("\n In ", cl, ": ",
                 "The second parameter of the rotated BB7 copula has to be in the interval [-75,0]",
                 call. = FALSE)
    } else if ((family == 30 || family == 40)) {
        if (par > -1 || par < -8)
            stop("\n In ", cl, ": ",
                 "The first parameter of the rotated BB8 copula has to be in the interval [-8,-1].",
                 call. = FALSE)
        if (par2 > -1e-4 || par2 < -1)
            stop("\n In ", cl, ": ",
                 "The second parameter of the rotated BB8 copula has to be in the interval [-1,-1e-4].",
                 call. = FALSE)
    } else if ((family == 41 || family == 51) && par <= 0) {
        stop("\n In ", cl, ": ",
             "The parameter of the reflection asymmetric copula has to be positive.",
             call. = FALSE)
    } else if ((family == 61 || family == 71) && par >= 0) {
        stop("\n In ", cl, ": ",
             "The parameter of the rotated reflection asymmetric copula has to be negative.",
             call. = FALSE)
    } else if (family == 42) {
        a <- par
        b <- par2
        limA <- (b - 3 - sqrt(9 + 6 * b - 3 * b^2))/2
        if (abs(b) > 1)
            stop("\n In ", cl, ": ",
                 "The second parameter of the two-parametric asymmetric copulas has to be in the interval [-1,1]",
                 call. = FALSE)
        if (a > 1 || a < limA)
            stop("\n In ", cl, ": ",
                 "The first parameter of the two-parametric asymmetric copula has to be in the interval [limA(par2),1]",
                 call. = FALSE)
    } else if (family == 104 || family == 114 || family == 204 || family == 214) {
        if (par < 1)
            stop("\n In ", cl, ": ",
                 "Please choose 'par' of the Tawn copula in [1,oo).",
                 call. = FALSE)
        if (par2 < 0 || par2 > 1)
            stop("\n In ", cl, ": ",
                 "Please choose 'par2' of the Tawn copula in [0,1].",
                 call. = FALSE)
    } else if ((family == 124 || family == 134 || family == 224 || family == 234)) {
        if (par > -1)
            stop("\n In ", cl, ": ",
                 "Please choose 'par' of the rotated Tawn copula in (-oo,-1].",
                 call. = FALSE)
        if (par2 < 0 || par2 > 1)
            stop("\n In ", cl, ": ",
                 "Please choose 'par2' of the rotated Tawn copula in [0,1].",
                 call. = FALSE)
    }
}


## check for family/parameter consistency
adjustPars <- function(family, par, par2) {
    if ((family %in% c(3, 13, 23, 33)) && (abs(par) > 100)) {
        par <- sign(par) * 100
    } else if ((family %in% c(4, 14, 24, 34)) && (abs(par) > 100)) {
        par <- sign(par) * 100
    } else if ((family == 5) && (abs(par) > 100)) {
        par <- sign(par) * 100
    } else if ((family %in% c(6, 16, 26, 36)) && (abs(par) > 50)) {
        par <- sign(par) * 50
    } else if (family %in% c(7, 17, 27, 37)) {
        if (abs(par) > 7)
            par <- sign(par) * 7
        if (abs(par2) > 7)
            par2 <- sign(par) * 7
    } else if (family %in% c(8, 18, 28, 38)) {
        if (abs(par) > 6)
            par <- sign(par) * 6
        if (abs(par2) > 8)
            par2 <- sign(par) * 8
    } else if (family %in% c(9, 19, 29, 39)) {
        if (abs(par) > 6)
            par <- sign(par) * 6
        if (abs(par2) > 75)
            par2 <- sign(par) * 75
    } else if (family %in% c(10, 20, 30, 40)) {
        if (abs(par) > 8)
            par <- sign(par) * 8
    }
    c(par, par2)
}

BiCopCheckTaus <- function(family, tau) {
    cl <- match.call()[1]
    ## check for family/tau consistency
    checkTaus<- function(x) {
        if (family %in% c(3, 13) && tau <= 0)
            stop("\n In ", cl, ": ",
                 "Clayton copula cannot be used for tau<=0.",
                 call. = FALSE)
        if (family %in% c(4, 14) && tau < 0)
            stop("\n In ", cl, ": ",
                 "Gumbel copula cannot be used for tau<0.",
                 call. = FALSE)
        if (family == 5 && tau == 0)
            stop("\n In ", cl, ": ",
                 "Frank copula cannot be used for tau=0",
                 call. = FALSE)
        if (family %in% c(6, 16) && tau < 0)
            stop("\n In ", cl, ": ",
                 "Joe copula cannot be used for tau<0.",
                 call. = FALSE)
        if (family %in% c(23, 33) && tau >= 0)
            stop("\n In ", cl, ": ",
                 "Rotated Clayton copula cannot be used for tau>=0.",
                 call. = FALSE)
        if (family %in% c(24, 34) && tau > 0)
            stop("\n In ", cl, ": ",
                 "Rotated Gumbel copula cannot be used for tau>0.",
                 call. = FALSE)
        if (family %in% c(26, 36) && tau > 0)
            stop("\n In ", cl, ": ",
                 "Rotated Joe copula cannot be used for tau>0.",
                 call. = FALSE)
    }
    apply(cbind(family, tau), 1, checkTaus)

    ## return TRUE if all checks pass
    TRUE
}
