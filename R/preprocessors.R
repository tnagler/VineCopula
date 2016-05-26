preproc <- function(args, ..., na.txt = NULL) {
    # augment arguments
    args$na.txt <- na.txt
    args$n <- length(args$u1)

    # what has to be done?
    tasks <- list(...)

    # perform all tasks sequentially
    for (i in seq_along(tasks)) {
        stopifnot(is.function(tasks[[i]]))
        args <- tasks[[i]](args)
    }

    # return preprocessed arguments
    args
}

#' check if all data have been provided and have same length
check_u <- function(args) {
    if (is.symbol(args$u1) | is.symbol(args$u2))
        stop("\n In ", args$call[1], ": ",
             "u1 and/or u2 are missing.",
             call. = FALSE)
    if (is.null(args$u1) == TRUE || is.null(args$u2) == TRUE)
        stop("\n In ", args$call[1], ": ",
             "u1 and/or u2 are not set or have length zero.",
             call. = FALSE)
    if (length(args$u1) != length(args$u2))
        stop("\n In ", args$call[1], ": ",
             "Lengths of 'u1' and 'u2' do not match.",
             call. = FALSE)

    args
}

#' set all NA values to 0.5, but store the index (will be reset to NA)
fix_nas <- function(args) {
    if (any(is.na(args$u1 + args$u2))) {
        # set warning message
        args$msg <- paste0("In ",
                           args$call[1],
                           ": ",
                           "Some of the data are NA.",
                           args$na.txt)
        # set NA values to 0.5 so that C code can operate (will be reset to NA)
        args$na.ind <- which(is.na(args$u1 + args$u2))
        args$u1[args$na.ind] <- 0.5
        args$u2[args$na.ind] <- 0.5
    } else {
        args$msg <- na.ind <- NULL
    }

    args
}

#' reset output to NA if input was
reset_nas <- function(out, args) {
    if (!is.null(args$na.ind))
        # set output to NA if input was
        out[args$na.ind] <- NA
    # print warning if necessary
    if (!is.null(args$msg))
        warning(args$msg, call. = FALSE)

    out
}

#' remove all NA values from the data
remove_nas <- function(args) {
    if (any(is.na(args$u1 + args$u2))) {
        # set warning message
        msg <- paste0("In ",
                      args$call[1],
                      ": ",
                      "Some of the data are NA.",
                      args$na.txt)
        warning(msg, call. = FALSE)
        # remove NAs
        args$na.ind <- which(is.na(args$u1 + args$u2))
        args$u1 <- args$u1[-args$na.ind]
        args$u2 <- args$u2[-args$na.ind]
    } else {
        args$msg <- na.ind <- NULL
    }

    args$n <- length(args$u1)
    args
}

#' check if all data are in (0, 1)^2
check_if_01 <- function(args) {
    if (any(args$u1 > 1) || any(args$u1 < 0))
        stop("\n In ", args$call[1], ": ",
             "Data has be in the interval [0,1].",
             call. = FALSE)
    if (any(args$u2 > 1) || any(args$u2 < 0))
        stop("\n In ", args$call[1], ": ",
             "Data has be in the interval [0,1].",
             call. = FALSE)

    args
}

#' make sure that family, par, par2 have the same length
match_spec_lengths <- function(args) {
    n <- length(args$u1)
    # if one vector is size n, expand all vectors to size n
    if (any(c(length(args$family), length(args$par), length(args$par2)) == n)) {
        if (length(args$family) == 1)
            args$family <- rep(args$family, n)
        if (length(args$par) == 1)
            args$par <- rep(args$par, n)
        if (length(args$par2) == 1)
            args$par2 <- rep(args$par2, n)
    }

    # check if input size was ok
    if (!(length(args$family) %in% c(1, n)))
        stop("\n In ", args$call[1], ": ",
             "'family' has to be a single number or a size n vector",
             call. = FALSE)
    if (!(length(args$par) %in% c(1, n)))
        stop("\n In ", args$call[1], ": ",
             "'par' has to be a single number or a size n vector",
             call. = FALSE)
    if (!(length(args$par2) %in% c(1, n)))
        stop("\n In ", args$call[1], ": ",
             "'par2' has to be a single number or a size n vector",
             call. = FALSE)

    args
}

#' extract family and parameters if BiCop object is provided
extract_from_BiCop <- function(args) {
    # set dummys if family and par are missing (-> when obj is provided)
    if (is.symbol(args$family))
        args$family <- NA
    if (is.symbol(args$par))
        args$par <- NA
    # has the short version BiCop...(u1, u2, obj) been used?
    if (class(args$family) == "BiCop")
        args$obj <- args$family
    # store info from obj into family par and par2
    if (!is.null(args$obj)) {
        stopifnot(class(args$obj) == "BiCop")
        args$family <- args$obj$family
        args$par <- args$obj$par
        args$par2 <- args$obj$par2
    }
    # set parameter if independence copula was specified without
    if (is.null(args$par) & (all(args$family == 0)))
        args$par <- 0
    # stop if insufficient info was procided
    if (any(is.na(args$family)) | any(is.na(args$par)))
        stop("\n In ", args$call[1], ": ",
             "Provide either 'family' and 'par' or 'obj'",
             call. = FALSE)

    args
}

#' sanity checks for family and parameters
check_fam_par <- function(args) {
    if (args$check.pars) {
        # check for family/parameter consistency (if not disabled)
        BiCopCheck(args$family, args$par, args$par2, call = args$call)
    } else {
        # allow zero parameter for Clayton an Frank otherwise
        args$family[(args$family %in% c(3, 13, 23, 33)) & (args$par == 0)] <- 0
        args$family[(args$family == 5) & (args$par == 0)] <- 0
    }

    args
}

#' check if more than one observation has been provided
check_nobs <- function(args) {
    if (length(args$u1) < 2)
        stop("\n In ", args$call[1], ": ",
             "Number of observations has to be at least 2.",
             call. = FALSE)
    args
}

#' add or remove families (rotations or negative index)
prep_familyset <- function(args) {
    if (is.na(args$familyset[1]))
        args$familyset <- allfams
    if (args$rotations)
        args$familyset <- with_rotations(args$familyset)
    if (any(args$familyset < 0)) {
        if (length(unique(sign(args$familyset))) != 1)
            stop("'familyset' must not contain positive AND negative numbers")
        args$familyset <- setdiff(allfams, -args$familyset)
    }
    args
}
