fix_nas <- function(args, add.txt = NULL) {
    if (any(is.na(args$u1 + args$u2))) {
        # set warning message
        args$msg <- paste0("In ",
                           args$call[1],
                           ": ",
                           "Some of the data are NA.",
                           add.txt)
        # set NA values to 0.5 so that C code can operate (will be reset to NA)
        args$na.ind <- which(is.na(args$u1 + args$u2))
        args$u1[args$na.ind] <- 0.5
        args$u2[args$na.ind] <- 0.5
    } else {
        args$msg <- na.ind <- NULL
    }

    args
}

reset_nas <- function(out, args) {
    if (!is.null(args$na.ind))
        # set output to NA if input was
        out[args$na.ind] <- NA
    # print warning if necessary
    if (!is.null(args$msg))
        warning(args$msg, call. = FALSE)

    out
}

check_if_01 <- function(args) {
    if (any(args$u1 > 1) || any(args$u1 < 0))
        stop("\n In ", args$call[1], ": ",
             "Data has be in the interval [0,1].",
             call. = FALSE)
    if (any(args$u2 > 1) || any(args$u2 < 0))
        stop("\n In ", args$call[1], ": ",
             "Data has be in the interval [0,1].",
             call. = FALSE)
}

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


extract_from_BiCop <- function(args) {
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
