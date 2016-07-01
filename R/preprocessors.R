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

## check if all data have been provided and have same length
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
             "Lengths of u1 and u2 do not match.",
             call. = FALSE)
    if ((!is.numeric(args$u1) & !all(is.na(args$u1))))
        stop("\n In ", args$call[1], ": ",
             "u1 has to be a numeric vector.",
             call. = FALSE)
    if ((!is.numeric(args$u2) & !all(is.na(args$u2))))
        stop("\n In ", args$call[1], ": ",
             "u2 has to be a numeric vector.",
             call. = FALSE)

    args
}



## set all NA values to 0.5, but store the index (will be reset to NA)
fix_nas <- function(args) {
    if (!is.null(args$data)) {
        if (any(is.na(args$data))) {
            num.na <- sum(!complete.cases(args$data))
            freq.na <- round(num.na / nrow(args$data) * 100, 1)
            # set warning message
            args$msg <- paste0(" In ",
                               args$call[1],
                               ": ",
                               num.na, " of the evaluation points",,
                               " (", freq.na, "%) contain",
                               ifelse(num.na == 1, "s", ""), " NAs.",
                               args$na.txt)
            # set NA values to 0.5 so that C code can operate (will be reset to NA)
            args$na.ind <- which(!complete.cases(args$data))
            args$data[args$na.ind, ] <- 0.5
        }
    } else if (any(is.na(args$u1 + args$u2))) {
        num.na <- sum(!complete.cases(args$u1 + args$u2))
        freq.na <- round(num.na / length(args$u1) * 100, 1)
        # set warning message
        args$msg <- paste0(" In ",
                           args$call[1],
                           ": ",
                           num.na, " of the evaluation points",
                           " (", freq.na, "%) contain",
                           ifelse(num.na == 1, "s", ""), " NAs.",
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

## reset output to NA if input was
reset_nas <- function(out, args) {
    # set output to NA if input was
    if (is.vector(out)) {
        if (!is.null(args$na.ind))
            out[args$na.ind] <- NA
    } else {
        if (length(dim(out)) == 2)
            out[args$na.ind, ] <- NA
        if (length(dim(out)) == 3)
            out[, , args$na.ind] <- NA
    }
    # print warning if necessary
    if (!is.null(args$msg))
        warning(args$msg, call. = FALSE)

    out
}

## remove all NA values from the data
remove_nas <- function(args) {
    if (!is.null(args$data)) {
        if (any(is.na(args$data))) {
            # set warning message
            num.na <- sum(!complete.cases(args$data))
            freq.na <- round(num.na / nrow(args$data) * 100, 1)
            warning(" In ",
                    args$call[1],
                    ": ",
                    num.na, " observation", ifelse(num.na > 1, "s", ""),
                    " (", freq.na, "%) contain",
                    ifelse(num.na == 1, "s", ""), " NAs.",
                    args$na.txt,
                    call. = FALSE)
            # remove NAs
            args$na.ind <- which(!complete.cases(args$data))
            args$data <- args$data[complete.cases(args$data), , drop = FALSE]
            args$n <- nrow(args$data)
        }
    } else {
        if (any(is.na(args$u1 + args$u2))) {
            # set warning message
            num.na <- sum(!complete.cases(args$u1 + args$u2))
            freq.na <- round(num.na / length(args$u1) * 100, 1)
            warning(" In ",
                    args$call[1],
                    ": ",
                    num.na, " observation", ifelse(num.na > 1, "s", ""),
                    " (", freq.na, "%) contain",
                    ifelse(num.na == 1, "s", ""), " NAs.",
                    args$na.txt,
                    call. = FALSE)
            # remove NAs
            args$na.ind <- which(is.na(args$u1 + args$u2))
            args$u1 <- args$u1[-args$na.ind]
            args$u2 <- args$u2[-args$na.ind]
        } else {
            args$msg <- na.ind <- NULL
        }
        args$n <- length(args$u1)
    }

    args
}

## check if all data are in (0, 1)^2
check_if_01 <- function(args) {
    u <- args$data
    if (!is.null(u)) {
        i <- complete.cases(u)
        if (sum(i) > 0) {
            if (any(u[i] > 1) || any(u[i] < 0))
                stop("\n In ", args$call[1], ": ",
                     "Data have to be in the interval [0,1]^d.",
                     call. = FALSE)
        }
    } else {
        if (any(args$u1 > 1) || any(args$u1 < 0))
            stop("\n In ", args$call[1], ": ",
                 "Data have to be in the interval [0,1].",
                 call. = FALSE)
        if (any(args$u2 > 1) || any(args$u2 < 0))
            stop("\n In ", args$call[1], ": ",
                 "Data have to be in the interval [0,1].",
                 call. = FALSE)
    }

    args
}

## make sure that family, par, par2 have the same length
match_spec_lengths <- function(args) {
    n <- ifelse(!is.null(args$u1),
                length(args$u1),
                max(length(args$family), length(args$par), length(args$par2)))
    args$family <- c(args$family)
    args$par <- c(args$par)
    args$par2 <- c(args$par2)

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

expand_lengths <- function(args) {
    n <- ifelse(!is.null(args$u1),
                length(args$u1),
                max(length(args$u1,
                           length(args$family),
                           length(args$par),
                           length(args$par2))))
    args$family <- c(args$family)
    args$par <- c(args$par)
    args$par2 <- c(args$par2)

    # if one vector is size n, expand all vectors to size n
    if (length(args$family) == 1)
        args$family <- rep(args$family, n)
    if (length(args$par) == 1)
        args$par <- rep(args$par, n)
    if (length(args$par2) == 1)
        args$par2 <- rep(args$par2, n)

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


## extract family and parameters if BiCop object is provided
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

## sanity checks for family and parameters
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

## check if more than one observation has been provided
check_nobs <- function(args) {
    if (!is.null(args$data)) {
        if (nrow(args$data) < 2)
            stop("\n In ", args$call[1], ": ",
                 "Number of observations has to be at least 2.",
                 call. = FALSE)
    } else {
        if (length(args$u1) < 2)
            stop("\n In ", args$call[1], ": ",
                 "Number of observations has to be at least 2.",
                 call. = FALSE)
    }

    args
}

## add or remove families (rotations or negative index)
prep_familyset <- function(args) {
    if (is.na(args$familyset[1]))
        args$familyset <- allfams
    if (args$rotations)
        args$familyset <- with_rotations(args$familyset)
    if (length(unique(sign(args$familyset[args$familyset != 0]))) > 1) {
        stop("\n In ", args$call[1], ": ",
             "'familyset' must not contain positive AND negative numbers.",
             call. = FALSE)
        args$familyset <- setdiff(allfams, -args$familyset)
    }
    if (!all(abs(args$familyset) %in% allfams))
        stop("\n In ", args$call[1], ": ",
             "Copula family not implemented.",
             call. = FALSE)
    args
}


## check if familyset has at least one family corresponding to empirical tau
check_fam_tau <- function(args) {
    if (is.null(args$weights))
        args$weights <- NA
    args$emp_tau <- fasttau(args$u1, args$u2, args$weights)
    if ((args$emp_tau > 0) & !any(args$familyset %in% c(0, posfams))) {
        stop("\n In ", args$call[1], ": ",
             "Empirical Kendall's tau is positive, but familyset only ",
             "contains families for negative dependence.",
             call. = FALSE)
    } else if ((args$emp_tau < 0) & !any(args$familyset %in% c(0, negfams))){
        stop("\n In ", args$call[1], ": ",
             "Empirical Kendall's tau is negative, but familyset only ",
             "contains families for positive dependence.",
             call. = FALSE)
    }

    args
}

todo_fams <- function(args) {
    if (args$emp_tau < 0) {
        todo <- c(0, negfams)
    } else if (args$emp_tau > 0) {
        todo <- c(0, posfams)
    } else {
        todo <- allfams
    }
    args$familyset <- todo[which(todo %in% args$familyset)]
    args
}

## check max.BB and max.df specifications
check_est_pars <- function(args) {
    if (!is.null(args$max.df)) {
        if (args$max.df <= 2)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the degrees of freedom parameter has to be",
                 "larger than 2.",
                 call. = FALSE)
    } else {
        args$max.df <- 30
    }

    if (!is.null(args$max.BB)) {
        if (!is.list(args$max.BB))
            stop("\n In ", args$call[1], ": ",
                 "max.BB has to be a list.",
                 call. = FALSE)
        if (!all(names(args$max.BB) == c("BB1", "BB6", "BB7", "BB8")))
            stop("\n In ", args$call[1], ": ",
                 'max.BB has to be a named list with entries "BB1", "BB6", "BB7", "BB8".',
                 call. = FALSE)
        if (any(sapply(args$max.BB, length) != 2))
            stop("\n In ", args$call[1], ": ",
                 'All components of max.BB have to be two-dimensional numeric vectors.',
                 call. = FALSE)
        if (args$max.BB$BB1[1] < 0.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the first parameter of the BB1 copula should ",
                 "be greater than 0.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB1[2] < 1.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the second parameter of the BB1 copula should ",
                 "be greater than 1.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB6[1] < 1.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the first parameter of the BB6 copula should ",
                 "be greater than 1.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB6[2] < 1.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the second parameter of the BB6 copula should ",
                 "be greater than 1.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB7[1] < 1.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the first parameter of the BB7 copula should ",
                 "be greater than 1.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB7[2] < 0.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the second parameter of the BB7 copula should ",
                 "be greater than 0.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB8[1] < 1.001)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the first parameter of the BB8 copula should ",
                 "be greater than 1.001 (lower bound for estimation).",
                 call. = FALSE)
        if (args$max.BB$BB8[2] < 0.001 || args$max.BB$BB8[2] > 1)
            stop("\n In ", args$call[1], ": ",
                 "The upper bound for the second parameter of the BB8 copula should ",
                 "be in the interval [0.001,1].",
                 call. = FALSE)
    } else {
        args$max.BB <- list(BB1 = c(5, 6),
                            BB6 = c(6, 6),
                            BB7 = c(5, 6),
                            BB8 = c(6, 1))
    }

    if (!is.null(args$family)) {
        if (!(all(args$family %in% allfams)))
            stop("\n In ", args$call[1], ": ",
                 "Copula family not implemented.",
                 call. = FALSE)
    }

    if (!is.null(args$familyset)) {
        if (!(all(abs(args$familyset) %in% allfams)))
            stop("\n In ", args$call[1], ": ",
                 "Copula family not implemented.",
                 call. = FALSE)
    }

    if (!is.null(args$method)) {
        if (args$method != "mle" && args$method != "itau")
            stop("\n In ", args$call[1], ": ",
                 "Estimation method has to be either 'mle' or 'itau'.",
                 call. = FALSE)
        if (!is.null(args$family)) {
            if ((args$method == "itau") && (!(args$family %in% c(0, allfams[onepar])))) {
                warning(" In ", args$call[1], ": ",
                        "For two parameter copulas the estimation method 'itau' cannot ",
                        "be used. The method is automatically set to 'mle'.",
                        call. = FALSE)
                args$method <- "mle"
            }
        }
        if (!is.null(args$familyset)) {
            if ((args$method == "itau") && (!all(args$familyset %in% c(0, allfams[onepar])))) {
                warning(" In ", args$call[1], ": ",
                        "For two parameter copulas the estimation method 'itau' cannot ",
                        "be used. The method is automatically set to 'mle'.",
                        call. = FALSE)
                args$method <- "mle"
            }
        }
    } else {
        args$method <- "mle"
    }

    if (!is.null(args$se)) {
        if (is.logical(args$se) == FALSE)
            stop("\n In ", args$call[1], ": ",
                 "se has to be a logical variable (TRUE or FALSE).",
                 call. = FALSE)
    } else {
        args$se = TRUE
    }

    args$weights <- ifelse(is.null(args$weights), NA, args$weights)

    args
}

check_data <- function(args) {
    if (is.symbol(args$data))
        stop("\n In ", args$call[1], ": ",
             "Argument data is missing.",
             call. = FALSE)
    if (is.null(args$data))
        stop("\n In ", args$call[1], ": ",
             "Argument data is missing.",
             call. = FALSE)
    if (is.vector(args$data)) {
        args$data <- t(as.matrix(args$data))
    } else {
        args$data <- as.matrix(args$data)
    }
    if (!is.numeric(args$data) & !all(is.na(args$data)))
        stop("\n In ", args$call[1], ": ",
             "Data have to be numeric.",
             call. = FALSE)
    if (ncol(args$data) < 2)
        stop("\n In ", args$call[1], ": ",
             "Dimension has to be at least 2.",
             call. = FALSE)
    if (is.null(colnames(args$data)))
        colnames(args$data) <- paste("V", seq.int(ncol(args$data)), sep = "")
    args$n <- nrow(args$data)
    args$d <- ncol(args$data)

    args
}

check_RVMs <- function(args) {
    if (!is.null(args$RVM))
        args$RVM <- check_RVM(args$RVM, "RVM",
                              args$call, args$d, args$check.pars)
    if (!is.null(args$RVM1))
        args$RVM1 <- check_RVM(args$RVM1, "RVM1",
                               args$call, args$d, args$check.pars)
    if (!is.null(args$RVM2))
        args$RVM2 <- check_RVM(args$RVM2, "RVM2",
                               args$call, args$d, args$check.pars)

    args
}

check_RVM <- function(RVM, name, call, d, check.pars) {
    if (!is(RVM, "RVineMatrix"))
        stop("\n In ", call[1], ": ",
             name,  " has to be an RVineMatrix object.",
             call. = FALSE)
    if (!is.null(d)) {
        if (d != dim(RVM))
            stop("\n In ", call[1], ": ",
                 "Dimensions of data and ", name, " do not match.",
                 call. = FALSE)
    }

    if (any(!is.na(RVM$par)) & (nrow(RVM$par) != ncol(RVM$par)))
        stop("\n In ", call[1], ": ",
             "Parameter matrix has to be quadratic.",
             call. = FALSE)
    if (any(!is.na(RVM$par2)) & (nrow(RVM$par2) != ncol(RVM$par2)))
        stop("\n In ", call[1], ": ",
             "Second parameter matrix has to be quadratic.",
             call. = FALSE)

    if (is.null(check.pars))
        check.pars <- TRUE
    if (!all(RVM$par %in% c(0, NA))) {
        for (i in 2:dim(RVM$Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if (check.pars) {
                    BiCopCheck(RVM$family[i, j],
                               RVM$par[i, j],
                               RVM$par2[i, j],
                               call = call)
                } else {
                    indep <- (RVM$family[i, j] %in% c(5, 3, 13, 23, 33)) & (RVM$par[i, j] == 0)
                    if (indep)
                        RVM$family[i, j] <- 0
                }
            }
        }
    }

    RVM
}

prep_RVMs <- function(args) {
    if (!is.null(args$RVM)) {
        args$RVM$par[is.na(args$RVM$par)] <- 0
        args$RVM$par[upper.tri(args$RVM$par, diag = TRUE)] <- 0
        args$RVM$par2[is.na(args$RVM$par2)] <- 0
        args$RVM$par2[upper.tri(args$RVM$par2, diag = TRUE)] <- 0
    } else {
        args$RVM1$par[is.na(args$RVM1$par)] <- 0
        args$RVM1$par[upper.tri(args$RVM1$par, diag = TRUE)] <- 0
        args$RVM1$par2[is.na(args$RVM1$par2)] <- 0
        args$RVM1$par2[upper.tri(args$RVM1$par2, diag = TRUE)] <- 0

        args$RVM2$par[is.na(args$RVM2$par)] <- 0
        args$RVM2$par[upper.tri(args$RVM2$par, diag = TRUE)] <- 0
        args$RVM2$par2[is.na(args$RVM2$par2)] <- 0
        args$RVM2$par2[upper.tri(args$RVM2$par2, diag = TRUE)] <- 0
    }

    args
}

check_matrix <- function(args) {
    if (is.symbol(args$Matrix))
        stop("\n In ", args$call[1], ": ",
             "Matrix is missing.",
             call. = FALSE)
    if (nrow(args$Matrix) != ncol(args$Matrix))
        stop("\n In ", args$call[1], ": ",
             "Structure matrix has to be quadratic.",
             call. = FALSE)
    if (max(args$Matrix) > nrow(args$Matrix))
        stop("\n In ", args$call[1], ": ",
             "Structure matrix can only contain numbers 0:",
             nrow(args$Matrix), ".",
             call. = FALSE)

    args$Matrix[is.na(args$Matrix)] <- 0
    args$Matrix <- ToLowerTri(args$Matrix)

    if (RVineMatrixCheck(args$Matrix) != 1)
        stop("\n In ", args$call[1], ": ",
             "Matrix is not a valid R-vine matrix.",
             call. = FALSE)

    args
}

check_parmat <- function(args) {
    args$family <- as.matrix(args$family)
    if (nrow(args$family) != ncol(args$family))
        stop("\n In ", args$call[1], ": ",
             "family has to be quadratic.",
             call. = FALSE)

    args$family[is.na(args$family)] <- 0
    args$family <- ToLowerTri(args$family)
    args$family[upper.tri(args$family, diag = T)] <- 0

    args
}

check_fammat <- function(args) {
    args$par <- as.matrix(args$par)
    if (nrow(args$par) != ncol(args$par))
        stop("\n In ", args$call[1], ": ",
             "par has to be quadratic.",
             call. = FALSE)

    args$par[is.na(args$par)] <- 0
    args$par <- ToLowerTri(args$par)
    args$par[upper.tri(args$par, diag = T)] <- 0

    args
}

check_par2mat <- function(args) {
    args$par2 <- as.matrix(args$par2)
    if (nrow(args$par2) != ncol(args$par2))
        stop("\n In ", args$call[1], ": ",
             "par2 has to be quadratic.",
             call. = FALSE)

    args$par2[is.na(args$par2)] <- 0
    args$par2 <- ToLowerTri(args$par2)
    args$par2[upper.tri(args$par2, diag = T)] <- 0

    args
}

set_dims <- function(family, par = 0, par2 = 0, tau = 0) {
    dims <- max(length(family), length(par), length(par2), length(tau))
    if (!is.null(dim(family)))
        dims <- dim(family)
    if (!is.null(dim(par)))
        dims <- dim(par)
    if (!is.null(dim(par2)))
        dims <- dim(par2)
    if (!is.null(dim(tau)))
        dims <- dim(tau)
    dims
}


