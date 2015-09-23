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
      if (family == 134) {
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
      if (family == 234) {
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
      } else {
        res <- rep(NA, length(u1))
      }

  ## return results
  res
}
