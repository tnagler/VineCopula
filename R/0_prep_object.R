copulaFromFamilyIndex <- function(family, par, par2 = 0) {
    constr <- switch(paste("fam", family, sep = ""),
                     fam0 = function(par) indepCopula(), 
                     fam1 = function(par) normalCopula(par[1]),
                     fam2 = function(par) tCopula(par[1], df = par[2]),
                     fam3 = function(par) claytonCopula(par[1]),
                     fam4 = function(par) gumbelCopula(par[1]), 
                     fam5 = function(par) frankCopula(par[1]),
                     fam6 = function(par) joeBiCopula(par[1]), 
                     fam7 = BB1Copula,
                     fam8 = BB6Copula,
                     fam9 = BB7Copula,
                     fam10 = BB8Copula, 
                     fam13 = function(par) surClaytonCopula(par[1]),
                     fam14 = function(par) surGumbelCopula(par[1]), 
                     fam16 = function(par) surJoeBiCopula(par[1]),
                     fam17 = surBB1Copula,
                     fam18 = surBB6Copula, 
                     fam19 = surBB7Copula,
                     fam20 = surBB8Copula,
                     fam23 = function(par) r90ClaytonCopula(par[1]), 
                     fam24 = function(par) r90GumbelCopula(par[1]),
                     fam26 = function(par) r90JoeBiCopula(par[1]), 
                     fam27 = r90BB1Copula,
                     fam28 = r90BB6Copula,
                     fam29 = r90BB7Copula,
                     fam30 = r90BB8Copula, 
                     fam33 = function(par) r270ClaytonCopula(par[1]),
                     fam34 = function(par) r270GumbelCopula(par[1]), 
                     fam36 = function(par) r270JoeBiCopula(par[1]),
                     fam37 = r270BB1Copula,
                     fam38 = r270BB6Copula, 
                     fam39 = r270BB7Copula,
                     fam40 = r270BB8Copula,
                     fam104 = tawnT1Copula,
                     fam114 = surTawnT1Copula, 
                     fam124 = r90TawnT1Copula,
                     fam134 = r270TawnT1Copula,
                     fam204 = tawnT2Copula, 
                     fam214 = surTawnT2Copula,
                     fam224 = r90TawnT2Copula,
                     fam234 = r270TawnT2Copula)
    constr(c(par, par2))
}

# generic fitting make fitCopula from copula generic
setGeneric("fitCopula", fitCopula)

####################### generic wrapper functions to the VineCopula package ##

# density from BiCopPDF
linkVineCop.PDF <- function(u, copula, log = FALSE) {
    param <- copula@parameters
    
    if (length(param) == 1) 
        param <- c(param, 0)
    n <- nrow(u)
    fam <- copula@family
    
    # coplik = RLL_mod_separate(fam, n, u, param)[[7]]
    coplik <- .C("LL_mod_seperate",
                 as.integer(fam),
                 as.integer(n),
                 as.double(u[, 1]),
                 as.double(u[, 2]),
                 as.double(param[1]),
                 as.double(param[2]),
                 as.double(rep(0, n)),
                 PACKAGE = "VineCopula")[[7]]
    if (log) return(coplik) else return(exp(coplik))
}

# cdf from BiCopCDF

# for 'standard' copulas: family %in% c(3:10)
linkVineCop.CDF <- function(u, copula) {
    param <- copula@parameters
    if (!is.matrix(u)) u <- matrix(u, ncol = 2)
    n <- nrow(u)
    fam <- copula@family
    
    res <- .C("archCDF",
              as.double(u[, 1]),
              as.double(u[, 2]),
              as.integer(n),
              as.double(param), 
              as.integer(fam),
              as.double(rep(0, n)),
              PACKAGE = "VineCopula")[[6]]
    return(res)
}

# for survival copulas: family %in% c(13, 14, 16:20)
linkVineCop.surCDF <- function(u, copula) {
    param <- copula@parameters
    if (!is.matrix(u)) u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    n <- nrow(u)
    fam <- copula@family
    
    res <- u1 + u2 - 1 + .C("archCDF",
                            as.double(1 - u1),
                            as.double(1 - u2),
                            as.integer(n), 
                            as.double(param),
                            as.integer(fam - 10),
                            as.double(rep(0, n)),
                            PACKAGE = "VineCopula")[[6]]
    return(res)
}

# for 90 deg rotated copulas: family %in% c(23, 24, 26:30)
linkVineCop.r90CDF <- function(u, copula) {
    param <- copula@parameters
    if (!is.matrix(u)) u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    n <- nrow(u)
    fam <- copula@family
    
    u2 - .C("archCDF",
            as.double(1 - u1),
            as.double(u2),
            as.integer(n),
            as.double(-param), 
            as.integer(fam - 20),
            as.double(rep(0, n)),
            PACKAGE = "VineCopula")[[6]]
}

# for 270 deg rotated copulas: family %in% c(33, 34, 36:40)
linkVineCop.r270CDF <- function(u, copula) {
    param <- copula@parameters
    if (!is.matrix(u)) u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    n <- nrow(u)
    fam <- copula@family
    
    u1 - .C("archCDF",
            as.double(u1),
            as.double(1 - u2),
            as.integer(n),
            as.double(-param), 
            as.integer(fam - 30),
            as.double(rep(0, n)),
            PACKAGE = "VineCopula")[[6]]
}

## for Tawn
# TawnC(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
linkVineCop.CDFtawn <- function(u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[, 1]
  u2 <- u[, 2]
  n <- nrow(u)
  fam <- copula@family
  
  if (fam == 104) {
    par3 <- 1
    res <- .C("TawnC",
              as.double(u1),
              as.double(u2),
              as.integer(n),
              as.double(param[1]),
              as.double(param[2]), 
              as.double(par3), 
              as.double(rep(0, n)),
              PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 114) {
    par3 <- 1
    res <- u1 + u2 - 1 + .C("TawnC",
                            as.double(1-u1),
                            as.double(1-u2),
                            as.integer(n),
                            as.double(param[1]),
                            as.double(param[2]), 
                            as.double(par3), 
                            as.double(rep(0, n)),
                            PACKAGE = "VineCopula")[[7]]
  }
  if (fam == 124) {
    par3 <- 1
    res <- u2 - .C("TawnC",
                   as.double(1-u1),
                   as.double(u2),
                   as.integer(n),
                   as.double(-param[1]),
                   as.double(param[2]), 
                   as.double(par3), 
                   as.double(rep(0, n)),
                   PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 134) {
    par3 <- 1
    res <- u1 - .C("TawnC",
                   as.double(u1),
                   as.double(1-u2),
                   as.integer(n),
                   as.double(-param[1]),
                   as.double(param[2]), 
                   as.double(par3), 
                   as.double(rep(0, n)),
                   PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 204) {
    par2 <- 1
    res <- .C("TawnC",
              as.double(u1),
              as.double(u2),
              as.integer(n),
              as.double(param[1]),
              as.double(par2), 
              as.double(param[2]), 
              as.double(rep(0, n)),
              PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 214) {
    par2 <- 1
    res <- u1 + u2 - 1 + .C("TawnC",
                            as.double(1-u1),
                            as.double(1-u2),
                            as.integer(n),
                            as.double(param[1]),
                            as.double(par2), 
                            as.double(param[2]), 
                            as.double(rep(0, n)),
                            PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 224) {
    par2 <- 1
    res <- u2 - .C("TawnC",
                   as.double(1-u1),
                   as.double(u2),
                   as.integer(n),
                   as.double(-param[1]),
                   as.double(par2), 
                   as.double(param[2]), 
                   as.double(rep(0, n)),
                   PACKAGE = "VineCopula")[[7]]
  } 
  if (fam == 234) {
    par2 <- 1
    res <- u1 - .C("TawnC",
                   as.double(u1),
                   as.double(1-u2),
                   as.integer(n),
                   as.double(-param[1]),
                   as.double(par2), 
                   as.double(param[2]), 
                   as.double(rep(0, n)),
                   PACKAGE = "VineCopula")[[7]]
  }
  return(res)
}

## derivtives/h-function from BiCopHfunc ddu
linkVineCop.ddu <- function(u, copula) {
    param <- copula@parameters
    
    if (length(param) == 1) param <- c(param, 0)
    
    u <- matrix(u, ncol = 2)
    n <- nrow(u)
    fam <- copula@family
    
    .C("Hfunc1",
       as.integer(fam),
       as.integer(n),
       as.double(u[, 2]),
       as.double(u[, 1]),
       as.double(param[1]),
       as.double(param[2]),
       as.double(rep(0, n)),
       PACKAGE = "VineCopula")[[7]]
}

# ddv
linkVineCop.ddv <- function(u, copula) {
    param <- copula@parameters
    
    if (length(param) == 1) param <- c(param, 0)
    
    u <- matrix(u, ncol = 2)
    n <- nrow(u)
    fam <- copula@family
    
    .C("Hfunc2",
       as.integer(fam),
       as.integer(n),
       as.double(u[, 1]),
       as.double(u[, 2]),
       as.double(param[1]),
       as.double(param[2]),
       as.double(rep(0, n)),
       PACKAGE = "VineCopula")[[7]]
}


## random numbers from VineCopulaSim
linkVineCop.r <- function(n, copula) {
    param <- copula@parameters
    
    if (length(param) == 1) param <- c(param, 0)
    
    fam <- copula@family
    if (is.na(param[2])) param <- c(param, 0)
    
    res <- .C("pcc",
              as.integer(n),
              as.integer(2),
              as.integer(fam),
              as.integer(1), 
              as.double(param[1]),
              as.double(param[2]),
              as.double(rep(0, n * 2)),
              PACKAGE = "VineCopula")[[7]]
    
    return(matrix(res, ncol = 2))
}

## Kendall's tau
linkVineCop.tau <- function(copula) {
    param <- copula@parameters
    if (length(param) == 1) param <- c(param, 0)
    
    BiCopPar2Tau(copula@family, param[1], param[2])
}

## get parameter from Kendall's tau (only for one parameter families)
linkVineCop.iTau <- function(copula, tau) {
    BiCopTau2Par(copula@family, tau)
}

## tailIndex
linkVineCop.tailIndex <- function(copula) {
    param <- copula@parameters
    if (length(param) == 1) param <- c(param, 0)
    
    unlist(BiCopPar2TailDep(copula@family, param[1], param[2]))
}

setGeneric("dduCopula", function(u, copula, ...) standardGeneric("dduCopula"))
setGeneric("ddvCopula", function(u, copula, ...) standardGeneric("ddvCopula"))