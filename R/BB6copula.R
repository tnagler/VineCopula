#####################
##                 ##
## the BB6 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param < lower))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB6Copula,
  contains = list("copula")
)

# constructor
BB6Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf, Inf) | param < c(1,1)))
    stop("Parameter value(s) out of bound(s): theta: [1,Inf), delta: [1,Inf).")
  new("BB6Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(1, 1), param.upbnd = c(Inf, Inf),
      family=8, fullname = "BB6 copula family. Number 8 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","BB6Copula"), 
          function(u, copula, log, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log, ...)
          })
setMethod("dCopula", signature("matrix","BB6Copula"), function(u, copula, log, ...) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","BB6Copula"), linkVineCop.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","BB6Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","BB6Copula"), linkVineCop.ddv)

## random number generater ??
setMethod("rCopula", signature("numeric","BB6Copula"), linkVineCop.r)

setMethod("tau",signature("BB6Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("BB6Copula"),linkVineCop.tailIndex)

#########################
## BB6 survival copula ##
#########################

setClass("surBB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB6Copula,
  contains = list("copula")
)

# constructor
surBB6Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf, Inf) | param < c(1,1)))
    stop("Parameter value(s) out of bound(s): theta: [1,Inf), delta: [1,Inf).")
  new("surBB6Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(1, 1), param.upbnd = c(Inf, Inf), 
      family=18, fullname = "Survival BB6 copula family. Number 18 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surBB6Copula"), 
          function(u, copula, log, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surBB6Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surBB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surBB6Copula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surBB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surBB6Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surBB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surBB6Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surBB6Copula"), linkVineCop.r)

setMethod("tau",signature("surBB6Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surBB6Copula"),linkVineCop.tailIndex)

#######################
## BB6 copula 90 deg ##
#######################

validRotBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("r90BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB6Copula,
  contains = list("copula")
)

# constructor
r90BB6Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param > c(-1,-1) | param <= c(-Inf,-Inf)))
    stop("Parameter value out of bound: theta: (-Inf,1], delta: (-Inf,1].")
  new("r90BB6Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, -1),
      family=28, fullname = "90 deg rotated BB6 copula family. Number 28 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90BB6Copula"), 
          function(u, copula, log, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula,log)
          })
setMethod("dCopula", signature("matrix","r90BB6Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90BB6Copula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90BB6Copula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90BB6Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90BB6Copula"), linkVineCop.r)

setMethod("tau",signature("r90BB6Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90BB6Copula"),linkVineCop.tailIndex)

###########################
## BB6 copula 270 degree ##
###########################

setClass("r270BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB6Copula,
  contains = list("copula")
)

# constructor
r270BB6Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param > c(-1,-1) | param <= c(-Inf,-Inf)))
    stop("Parameter value out of bound: theta: (-Inf,1], delta: (-Inf,1].")
  new("r270BB6Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, -1),
      family=38, fullname = "270 deg rotated BB6 copula family. Number 38 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270BB6Copula"), 
          function(u, copula, log, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension, log),copula)
          })
setMethod("dCopula", signature("matrix","r270BB6Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270BB6Copula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270BB6Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270BB6Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270BB6Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270BB6Copula"), linkVineCop.r)

setMethod("tau",signature("r270BB6Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270BB6Copula"),linkVineCop.tailIndex)