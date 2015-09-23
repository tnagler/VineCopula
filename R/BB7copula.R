#####################
##                 ##
## the BB7 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB7Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB7 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB7Copula,
  contains = list("copula")
)

# constructor
BB7Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf, Inf) | param[1] < 1 | param[2] <= 0))
    stop(paste("Parameter values out of bounds: theta: [1,Inf), delta: (0,Inf)."))
  new("BB7Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, Inf),
      family=9, fullname = "BB7 copula family. Number 9 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","BB7Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","BB7Copula"), function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","BB7Copula"), linkVineCop.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","BB7Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","BB7Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","BB7Copula"), linkVineCop.r)

setMethod("tau",signature("BB7Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("BB7Copula"),linkVineCop.tailIndex)


#########################
## BB7 survival copula ##
#########################

setClass("surBB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB7Copula,
  contains = list("copula")
)

# constructor
surBB7Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf, Inf) | param[1] < 1 | param[2] <= 0))
    stop(paste("Parameter values out of bounds: theta: [1,Inf), delta: (0,Inf)."))
  new("surBB7Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, Inf),
      family= 19, fullname = "Survival BB7 copula family. Number 19 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surBB7Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension,),copula,log=log)
          })
setMethod("dCopula", signature("matrix","surBB7Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surBB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surBB7Copula"), linkVineCop.surCDF)
  
## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surBB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surBB7Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surBB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surBB7Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surBB7Copula"), linkVineCop.r)

setMethod("tau",signature("surBB7Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surBB7Copula"),linkVineCop.tailIndex)

###################
## BB7 copula 90 ##
###################

validRotBB7Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB7 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param[1] > upper[1] | param[2] >= upper[2] | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB7Copula,
  contains = list("copula")
)

# constructor
r90BB7Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param[1] > -1 | param[2] >= 0 | param <= c(-Inf,-Inf)))
    stop(paste("Parameter values out of bounds: theta: (-Inf,-1], delta: (-Inf,0)."))
  new("r90BB7Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, 0),
      family=29, fullname = "90 deg rotated BB7 copula family. Number 29 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90BB7Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90BB7Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90BB7Copula"), linkVineCop.r90CDF)
  
## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90BB7Copula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90BB7Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90BB7Copula"), linkVineCop.r)

setMethod("tau",signature("r90BB7Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90BB7Copula"),linkVineCop.tailIndex)

########################
## BB7 copula 270 deg ##
########################

setClass("r270BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB7Copula,
  contains = list("copula")
)

# constructor
r270BB7Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param[1] > -1 | param[2] >= 0 | param <= c(-Inf,-Inf)))
    stop(paste("Parameter values out of bounds: theta: (-Inf,-1], delta: (-Inf,0)."))
  new("r270BB7Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, 0), 
      family=39, fullname = "270 deg rotated BB7 copula family. Number 39 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270BB7Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270BB7Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270BB7Copula"), linkVineCop.r270CDF)
  
## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270BB7Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270BB7Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270BB7Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270BB7Copula"), linkVineCop.r)

setMethod("tau",signature("r270BB7Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270BB7Copula"),linkVineCop.tailIndex)