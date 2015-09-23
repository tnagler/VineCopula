####################
##                ##
## the Joe copula ##
##                ##
####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validJoeBiCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("joeBiCopula",
  representation = representation("copula", family="numeric"),
  validity = validJoeBiCopula,
  contains = list("copula")
)

# constructor
joeBiCopula <- function (param=2) {
  if (any(is.na(param) | param >= Inf | param <= 1 ))
    stop("Parameter is outside of the allowed interval (1,Inf).")
  new("joeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 1, param.upbnd = Inf, family=6, 
      fullname = "Joe copula family. Number 6 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","joeBiCopula"), 
          function(u, copula, log, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","joeBiCopula"), function(u, copula, log, ...) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","joeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","joeBiCopula"), linkVineCop.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","joeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","joeBiCopula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","joeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","joeBiCopula"), linkVineCop.ddv)

## random number generater
setMethod("rCopula", signature("numeric","joeBiCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("joeBiCopula"), 
          function(copula, tau) {
            if(tau <= 0) 
              return(NA)
            linkVineCop.iTau(copula, max(1e-6,abs(tau)))
          })

setMethod("tau",signature("joeBiCopula"),linkVineCop.tau)
setMethod("tailIndex",signature("joeBiCopula"),linkVineCop.tailIndex)


#########################
## Joe survival copula ##
#########################

setClass("surJoeBiCopula",
  representation = representation("copula", family="numeric"),
  validity = validJoeBiCopula,
  contains = list("copula")
)

# constructor
surJoeBiCopula <- function (param=2) {
  if (any(is.na(param) | param >= Inf | param <= 1 ))
    stop("Parameter is outside of the allowed interval (1,Inf).")
  new("surJoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = 1, param.upbnd = Inf, family=16, 
      fullname = "Survival Joe copula family. Number 16 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surJoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension), copula, ...)
          })
setMethod("dCopula", signature("matrix","surJoeBiCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surJoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surJoeBiCopula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surJoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surJoeBiCopula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surJoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surJoeBiCopula"), linkVineCop.ddv)

## random number generater
setMethod("rCopula", signature("numeric","surJoeBiCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("surJoeBiCopula"), 
          function(copula, tau) {
            if(tau <= 0) 
              return(NA)
            linkVineCop.iTau(copula, max(1e-6,abs(tau)))
          })

setMethod("tau",signature("surJoeBiCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("surJoeBiCopula"),linkVineCop.tailIndex)

###################
## Joe copula 90 ##
###################

validRotJoeBiCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param <= lower))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("r90JoeBiCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeBiCopula,
  contains = list("copula")
)

# constructor
r90JoeBiCopula <- function (param=-2) {
  if (any(is.na(param) | param >= -1 | param <= -Inf ))
    stop("Parameter is outside of the allowed interval (-Inf,-1).")
  new("r90JoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = -1, family=26, 
      fullname = "90 deg rotated Joe copula family. Number 26 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","r90JoeBiCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90JoeBiCopula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90JoeBiCopula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90JoeBiCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90JoeBiCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r90JoeBiCopula"),
          function(copula, tau) {
            if(tau >= 0) 
              return(NA)
            linkVineCop.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r90JoeBiCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("r90JoeBiCopula"),linkVineCop.tailIndex)

####################
## Joe copula 270 ##
####################

setClass("r270JoeBiCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeBiCopula,
  contains = list("copula")
)

# constructor
r270JoeBiCopula <- function (param=-2) {
  if (any(is.na(param) | param >= -1 | param <= -Inf ))
    stop("Parameter is outside of the allowed interval (-Inf,-1).")
  new("r270JoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = -1, family=36, 
      fullname = "270 deg rotated Joe copula family. Number 36 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","r270JoeBiCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270JoeBiCopula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270JoeBiCopula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270JoeBiCopula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270JoeBiCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270JoeBiCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r270JoeBiCopula"), 
          function(copula, tau) {
            if(tau >= 0) 
              return(NA)
            linkVineCop.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r270JoeBiCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("r270JoeBiCopula"),linkVineCop.tailIndex)