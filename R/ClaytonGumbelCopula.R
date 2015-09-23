#####################################
##                                 ##
## additions to the Clayton copula ##
##                                 ##
#####################################

#############################
## Clayton survival copula ##
#############################

validClaytonCopula = function(object) {
  if (object@dimension != 2)
    return("Only Clayton copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param <= lower ))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("surClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validClaytonCopula,
  contains = list("copula")
)

# constructor
surClaytonCopula <- function (param=1) {
  new("surClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 0, param.upbnd = Inf, family=13, 
      fullname = "Survival Clayton copula family. Number 13 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surClaytonCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula) {
            linkVineCop.surCDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surClaytonCopula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surClaytonCopula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surClaytonCopula"), linkVineCop.ddv)

## random number generater ??
setMethod("rCopula", signature("numeric","surClaytonCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("surClaytonCopula"), 
          function(copula, tau) {
            if(tau <= 0) 
              return(NA)
            linkVineCop.iTau(copula, max(1e-6,abs(tau)))
          })

setMethod("tau",signature("surClaytonCopula"),linkVineCop.tau)
setMethod("tailIndex",signature("surClaytonCopula"),linkVineCop.tailIndex)

#######################
## Clayton copula 90 ##
#######################

validRotClaytonCopula = function(object) {
  if (object@dimension != 2)
    return("Only Clayton copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90ClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotClaytonCopula,
  contains = list("copula")
)

# constructor
r90ClaytonCopula <- function (param=-1) {
  new("r90ClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = 0, family=23, 
      fullname = "90 deg rotated Clayton copula family. Number 23 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90ClaytonCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.r90CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90ClaytonCopula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90ClaytonCopula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90ClaytonCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90ClaytonCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r90ClaytonCopula"),
          function(copula, tau) {
            if(tau >= 0) 
              return(NA)
            linkVineCop.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r90ClaytonCopula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90ClaytonCopula"),linkVineCop.tailIndex)

########################
## Clayton copula 270 ##
########################

setClass("r270ClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotClaytonCopula,
  contains = list("copula")
)

# constructor
r270ClaytonCopula <- function (param=-1) {
  new("r270ClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = 0, family=33, 
      fullname = "270 deg rotated Clayton copula family. Number 33 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270ClaytonCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.r270CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270ClaytonCopula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270ClaytonCopula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270ClaytonCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270ClaytonCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r270ClaytonCopula"), 
          function(copula, tau) {
            if(tau >= 0) 
              return(NA)
            linkVineCop.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r270ClaytonCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("r270ClaytonCopula"),linkVineCop.tailIndex)

####################################
##                                ##
## additions to the Gumbel copula ##
##                                ##
####################################

validGumbelCopula = function(object) {
  if (object@dimension != 2)
    return("Only Gumbel copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param < lower ))
    return("Parameter value out of bound.")
  else return (TRUE)
}

############################
## Gumbel survival copula ##
############################

setClass("surGumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validGumbelCopula,
         contains = list("copula")
         )

# constructor
surGumbelCopula <- function (param=1) {
  new("surGumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 1, param.upbnd = Inf, family=14, 
      fullname = "Survival Gumbel copula family. Number 14 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surGumbelCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkVineCop.surCDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surGumbelCopula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surGumbelCopula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surGumbelCopula"), linkVineCop.ddv)

## random number generater ??
setMethod("rCopula", signature("numeric","surGumbelCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("surGumbelCopula"), 
          function(copula, tau) {
            if(tau < 0) 
              return(NA)
            linkVineCop.iTau(copula, max(0,abs(tau)))
          })

setMethod("tau",signature("surGumbelCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("surGumbelCopula"),linkVineCop.tailIndex)

#######################
## Gumbel copula 90 ##
#######################

validRotGumbelCopula = function(object) {
  if (object@dimension != 2)
    return("Only Gumbel copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param > upper | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90GumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validRotGumbelCopula,
         contains = list("copula")
         )

# constructor
r90GumbelCopula <- function (param=-1) {
  new("r90GumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = -1, family=24, 
      fullname = "90 deg rotated Gumbel copula family. Number 24 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90GumbelCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkVineCop.r90CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90GumbelCopula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90GumbelCopula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90GumbelCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90GumbelCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r90GumbelCopula"),
          function(copula, tau) {
            if(tau > 0) 
              return(NA)
            linkVineCop.iTau(copula, min(0,-abs(tau)))
          })

setMethod("tau",signature("r90GumbelCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("r90GumbelCopula"),linkVineCop.tailIndex)

########################
## Gumbel copula 270 ##
########################

setClass("r270GumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validRotGumbelCopula,
         contains = list("copula")
         )

# constructor
r270GumbelCopula <- function (param=-1) {
  new("r270GumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = -1, family=34, 
      fullname = "270 deg rotated Gumbel copula family. Number 34 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270GumbelCopula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkVineCop.r270CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270GumbelCopula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270GumbelCopula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270GumbelCopula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270GumbelCopula"), linkVineCop.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r270GumbelCopula"), 
          function(copula, tau) {
            if(tau >= 0)
              return(NA)
            linkVineCop.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r270GumbelCopula"),linkVineCop.tau)

setMethod("tailIndex",signature("r270GumbelCopula"),linkVineCop.tailIndex)
