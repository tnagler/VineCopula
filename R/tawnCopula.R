#########################
##                     ##
## Tawn type 1 copulas ##
##                     ##
#########################

validTawnCopula = function(object) {
  if (object@dimension != 2)
    return("Only Tawn copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("tawnT1Copula",
  representation = representation("copula", family="numeric"),
  validity = validTawnCopula,
  contains = list("copula")
)

# constructor
tawnT1Copula <- function (param=c(2,0.5)) {
  if (any(is.na(param) | param < c(1,0) | param[1] == Inf | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("tawnT1Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=104, fullname = "Tawn type 1 copula family. Number 104 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","tawnT1Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","tawnT1Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","tawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","tawnT1Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","tawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","tawnT1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","tawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","tawnT1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","tawnT1Copula"), linkVineCop.r)

setMethod("tau",signature("tawnT1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("tawnT1Copula"),linkVineCop.tailIndex)

# Pickand's A
# c-code: Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out)
setMethod("A", signature("tawnT1Copula"), function(copula, w) {
  .C("Tawn2",as.double(w), as.integer(length(w)), 
     as.double(copula@parameters[1]), as.double(copula@parameters[2]),
     as.double(1), as.double(rep(0,length(w))), PACKAGE = "VineCopula")[[6]]
})

#################################
## Tawn type 1 survival copula ##
#################################

setClass("surTawnT1Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
surTawnT1Copula <- function (param=c(2,0.5)) {
  if (any(is.na(param) | param < c(1,0) | param[1] == Inf | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("surTawnT1Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=114, fullname = "Survival Tawn type 1 copula family. Number 114 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surTawnT1Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surTawnT1Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","surTawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surTawnT1Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surTawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surTawnT1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surTawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surTawnT1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surTawnT1Copula"), linkVineCop.r)

setMethod("tau",signature("surTawnT1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surTawnT1Copula"),linkVineCop.tailIndex)

# Pickand's A
# c-code: Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out)
setMethod("A", signature("surTawnT1Copula"), function(copula, w) {
  u <- -expm1(-1+w)
  v <- -expm1(-w)
  
  surA <- .C("Tawn2",as.double(log(v)/log(u*v)), as.integer(length(w)), 
             as.double(copula@parameters[1]), as.double(copula@parameters[2]),
             as.double(1), as.double(rep(0,length(w))), PACKAGE = "VineCopula")[[6]]
  -log(1-u + 1-v - 1 + (u*v)^surA)
  # -log(1-u + 1-v - 1 + VineCopula:::linkVineCop.CDFtawn(cbind(u,v), tawnT1Copula(copula@parameters)))
})

#######################################
## Tawn type 1 90 deg. rotate copula ##
#######################################

setClass("r90TawnT1Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
r90TawnT1Copula <- function (param=c(-2, 0.5)) {
  if (any(is.na(param) | param[1] == -Inf | param[1] > -1 | param[2] < 0 | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("r90TawnT1Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(-Inf, 0), param.upbnd = c(-1, 1),
      family=124, fullname = "90 deg rotated Tawn type 1 copula family. Number 124 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90TawnT1Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90TawnT1Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","r90TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90TawnT1Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90TawnT1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r90TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90TawnT1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90TawnT1Copula"), linkVineCop.r)

setMethod("tau",signature("r90TawnT1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90TawnT1Copula"),linkVineCop.tailIndex)

########################################
## Tawn type 1 270 deg. rotate copula ##
########################################

setClass("r270TawnT1Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
r270TawnT1Copula <- function (param=c(-2, 0.5)) {
  if (any(is.na(param) | param[1] == -Inf | param[1] > -1 | param[2] < 0 | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("r270TawnT1Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(-Inf, 0), param.upbnd = c(-1, 1),
      family=134, fullname = "270 deg rotated Tawn type 1 copula family. Number 134 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270TawnT1Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270TawnT1Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","r270TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270TawnT1Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270TawnT1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270TawnT1Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270TawnT1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270TawnT1Copula"), linkVineCop.r)

setMethod("tau",signature("r270TawnT1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270TawnT1Copula"),linkVineCop.tailIndex)

#########################
##                     ##
## Tawn type 2 copulas ##
##                     ##
#########################

setClass("tawnT2Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
tawnT2Copula <- function (param=c(2,0.5)) {
  if (any(is.na(param) | param < c(1,0) | param[1] == Inf | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("tawnT2Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=204, fullname = "Tawn type 2 copula family. Number 204 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","tawnT2Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","tawnT2Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","tawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","tawnT2Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","tawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","tawnT2Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","tawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","tawnT2Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","tawnT2Copula"), linkVineCop.r)

setMethod("tau",signature("tawnT2Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("tawnT2Copula"),linkVineCop.tailIndex)

# Pickand's A
# c-code: Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out)
setMethod("A", signature("tawnT2Copula"), function(copula, w) {
  .C("Tawn2",as.double(w), as.integer(length(w)), 
     as.double(copula@parameters[1]), as.double(1),
     as.double(copula@parameters[2]),
     as.double(rep(0,length(w))), PACKAGE = "VineCopula")[[6]]
})

#################################
## Tawn type 2 survival copula ##
#################################

setClass("surTawnT2Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
surTawnT2Copula <- function (param=c(2,0.5)) {
  if (any(is.na(param) | param < c(1,0) | param[1] == Inf | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("surTawnT2Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=214, fullname = "Survival Tawn type 2 copula family. Number 214 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surTawnT2Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surTawnT2Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","surTawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surTawnT2Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surTawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surTawnT2Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surTawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surTawnT2Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surTawnT2Copula"), linkVineCop.r)

setMethod("tau",signature("surTawnT2Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surTawnT2Copula"),linkVineCop.tailIndex)

# Pickand's A
# c-code: Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out)
setMethod("A", signature("surTawnT2Copula"), function(copula, w) {
  u <- -expm1(-1+w) # 1-u
  v <- -expm1(-w)   # 1-v
  
  surA <- .C("Tawn2",as.double(log(v)/log(u*v)), as.integer(length(w)), 
             as.double(copula@parameters[1]), as.double(1),
             as.double(copula@parameters[2]), 
             as.double(rep(0,length(w))), PACKAGE = "VineCopula")[[6]]
  -log(1-u + 1-v - 1 + (u*v)^surA)
})

#######################################
## Tawn type 2 90 deg. rotate copula ##
#######################################

setClass("r90TawnT2Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
r90TawnT2Copula <- function (param=c(-2, 0.5)) {
  if (any(is.na(param) | param[1] == -Inf | param[1] > -1 | param[2] < 0 | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("r90TawnT2Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(-Inf, 0), param.upbnd = c(-1, 1),
      family=224, fullname = "90 deg rotated Tawn type 2 copula family. Number 224 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90TawnT2Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90TawnT2Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","r90TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90TawnT2Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90TawnT2Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r90TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90TawnT2Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90TawnT2Copula"), linkVineCop.r)

setMethod("tau",signature("r90TawnT2Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90TawnT2Copula"),linkVineCop.tailIndex)

########################################
## Tawn type 2 270 deg. rotate copula ##
########################################

setClass("r270TawnT2Copula",
         representation = representation("copula", family="numeric"),
         validity = validTawnCopula,
         contains = list("copula")
)

# constructor
r270TawnT2Copula <- function (param=c(-2, 0.5)) {
  if (any(is.na(param) | param[1] == -Inf | param[1] > -1 | param[2] < 0 | param[2] > 1))
    stop(paste("Parameter values out of bounds: param1: [1,Inf), param2: [0,1]."))
  new("r270TawnT2Copula", dimension = as.integer(2), parameters = param, 
      param.names = c("param1", "param2"), param.lowbnd = c(-Inf, 0), param.upbnd = c(-1, 1),
      family=234, fullname = "270 deg rotated Tawn type 2 copula family. Number 234 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270TawnT2Copula"), 
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270TawnT2Copula"), 
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","r270TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.CDFtawn(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270TawnT2Copula"), linkVineCop.CDFtawn)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270TawnT2Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270TawnT2Copula"), 
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270TawnT2Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270TawnT2Copula"), linkVineCop.r)

setMethod("tau",signature("r270TawnT2Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270TawnT2Copula"),linkVineCop.tailIndex)
