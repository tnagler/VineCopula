#' Classes \code{"BB1Copula"}, \code{"surBB1Copula"}, \code{"r90BB1Copula"} and
#' \code{"r270BB1Copula"}
#'
#' Wrapper classes representing the BB1, survival BB1, 90 degree and 270 degree
#' rotated BB1 copula families (Joe 1997) from
#' \code{\link{VineCopula-package}}.
#'
#'
#' @name BB1Copula-class
#' @aliases BB1Copula-class dduCopula,numeric,BB1Copula-method
#' ddvCopula,numeric,BB1Copula-method dduCopula,matrix,BB1Copula-method
#' ddvCopula,matrix,BB1Copula-method getKendallDistr,BB1Copula-method
#' kendallDistribution,BB1Copula-method surBB1Copula-class
#' dduCopula,numeric,surBB1Copula-method ddvCopula,numeric,surBB1Copula-method
#' dduCopula,matrix,surBB1Copula-method ddvCopula,matrix,surBB1Copula-method
#' r90BB1Copula-class dduCopula,numeric,r90BB1Copula-method
#' ddvCopula,numeric,r90BB1Copula-method dduCopula,matrix,r90BB1Copula-method
#' ddvCopula,matrix,r90BB1Copula-method r270BB1Copula-class
#' dduCopula,numeric,r270BB1Copula-method
#' ddvCopula,numeric,r270BB1Copula-method dduCopula,matrix,r270BB1Copula-method
#' ddvCopula,matrix,r270BB1Copula-method
#' @docType class
#' @section Objects from the Classes: Objects can be created by calls of the
#' form \code{new("BB1Copula", ...)}, \code{new("surBB1Copula", ...)},
#' \code{new("r90BB1Copula", ...)} and \code{new("r270BB1Copula", ...)} or by
#' the functions \code{\link{BB1Copula}}, \code{\link{surBB1Copula}},
#' \code{\link{r90BB1Copula}} and \code{\link{r270BB1Copula}}.
#' @author Benedikt Graeler
#' @seealso See also \code{\linkS4class{BB6Copula}},
#' \code{\linkS4class{BB7Copula}}, \code{\linkS4class{BB8Copula}} and
#' \code{\linkS4class{joeCopula}} for further wrapper classes to the
#' \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords classes
#' @examples
#'
#' showClass("BB1Copula")
#'
NULL


validBB1Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB1 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB1Copula,
  contains = list("copula")
)

# constructor


#' Constructor of the BB1 Family and Rotated Versions thereof
#'
#' Constructs an object of the \code{\linkS4class{BB1Copula}} (survival
#' \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270})
#' family for given parameters.
#'
#'
#' @aliases BB1Copula surBB1Copula r90BB1Copula r270BB1Copula
#' @param param The parameter \code{param} defines the copula through
#' \code{theta} and \code{delta}.
#' @return One of the respective BB1 copula classes
#' (\code{\linkS4class{BB1Copula}}, \code{\linkS4class{surBB1Copula}},
#' \code{\linkS4class{r90BB1Copula}}, \code{\linkS4class{r270BB1Copula}}).
#' @author Benedikt Graeler
#' @seealso See also \code{\link{BB6Copula}}, \code{\link{BB7Copula}},
#' \code{\link{BB8Copula}} and \code{\link{joeCopula}} for further wrapper
#' functions to the \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords distribution copula
#' @examples
#'
#' library(copula)
#'
#' persp(BB1Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(surBB1Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(r90BB1Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#' persp(r270BB1Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#'
BB1Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf,Inf) | param[1] <= 0 | param[2] < 1))
    stop(paste("Parameter values out of bounds: theta: (0,Inf), delta: [1,Inf)."))
  new("BB1Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(0, 1), param.upbnd = c(Inf, Inf),
      family=7, fullname = "BB1 copula family. Number 7 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","BB1Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","BB1Copula"),
          function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","BB1Copula"), linkVineCop.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","BB1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","BB1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","BB1Copula"), linkVineCop.r)

setMethod("tau",signature("BB1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("BB1Copula"),linkVineCop.tailIndex)

#########################
## BB1 survival copula ##
#########################

setClass("surBB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB1Copula,
  contains = list("copula")
)

# constructor
surBB1Copula <- function (param=c(1,1)) {
  if (any(is.na(param) | param >= c(Inf,Inf) | param[1] <= 0 | param[2] < 1))
    stop(paste("Parameter values out of bounds: theta: (0,Inf), delta: [1,Inf)."))
  new("surBB1Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(0, 1), param.upbnd = c(Inf, Inf),
      family=17, fullname = "Survival BB1 copula family. Number 17 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surBB1Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surBB1Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surBB1Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surBB1Copula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surBB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surBB1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surBB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surBB1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surBB1Copula"), linkVineCop.r)

setMethod("tau",signature("surBB1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surBB1Copula"),linkVineCop.tailIndex)

#######################
## BB1 copula 90 deg ##
#######################

validRotBB1Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB1 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("r90BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB1Copula,
  contains = list("copula")
)

# constructor
r90BB1Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param[1] >= 0 | param[2] > -1 | param <= c(-Inf,-Inf)))
    stop(paste("Parameter values out of bounds: theta: (-Inf,0), delta: (-Inf,-1]."))
  new("r90BB1Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(0, -1),
      family=27, fullname = "90 deg rotated BB1 copula family. Number 27 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90BB1Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90BB1Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90BB1Copula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90BB1Copula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90BB1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90BB1Copula"), linkVineCop.r)

setMethod("tau",signature("r90BB1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90BB1Copula"),linkVineCop.tailIndex)

########################
## BB1 copula 270 deg ##
########################

setClass("r270BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB1Copula,
  contains = list("copula")
)

# constructor
r270BB1Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param[1] >= 0 | param[2] > -1 | param <= c(-Inf,-Inf)))
    stop(paste("Parameter values out of bounds: theta: (-Inf,0), delta: (-Inf,-1]."))
  new("r270BB1Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(0, -1),
      family=37, fullname = "270 deg rotated BB1 copula family. Number 37 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r270BB1Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270BB1Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270BB1Copula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270BB1Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270BB1Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270BB1Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270BB1Copula"), linkVineCop.r)

setMethod("tau",signature("r270BB1Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270BB1Copula"),linkVineCop.tailIndex)
