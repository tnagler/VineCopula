
#' Classes \code{"BB7Copula"}, \code{"surBB7Copula"}, \code{"r90BB7Copula"} and
#' \code{"r270BB7Copula"}
#'
#' Wrapper classes representing the BB7, survival BB7, 90 degree and 270 degree
#' rotated BB7 copula families (Joe 1997) from the
#' \code{\link{VineCopula-package}} package.
#'
#'
#' @name BB7Copula-class
#' @aliases BB7Copula-class dduCopula,numeric,BB7Copula-method
#' ddvCopula,numeric,BB7Copula-method dduCopula,matrix,BB7Copula-method
#' ddvCopula,matrix,BB7Copula-method getKendallDistr,BB7Copula-method
#' kendallDistribution,BB7Copula-method surBB7Copula-class
#' dduCopula,numeric,surBB7Copula-method ddvCopula,numeric,surBB7Copula-method
#' dduCopula,matrix,surBB7Copula-method ddvCopula,matrix,surBB7Copula-method
#' r90BB7Copula-class dduCopula,numeric,r90BB7Copula-method
#' ddvCopula,numeric,r90BB7Copula-method dduCopula,matrix,r90BB7Copula-method
#' ddvCopula,matrix,r90BB7Copula-method r270BB7Copula-class
#' dduCopula,numeric,r270BB7Copula-method
#' ddvCopula,numeric,r270BB7Copula-method dduCopula,matrix,r270BB7Copula-method
#' ddvCopula,matrix,r270BB7Copula-method
#' @docType class
#' @section Objects from the Classes: Objects can be created by calls of the
#' form \code{new("BB7Copula", ...)}, \code{new("surBB7Copula", ...)},
#' \code{new("r90BB7Copula", ...)} and \code{new("r270BB7Copula", ...)} or by
#' the functions \code{\link{BB7Copula}}, \code{\link{surBB7Copula}},
#' \code{\link{r90BB7Copula}} and \code{\link{r270BB7Copula}}.
#' @author Benedikt Graeler
#' @seealso See also \code{\linkS4class{BB1Copula}},
#' \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB8Copula}} and
#' \code{\linkS4class{joeCopula}} for further wrapper classes to the
#' \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords classes
#' @examples
#'
#' showClass("BB7Copula")
#'
NULL



validBB7Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB7 copulas of dimension 2 are supported.")
    p.n <- length(object@parameters)
    if (p.n != length(object@param.upbnd))
        return("Parameter and upper bound have non-equal length.")
    if (p.n != length(object@param.lowbnd))
        return("Parameter and lower bound have non-equal length.")
    if (p.n != length(object@param.names))
        return("Parameter and parameter names have non-equal length.")
  else return (TRUE)
}

setClass("BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB7Copula,
  contains = list("copula")
)

# constructor


#' Constructor of the BB7 Family and Rotated Versions thereof
#'
#' Constructs an object of the \code{\linkS4class{BB7Copula}} (survival
#' \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270})
#' family for given parameters.
#'
#'
#' @aliases BB7Copula surBB7Copula r90BB7Copula r270BB7Copula
#' @param param The parameter \code{param} defines the copula through
#' \code{theta} and \code{delta}.
#' @return One of the respective BB7 copula classes
#' (\code{\linkS4class{BB7Copula}}, \code{\linkS4class{surBB7Copula}},
#' \code{\linkS4class{r90BB7Copula}}, \code{\linkS4class{r270BB7Copula}}).
#' @author Benedikt Graeler
#' @seealso See also \code{\link{BB6Copula}}, \code{\link{BB7Copula}},
#' \code{\link{BB8Copula}} and \code{\link{joeCopula}} for further wrapper
#' functions to the \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @examples
#'
#' library(copula)
#'
#' persp(BB7Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(surBB7Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(r90BB7Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#' persp(r270BB7Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#'
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
setMethod("lambda",signature("BB7Copula"),linkVineCop.tailIndex)


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
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula,log=log)
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
setMethod("lambda",signature("surBB7Copula"),linkVineCop.tailIndex)

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
    return("Parameter and upper bound have non-equal length.")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length.")
  if (length(param) != length(object@param.names))
      return("Parameter and parameter names have non-equal length.")
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
setMethod("lambda",signature("r90BB7Copula"),linkVineCop.tailIndex)

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
setMethod("lambda",signature("r270BB7Copula"),linkVineCop.tailIndex)
