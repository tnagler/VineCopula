#' Classes \code{"BB6Copula"}, \code{"surBB6Copula"}, \code{"r90BB6Copula"} and
#' \code{"r270BB6Copula"}
#'
#' Wrapper classes representing the BB6, survival BB6, 90 degree and 270 degree
#' rotated BB6 copula families (Joe 1997) from the
#' \code{\link{VineCopula-package}}.
#'
#'
#' @name BB6Copula-class
#' @aliases BB6Copula-class dduCopula,numeric,BB6Copula-method
#' ddvCopula,numeric,BB6Copula-method dduCopula,matrix,BB6Copula-method
#' ddvCopula,matrix,BB6Copula-method getKendallDistr,BB6Copula-method
#' kendallDistribution,BB6Copula-method surBB6Copula-class
#' dduCopula,numeric,surBB6Copula-method ddvCopula,numeric,surBB6Copula-method
#' dduCopula,matrix,surBB6Copula-method ddvCopula,matrix,surBB6Copula-method
#' r90BB6Copula-class dduCopula,numeric,r90BB6Copula-method
#' ddvCopula,numeric,r90BB6Copula-method dduCopula,matrix,r90BB6Copula-method
#' ddvCopula,matrix,r90BB6Copula-method r270BB6Copula-class
#' dduCopula,numeric,r270BB6Copula-method
#' ddvCopula,numeric,r270BB6Copula-method dduCopula,matrix,r270BB6Copula-method
#' ddvCopula,matrix,r270BB6Copula-method
#' @docType class
#' @section Objects from the Classes: Objects can be created by calls of the
#' form \code{new("BB6Copula", ...)}, \code{new("surBB6Copula", ...)},
#' \code{new("r90BB6Copula", ...)} and \code{new("r270BB6Copula", ...)} or by
#' the functions \code{\link{BB6Copula}}, \code{\link{surBB6Copula}},
#' \code{\link{r90BB6Copula}} and \code{\link{r270BB6Copula}}.
#' @author Benedikt Graeler
#' @seealso See also \code{\linkS4class{BB1Copula}},
#' \code{\linkS4class{BB7Copula}}, \code{\linkS4class{BB8Copula}} and
#' \code{\linkS4class{joeCopula}} for further wrapper classes to the
#' \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords classes
#' @examples
#'
#' showClass("BB6Copula")
#'
NULL


validBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
  param <- object@parameters
  p.n <- length(param)
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (p.n != length(upper))
      return("Parameter and upper bound have non-equal length.")
  if (p.n != length(lower))
      return("Parameter and lower bound have non-equal length.")
  if (p.n != length(object@param.names))
      return("Parameter and parameter names have non-equal length.")
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


#' Constructor of the BB6 Family and Rotated Versions thereof
#'
#' Constructs an object of the \code{\linkS4class{BB6Copula}} (survival
#' \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270})
#' family for given parameters.
#'
#'
#' @aliases BB6Copula surBB6Copula r90BB6Copula r270BB6Copula
#' @param param The parameter \code{param} defines the copula through
#' \code{theta} and \code{delta}.
#' @return One of the respective BB6 copula classes
#' (\code{\linkS4class{BB6Copula}}, \code{\linkS4class{surBB6Copula}},
#' \code{\linkS4class{r90BB6Copula}}, \code{\linkS4class{r270BB6Copula}}).
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
#' persp(BB6Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(surBB6Copula(c(1,1.5)), dCopula, zlim = c(0,10))
#' persp(r90BB6Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#' persp(r270BB6Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
#'
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
setMethod("lambda",signature("BB6Copula"),linkVineCop.tailIndex)

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
setMethod("lambda",signature("surBB6Copula"),linkVineCop.tailIndex)

#######################
## BB6 copula 90 deg ##
#######################

validRotBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
    param <- object@parameters
    p.n <- length(param)
    upper <- object@param.upbnd
    lower <- object@param.lowbnd
    if (p.n != length(upper))
        return("Parameter and upper bound have non-equal length.")
    if (p.n != length(lower))
        return("Parameter and lower bound have non-equal length.")
    if (p.n != length(object@param.names))
        return("Parameter and parameter names have non-equal length.")
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
setMethod("lambda",signature("r90BB6Copula"),linkVineCop.tailIndex)

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
setMethod("lambda",signature("r270BB6Copula"),linkVineCop.tailIndex)
