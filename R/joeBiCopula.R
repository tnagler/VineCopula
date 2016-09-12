
#' Classes \code{"joeBiCopula"}, \code{"surJoeBiCopula"},
#' \code{"r90JoeBiCopula"} and \code{"r270JoeBiCopula"}
#'
#' Wrapper classes representing the bivariate Joe, survival Joe, 90 degree and
#' 270 degree rotated Joe copula families (Joe 1997) from
#' \code{\link{VineCopula-package}}. Note that package
#' \code{\link{copula-package}} provides a class \code{\linkS4class{joeCopula}}
#' as well.
#'
#'
#' @name joeBiCopula-class
#' @aliases joeBiCopula-class dduCopula,numeric,joeBiCopula-method
#' ddvCopula,numeric,joeBiCopula-method dduCopula,matrix,joeBiCopula-method
#' ddvCopula,matrix,joeBiCopula-method getKendallDistr,joeBiCopula-method
#' kendallDistribution,joeBiCopula-method surJoeBiCopula-class
#' dduCopula,numeric,surJoeBiCopula-method
#' ddvCopula,numeric,surJoeBiCopula-method
#' dduCopula,matrix,surJoeBiCopula-method
#' ddvCopula,matrix,surJoeBiCopula-method r90JoeBiCopula-class
#' dduCopula,numeric,r90JoeBiCopula-method
#' ddvCopula,numeric,r90JoeBiCopula-method
#' dduCopula,matrix,r90JoeBiCopula-method
#' ddvCopula,matrix,r90JoeBiCopula-method r270JoeBiCopula-class
#' dduCopula,numeric,r270JoeBiCopula-method
#' ddvCopula,numeric,r270JoeBiCopula-method
#' dduCopula,matrix,r270JoeBiCopula-method
#' ddvCopula,matrix,r270JoeBiCopula-method
#' @docType class
#' @section Objects from the Classes: Objects can be created by calls of the
#' form \code{new("joeBiCopula", ...)}, \code{new("surJoeBiCopula", ...)},
#' \code{new("r90JoeBiCopula", ...)} and \code{new("r270JoeBiCopula", ...)} or
#' by the functions \code{\link{joeBiCopula}}, \code{\link{surJoeBiCopula}},
#' \code{\link{r90JoeBiCopula}} and \code{\link{r270JoeBiCopula}}.
#' @author Benedikt Graeler
#' @seealso See also \code{\linkS4class{BB1Copula}},
#' \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}} and
#' \code{\linkS4class{BB8Copula}} for further wrapper classes to the
#' \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords classes
#' @examples
#'
#' showClass("surJoeBiCopula")
#'
NULL


validJoeBiCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
    p.n <- length(object@parameters)
    if (p.n != length(object@param.upbnd))
        return("Parameter and upper bound have non-equal length.")
    if (p.n != length(object@param.lowbnd))
        return("Parameter and lower bound have non-equal length.")
    if (p.n != length(object@param.names))
        return("Parameter and parameter names have non-equal length.")
  else return (TRUE)
}

setClass("joeBiCopula",
  representation = representation("copula", family="numeric"),
  validity = validJoeBiCopula,
  contains = list("copula")
)

# constructor


#' Constructor of the Joe Family and Rotated Versions thereof
#'
#' Constructs an object of the (survival \code{surJoeBiCopula}, 90 degree
#' rotated \code{r90JoeBiCopula} and 270 degree rotated \code{r270JoeBiCopula})
#' family for a given parameter. Note that package \code{\link{copula-package}}
#' provides a class \code{\linkS4class{joeCopula}} as well.
#'
#'
#' @aliases joeBiCopula surJoeBiCopula r90JoeBiCopula r270JoeBiCopula
#' @param param The parameter \code{param} defines the copula through
#' \code{theta}.
#' @return One of the respective Joe copula classes
#' (\code{\linkS4class{joeBiCopula}}, \code{\linkS4class{surJoeBiCopula}},
#' \code{\linkS4class{r90JoeBiCopula}}, \code{\linkS4class{r270JoeBiCopula}}).
#' @author Benedikt Graeler
#' @seealso See also \code{\link{BB1Copula}}, \code{\link{BB6Copula}},
#' \code{\link{BB7Copula}} and \code{\link{BB8Copula}} for further wrapper
#' functions to the \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @examples
#'
#' library(copula)
#'
#' persp(surJoeBiCopula(1.5), dCopula, zlim = c(0,10))
#' persp(r90JoeBiCopula(-1.5), dCopula, zlim = c(0,10))
#' persp(r270JoeBiCopula(-1.5), dCopula, zlim = c(0,10))
#'
#' @export joeBiCopula
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
setMethod("lambda",signature("joeBiCopula"),linkVineCop.tailIndex)


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
  new("surJoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
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

setMethod("lambda",signature("surJoeBiCopula"),linkVineCop.tailIndex)

###################
## Joe copula 90 ##
###################

validRotJoeBiCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
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
  new("r90JoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
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

setMethod("lambda",signature("r90JoeBiCopula"),linkVineCop.tailIndex)

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
  new("r270JoeBiCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
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

setMethod("lambda",signature("r270JoeBiCopula"),linkVineCop.tailIndex)
