#' Classes \code{"BB8Copula"}, \code{"surBB8Copula"}, \code{"r90BB8Copula"} and
#' \code{"r270BB8Copula"}
#'
#' Wrapper classes representing the BB8, survival BB8, 90 degree and 270 degree
#' rotated BB8 copula families (Joe 1997) from the
#' \code{\link{VineCopula-package}} package.
#'
#'
#' @name BB8Copula-class
#' @aliases BB8Copula-class dduCopula,numeric,BB8Copula-method
#' ddvCopula,numeric,BB8Copula-method dduCopula,matrix,BB8Copula-method
#' ddvCopula,matrix,BB8Copula-method getKendallDistr,BB8Copula-method
#' kendallDistribution,BB8Copula-method surBB8Copula-class
#' dduCopula,numeric,surBB8Copula-method ddvCopula,numeric,surBB8Copula-method
#' dduCopula,matrix,surBB8Copula-method ddvCopula,matrix,surBB8Copula-method
#' r90BB8Copula-class dduCopula,numeric,r90BB8Copula-method
#' ddvCopula,numeric,r90BB8Copula-method dduCopula,matrix,r90BB8Copula-method
#' ddvCopula,matrix,r90BB8Copula-method r270BB8Copula-class
#' dduCopula,numeric,r270BB8Copula-method
#' ddvCopula,numeric,r270BB8Copula-method dduCopula,matrix,r270BB8Copula-method
#' ddvCopula,matrix,r270BB8Copula-method fitCopula,twoParamBiCop-method
#' @docType class
#' @section Objects from the Classes: Objects can be created by calls of the
#' form \code{new("BB8Copula", ...)}, \code{new("surBB8Copula", ...)},
#' \code{new("r90BB8Copula", ...)} and \code{new("r270BB8Copula", ...)} or by
#' the functions \code{\link{BB8Copula}}, \code{\link{surBB8Copula}},
#' \code{\link{r90BB8Copula}} and \code{\link{r270BB8Copula}}.
#' @author Benedikt Graeler
#' @seealso See also \code{\linkS4class{BB1Copula}},
#' \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}} and
#' \code{\linkS4class{joeCopula}} for further wrapper classes to the
#' \code{\link{VineCopula-package}}.
#' @references Joe, H., (1997). Multivariate Models and Dependence Concepts.
#' Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
#' @keywords classes
#' @examples
#'
#' showClass("BB8Copula")
#'
NULL


validBB8Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB8 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param)) | param[1] >= upper[1] | param[2] > upper[2] | param[1] < lower[1] | param[2] <= lower[2])
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("BB8Copula",
  representation = representation("copula",family="numeric"),
  validity = validBB8Copula,
  contains = list("copula")
)

# constructor


#' Constructor of the BB8 Family and Rotated Versions thereof
#'
#' Constructs an object of the \code{\linkS4class{BB8Copula}} (survival
#' \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270})
#' family for given parameters.
#'
#'
#' @aliases BB8Copula surBB8Copula r90BB8Copula r270BB8Copula
#' @param param The parameter \code{param} defines the copula through
#' \code{theta} and \code{delta}.
#' @return One of the respective BB8 copula classes
#' (\code{\linkS4class{BB8Copula}}, \code{\linkS4class{surBB8Copula}},
#' \code{\linkS4class{r90BB8Copula}}, \code{\linkS4class{r270BB8Copula}}).
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
#' persp(BB8Copula(c(2,0.9)), dCopula, zlim = c(0,10))
#' persp(surBB8Copula(c(2,0.9)), dCopula, zlim = c(0,10))
#' persp(r90BB8Copula(c(-2,-0.9)), dCopula, zlim = c(0,10))
#' persp(r270BB8Copula(c(-2,-0.9)), dCopula, zlim = c(0,10))
#'
BB8Copula <- function (param=c(1,1)) {
  if (any(is.na(param)) | param[1] >= Inf | param[2] > 1 | param[1] < 1 | param[2] <= 0)
    stop("Parameter value out of bound: theta: [1,Inf), delta: (0,1].")
  new("BB8Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=10, fullname = "BB8 copula family. Number 10 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","BB8Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","BB8Copula"), function(u, copula, log) linkVineCop.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","BB8Copula"), linkVineCop.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","BB8Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","BB8Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","BB8Copula"),linkVineCop.r)

setMethod("tau",signature("BB8Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("BB8Copula"),linkVineCop.tailIndex)

#########################
## BB8 survival copula ##
#########################

setClass("surBB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB8Copula,
  contains = list("copula")
)

# constructor
surBB8Copula <- function (param=c(1,1)) {
  if (any(is.na(param)) | param[1] >= Inf | param[2] > 1 | param[1] < 1 | param[2] <= 0)
    stop("Parameter value out of bound: theta: [1,Inf), delta: (0,1].")
  new("surBB8Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1),
      family=20, fullname = "Survival BB8 copula family. Number 20 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","surBB8Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surBB8Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surBB8Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surBB8Copula"), linkVineCop.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surBB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surBB8Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surBB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surBB8Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surBB8Copula"), linkVineCop.r)

setMethod("tau",signature("surBB8Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("surBB8Copula"),linkVineCop.tailIndex)

#######################
## BB8 copula 90 deg ##
#######################

validRotBB8Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB8 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  else return (TRUE)
}

setClass("r90BB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB8Copula,
  contains = list("copula")
)

# constructor
r90BB8Copula <- function (param=c(-1,-1)) {
  if (any(is.na(param) | param[1] > -1 | param[2] >= 0 | param[1] <= -Inf | param[2] < -1))
    stop("Parameter value out of bound: theta: (-Inf,-1], delta: [-1,0).")
  new("r90BB8Copula", dimension = as.integer(2), parameters = param,
      param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0),
      family=30, fullname = "90 deg rotated BB8 copula family. Number 30 in VineCopula.")
}

## density ##
setMethod("dCopula", signature("numeric","r90BB8Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90BB8Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90BB8Copula"), linkVineCop.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90BB8Copula"), linkVineCop.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90BB8Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90BB8Copula"), linkVineCop.r)

setMethod("tau",signature("r90BB8Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r90BB8Copula"),linkVineCop.tailIndex)

###########################
## BB8 copula 270 degree ##
###########################

setClass("r270BB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB8Copula,
  contains = list("copula")
)

# constructor
r270BB8Copula <- function (param=c(-1,-1)) {
  val <- new("r270BB8Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0), family=40, fullname = "270 deg rotated BB8 copula family. Number 40 in VineCopula.")
  val
}

## density ##
setMethod("dCopula", signature("numeric","r270BB8Copula"),
          function(u, copula, log) {
            linkVineCop.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270BB8Copula"), linkVineCop.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270BB8Copula"), linkVineCop.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270BB8Copula"), linkVineCop.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270BB8Copula"),
          function(u, copula, ...) {
            linkVineCop.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270BB8Copula"), linkVineCop.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270BB8Copula"), linkVineCop.r)

setMethod("tau",signature("r270BB8Copula"),linkVineCop.tau)
setMethod("tailIndex",signature("r270BB8Copula"),linkVineCop.tailIndex)

### set union

setClassUnion("twoParamBiCop",c("BB1Copula","BB6Copula","BB7Copula","BB8Copula",
                                "surBB1Copula","surBB6Copula","surBB7Copula","surBB8Copula",
                                "r90BB1Copula","r90BB6Copula","r90BB7Copula","r90BB8Copula",
                                "r270BB1Copula","r270BB6Copula","r270BB7Copula","r270BB8Copula"))

fitCopula.twoParamBiCop <- function(copula, data, method = "mpl",
                                    estimate.variance = FALSE) {
  stopifnot(method=="mpl")
  fit <- BiCopEst(data[,1], data[,2], copula@family, "mle",
                  se=estimate.variance)

  if(!estimate.variance) {
    fit$se <- NA
    fit$se2 <- NA
  }

  copFit <- copulaFromFamilyIndex(copula@family, fit$par, fit$par2)
  new("fitCopula", estimate = c(fit$par, fit$par2), var.est = cbind(fit$se, fit$se2),
      method = "maximum pseudo-likelihood via BiCopEst",
      loglik = sum(dCopula(data, copFit, log=T)),
      fitting.stats=list(convergence = as.integer(NA)), nsample = nrow(data),
      copula=copFit)
}

setMethod("fitCopula", signature("twoParamBiCop"), fitCopula.twoParamBiCop)
