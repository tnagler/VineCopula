\name{surClaytonCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{surClaytonCopula-class}
\alias{dduCopula,matrix,surClaytonCopula-method}
\alias{dduCopula,numeric,surClaytonCopula-method}
\alias{ddvCopula,matrix,surClaytonCopula-method}
\alias{ddvCopula,numeric,surClaytonCopula-method}

\alias{r90ClaytonCopula-class}
\alias{dduCopula,matrix,r90ClaytonCopula-method}
\alias{dduCopula,numeric,r90ClaytonCopula-method}
\alias{ddvCopula,matrix,r90ClaytonCopula-method}
\alias{ddvCopula,numeric,r90ClaytonCopula-method}

\alias{r270ClaytonCopula-class}
\alias{dduCopula,matrix,r270ClaytonCopula-method}
\alias{dduCopula,numeric,r270ClaytonCopula-method}
\alias{ddvCopula,matrix,r270ClaytonCopula-method}
\alias{ddvCopula,numeric,r270ClaytonCopula-method}

\title{Classes \code{"surClaytonCopula"}, \code{"r90ClaytonCopula"} and \code{"r270ClaytonCopula"}}
\description{
A class representing rotated versions of the Clayton copula family (survival, 90 and 270 degree rotated).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("surClaytonCopula", ...)}, \code{new("r90ClaytonCopula", ...)} and \code{new("r270ClaytonCopula", ...)} or by the function \code{\link{surClaytonCopula}}, \code{\link{r90ClaytonCopula}} and \code{\link{r270ClaytonCopula}} respectively.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"numeric"} The family number in \code{\link{VineCopula-package}}}
    \item{\code{dimension}:}{Object of class \code{"integer"} The dimension of the copula (2). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} The parameter }
    \item{\code{param.names}:}{Object of class \code{"character"} name of the parameter }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} lower bound of the parameter }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} upper bound of the parameter }
    \item{\code{fullname}:}{Object of class \code{"character"} descriptive name of the family }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "surClaytonCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "surClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "surClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "surClaytonCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "matrix", copula = "r90ClaytonCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "r90ClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "r90ClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "r90ClaytonCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "matrix", copula = "r270ClaytonCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "r270ClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "r270ClaytonCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "r270ClaytonCopula")}: ... }
	 }
}
\author{
Benedikt Graeler
}
\seealso{
\code{\link{VineCopula-package}}
}
\examples{
library(copula)

persp(surClaytonCopula(.5),dCopula,zlim=c(0,10))
persp(r90ClaytonCopula(-.5),dCopula,zlim=c(0,10))
persp(r270ClaytonCopula(-.5),dCopula,zlim=c(0,10))
}
\keyword{classes}
