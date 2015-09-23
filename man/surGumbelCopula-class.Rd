\name{surGumbelCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{surGumbelCopula-class}
\alias{dduCopula,matrix,surGumbelCopula-method}
\alias{dduCopula,numeric,surGumbelCopula-method}
\alias{ddvCopula,matrix,surGumbelCopula-method}
\alias{ddvCopula,numeric,surGumbelCopula-method}

\alias{r90GumbelCopula-class}
\alias{dduCopula,matrix,r90GumbelCopula-method}
\alias{dduCopula,numeric,r90GumbelCopula-method}
\alias{ddvCopula,matrix,r90GumbelCopula-method}
\alias{ddvCopula,numeric,r90GumbelCopula-method}

\alias{r270GumbelCopula-class}
\alias{dduCopula,matrix,r270GumbelCopula-method}
\alias{dduCopula,numeric,r270GumbelCopula-method}
\alias{ddvCopula,matrix,r270GumbelCopula-method}
\alias{ddvCopula,numeric,r270GumbelCopula-method}

\title{Classes \code{"surGumbelCopula"}, \code{"r90GumbelCopula"} and \code{"r270GumbelCopula"}}
\description{
A class representing rotated versions of the Gumbel copula family (survival, 90 and 270 degree rotated).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("surGumbelCopula", ...)}, \code{new("r90GumbelCopula", ...)} and \code{new("r270GumbelCopula", ...)} or by the function \code{\link{surGumbelCopula}}, \code{\link{r90GumbelCopula}} and \code{\link{r270GumbelCopula}} respectively.
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
    \item{dduCopula}{\code{signature(u = "matrix", copula = "surGumbelCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "surGumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "surGumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "surGumbelCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "matrix", copula = "r90GumbelCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "r90GumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "r90GumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "r90GumbelCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "matrix", copula = "r270GumbelCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "r270GumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "r270GumbelCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "r270GumbelCopula")}: ... }
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

persp(surGumbelCopula(1.5),dCopula,zlim=c(0,10))
persp(r90GumbelCopula(-1.5),dCopula,zlim=c(0,10))
persp(r270GumbelCopula(-1.5),dCopula,zlim=c(0,10))
}
\keyword{classes}
