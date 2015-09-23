\name{tawnT2Copula-class}
\Rdversion{1.1}
\docType{class}
\alias{tawnT2Copula-class}
\alias{dduCopula,matrix,tawnT2Copula-method}
\alias{dduCopula,numeric,tawnT2Copula-method}
\alias{ddvCopula,matrix,tawnT2Copula-method}
\alias{ddvCopula,numeric,tawnT2Copula-method}

\alias{surTawnT2Copula-class}
\alias{dduCopula,matrix,surTawnT2Copula-method}
\alias{dduCopula,numeric,surTawnT2Copula-method}
\alias{ddvCopula,matrix,surTawnT2Copula-method}
\alias{ddvCopula,numeric,surTawnT2Copula-method}

\alias{r90TawnT2Copula-class}
\alias{dduCopula,matrix,r90TawnT2Copula-method}
\alias{dduCopula,numeric,r90TawnT2Copula-method}
\alias{ddvCopula,matrix,r90TawnT2Copula-method}
\alias{ddvCopula,numeric,r90TawnT2Copula-method}

\alias{r270TawnT2Copula-class}
\alias{dduCopula,matrix,r270TawnT2Copula-method}
\alias{dduCopula,numeric,r270TawnT2Copula-method}
\alias{ddvCopula,matrix,r270TawnT2Copula-method}
\alias{ddvCopula,numeric,r270TawnT2Copula-method}

\title{Class \code{"tawnT2Copula"}}
\description{
S4-class representation of the Tawn Copula family of type 2 and rotated versions there of.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tawnT2Copula", ...)}, or through the explicit constructors \code{\link{tawnT2Copula}}, \code{\link{surTawnT2Copula}}, \code{\link{r90TawnT2Copula}} and \code{\link{r270TawnT2Copula}} respectively.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"numeric"} providing the unique number in VineCopula.}
    \item{\code{dimension}:}{Object of class \code{"integer"} and fixed to \code{2L}.}
    \item{\code{parameters}:}{Object of class \code{"numeric"} representing the two parameters.}
    \item{\code{param.names}:}{Object of class \code{"character"} providing the names of the parameters.}
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} providing the lower bounds of the parameters.}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} providing the upper bounds of the parameters.}
    \item{\code{fullname}:}{Object of class \code{"character"} providing a textual summary of the copula class.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(u = "matrix", copula = "tawnT2Copula")}: ... }
    \item{dCopula}{\code{signature(u = "numeric", copula = "tawnT2Copula")}: ... }
    \item{dduCopula}{\code{signature(u = "matrix", copula = "tawnT2Copula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "tawnT2Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "tawnT2Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "tawnT2Copula")}: ... }
    \item{pCopula}{\code{signature(u = "matrix", copula = "tawnT2Copula")}: ... }
    \item{pCopula}{\code{signature(u = "numeric", copula = "tawnT2Copula")}: ... }
    \item{rCopula}{\code{signature(n = "numeric", copula = "tawnT2Copula")}: ... }
    \item{tailIndex}{\code{signature(copula = "tawnT2Copula")}: ... }
    \item{tau}{\code{signature(copula = "tawnT2Copula")}: ... }
	 }
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{tawnT1Copula}} and the package \code{\link{VineCopula-package}} for implementation details.
}
\examples{
showClass("tawnT2Copula")
}
\keyword{classes}
