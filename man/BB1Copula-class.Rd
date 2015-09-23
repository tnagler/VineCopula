\name{BB1Copula-class}
\Rdversion{1.1}
\docType{class}
\alias{BB1Copula-class}
\alias{dduCopula,numeric,BB1Copula-method}
\alias{ddvCopula,numeric,BB1Copula-method}
\alias{dduCopula,matrix,BB1Copula-method}
\alias{ddvCopula,matrix,BB1Copula-method}

\alias{getKendallDistr,BB1Copula-method}
\alias{kendallDistribution,BB1Copula-method}

\alias{surBB1Copula-class}
\alias{dduCopula,numeric,surBB1Copula-method}
\alias{ddvCopula,numeric,surBB1Copula-method}
\alias{dduCopula,matrix,surBB1Copula-method}
\alias{ddvCopula,matrix,surBB1Copula-method}

\alias{r90BB1Copula-class}
\alias{dduCopula,numeric,r90BB1Copula-method}
\alias{ddvCopula,numeric,r90BB1Copula-method}
\alias{dduCopula,matrix,r90BB1Copula-method}
\alias{ddvCopula,matrix,r90BB1Copula-method}

\alias{r270BB1Copula-class}
\alias{dduCopula,numeric,r270BB1Copula-method}
\alias{ddvCopula,numeric,r270BB1Copula-method}
\alias{dduCopula,matrix,r270BB1Copula-method}
\alias{ddvCopula,matrix,r270BB1Copula-method}

\title{Classes \code{"BB1Copula"}, \code{"surBB1Copula"}, \code{"r90BB1Copula"} and \code{"r270BB1Copula"}}
\description{
Wrapper classes representing the BB1, survival BB1, 90 degree and 270 degree rotated BB1 copula families (Joe 1997) from \code{\link{VineCopula-package}}.
}
\section{Objects from the Classes}{
Objects can be created by calls of the form \code{new("BB1Copula", ...)}, \code{new("surBB1Copula", ...)}, \code{new("r90BB1Copula", ...)} and \code{new("r270BB1Copula", ...)} or by the functions \code{\link{BB1Copula}}, \code{\link{surBB1Copula}}, \code{\link{r90BB1Copula}} and \code{\link{r270BB1Copula}}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"numeric"} defining the family number in \code{\link{VineCopula-package}} }
    \item{\code{dimension}:}{Object of class \code{"integer"} defining the dimension of the copula }
    \item{\code{parameters}:}{Object of class \code{"numeric"} the two-place parameter vector}
    \item{\code{param.names}:}{Object of class \code{"character"}, parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}, lower bounds of the copula parameters}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}, upper bounds of the copula parameters}
    \item{\code{fullname}:}{Object of class \code{"character"}, family name of the copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "BB1Copula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "BB1Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "BB1Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "BB1Copula")}: ... }
    \item{getKendallDistr}{\code{signature(copula = "BB1Copula")}: ... }
    \item{kendallDistribution}{\code{signature(copula = "BB1Copula")}: ... }
	 }
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}}, \code{\linkS4class{BB8Copula}} and \code{\linkS4class{joeCopula}} for further wrapper classes to the \code{\link{VineCopula-package}}.
}
\examples{
showClass("BB1Copula")
}
\keyword{classes}