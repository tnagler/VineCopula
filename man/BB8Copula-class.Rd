\name{BB8Copula-class}
\Rdversion{1.1}
\docType{class}
\alias{BB8Copula-class}
\alias{dduCopula,numeric,BB8Copula-method}
\alias{ddvCopula,numeric,BB8Copula-method}
\alias{dduCopula,matrix,BB8Copula-method}
\alias{ddvCopula,matrix,BB8Copula-method}

\alias{getKendallDistr,BB8Copula-method}
\alias{kendallDistribution,BB8Copula-method}

\alias{surBB8Copula-class}
\alias{dduCopula,numeric,surBB8Copula-method}
\alias{ddvCopula,numeric,surBB8Copula-method}
\alias{dduCopula,matrix,surBB8Copula-method}
\alias{ddvCopula,matrix,surBB8Copula-method}

\alias{r90BB8Copula-class}
\alias{dduCopula,numeric,r90BB8Copula-method}
\alias{ddvCopula,numeric,r90BB8Copula-method}
\alias{dduCopula,matrix,r90BB8Copula-method}
\alias{ddvCopula,matrix,r90BB8Copula-method}

\alias{r270BB8Copula-class}
\alias{dduCopula,numeric,r270BB8Copula-method}
\alias{ddvCopula,numeric,r270BB8Copula-method}
\alias{dduCopula,matrix,r270BB8Copula-method}
\alias{ddvCopula,matrix,r270BB8Copula-method}

\alias{fitCopula,twoParamBiCop-method}

\title{Classes \code{"BB8Copula"}, \code{"surBB8Copula"}, \code{"r90BB8Copula"} and \code{"r270BB8Copula"}}
\description{
Wrapper classes representing the BB8, survival BB8, 90 degree and 270 degree rotated BB8 copula families (Joe 1997) from the \code{\link{VineCopula-package}} package.
}
\section{Objects from the Classes}{
Objects can be created by calls of the form \code{new("BB8Copula", ...)}, \code{new("surBB8Copula", ...)}, \code{new("r90BB8Copula", ...)} and \code{new("r270BB8Copula", ...)} or by the functions \code{\link{BB8Copula}}, \code{\link{surBB8Copula}}, \code{\link{r90BB8Copula}} and \code{\link{r270BB8Copula}}.
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
    \item{dduCopula}{\code{signature(u = "matrix", copula = "BB8Copula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "BB8Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "BB8Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "BB8Copula")}: ... }
    \item{getKendallDistr}{\code{signature(copula = "BB8Copula")}: ... }
    \item{kendallDistribution}{\code{signature(copula = "BB8Copula")}: ... }
   }
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\linkS4class{BB1Copula}}, \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}} and \code{\linkS4class{joeCopula}} for further wrapper classes to the \code{\link{VineCopula-package}}.
}
\examples{
showClass("BB8Copula")
}
\keyword{classes}