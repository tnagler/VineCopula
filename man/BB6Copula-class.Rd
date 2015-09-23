\name{BB6Copula-class}
\Rdversion{1.1}
\docType{class}
\alias{BB6Copula-class}
\alias{dduCopula,numeric,BB6Copula-method}
\alias{ddvCopula,numeric,BB6Copula-method}
\alias{dduCopula,matrix,BB6Copula-method}
\alias{ddvCopula,matrix,BB6Copula-method}

\alias{getKendallDistr,BB6Copula-method}
\alias{kendallDistribution,BB6Copula-method}

\alias{surBB6Copula-class}
\alias{dduCopula,numeric,surBB6Copula-method}
\alias{ddvCopula,numeric,surBB6Copula-method}
\alias{dduCopula,matrix,surBB6Copula-method}
\alias{ddvCopula,matrix,surBB6Copula-method}

\alias{r90BB6Copula-class}
\alias{dduCopula,numeric,r90BB6Copula-method}
\alias{ddvCopula,numeric,r90BB6Copula-method}
\alias{dduCopula,matrix,r90BB6Copula-method}
\alias{ddvCopula,matrix,r90BB6Copula-method}

\alias{r270BB6Copula-class}
\alias{dduCopula,numeric,r270BB6Copula-method}
\alias{ddvCopula,numeric,r270BB6Copula-method}
\alias{dduCopula,matrix,r270BB6Copula-method}
\alias{ddvCopula,matrix,r270BB6Copula-method}

\title{Classes \code{"BB6Copula"}, \code{"surBB6Copula"}, \code{"r90BB6Copula"} and \code{"r270BB6Copula"}}
\description{
Wrapper classes representing the BB6, survival BB6, 90 degree and 270 degree rotated BB6 copula families (Joe 1997) from the \code{\link{VineCopula-package}}.
}
\section{Objects from the Classes}{
Objects can be created by calls of the form \code{new("BB6Copula", ...)}, \code{new("surBB6Copula", ...)}, \code{new("r90BB6Copula", ...)} and \code{new("r270BB6Copula", ...)} or by the functions \code{\link{BB6Copula}}, \code{\link{surBB6Copula}}, \code{\link{r90BB6Copula}} and \code{\link{r270BB6Copula}}.
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
    \item{dduCopula}{\code{signature(u = "matrix", copula = "BB6Copula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "BB6Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "BB6Copula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "BB6Copula")}: ... }
    \item{getKendallDistr}{\code{signature(copula = "BB6Copula")}: ... }
    \item{kendallDistribution}{\code{signature(copula = "BB6Copula")}: ... }
   }
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\linkS4class{BB1Copula}}, \code{\linkS4class{BB7Copula}}, \code{\linkS4class{BB8Copula}} and \code{\linkS4class{joeCopula}} for further wrapper classes to the \code{\link{VineCopula-package}}.
}
\examples{
showClass("BB6Copula")
}
\keyword{classes}