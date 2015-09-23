\name{joeBiCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{joeBiCopula-class}
\alias{dduCopula,numeric,joeBiCopula-method}
\alias{ddvCopula,numeric,joeBiCopula-method}
\alias{dduCopula,matrix,joeBiCopula-method}
\alias{ddvCopula,matrix,joeBiCopula-method}

\alias{getKendallDistr,joeBiCopula-method}
\alias{kendallDistribution,joeBiCopula-method}

\alias{surJoeBiCopula-class}
\alias{dduCopula,numeric,surJoeBiCopula-method}
\alias{ddvCopula,numeric,surJoeBiCopula-method}
\alias{dduCopula,matrix,surJoeBiCopula-method}
\alias{ddvCopula,matrix,surJoeBiCopula-method}


\alias{r90JoeBiCopula-class}
\alias{dduCopula,numeric,r90JoeBiCopula-method}
\alias{ddvCopula,numeric,r90JoeBiCopula-method}
\alias{dduCopula,matrix,r90JoeBiCopula-method}
\alias{ddvCopula,matrix,r90JoeBiCopula-method}

\alias{r270JoeBiCopula-class}
\alias{dduCopula,numeric,r270JoeBiCopula-method}
\alias{ddvCopula,numeric,r270JoeBiCopula-method}
\alias{dduCopula,matrix,r270JoeBiCopula-method}
\alias{ddvCopula,matrix,r270JoeBiCopula-method}

\title{Classes \code{"joeBiCopula"}, \code{"surJoeBiCopula"}, \code{"r90JoeBiCopula"} and \code{"r270JoeBiCopula"}}
\description{
Wrapper classes representing the bivariate Joe, survival Joe, 90 degree and 270 degree rotated Joe copula families (Joe 1997) from \code{\link{VineCopula-package}}. Note that package \code{\link{copula-package}} provides a class \code{\linkS4class{joeCopula}} as well.
}
\section{Objects from the Classes}{
Objects can be created by calls of the form \code{new("joeBiCopula", ...)}, \code{new("surJoeBiCopula", ...)}, \code{new("r90JoeBiCopula", ...)} and \code{new("r270JoeBiCopula", ...)} or by the functions \code{\link{joeBiCopula}}, \code{\link{surJoeBiCopula}}, \code{\link{r90JoeBiCopula}} and \code{\link{r270JoeBiCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"numeric"} defining the family number in \code{\link{VineCopula-package}} }
    \item{\code{dimension}:}{Object of class \code{"integer"} defining the dimension of the copula }
    \item{\code{parameters}:}{Object of class \code{"numeric"} the single parameter}
    \item{\code{param.names}:}{Object of class \code{"character"}, parameter name.}
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}, lower bound of the copula parameter}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}, upper bound of the copula parameter}
    \item{\code{fullname}:}{Object of class \code{"character"}, family name of the copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "joeBiCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "joeBiCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "joeBiCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "joeBiCopula")}: ... }
    \item{getKendallDistr}{\code{signature(copula = "joeBiCopula")}: ... }
    \item{kendallDistribution}{\code{signature(copula = "joeBiCopula")}: ... }
   }
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\linkS4class{BB1Copula}}, \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}} and \code{\linkS4class{BB8Copula}} for further wrapper classes to the \code{\link{VineCopula-package}}.
}
\examples{
showClass("surJoeBiCopula")
}
\keyword{classes}