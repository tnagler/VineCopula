\name{vineCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{vineCopula-class}
\alias{fitCopula,vineCopula-method}

\title{Class \code{"vineCopula"}}
\description{
A class representing vine copulas in a object oriented implementations. Many functions go back to the package \code{\link{VineCopula-package}}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("vineCopula", ...)} or through the function \code{\link{vineCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{RVM}:}{An \code{RVineMatrix} object from \code{\link{RVineMatrix}} describing the R-Vine structure.}
    \item{\code{copulas}:}{Object of class \code{"list"} holding all copulas. }
    \item{\code{dimension}:}{Object of class \code{"integer"}; the vines dimension. }
    \item{\code{parameters}:}{Object of class \code{"numeric"}: empty }
    \item{\code{param.names}:}{Object of class \code{"character"}: empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}: empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}: empty }
    \item{\code{fullname}:}{Object of class \code{"character"} providing a descriptive name of the vine copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
No additional methods yet, but uses e.g. \code{\link{dCopula}}, \code{\link{pCopula}}, \code{\link{rCopula}} and \code{\link{rCopula}} as any other copula. Via the \code{method} argument in \code{fitCopula}, control over the  fit of the RVine can be taken via entries \code{StructureSelect}, \code{indeptest} and \code{familyset}. See \code{\link{RVineCopSelect}} and \code{\link{RVineStructureSelect}} for further details on the underlying functions. Missing entries are treated as default values, i.e. \code{StructureSelect=FALSE}, \code{indeptest=FALSE} and \code{familyset=NA}
}

\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence Insurance: Mathematics and Economics 44 (2), 182-198. 
}
\seealso{
\code{\link{RVineMatrix}} from package \code{\link{VineCopula-package}}
}
\author{
Benedikt Graeler
}

\examples{
showClass("vineCopula")
}
\keyword{classes}
