\name{tawnT1Copula}
\alias{tawnT1Copula}
\alias{surTawnT1Copula}
\alias{r90TawnT1Copula}
\alias{r270TawnT1Copula}

\title{
Constructor of the Tawn Type 1 Family and Rotated Versions thereof
}
\description{
Constructs an object of the \code{\linkS4class{tawnT1Copula}} (survival \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270}) family for given parameters.
}
\usage{
tawnT1Copula(param = c(2, 0.5))
surTawnT1Copula(param = c(2, 0.5))
r90TawnT1Copula(param = c(-2, 0.5))
r270TawnT1Copula(param = c(-2, 0.5))
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{param1} and \code{param2}.
}
}

\value{
One of the Tawn type 1 copula classes (\code{\linkS4class{tawnT1Copula}}, \code{\linkS4class{surTawnT1Copula}}, \code{\linkS4class{r90TawnT1Copula}}, \code{\linkS4class{r270TawnT1Copula}}).
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{tawnT2Copula}} and the package \code{\link{VineCopula-package}} for implementation details.
}
\examples{
library(copula)

persp(tawnT1Copula(), dCopula, zlim = c(0,10))
persp(surTawnT1Copula(), dCopula, zlim = c(0,10))
persp(r90TawnT1Copula(), dCopula, zlim = c(0,10))
persp(r270TawnT1Copula(), dCopula, zlim = c(0,10))
}
\keyword{ distribution }
\keyword{ copula } 
