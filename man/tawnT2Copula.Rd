\name{tawnT2Copula}
\alias{tawnT2Copula}
\alias{surTawnT2Copula}
\alias{r90TawnT2Copula}
\alias{r270TawnT2Copula}

\title{
Constructor of the Tawn Type 2 Family and Rotated Versions thereof
}
\description{
Constructs an object of the \code{\linkS4class{tawnT2Copula}} (survival \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270}) family for given parameters.
}
\usage{
tawnT2Copula(param = c(2, 0.5))
surTawnT2Copula(param = c(2, 0.5))
r90TawnT2Copula(param = c(-2, 0.5))
r270TawnT2Copula(param = c(-2, 0.5))
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{param1} and \code{param2}.
}
}

\value{
One of the Tawn type 2 copula classes (\code{\linkS4class{tawnT2Copula}}, \code{\linkS4class{surTawnT2Copula}}, \code{\linkS4class{r90TawnT2Copula}}, \code{\linkS4class{r270TawnT2Copula}}).
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{tawnT2Copula}} and the package \code{\link{VineCopula-package}} for implementation details.
}
\examples{
library(copula)

persp(tawnT2Copula(), dCopula, zlim = c(0,10))
persp(surTawnT2Copula(), dCopula, zlim = c(0,10))
persp(r90TawnT2Copula(), dCopula, zlim = c(0,10))
persp(r270TawnT2Copula(), dCopula, zlim = c(0,10))
}
\keyword{ distribution }
\keyword{ copula } 
