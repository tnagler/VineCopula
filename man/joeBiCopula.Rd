\name{joeBiCopula}
\alias{joeBiCopula}
\alias{surJoeBiCopula}
\alias{r90JoeBiCopula}
\alias{r270JoeBiCopula}

\title{
Constructor of the Joe Family and Rotated Versions thereof
}
\description{
Constructs an object of the (survival \code{surJoeBiCopula}, 90 degree rotated \code{r90JoeBiCopula} and 270 degree rotated \code{r270JoeBiCopula}) family for a given parameter. Note that package \code{\link{copula-package}} provides a class \code{\linkS4class{joeCopula}} as well.
}
\usage{
surJoeBiCopula(param)
r90JoeBiCopula(param)
r270JoeBiCopula(param)
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{theta} and \code{delta}.
}
}
\value{
One of the respective Joe copula classes (\code{\linkS4class{joeBiCopula}}, \code{\linkS4class{surJoeBiCopula}}, \code{\linkS4class{r90JoeBiCopula}}, \code{\linkS4class{r270JoeBiCopula}}).
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}
\seealso{
See also \code{\link{BB1Copula}}, \code{\link{BB6Copula}}, \code{\link{BB7Copula}} and \code{\link{BB8Copula}} for further wrapper functions to the \code{\link{VineCopula-package}}. 
}
\examples{
library(copula)

persp(surJoeBiCopula(1.5), dCopula, zlim = c(0,10))
persp(r90JoeBiCopula(-1.5), dCopula, zlim = c(0,10))
persp(r270JoeBiCopula(-1.5), dCopula, zlim = c(0,10))
}