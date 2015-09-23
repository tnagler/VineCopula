\name{surClaytonCopula}
\alias{surClaytonCopula}
\alias{r90ClaytonCopula}
\alias{r270ClaytonCopula}

\title{
Survival and Rotated Clayton Copulas
}
\description{
These are wrappers to functions from \code{\link{VineCopula-package}}
}
\usage{
surClaytonCopula(param)
r90ClaytonCopula(param)
r270ClaytonCopula(param)
}
\arguments{
  \item{param}{
A single parameter defining the Copula.
}
}
\value{
An object of class \code{\linkS4class{surClaytonCopula}}, \code{\linkS4class{r90ClaytonCopula}} or \code{\linkS4class{r270ClaytonCopula}} respectively.
}
\author{
Benedikt Graeler
}
\examples{
library(copula)

persp(surClaytonCopula(1.5), dCopula, zlim = c(0,10))
persp(r90ClaytonCopula(-1.5), dCopula, zlim = c(0,10))
persp(r270ClaytonCopula(-1.5), dCopula, zlim = c(0,10))
}
\keyword{copula}
