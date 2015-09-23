\name{BB6Copula}
\alias{BB6Copula}
\alias{surBB6Copula}
\alias{r90BB6Copula}
\alias{r270BB6Copula}

\title{
Constructor of the BB6 Family and Rotated Versions thereof
}
\description{
Constructs an object of the \code{\linkS4class{BB6Copula}} (survival \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270}) family for given parameters.
}
\usage{
BB6Copula(param)
surBB6Copula(param)
r90BB6Copula(param)
r270BB6Copula(param)
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{theta} and \code{delta}.
}
}
\value{
One of the respective BB6 copula classes (\code{\linkS4class{BB6Copula}}, \code{\linkS4class{surBB6Copula}}, \code{\linkS4class{r90BB6Copula}}, \code{\linkS4class{r270BB6Copula}}).
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}
\seealso{
See also \code{\link{BB6Copula}}, \code{\link{BB7Copula}}, \code{\link{BB8Copula}} and \code{\link{joeCopula}} for further wrapper functions to the \code{\link{VineCopula-package}}. 
}
\examples{
library(copula)

persp(BB6Copula(c(1,1.5)), dCopula, zlim = c(0,10))
persp(surBB6Copula(c(1,1.5)), dCopula, zlim = c(0,10))
persp(r90BB6Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
persp(r270BB6Copula(c(-1,-1.5)), dCopula, zlim = c(0,10))
}