\name{BB8Copula}
\alias{BB8Copula}
\alias{surBB8Copula}
\alias{r90BB8Copula}
\alias{r270BB8Copula}

\title{
Constructor of the BB8 Family and Rotated Versions thereof
}
\description{
Constructs an object of the \code{\linkS4class{BB8Copula}} (survival \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270}) family for given parameters.
}
\usage{
BB8Copula(param)
surBB8Copula(param)
r90BB8Copula(param)
r270BB8Copula(param)
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{theta} and \code{delta}.
}
}
\value{
One of the respective BB8 copula classes (\code{\linkS4class{BB8Copula}}, \code{\linkS4class{surBB8Copula}}, \code{\linkS4class{r90BB8Copula}}, \code{\linkS4class{r270BB8Copula}}).
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

persp(BB8Copula(c(2,0.9)), dCopula, zlim = c(0,10))
persp(surBB8Copula(c(2,0.9)), dCopula, zlim = c(0,10))
persp(r90BB8Copula(c(-2,-0.9)), dCopula, zlim = c(0,10))
persp(r270BB8Copula(c(-2,-0.9)), dCopula, zlim = c(0,10))
}