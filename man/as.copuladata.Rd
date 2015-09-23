\name{as.copuladata}
\alias{as.copuladata}
\alias{as.copuladata.data.frame}
\alias{as.copuladata.matrix}
\alias{as.copuladata.list}

\title{
Copula Data Objects
}

\description{
The function \code{as.copuladata} coerces an object (\code{data.frame}, \code{matrix}, \code{list}) to a \code{copuladata} object.
}

\usage{
as.copuladata(data)
}

\arguments{
  \item{data}{Either a \code{data.frame}, a \code{matrix} or a \code{list} containing copula data (i.e. data with uniform margins on [0,1]). The \code{list} elements have to be vectors of identical length.}
}

\author{Tobias Erhardt}

\seealso{\code{\link{pobs}}, \code{\link{pairs.copuladata}}}

\examples{
data(daxreturns)

data <- as(daxreturns, "matrix")
class(as.copuladata(data))

data <- as(daxreturns, "data.frame")
class(as.copuladata(data))

data <- as(daxreturns, "list")
names(data) <- names(daxreturns)
class(as.copuladata(data))
}
