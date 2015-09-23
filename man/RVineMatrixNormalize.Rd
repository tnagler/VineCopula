\name{RVineMatrixNormalize}
\Rdversion{1.1}
\alias{RVineMatrixNormalize}

\title{Normalization of R-Vine Matrix}
\description{
An \code{\link{RVineMatrix}} is permuted to achieve a natural ordering (i.e. \code{diag(RVM$Matrix) == d:1})
}

\usage{
RVineMatrixNormalize(RVM)
}
\arguments{
  \item{RVM}{\code{\link{RVineMatrix}} defining the R-vine structure}

}
\value{
  \item{RVM}{An \code{\link{RVineMatrix}} in natural ordering with entries in \code{RVM$names} keeping track of the reordering.}
}

\examples{
Matrix <- matrix(c(5, 2, 3, 1, 4,
                   0, 2, 3, 4, 1,
                   0, 0, 3, 4, 1,
                   0, 0, 0, 4, 1,
                   0, 0, 0, 0, 1), 5, 5)
family <- matrix(1,5,5)

par <- matrix(c(0, 0.2, 0.9, 0.5, 0.8,
                0,   0, 0.1, 0.6, 0.9,
                0,   0,   0, 0.7, 0.5,
                0,   0,   0,   0, 0.8,
                0,   0,   0,   0,   0), 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix, family, par)

# normalise the RVine
RVineMatrixNormalize(RVM)
}

\keyword{vine}
