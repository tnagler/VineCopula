\name{TauMatrix}      
\alias{TauMatrix}

\title{Matrix of Empirical Kendall's Tau Values}

\description{
This function computes the empirical Kendall's tau using the algorithm by Knight (1966).
}

\usage{
TauMatrix(data, weights = NA)
}

\arguments{
  \item{data}{An N x d data matrix.}
  \item{weights}{Numerical; weights for each observation (opitional).}
}

\value{
Matrix of the empirical Kendall's taus. 
}

\references{
Knight, W. R. (1966).
A computer method for calculating Kendall's tau with ungrouped data.
Journal of the American Statistical Association 61 (314), 436-439.
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{BiCopTau2Par}}, \code{\link{BiCopPar2Tau}}, \code{\link{BiCopEst}}}

\examples{
data(daxreturns)
Data <- as.matrix(daxreturns)

# compute the empirical Kendall's taus
TauMatrix(Data)
}

