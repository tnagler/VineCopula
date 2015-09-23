\name{BetaMatrix}      
\alias{BetaMatrix}

\title{Matrix of Empirical Blomqvist's Beta Values}

\description{
This function computes the empirical Blomqvist's beta.
}

\usage{
BetaMatrix(data)
}

\arguments{
  \item{data}{An N x d data matrix.}
}

\value{
Matrix of the empirical Blomqvist's betas. 
}

\references{
Blomqvist, N. (1950). 
On a measure of dependence between two random variables.
The Annals of Mathematical Statistics, 21(4), 593-600.

Nelsen, R. (2006).
An introduction to copulas. 
Springer
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{TauMatrix}}, \code{\link{BiCopPar2Beta}}, \code{\link{RVinePar2Beta}}}

\examples{
data(daxreturns)
Data <- as.matrix(daxreturns)

# compute the empirical Blomqvist's betas
BetaMatrix(Data)
}

