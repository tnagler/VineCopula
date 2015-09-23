\name{daxreturns}
\docType{data}
\alias{daxreturns}

\title{Major German Stocks}

\description{
This data set contains transformed standardized residuals of daily log returns of 15 major German stocks represented in the index DAX observed from January 2005 to August 2009.
Each time series is filtered using a GARCH(1,1) model with Student t innovations.
}

\format{
A data frame with 1158 observations on 15 variables. Column names correspond to ticker symbols of the stocks.
}

\source{Yahoo! Finance}

\seealso{\code{\link{RVineStructureSelect}}}

\examples{
# load the data set
data(daxreturns)

# compute the empirical Kendall's tau matrix
TauMatrix(daxreturns)
}
