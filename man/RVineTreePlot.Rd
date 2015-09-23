\name{RVineTreePlot}      
\alias{RVineTreePlot}

\title{Visualisation of R-Vine Tree Structure}

\description{
Function is deprecated since \code{VineCopula 2.0}. Use \code{\link[VineCopula:plot.RVineMatrix]{plot.RVineMatrix}} instead.
}

\usage{
RVineTreePlot(x, tree = "ALL", type = 0, edge.labels = NULL, legend.pos = "bottomleft",
              interactive = FALSE, ...)
}

\arguments{
  \item{x}{\code{RVineMatrix} object.}
  \item{tree}{\code{"ALL"} or integer vector; specifies which trees are plotted.}
  \item{type}{integer; specifies how to make use of variable names: \cr
  \code{0} = variable names are ignored, \cr 
  \code{1} = variable names are used to annotate vertices, \cr
  \code{2} = uses numbers in plot and adds a legend for variable names.}
  \item{edge.labels}{character; either a vector of edge labels
  or one of the following: \cr
  \code{"family"} = pair-copula family abbreviation (see \code{\link[VineCopula:BiCopName]{BiCopName}}), \cr
  \code{"par"} = pair-copula parameters, \cr
  \code{"tau"} = pair-copula Kendall's tau (by conversion of parameters) \cr
  \code{"family-par"} = pair-copula family and parameters \cr
  \code{"family-tau"} = pair-copula family and Kendall's tau.
  }
  \item{legend.pos}{the \code{x} argument for \code{\link[graphics:legend]{legend}}.}
  \item{interactive}{logical; if TRUE, the user is asked to adjust the positioning of 
  vertices with his mouse.}
  \item{\dots}{
Arguments passed to \code{\link[network:plot.network]{plot.network}}.}
}


\author{Thomas Nagler}

\seealso{\code{\link[VineCopula:plot.RVineMatrix]{plot.RVineMatrix}}}
