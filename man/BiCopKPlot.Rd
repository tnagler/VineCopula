\name{BiCopKPlot}           
\alias{BiCopKPlot}         

\title{Kendall's Plot for Bivariate Copula Data}

\description{
This function creates a Kendall's plot (K-plot) of given bivariate copula data.
}

\usage{
BiCopKPlot(u1, u2, PLOT = TRUE, ...)
}

\arguments{
  \item{u1, u2}{Data vectors of equal length with values in [0,1].}
  \item{PLOT}{Logical; whether the results are plotted. If \code{PLOT = FALSE}, the values \code{W.in} and \code{Hi.sort} are returned (see below; default: \code{PLOT = TRUE}).}
  \item{...}{Additional plot arguments.}
}

\value{
  \item{W.in}{W-statistics (x-axis).}
  \item{Hi.sort}{H-statistics (y-axis).}
}

\details{
For observations \eqn{u_{i,j},\ i=1,...,N,\ j=1,2,}{u_{i,j}, i=1,...,N, j=1,2,} the K-plot considers two quantities:
First, the ordered values of the empirical bivariate distribution function \eqn{H_i:=\hat{F}_{U_1U_2}(u_{i,1},u_{i,2})}
and, second, \eqn{W_{i:N}}, which are the expected values of the order statistics
from a random sample of size \eqn{N} of the random variable \eqn{W=C(U_1,U_2)} under the null hypothesis of independence between \eqn{U_1} and \eqn{U_2}.
\eqn{W_{i:N}} can be calculated as follows
\deqn{
W_{i:n}= N {N-1 \choose i-1} \int\limits_{0}^1 \omega k_0(\omega) ( K_0(\omega) )^{i-1} ( 1-K_0(\omega) )^{N-i} d\omega,
}
where
\deqn{
K_0(\omega)=\omega - \omega \log(\omega), 
\nonumber
}{
K_=(\omega)=\omega - \omega log(\omega)
}
and \eqn{k_0(\cdot)}{k_0()} is the corresponding density.

K-plots can be seen as the bivariate copula equivalent to QQ-plots.
If the points of a K-plot lie approximately on the diagonal \eqn{y=x}, then \eqn{U_1} and \eqn{U_2} are approximately independent.
Any deviation from the diagonal line points towards dependence.
In case of positive dependence, the points of the K-plot should be located above the diagonal line, and vice versa for negative dependence.
The larger the deviation from the diagonal, the stronger is the degree of dependency.
There is a perfect positive dependence if points \eqn{\left(W_{i:N},H_i\right)} lie on the curve \eqn{K_0(\omega)} located above the main diagonal.
If points \eqn{\left(W_{i:N},H_i\right)}{(W_{i:N},H_i)} however lie on the x-axis, this indicates a perfect negative dependence between \eqn{U_1} and \eqn{U_2}. 
}

\author{Natalia Belgorodski, Ulf Schepsmeier}

\references{
Genest, C. and A. C. Favre (2007).
Everything you always wanted to know about copula modeling but were afraid to ask.
Journal of Hydrologic Engineering, 12 (4), 347-368.
}

\seealso{\code{\link{BiCopMetaContour}}, \code{\link{BiCopChiPlot}}, \code{\link{BiCopLambda}}, \code{\link{BiCopGofTest}}}

\examples{
# Gaussian and Clayton copulas
n <- 500
tau <- 0.5

# simulate from Gaussian copula
fam1 <- 1  
theta1 <- BiCopTau2Par(fam1, tau)
set.seed(123)
dat1 <- BiCopSim(n, fam1, theta1)  

# simulate from Clayton copula
fam2 <- 3
theta2 <- BiCopTau2Par(fam2, tau)
set.seed(123)
dat2 <- BiCopSim(n, fam2, theta2)

# create K-plots
par(mfrow=c(1,2))
BiCopKPlot(dat1[,1], dat1[,2], main = "Gaussian copula")
BiCopKPlot(dat2[,1], dat2[,2], main = "Clayton copula")
}
