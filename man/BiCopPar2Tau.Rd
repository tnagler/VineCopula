\name{BiCopPar2Tau}      
\alias{BiCopPar2Tau} 

\title{Kendall's Tau Value of a Bivariate Copula}

\description{
This function computes the theoretical Kendall's tau value of a bivariate copula for given parameter values. 
}

\usage{
BiCopPar2Tau(family, par, par2 = 0, obj = NULL, check.pars = TRUE)
}

\arguments{
  \item{family}{integer; single number or vector of size \code{m}; defines the bivariate copula family: \cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
	        \code{3} = Clayton copula \cr
	        \code{4} = Gumbel copula \cr
	        \code{5} = Frank copula \cr
	        \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
		\code{8} = BB6 copula \cr
		\code{9} = BB7 copula \cr
		\code{10} = BB8 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees) \cr
		\code{28} = rotated BB6 copula (90 degrees) \cr
		\code{29} = rotated BB7 copula (90 degrees) \cr
		\code{30} = rotated BB8 copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees) \cr
		\code{38} = rotated BB6 copula (270 degrees) \cr
		\code{39} = rotated BB7 copula (270 degrees) \cr
		\code{40} = rotated BB8 copula (270 degrees) \cr
    \code{104} = Tawn type 1 copula \cr
    \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
    \code{124} = rotated Tawn type 1 copula (90 degrees)  \cr
    \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
    \code{204} = Tawn type 2 copula  \cr
    \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
    \code{224} = rotated Tawn type 2 copula (90 degrees)  \cr
    \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
		}
  \item{par}{numeric; single number or vector of size \code{m}; copula parameter.}
  \item{par2}{numeric; single number or vector of size \code{m}; second parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default: \code{par2 = 0}).
    Note that the degrees of freedom parameter of the t-copula does not need to be set,
    because the theoretical Kendall's tau value of the t-copula is independent of this choice.}
  \item{obj}{\code{BiCop} object containing the family and parameter specification.}
  \item{check.pars}{logical; default is \code{TRUE}; if \code{FALSE}, checks for family/parameter-consistency are ommited (should only be used with care).}
}                                                                                  

\details{
If the family and parameter specification is stored in a \code{\link{BiCop}} object \code{obj}, the alternative version \cr
\preformatted{BiCopPar2Tau(obj)} 
can be used.
}

\note{
The number \code{m} can be chosen arbitrarily. 
}

\value{
Theoretical value of Kendall's tau (vector) corresponding to the bivariate copula \code{family} and parameter vector \eqn{(\theta, \delta) =} \code{(par, par2)}.
\tabular{ll}{
No. (\code{family}) \tab Kendall's tau (\code{tau}) \cr
\code{1, 2} \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
\code{3, 13} \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
\code{4, 14} \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
\code{5} \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta + 4 D_1(\theta)/\theta} \cr
	 \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{x/\theta}{\exp(x)-1}dx}{D_1(\theta)=\int_0^\theta (x/\theta)/(exp(x)-1)dx} (Debye function) \cr
\code{6, 16} \tab \eqn{1+\frac{4}{\theta^2}\int_0^1 x\log(x)(1-x)^{2(1-\theta)/\theta}dx}{1+4/\theta^2\int_0^1 x\log(x)(1-x)^{2(1-\theta)/\theta}dx} \cr
\code{7, 17} \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
\code{8, 18} \tab \eqn{1+4\int_0^1 -\log(-(1-t)^\theta+1)(1-t-(1-t)^{-\theta}+(1-t)^{-\theta}t)/(\delta\theta) dt} \cr
\code{9, 19} \tab \eqn{1+4\int_0^1 ( (1-(1-t)^{\theta})^{-\delta} - )/( -\theta\delta(1-t)^{\theta-1}(1-(1-t)^{\theta})^{-\delta-1} ) dt} \cr
\code{10, 20} \tab \eqn{1+4\int_0^1 -\log \left( ((1-t\delta)^\theta-1)/((1-\delta)^\theta-1) \right) } \cr
	 \tab \eqn{* (1-t\delta-(1-t\delta)^{-\theta}+(1-t\delta)^{-\theta}t\delta)/(\theta\delta) dt} \cr
\code{23, 33} \tab \eqn{\frac{\theta}{2-\theta}}{\theta/(2-\theta)} \cr
\code{24, 34} \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
\code{26, 36} \tab \eqn{-1-\frac{4}{\theta^2}\int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx}{-1-4/\theta^2\int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx} \cr
\code{27, 37} \tab \eqn{-1-\frac{2}{\delta(2-\theta)}}{1-2/(\delta(\theta+2))} \cr
\code{28, 38} \tab \eqn{-1-4\int_0^1 -\log(-(1-t)^{-\theta}+1)(1-t-(1-t)^{\theta}+(1-t)^{\theta}t)/(\delta\theta) dt} \cr
\code{29, 39} \tab \eqn{-1-4\int_0^1 ( (1-(1-t)^{-\theta})^{\delta} - )/( -\theta\delta(1-t)^{-\theta-1}(1-(1-t)^{-\theta})^{\delta-1} ) dt} \cr
\code{30, 40} \tab \eqn{-1-4\int_0^1 -\log \left( ((1+t\delta)^{-\theta}-1)/((1+\delta)^{-\theta}-1) \right)} \cr 
	 \tab \eqn{* (1+t\delta-(1+t\delta)^{\theta}-(1+t\delta)^{\theta}t\delta)/(\theta\delta) dt} \cr
   \code{104,114} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{\theta}+t^{\theta}]^{1/\theta}} \cr
   \code{204,214} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
   \code{124,134} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{-\theta}+t^{-\theta}]^{-1/\theta}} \cr
   \code{224,234} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
   
}
}

\author{Ulf Schepsmeier, Tobias Erhardt}

\references{
Joe, H. (1997).
Multivariate Models and Dependence Concepts.
Chapman and Hall, London.

Czado, C., U. Schepsmeier, and A. Min (2012).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Statistical Modelling, 12(3), 229-255.
}

\seealso{\code{\link{BiCopTau2Par}}, \code{\link{BiCop}}}

\examples{
## Example 1: Gaussian copula
tau0 <- 0.5
rho <- BiCopTau2Par(family = 1, tau = tau0)

# transform back
tau <- BiCopPar2Tau(family = 1, par = rho)
tau - 2/pi*asin(rho)


## Example 2: Clayton copula
theta <- BiCopTau2Par(family = 3, tau = c(0.4, 0.5, 0.6))
BiCopPar2Tau(family = 3, par = theta)


## Example 3:
vpar <- seq(from = 1.1, to = 10, length.out = 100)
tauC <- BiCopPar2Tau(family = 3, par = vpar)
tauG <- BiCopPar2Tau(family = 4, par = vpar)
tauF <- BiCopPar2Tau(family = 5, par = vpar)
tauJ <- BiCopPar2Tau(family = 6, par = vpar)
plot(tauC ~ vpar, type = "l", ylim = c(0,1))
lines(tauG ~ vpar, col = 2)
lines(tauF ~ vpar, col = 3)
lines(tauJ ~ vpar, col = 4)

\dontshow{
# Test BiCopPar2Tau (one parametric families)
theta <- BiCopTau2Par(family = 0, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 0, par = theta)
theta <- BiCopTau2Par(family = 1, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 1, par = theta)
theta <- BiCopTau2Par(family = 3, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 3, par = theta)
theta <- BiCopTau2Par(family = 4, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 4, par = theta)
theta <- BiCopTau2Par(family = 5, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 5, par = theta)
theta <- BiCopTau2Par(family = 6, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 6, par = theta)
theta <- BiCopTau2Par(family = 13, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 13, par = theta)
theta <- BiCopTau2Par(family = 14, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 14, par = theta)
theta <- BiCopTau2Par(family = 16, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 16, par = theta)
theta <- BiCopTau2Par(family = 23, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 23, par = theta)
theta <- BiCopTau2Par(family = 24, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 24, par = theta)
theta <- BiCopTau2Par(family = 26, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 26, par = theta)
theta <- BiCopTau2Par(family = 33, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 33, par = theta)
theta <- BiCopTau2Par(family = 34, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 34, par = theta)
theta <- BiCopTau2Par(family = 36, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 36, par = theta)
theta <- BiCopTau2Par(family = 41, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 41, par = theta)
theta <- BiCopTau2Par(family = 51, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 51, par = theta)
theta <- BiCopTau2Par(family = 61, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 61, par = theta)
theta <- BiCopTau2Par(family = 71, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 71, par = theta)
theta <- BiCopTau2Par(family = 41, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 41, par = theta)
theta <- BiCopTau2Par(family = 51, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 51, par = theta)
theta <- BiCopTau2Par(family = 61, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 61, par = theta)
theta <- BiCopTau2Par(family = 71, tau = -c(0.4,0.5,0.6))
BiCopPar2Tau(family = 71, par = theta)

# Test BiCopPar2Tau (two parametric families)
theta <- BiCopTau2Par(family = 2, tau = c(0.4,0.5,0.6))
BiCopPar2Tau(family = 2, par = theta)
theta <- 1:3
delta <- 1:3
BiCopPar2Tau(family = 7, par = theta, par2 = delta)
BiCopPar2Tau(family = 17, par = theta, par2 = delta)
theta <- -(1:3)
delta <- -(1:3)
BiCopPar2Tau(family = 27, par = theta, par2 = delta)
BiCopPar2Tau(family = 37, par = theta, par2 = delta)
theta <- 2:4
delta <- 1:3
BiCopPar2Tau(family = 8, par = theta, par2 = delta)
BiCopPar2Tau(family = 18, par = theta, par2 = delta)
theta <- -(2:4)
delta <- -(1:3)
BiCopPar2Tau(family = 28, par = theta, par2 = delta)
BiCopPar2Tau(family = 38, par = theta, par2 = delta)
theta <- 1:3
delta <- 1:3
BiCopPar2Tau(family = 9, par = theta, par2 = delta)
BiCopPar2Tau(family = 19, par = theta, par2 = delta)
theta <- -(1:3)
delta <- -(1:3)
BiCopPar2Tau(family = 29, par = theta, par2 = delta)
BiCopPar2Tau(family = 39, par = theta, par2 = delta)
theta <- 2:4
delta <- c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 10, par = theta, par2 = delta)
BiCopPar2Tau(family = 20, par = theta, par2 = delta)
theta <- -(2:4)
delta <- -c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 30, par = theta, par2 = delta)
BiCopPar2Tau(family = 40, par = theta, par2 = delta)

theta <- 2:4
delta <- c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 104, par = theta, par2 = delta)
BiCopPar2Tau(family = 114, par = theta, par2 = delta)
theta <- -(2:4)
delta <- c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 124, par = theta, par2 = delta)
BiCopPar2Tau(family = 134, par = theta, par2 = delta)

theta <- 2:4
delta <- c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 204, par = theta, par2 = delta)
BiCopPar2Tau(family = 214, par = theta, par2 = delta)
theta <- -(2:4)
delta <- c(0.1, 0.5, 0.9)
BiCopPar2Tau(family = 224, par = theta, par2 = delta)
BiCopPar2Tau(family = 234, par = theta, par2 = delta)
}
}

