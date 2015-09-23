\name{RVineMatrix}
\alias{RVineMatrix}

\title{R-Vine Copula Model in Matrix Notation}

\description{
This function creates an \code{\link{RVineMatrix}} object which encodes an R-vine copula model.
It contains the matrix identifying the R-vine tree structure, the matrix identifying the copula families
utilized and two matrices for corresponding parameter values.
}

\usage{
RVineMatrix(Matrix, family = array(0, dim = dim(Matrix)),
            par = array(NA, dim = dim(Matrix)),
            par2 = array(NA, dim = dim(Matrix)), names=NULL, check.pars = TRUE)
}

\arguments{
  \item{Matrix}{Lower (or upper) triangular d x d matrix that defines the R-vine tree structure.}
  \item{family}{Lower  (or upper) triangular d x d matrix with zero diagonal entries that assigns the pair-copula families
    to each (conditional) pair defined by \code{Matrix} (default: \code{family = array(0,dim=dim(Matrix))}).
    The bivariate copula families are defined as follows:\cr
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
  \item{par}{Lower (or upper) triangular d x d matrix with zero diagonal entries that assigns the (first) pair-copula parameter
    to each (conditional) pair defined by \code{Matrix} \cr
    (default: \code{par = array(NA, dim = dim(Matrix))}).}
  \item{par2}{Lower (or upper) triangular d x d matrix with zero diagonal entries that assigns the second parameter
    for pair-copula families with two parameters to each (conditional) pair defined by \code{Matrix} (default: \code{par2 = array(NA, dim = dim(Matrix))}).}
  \item{names}{A vector of names for the d variables; default: \code{names = NULL}.}
  \item{check.pars}{logical; default is \code{TRUE}; if \code{FALSE}, checks for family/parameter-consistency are ommited (should only be used with care).}

}

\value{
	An \code{\link{RVineMatrix}} object with the following matrix components:
	\item{Matrix}{R-vine tree structure matrix.}
	\item{family}{Pair-copula family matrix with values as above.}
	\item{par}{Pair-copula parameter matrix.}
	\item{par2}{Second pair-copula parameter matrix with parameters necessary for pair-copula families with two parameters.}
}

\note{
The \code{print} function writes the R-vine matrix defined by \code{Matrix}. A detailed output is given by \code{print(RVM, detail=TRUE)},
where \code{RVM} is the \code{\link{RVineMatrix}} object. \cr
The \code{\link{RVineMatrix}} function automatically checks if the given matrix is a valid R-vine matrix (see \code{\link{RVineMatrixCheck}}). \cr
Although the function allows upper triangular matrices as its input, it will always store them as lower triangular matrices.
}

\references{
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.
}

\author{Jeffrey Dissmann}

\seealso{\code{\link{RVineMatrixCheck}}, \code{\link{RVineMLE}}, \code{\link{RVineSim}}, \code{\link{C2RVine}}, \code{\link{D2RVine}}}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))

# Print detailed information
print(RVM, detail = TRUE)
}
