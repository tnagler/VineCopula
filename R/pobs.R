#' Pseudo-Observations
#'
#' Compute the pseudo-observations for the given data matrix.
#'
#' Given \eqn{n} realizations
#' \eqn{\bm{x}_i=(x_{i1},\dots,x_{id})^T}{x_i=(x_{i1},...,x_{id})},
#' \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} of a random vector \eqn{\bm{X}}{X},
#' the pseudo-observations are defined via \eqn{u_{ij}=r_{ij}/(n+1)} for
#' \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} and \eqn{j\in\{1,\dots,d\}}{j in
#' {1,...,d}}, where \eqn{r_{ij}} denotes the rank of \eqn{x_{ij}} among all
#' \eqn{x_{kj}}, \eqn{k\in\{1,\dots,n\}}{k in {1,...,n}}.  The
#' pseudo-observations can thus also be computed by component-wise applying the
#' empirical distribution functions to the data and scaling the result by
#' \eqn{n/(n+1)}.  This asymptotically negligible scaling factor is used to
#' force the variates to fall inside the open unit hypercube, for example, to
#' avoid problems with density evaluation at the boundaries. Note that
#' \code{pobs(, lower.tail=FALSE)} simply returns \code{1-pobs()}.
#'
#' @name pobs
#' @param x \eqn{n\times d}{n x d}-matrix of random variates to be converted to
#' pseudo-observations.
#' @param na.last,ties.method are passed to \code{\link{rank}}; see there.
#' @param lower.tail \code{\link{logical}} which, if \code{FALSE}, returns the
#' pseudo-observations when applying the empirical marginal survival functions.
#' @return matrix of the same dimensions as \code{x} containing the
#' pseudo-observations.
#' @note This function is borrowed from the
#' \code{\link[copula]{copula-package}}, see \code{\link[copula]{pobs}}.
#' @author Marius Hofert
#' @examples
#'
#' ## Simple definition of the function:
#' pobs
#'
#' ## simulate data from a multivariate normal distribution
#' library(mvtnorm)
#' set.seed(123)
#' Sigma <- matrix(c(2, 1, -0.2, 1, 1, 0.3, -0.2, 0.3, 0.5), 3, 3)
#' mu <- c(-3, 2, 1)
#' dat <- rmvnorm(500, sigma = Sigma)
#' pairs(dat)  # plot observations
#'
#' ## compute pseudo-observations for copula inference
#' udat <- pobs(dat)
#' pairs(udat)
#' # estimate vine copula model
#' fit <- RVineStructureSelect(udat, familyset = c(1, 2))
#'
NULL
