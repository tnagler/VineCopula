#' Kendall's Tau Value of a Bivariate Copula
#'
#' This function computes the theoretical Kendall's tau value of a bivariate
#' copula for given parameter values.
#'
#' If the family and parameter specification is stored in a \code{\link{BiCop}}
#' object \code{obj}, the alternative version \cr
#' \preformatted{BiCopPar2Tau(obj)} can be used.
#'
#' @param family integer; single number or vector of size \code{m}; defines the
#' bivariate copula family: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param par numeric; single number or vector of size \code{n}; copula
#' parameter.
#' @param par2 numeric; single number or vector of size \code{n}; second
#' parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8,
#' Tawn type 1 and type 2; default: \code{par2 = 0}).  Note that the degrees of
#' freedom parameter of the t-copula does not need to be set, because the
#' theoretical Kendall's tau value of the t-copula is independent of this
#' choice.
#' @param obj \code{BiCop} object containing the family and parameter
#' specification.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are omitted (should only be used with
#' care).
#' @return Theoretical value of Kendall's tau (vector) corresponding to the
#' bivariate copula \code{family} and parameter vector \eqn{(\theta, \delta) =}
#' \code{(par, par2)}.
#' \tabular{ll}{ No. (\code{family}) \tab Kendall's tau (\code{tau}) \cr
#' \code{1, 2} \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
#' \code{3, 13} \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
#' \code{4, 14} \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
#' \code{5} \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta +
#' 4 D_1(\theta)/\theta} \cr
#' \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{x/\theta}{\exp(x)-1}dx}{D_1(\theta)=
#' \int_0^\theta (x/\theta)/(exp(x)-1)dx} (Debye function) \cr
#' \code{6, 16} \tab \eqn{1+\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx}{1+4/\theta^2\int_0^1
#' x\log(x)(1-x)^{2(1-\theta)/\theta}dx} \cr
#' \code{7, 17} \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
#' \code{8, 18} \tab \eqn{1+4\int_0^1 -\log(-(1-t)^\theta+1)
#' (1-t-(1-t)^{-\theta}+(1-t)^{-\theta}t)/(\delta\theta) dt} \cr
#' \code{9, 19} \tab \eqn{1+4\int_0^1 ( (1-(1-t)^{\theta})^{-\delta} - 1)
#' /( -\theta\delta(1-t)^{\theta-1}(1-(1-t)^{\theta})^{-\delta-1} ) dt} \cr
#' \code{10, 20} \tab \eqn{1+4\int_0^1
#' -\log \left(((1-t\delta)^\theta-1)/((1-\delta)^\theta-1) \right) } \cr
#' \tab \eqn{* (1-t\delta-(1-t\delta)^{-\theta}+(1-t\delta)^{-\theta}t\delta)/(\theta\delta) dt} \cr
#' \code{23, 33} \tab \eqn{\frac{\theta}{2-\theta}}{\theta/(2-\theta)} \cr
#' \code{24, 34} \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
#' \code{26, 36} \tab \eqn{-1-\frac{4}{\theta^2}\int_0^1
#' x\log(x)(1-x)^{-2(1+\theta)/\theta}dx}{-1-4/\theta^2
#' \int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx} \cr
#' \code{27, 37} \tab \eqn{-1-\frac{2}{\delta(2-\theta)}}{1-2/(\delta(\theta+2))} \cr
#' \code{28, 38} \tab \eqn{-1-4\int_0^1 -\log(-(1-t)^{-\theta}+1)
#' (1-t-(1-t)^{\theta}+(1-t)^{\theta}t)/(\delta\theta) dt} \cr
#' \code{29, 39} \tab \eqn{-1-4\int_0^1 ( (1-(1-t)^{-\theta})^{\delta} - 1)
#' /( -\theta\delta(1-t)^{-\theta-1}(1-(1-t)^{-\theta})^{\delta-1} ) dt} \cr
#' \code{30, 40} \tab \eqn{-1-4\int_0^1 -\log
#' \left( ((1+t\delta)^{-\theta}-1)/((1+\delta)^{-\theta}-1) \right)} \cr
#' \tab \eqn{* (1+t\delta-(1+t\delta)^{\theta}-(1+t\delta)^{\theta}t\delta)/(\theta\delta) dt} \cr
#' \code{104,114} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{\theta}+t^{\theta}]^{1/\theta}} \cr
#' \code{204,214} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt}  \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#' \code{124,134} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{-\theta}+t^{-\theta}]^{-1/\theta}} \cr
#' \code{224,234} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
#' \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
#'
#' }
#'
#' @note The number \code{n} can be chosen arbitrarily, but must agree across
#' arguments.
#'
#' @author Ulf Schepsmeier, Tobias Erhardt
#'
#' @seealso \code{\link{BiCopTau2Par}}, \code{\link{BiCop}}
#'
#' @references Joe, H. (1997). Multivariate Models and Dependence Concepts.
#' Chapman and Hall, London.
#'
#' Czado, C., U. Schepsmeier, and A. Min (2012). Maximum likelihood estimation
#' of mixed C-vines with application to exchange rates. Statistical Modelling,
#' 12(3), 229-255.
#'
#' @examples
#' ## Example 1: Gaussian copula
#' tau0 <- 0.5
#' rho <- BiCopTau2Par(family = 1, tau = tau0)
#' # transform back
#' tau <- BiCopPar2Tau(family = 1, par = rho)
#' tau - 2/pi*asin(rho)
#'
#' ## Example 2:
#' vpar <- seq(from = 1.1, to = 10, length.out = 100)
#' tauC <- BiCopPar2Tau(family = 3, par = vpar)
#' tauG <- BiCopPar2Tau(family = 4, par = vpar)
#' tauF <- BiCopPar2Tau(family = 5, par = vpar)
#' tauJ <- BiCopPar2Tau(family = 6, par = vpar)
#' plot(tauC ~ vpar, type = "l", ylim = c(0,1))
#' lines(tauG ~ vpar, col = 2)
#' lines(tauF ~ vpar, col = 3)
#' lines(tauJ ~ vpar, col = 4)
#'
#' ## Example 3: different copula families
#' theta <- BiCopTau2Par(family = c(3,4,6), tau = c(0.4, 0.5, 0.6))
#' BiCopPar2Tau(family = c(3,4,6), par = theta)
#'
#' \dontshow{
#' # Test BiCopPar2Tau (one parametric families)
#' theta <- BiCopTau2Par(family = 0, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 0, par = theta)
#' theta <- BiCopTau2Par(family = 1, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 1, par = theta)
#' theta <- BiCopTau2Par(family = 3, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 3, par = theta)
#' theta <- BiCopTau2Par(family = 4, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 4, par = theta)
#' theta <- BiCopTau2Par(family = 5, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 5, par = theta)
#' theta <- BiCopTau2Par(family = 6, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 6, par = theta)
#' theta <- BiCopTau2Par(family = 13, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 13, par = theta)
#' theta <- BiCopTau2Par(family = 14, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 14, par = theta)
#' theta <- BiCopTau2Par(family = 16, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 16, par = theta)
#' theta <- BiCopTau2Par(family = 23, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 23, par = theta)
#' theta <- BiCopTau2Par(family = 24, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 24, par = theta)
#' theta <- BiCopTau2Par(family = 26, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 26, par = theta)
#' theta <- BiCopTau2Par(family = 33, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 33, par = theta)
#' theta <- BiCopTau2Par(family = 34, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 34, par = theta)
#' theta <- BiCopTau2Par(family = 36, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 36, par = theta)
#' theta <- BiCopTau2Par(family = 41, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 41, par = theta)
#' theta <- BiCopTau2Par(family = 51, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 51, par = theta)
#' theta <- BiCopTau2Par(family = 61, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 61, par = theta)
#' theta <- BiCopTau2Par(family = 71, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 71, par = theta)
#' theta <- BiCopTau2Par(family = 41, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 41, par = theta)
#' theta <- BiCopTau2Par(family = 51, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 51, par = theta)
#' theta <- BiCopTau2Par(family = 61, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 61, par = theta)
#' theta <- BiCopTau2Par(family = 71, tau = -c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 71, par = theta)
#'
#' # Test BiCopPar2Tau (two parametric families)
#' theta <- BiCopTau2Par(family = 2, tau = c(0.4,0.5,0.6))
#' BiCopPar2Tau(family = 2, par = theta)
#' theta <- 1:3
#' delta <- 1:3
#' BiCopPar2Tau(family = 7, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 17, par = theta, par2 = delta)
#' theta <- -(1:3)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 27, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 37, par = theta, par2 = delta)
#' theta <- 2:4
#' delta <- 1:3
#' BiCopPar2Tau(family = 8, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 18, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 28, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 38, par = theta, par2 = delta)
#' theta <- 1:3
#' delta <- 1:3
#' BiCopPar2Tau(family = 9, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 19, par = theta, par2 = delta)
#' theta <- -(1:3)
#' delta <- -(1:3)
#' BiCopPar2Tau(family = 29, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 39, par = theta, par2 = delta)
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 10, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 20, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- -c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 30, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 40, par = theta, par2 = delta)
#'
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 104, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 114, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 124, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 134, par = theta, par2 = delta)
#'
#' theta <- 2:4
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 204, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 214, par = theta, par2 = delta)
#' theta <- -(2:4)
#' delta <- c(0.1, 0.5, 0.9)
#' BiCopPar2Tau(family = 224, par = theta, par2 = delta)
#' BiCopPar2Tau(family = 234, par = theta, par2 = delta)
#' }
#'
BiCopPar2Tau <- function(family, par, par2 = 0, obj = NULL, check.pars = TRUE) {
    # fix for SemiParBIVProbit package
    dims <- set_dims(family, par, par2)
    # set arbitrary par2 for t-copula
    if (class(family) != "BiCop")
        par2[family == 2] <- par2[family == 2] + 4
    ## preprocessing of arguments
    args <- preproc(c(as.list(environment()), call = match.call()),
                    extract_from_BiCop,
                    match_spec_lengths,
                    check_fam_par)
    list2env(args, environment())

    ## calculate Kendall's tau
    out <- vapply(1:length(par),
                  function(i) calcTau(family[i], par[i], par2[i]),
                  numeric(1))

    ## return result
    if (length(dims) > 1)
        out <- array(out, dim = dims)
    out
}



# frankParGrid <- seq(-101, 101, l = 100)
frankParGrid <- c(-101, -98.959595959596, -96.9191919191919, -94.8787878787879,
                  -92.8383838383838, -90.7979797979798, -88.7575757575758, -86.7171717171717,
                  -84.6767676767677, -82.6363636363636, -80.5959595959596, -78.5555555555556,
                  -76.5151515151515, -74.4747474747475, -72.4343434343434, -70.3939393939394,
                  -68.3535353535354, -66.3131313131313, -64.2727272727273, -62.2323232323232,
                  -60.1919191919192, -58.1515151515152, -56.1111111111111, -54.0707070707071,
                  -52.030303030303, -49.989898989899, -47.949494949495, -45.9090909090909,
                  -43.8686868686869, -41.8282828282828, -39.7878787878788, -37.7474747474748,
                  -35.7070707070707, -33.6666666666667, -31.6262626262626, -29.5858585858586,
                  -27.5454545454545, -25.5050505050505, -23.4646464646465, -21.4242424242424,
                  -19.3838383838384, -17.3434343434344, -15.3030303030303, -13.2626262626263,
                  -11.2222222222222, -9.18181818181819, -7.14141414141415, -5.1010101010101,
                  -3.06060606060606, -1.02020202020203, 1.02020202020201, 3.06060606060605,
                  5.10101010101009, 7.14141414141413, 9.18181818181817, 11.2222222222222,
                  13.2626262626263, 15.3030303030303, 17.3434343434343, 19.3838383838384,
                  21.4242424242424, 23.4646464646464, 25.5050505050505, 27.5454545454545,
                  29.5858585858586, 31.6262626262626, 33.6666666666667, 35.7070707070707,
                  37.7474747474747, 39.7878787878788, 41.8282828282828, 43.8686868686869,
                  45.9090909090909, 47.9494949494949, 49.989898989899, 52.030303030303,
                  54.070707070707, 56.1111111111111, 58.1515151515151, 60.1919191919192,
                  62.2323232323232, 64.2727272727273, 66.3131313131313, 68.3535353535353,
                  70.3939393939394, 72.4343434343434, 74.4747474747475, 76.5151515151515,
                  78.5555555555555, 80.5959595959596, 82.6363636363636, 84.6767676767677,
                  86.7171717171717, 88.7575757575758, 90.7979797979798, 92.8383838383838,
                  94.8787878787879, 96.9191919191919, 98.9595959595959, 101)
# frankTauVals <- 1 - 4/frankParGrid + 4/frankParGrid * copula::debye1(frankParGrid)
frankTauVals <- c(-0.961041048550867, -0.960251344564296, -0.95942897342536,
                 -0.958571866033333, -0.957677774843704, -0.956744254215267, -0.955768638102463,
                 -0.954748014665184, -0.953679197287614, -0.952558691399602, -0.951382656374172,
                 -0.950146861627529, -0.948846635866317, -0.947476808201624, -0.946031639568512,
                 -0.944504742537972, -0.942888987164558, -0.941176389950323, -0.93935798228736,
                 -0.937423653818046, -0.935361964957013, -0.933159921260012, -0.930802700275106,
                 -0.928273318793382, -0.925552224778557, -0.922616793339209, -0.919440699396042,
                 -0.91599313043174, -0.912237789768897, -0.908131622511068, -0.903623170019409,
                 -0.898650420568253, -0.893137967278655, -0.886993199334039, -0.880101121990825,
                 -0.872317196660623, -0.863457265500544, -0.853283089132126, -0.84148112450054,
                 -0.827630609937665, -0.811154246217934, -0.791239672221123, -0.766710388198721,
                 -0.73580674780385, -0.69580488947752, -0.642352809081713, -0.568396118348008,
                 -0.462982587231658, -0.312511738501557, -0.112196429394003, 0.112196429394,
                 0.312511738501556, 0.462982587231657, 0.568396118348008, 0.642352809081712,
                 0.695804889477519, 0.73580674780385, 0.76671038819872, 0.791239672221123,
                 0.811154246217934, 0.827630609937665, 0.84148112450054, 0.853283089132126,
                 0.863457265500544, 0.872317196660623, 0.880101121990825, 0.886993199334039,
                 0.893137967278656, 0.898650420568252, 0.903623170019409, 0.908131622511068,
                 0.912237789768897, 0.91599313043174, 0.919440699396043, 0.92261679333921,
                 0.925552224778557, 0.928273318793382, 0.930802700275106, 0.933159921260012,
                 0.935361964957013, 0.937423653818045, 0.93935798228736, 0.941176389950323,
                 0.942888987164558, 0.944504742537972, 0.946031639568512, 0.947476808201624,
                 0.948846635866317, 0.950146861627529, 0.951382656374172, 0.952558691399602,
                 0.953679197287614, 0.954748014665184, 0.955768638102463, 0.956744254215266,
                 0.957677774843704, 0.958571866033333, 0.95942897342536, 0.960251344564296,
                 0.961041048550867)
frankTau <- function(par) {
    approx(x = frankParGrid, y = frankTauVals, xout = par)$y
}

calcTau <- function(family, par, par2) {
    ## calculation of tau(s) depending on pair-copula family
    if (family == 0) {
        tau <- rep(0, times = length(par))
    } else if (family == 1 | family == 2) {
        tau <- 2/pi * asin(par)
    } else if (family == 3 || family == 13) {
        tau <- par/(par + 2)
    } else if (family == 4 || family == 14) {
        tau <- 1 - 1/par
    } else if (family == 5) {
        tau <- frankTau(par)
    } else if (family == 6 || family == 16) {
        # tau = 1 + 4/par^2 * integrate(function(x) log(x)*x*(1-x)^(2*(1-par)/par), 0,
        # 1)$value
        param1 <- 2/par + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - par)
        tau[par == 2] <- 1 - trigamma(2)
    } else if (family == 7 || family == 17) {
        theta <- par
        delta <- par2
        tau <- 1 - 2/(delta * (theta + 2))
    } else if (family == 8 || family == 18) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 9 || family == 19) {
        theta <- par
        delta <- par2

        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^-de - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t)  kt(t, th = theta, de = delta), 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 10 || family == 20) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 23 || family == 33) {
        tau <- par/(-par + 2)
    } else if (family == 24 || family == 34) {
        tau <- -1 - 1/par
    } else if (family == 26 || family == 36) {
        theta <- -par
        param1 <- 2/theta + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - theta)
        tau[theta == 2] <- 1 - trigamma(2)
        tau <- -tau
    } else if (family == 27 || family == 37) {
        theta <- -par
        delta <- -par2
        tau <- 1 - 2/(delta * (theta + 2))
        tau <- -tau
    } else if (family == 28 || family == 38) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 29 || family == 39) {
        theta <- -par
        delta <- -par2

        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^(-de) - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 30 || family == 40) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- try(1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    } else if (family == 41 || family == 51) {
        de <- par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
    } else if (family == 61 || family == 71) {
        de <- -par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
        tau <- -tau
    } else if (family == 42) {
        tau <- (75 * par2 - par2^2 + par * (25 - par2))/450
    } else if (family == 104 || family == 114 || family == 204 || family == 214) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2",
                        as.double(t),
                        as.integer(length(t)),
                        as.double(th),
                        as.double(de),
                        as.double(1),
                        as.double(rep(0, length(t))),
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn",
                              as.double(t),
                              as.integer(length(t)),
                              as.double(th),
                              as.double(de),
                              as.double(1),
                              as.double(rep(0, length(t))),
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- try(mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
    } else if (family == 124 || family == 134 || family == 224 || family == 234) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2",
                        as.double(t),
                        as.integer(length(t)),
                        as.double(-th),
                        as.double(de),
                        as.double(1),
                        as.double(rep(0, length(t))),
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn",
                              as.double(t),
                              as.integer(length(t)),
                              as.double(-th),
                              as.double(de),
                              as.double(1),
                              as.double(rep(0, length(t))),
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- try(mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2), silent = TRUE)
        if (inherits(tau, "try-error"))
            tau <- NA
        tau <- -tau
    }

    ## return result
    tau
}
