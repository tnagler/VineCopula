"CvM" <- function(cdf = NULL) {
    # Cumulative distribution function test: Function that computes the
    # Cramer-von-Mises test statistic
    # --------------------------------------------------------------------------
    # INPUT: cdf CDF for which to compute CvM test OUTPUT: CvM
    # Cramer-von-Mises test statistic
    # --------------------------------------------------------------------------
    # Author: Daniel Berg <daniel at danielberg.no> Date: 11.05.2005
    # Version: 1.0.1
    # --------------------------------------------------------------------------
    n <- length(cdf)
    CvM <- .C("CvMtest", as.double(cdf), as.integer(n), as.double(0))[[3]]
    CvM
}
