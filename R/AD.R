"AD" <- function(cdf = NULL) {
    # Cumulative distribution function test: Function that computes the
    # Anderson-Darling test statistic
    # --------------------------------------------------------------------------
    # INPUT: cdf CDF for which to compute AD test OUTPUT:
    # AD Anderson-Darling test statistic
    # --------------------------------------------------------------------------
    # Author: Daniel Berg <daniel at danielberg.no> Date: 27.03.2006 Version: 1.0.1
    # --------------------------------------------------------------------------
    n <- length(cdf)
    AD <- .C("ADtest",
             as.double(cdf),
             as.integer(n),
             as.double(0),
             PACKAGE = "VineCopula")[[3]]
    AD
}
