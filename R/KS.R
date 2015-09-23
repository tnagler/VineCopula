"KS" <- function(cdf = NULL) {
    # Cumulative distribution function test:
    # Function that computes the Kolmogorov-Smirnov test statistic
    #--------------------------------------------------------------------------
    # INPUT:
    #   cdf        CDF for which to compute KS test
    # OUTPUT:
    #   KS         Kolmogorov-Smirnov test statistic
    #--------------------------------------------------------------------------
    # Author: Daniel Berg <daniel at danielberg.no>
    # Date: 11.05.2005
    # Version: 1.0.1
    #-------------------------------------------------------------------------- 
    n <- length(cdf)
    KS <- .C("KStest", 
             as.double(cdf), 
             as.integer(n),
             as.double(0),
             PACKAGE = "VineCopula")[[3]]
    KS
}
