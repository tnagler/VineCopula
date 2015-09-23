"Fhat" <- function(fj, m = NULL) {
    # Function that computes CDF of input vector fj:  P[ fj  <=  w ]
    #--------------------------------------------------------------------------
    # INPUT:
    #   fj                Data vector for which to compute CDF
    #   m                 Resolution of the cdf, that is at how many points to compute the cdf
    #                      If missing then m is taken to equal n, the length of fj. 
    #                      m must be smaller than n. 
    # OUTPUT:
    #   Fhat           CDF
    #--------------------------------------------------------------------------
    # Author: Henrik Bakken <henrikb at stud.ntnu.no>, Daniel Berg <daniel at danielberg.no>
    # Date: 06 September 2006
    # Version: 1.0.2
    #--------------------------------------------------------------------------
    n <- length(fj)
    if (missing(m) || m > n)  m <- n
    Fhat <- rep(0, m)
    Fhat <- .C("CumDist",
               as.double(fj),
               as.integer(n), 
               as.integer(m), 
               as.double(Fhat),
               PACKAGE = "VineCopula")[[4]]
    Fhat
}
