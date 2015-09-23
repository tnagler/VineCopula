BiCop <- function(family, par, par2 = 0) {
    ## family/parameter consistency checks
    BiCopCheck(family, par, par2)
    
    ## return BiCop object
    out <- list(family = family, par = par, par2 = par2)
    class(out) <- "BiCop"
    out
}