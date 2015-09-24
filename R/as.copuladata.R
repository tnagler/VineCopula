#' Copula Data Objects
#' 
#' The function \code{as.copuladata} coerces an object (\code{data.frame},
#' \code{matrix}, \code{list}) to a \code{copuladata} object.
#' 
#' 
#' @aliases as.copuladata as.copuladata.data.frame as.copuladata.matrix
#' as.copuladata.list
#' @param data Either a \code{data.frame}, a \code{matrix} or a \code{list}
#' containing copula data (i.e. data with uniform margins on [0,1]). The
#' \code{list} elements have to be vectors of identical length.
#' @author Tobias Erhardt
#' @seealso \code{\link{pobs}}, \code{\link{pairs.copuladata}}
#' @examples
#' 
#' data(daxreturns)
#' 
#' data <- as(daxreturns, "matrix")
#' class(as.copuladata(data))
#' 
#' data <- as(daxreturns, "data.frame")
#' class(as.copuladata(data))
#' 
#' data <- as(daxreturns, "list")
#' names(data) <- names(daxreturns)
#' class(as.copuladata(data))
#' 
#' @export as.copuladata
as.copuladata <- function(data) {
    ## generic function for coercion to 'copuladata'
    UseMethod("as.copuladata", data)
}

as.copuladata.data.frame <- function(data) {
    ## coercion of 'data.frame' to 'copuladata'
    if (any(sapply(data, mode) != "numeric")) 
        stop("Data has to be numeric.")
    if (any(data > 1 || data < 0)) 
        stop("Data has to be in the interval [0,1].")
    class(data) <- append("copuladata", class(data))
    return(data)
}

as.copuladata.matrix <- function(data) {
    ## coercion of 'matrix' to 'copuladata'
    if (mode(data) != "numeric") 
        stop("Data has to be numeric.")
    if (any(data > 1 || data < 0)) 
        stop("Data has to be in the interval [0,1].")
    data <- data.frame(data)
    class(data) <- append("copuladata", class(data))
    return(data)
}

as.copuladata.list <- function(data) {
    ## coercion of 'list' to 'copuladata'
    if (any(sapply(data, mode) != "numeric")) 
        stop("Data has to be numeric.")
    if (any(sapply(data, length) != length(data[[1]]))) 
        stop("All list entries have to be of same length.")
    data <- data.frame(data)
    if (any(data > 1 || data < 0)) 
        stop("Data has to be in the interval [0,1].")
    class(data) <- append("copuladata", class(data))
    return(data)
}
