#' Copula Data Objects
#'
#' The function `as.copuladata` coerces an object (`data.frame`,
#' `matrix`, `list`) to a `copuladata` object.
#'
#'
#' @aliases as.copuladata as.copuladata.data.frame as.copuladata.matrix
#' as.copuladata.list
#' @param data Either a `data.frame`, a `matrix` or a `list`
#' containing copula data (i.e. data with uniform margins on \eqn{[0,1]}). The
#' `list` elements have to be vectors of identical length.
#' @author Tobias Erhardt
#' @seealso [pobs()], [pairs.copuladata()]
#' @examples
#'
#' data(daxreturns)
#'
#' data <- as.matrix(daxreturns)
#' class(as.copuladata(data))
#'
#' data <- as.data.frame(daxreturns)
#' class(as.copuladata(data))
#'
#' data <- as.list(daxreturns)
#' names(data) <- names(daxreturns)
#' class(as.copuladata(data))
#'
as.copuladata <- function(data) {
    ## generic function for coercion to 'copuladata'
    UseMethod("as.copuladata", data)
}

as.copuladata.data.frame <- function(data) {
    ## coercion of 'data.frame' to 'copuladata'
    if (!all(sapply(data, is.numeric)))
        stop("Data has to be numeric.")
    if (any(data > 1 | data < 0))
        stop("Data has to be in the interval [0,1].")
    class(data) <- append("copuladata", class(data))
    return(data)
}

as.copuladata.matrix <- function(data) {
    ## coercion of 'matrix' to 'copuladata'
    if (!all(is.numeric(data)))
        stop("Data has to be numeric.")
    if (any(data > 1 | data < 0))
        stop("Data has to be in the interval [0,1].")
    data <- data.frame(data)
    class(data) <- append("copuladata", class(data))
    return(data)
}

as.copuladata.list <- function(data) {
    ## coercion of 'list' to 'copuladata'
    if (!all(sapply(data, is.numeric)))
        stop("Data has to be numeric.")
    if (any(sapply(data, length) != length(data[[1]])))
        stop("All list entries have to be of same length.")
    data <- data.frame(data)
    if (any(data > 1 | data < 0))
        stop("Data has to be in the interval [0,1].")
    class(data) <- append("copuladata", class(data))
    return(data)
}
