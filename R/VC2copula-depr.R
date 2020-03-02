#' Deprecated
#'
#' This functionality is deprecated in 'VineCopula'. Use the package 'VC2copula'
#' instead.
#' @name VC2copula-deprecated
#' @param family ..
#' @param par ...
#' @param par2 ...
#' @aliases copulaFromFamilyIndex
copulaFromFamilyIndex <- function(family, par, par2 = 0) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases surClaytonCopula surClaytonCopula-class dduCopula,matrix,surClaytonCopula-method
#' dduCopula,numeric,surClaytonCopula-method
#' ddvCopula,matrix,surClaytonCopula-method
#' ddvCopula,numeric,surClaytonCopula-method r90ClaytonCopula-class
#' dduCopula,matrix,r90ClaytonCopula-method
#' dduCopula,numeric,r90ClaytonCopula-method
#' ddvCopula,matrix,r90ClaytonCopula-method
#' ddvCopula,numeric,r90ClaytonCopula-method r270ClaytonCopula-class
#' dduCopula,matrix,r270ClaytonCopula-method
#' dduCopula,numeric,r270ClaytonCopula-method
#' ddvCopula,matrix,r270ClaytonCopula-method
#' ddvCopula,numeric,r270ClaytonCopula-method
setClass("surClaytonCopula", slots = c(param = "numeric"))
setClass("r90ClaytonCopula", slots = c(param = "numeric"))
setClass("r270ClaytonCopula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
surClaytonCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90ClaytonCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270ClaytonCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases surGumbelCopula surGumbelCopula-class dduCopula,matrix,surGumbelCopula-method
#' dduCopula,numeric,surGumbelCopula-method
#' ddvCopula,matrix,surGumbelCopula-method
#' ddvCopula,numeric,surGumbelCopula-method r90GumbelCopula-class
#' dduCopula,matrix,r90GumbelCopula-method
#' dduCopula,numeric,r90GumbelCopula-method
#' ddvCopula,matrix,r90GumbelCopula-method
#' ddvCopula,numeric,r90GumbelCopula-method r270GumbelCopula-class
#' dduCopula,matrix,r270GumbelCopula-method
#' dduCopula,numeric,r270GumbelCopula-method
#' ddvCopula,matrix,r270GumbelCopula-method
#' ddvCopula,numeric,r270GumbelCopula-method
setClass("surGumbelCopula", slots = c(param = "numeric"))
setClass("r90GumbelCopula", slots = c(param = "numeric"))
setClass("r270GumbelCopula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
surGumbelCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90GumbelCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270GumbelCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}


#' @rdname VC2copula-deprecated
#' @aliases joeBiCopula joeBiCopula-class dduCopula,numeric,joeBiCopula-method
#' ddvCopula,numeric,joeBiCopula-method dduCopula,matrix,joeBiCopula-method
#' ddvCopula,matrix,joeBiCopula-method getKendallDistr,joeBiCopula-method
#' kendallDistribution,joeBiCopula-method surJoeBiCopula-class
#' dduCopula,numeric,surJoeBiCopula-method
#' ddvCopula,numeric,surJoeBiCopula-method
#' dduCopula,matrix,surJoeBiCopula-method
#' ddvCopula,matrix,surJoeBiCopula-method r90JoeBiCopula-class
#' dduCopula,numeric,r90JoeBiCopula-method
#' ddvCopula,numeric,r90JoeBiCopula-method
#' dduCopula,matrix,r90JoeBiCopula-method
#' ddvCopula,matrix,r90JoeBiCopula-method r270JoeBiCopula-class
#' dduCopula,numeric,r270JoeBiCopula-method
#' ddvCopula,numeric,r270JoeBiCopula-method
#' dduCopula,matrix,r270JoeBiCopula-method
#' ddvCopula,matrix,r270JoeBiCopula-method
setClass("joeBiCopula", slots = c(param = "numeric"))
setClass("surJoeBiCopula", slots = c(param = "numeric"))
setClass("r90JoeBiCopula", slots = c(param = "numeric"))
setClass("r270JoeBiCopula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
joeBiCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
surJoeBiCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90JoeBiCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270JoeBiCopula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases BB1Copula BB1Copula-class dduCopula,numeric,BB1Copula-method
#' ddvCopula,numeric,BB1Copula-method dduCopula,matrix,BB1Copula-method
#' ddvCopula,matrix,BB1Copula-method getKendallDistr,BB1Copula-method
#' kendallDistribution,BB1Copula-method surBB1Copula-class
#' dduCopula,numeric,surBB1Copula-method ddvCopula,numeric,surBB1Copula-method
#' dduCopula,matrix,surBB1Copula-method ddvCopula,matrix,surBB1Copula-method
#' r90BB1Copula-class dduCopula,numeric,r90BB1Copula-method
#' ddvCopula,numeric,r90BB1Copula-method dduCopula,matrix,r90BB1Copula-method
#' ddvCopula,matrix,r90BB1Copula-method r270BB1Copula-class
#' dduCopula,numeric,r270BB1Copula-method
#' ddvCopula,numeric,r270BB1Copula-method dduCopula,matrix,r270BB1Copula-method
#' ddvCopula,matrix,r270BB1Copula-method
setClass("BB1Copula", slots = c(param = "numeric"))
setClass("surBB1Copula", slots = c(param = "numeric"))
setClass("r90BB1Copula", slots = c(param = "numeric"))
setClass("r270BB1Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
BB1Copula <- function(param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
surBB1Copula <- function(param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r90BB1Copula <- function(param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r270BB1Copula <- function(param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}


#' @rdname VC2copula-deprecated
#' @aliases BB6Copula BB6Copula-class dduCopula,numeric,BB6Copula-method
#' ddvCopula,numeric,BB6Copula-method dduCopula,matrix,BB6Copula-method
#' ddvCopula,matrix,BB6Copula-method getKendallDistr,BB6Copula-method
#' kendallDistribution,BB6Copula-method surBB6Copula-class
#' dduCopula,numeric,surBB6Copula-method ddvCopula,numeric,surBB6Copula-method
#' dduCopula,matrix,surBB6Copula-method ddvCopula,matrix,surBB6Copula-method
#' r90BB6Copula-class dduCopula,numeric,r90BB6Copula-method
#' ddvCopula,numeric,r90BB6Copula-method dduCopula,matrix,r90BB6Copula-method
#' ddvCopula,matrix,r90BB6Copula-method r270BB6Copula-class
#' dduCopula,numeric,r270BB6Copula-method
#' ddvCopula,numeric,r270BB6Copula-method dduCopula,matrix,r270BB6Copula-method
#' ddvCopula,matrix,r270BB6Copula-method
setClass("BB6Copula", slots = c(param = "numeric"))
setClass("surBB6Copula", slots = c(param = "numeric"))
setClass("r90BB6Copula", slots = c(param = "numeric"))
setClass("r270BB6Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
BB6Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
surBB6Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r90BB6Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r270BB6Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases BB7Copula BB7Copula-class dduCopula,numeric,BB7Copula-method
#' ddvCopula,numeric,BB7Copula-method dduCopula,matrix,BB7Copula-method
#' ddvCopula,matrix,BB7Copula-method getKendallDistr,BB7Copula-method
#' kendallDistribution,BB7Copula-method surBB7Copula-class
#' dduCopula,numeric,surBB7Copula-method ddvCopula,numeric,surBB7Copula-method
#' dduCopula,matrix,surBB7Copula-method ddvCopula,matrix,surBB7Copula-method
#' r90BB7Copula-class dduCopula,numeric,r90BB7Copula-method
#' ddvCopula,numeric,r90BB7Copula-method dduCopula,matrix,r90BB7Copula-method
#' ddvCopula,matrix,r90BB7Copula-method r270BB7Copula-class
#' dduCopula,numeric,r270BB7Copula-method
#' ddvCopula,numeric,r270BB7Copula-method dduCopula,matrix,r270BB7Copula-method
#' ddvCopula,matrix,r270BB7Copula-method
setClass("BB7Copula", slots = c(param = "numeric"))
setClass("surBB7Copula", slots = c(param = "numeric"))
setClass("r90BB7Copula", slots = c(param = "numeric"))
setClass("r270BB7Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
BB7Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
surBB7Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r90BB7Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}
#' @rdname VC2copula-deprecated
r270BB7Copula <- function (param=c(1,1)) {
  stop("Function is deprecated in package VineCopula. ",
       "Use the 'VC2copula' package instead.")
}


#' @rdname VC2copula-deprecated
#' @aliases BBB8Copula B8Copula-class dduCopula,numeric,BB8Copula-method
#' ddvCopula,numeric,BB8Copula-method dduCopula,matrix,BB8Copula-method
#' ddvCopula,matrix,BB8Copula-method getKendallDistr,BB8Copula-method
#' kendallDistribution,BB8Copula-method surBB8Copula-class
#' dduCopula,numeric,surBB8Copula-method ddvCopula,numeric,surBB8Copula-method
#' dduCopula,matrix,surBB8Copula-method ddvCopula,matrix,surBB8Copula-method
#' r90BB8Copula-class dduCopula,numeric,r90BB8Copula-method
#' ddvCopula,numeric,r90BB8Copula-method dduCopula,matrix,r90BB8Copula-method
#' ddvCopula,matrix,r90BB8Copula-method r270BB8Copula-class
#' dduCopula,numeric,r270BB8Copula-method
#' ddvCopula,numeric,r270BB8Copula-method dduCopula,matrix,r270BB8Copula-method
#' ddvCopula,matrix,r270BB8Copula-method
setClass("BB8Copula", slots = c(param = "numeric"))
setClass("surBB8Copula", slots = c(param = "numeric"))
setClass("r90BB8Copula", slots = c(param = "numeric"))
setClass("r270BB8Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
BB8Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
surBB8Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90BB8Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270BB8Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases tawnT1Copula tawnT1Copula-class dduCopula,matrix,tawnT1Copula-method
#' dduCopula,numeric,tawnT1Copula-method ddvCopula,matrix,tawnT1Copula-method
#' ddvCopula,numeric,tawnT1Copula-method surTawnT1Copula-class
#' dduCopula,matrix,surTawnT1Copula-method
#' dduCopula,numeric,surTawnT1Copula-method
#' ddvCopula,matrix,surTawnT1Copula-method
#' ddvCopula,numeric,surTawnT1Copula-method r90TawnT1Copula-class
#' dduCopula,matrix,r90TawnT1Copula-method
#' dduCopula,numeric,r90TawnT1Copula-method
#' ddvCopula,matrix,r90TawnT1Copula-method
#' ddvCopula,numeric,r90TawnT1Copula-method r270TawnT1Copula-class
#' dduCopula,matrix,r270TawnT1Copula-method
#' dduCopula,numeric,r270TawnT1Copula-method
#' ddvCopula,matrix,r270TawnT1Copula-method
#' ddvCopula,numeric,r270TawnT1Copula-method
setClass("tawnT1Copula", slots = c(param = "numeric"))
setClass("surTawnT1Copula", slots = c(param = "numeric"))
setClass("r90TawnT1Copula", slots = c(param = "numeric"))
setClass("r270TawnT1Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
tawnT1Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
surTawnT1Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90TawnT1Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270TawnT1Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

#' @rdname VC2copula-deprecated
#' @aliases tawnT2Copula tawnT2Copula-class dduCopula,matrix,tawnT2Copula-method
#' dduCopula,numeric,tawnT2Copula-method ddvCopula,matrix,tawnT2Copula-method
#' ddvCopula,numeric,tawnT2Copula-method surTawnT2Copula-class
#' dduCopula,matrix,surTawnT2Copula-method
#' dduCopula,numeric,surTawnT2Copula-method
#' ddvCopula,matrix,surTawnT2Copula-method
#' ddvCopula,numeric,surTawnT2Copula-method r90TawnT2Copula-class
#' dduCopula,matrix,r90TawnT2Copula-method
#' dduCopula,numeric,r90TawnT2Copula-method
#' ddvCopula,matrix,r90TawnT2Copula-method
#' ddvCopula,numeric,r90TawnT2Copula-method r270TawnT2Copula-class
#' dduCopula,matrix,r270TawnT2Copula-method
#' dduCopula,numeric,r270TawnT2Copula-method
#' ddvCopula,matrix,r270TawnT2Copula-method
#' ddvCopula,numeric,r270TawnT2Copula-method
setClass("tawnT2Copula", slots = c(param = "numeric"))
setClass("surTawnT2Copula", slots = c(param = "numeric"))
setClass("r90TawnT2Copula", slots = c(param = "numeric"))
setClass("r270TawnT2Copula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
tawnT2Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
surTawnT2Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r90TawnT2Copula <- function (param=c(1,1)) {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
#' @rdname VC2copula-deprecated
r270TawnT2Copula <- function (param=c(1,1)) {
  stop("This functionality is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}

setClass("vineCopula", slots = c(param = "numeric"))
#' @rdname VC2copula-deprecated
#' @aliases vineCopula-class vineCopula fitCopula vineCopula-method
#' @param RVM ...
#' @param type ...
#' @param param ...
vineCopula <- function(RVM, type="CVine") {
  stop("This functionality  is deprecated in 'VineCopula'. ",
       "Use the package 'VC2copula' instead.")
}
