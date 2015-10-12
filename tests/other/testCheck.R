##' Testsuite - Check
##'
##' Run several tests for the BiCop-functions of the VineCopula-package
##'
##' @author Dr. Ulf Schepsmeier
##' @param results list of results returned from testRun*

## sub function
isFiniteCheck <- function(res){
    if(any(!is.finite(res))){
        return(FALSE)
    }else{
        return(TRUE)
    }
}


## This check tests results of testRunBiCopPar
testCheck <- function(results){
  ## length of results
  n <- length(results)

  check <- sapply(results, FUN = isFiniteCheck)

  for(i in 1:n){
    ## Check 4: in range
    if(check[i]){
        if(names(results)[i] %in% c(1:10,13,14,16:20,104,114,204,214)){
          if(any( results[[i]] < 0 || results[[i]] > 1 ) ) check[i] <- FALSE
        } else {
          if(any( results[[i]] > 0 || results[[i]] < -1 ) ) check[i] <- FALSE
        }
    }
    ## check for jumps
    ## TODO

  }

  return(check)
}


## This check tests results of testRunBiCopTau
testCheck2 <- function(results){
  ## length of results
  n <- length(results)

  check <- sapply(results, FUN = isFiniteCheck)

  for(i in 1:n){
    if(check[i]){
        tmp <- VineCopula:::BiCopCheck(family=as.numeric(names(results)[i]),
                                   par=results[[i]], par2=rep(5,length(results[[i]])))
        if(!tmp) check[i] <- FALSE
    }
    ## check for jumps
    ## TODO
  }

  return(check)
}


## This check tests results of testRunBiCop
testCheck3 <- function(results){

  check <- sapply(results, FUN = isFiniteCheck)

  return(check)
}
