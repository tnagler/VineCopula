##' Testsuite
##'
##' Tests for the VineCopula package
##'
##' @author Dr. Ulf Schepsmeier
##'

## Main function

library(VineCopula)

source("tests/other/testRun.r")
source("tests/other/testCheck.r")

# BiCopPar2Tau
results_BiCopPar2Tau <- testRunBiCopPar("BiCopPar2Tau")
check_BiCopPar2Tau <- testCheck(results_BiCopPar2Tau)
if(!all(check_BiCopPar2Tau)){
  print(check_BiCopPar2Tau)
} else {
  rm(results_BiCopPar2Tau)
  gc()
}

# BiCopPar2Beta
results_BiCopPar2Beta <- testRunBiCopPar("BiCopPar2Beta")
check_BiCopPar2Beta <- testCheck(results_BiCopPar2Beta)
if(!all(check_BiCopPar2Beta)){
  print(check_BiCopPar2Beta)
} else {
  rm(results_BiCopPar2Beta)
  gc()
}

# BiCopPar2TailDep
results_BiCopPar2TailDep <- testRunBiCopPar("BiCopPar2TailDep")
check_BiCopPar2TailDep <- testCheck(results_BiCopPar2TailDep)
if(!all(check_BiCopPar2TailDep)){
  print(check_BiCopPar2TailDep)
} else {
  rm(results_BiCopPar2TailDep)
  gc()
}


# BiCopTau2Par
results_BiCopTau2Par <- testRunBiCopTau("BiCopTau2Par")
check_BiCopTau2Par <- testCheck2(results_BiCopTau2Par)
if(!all(check_BiCopTau2Par)){
  print(check_BiCopTau2Par)
} else {
  rm(results_BiCopTau2Par)
  gc()
}


# BiCopCDF
results_BiCopCDF <- testRunBiCop("BiCopCDF")
check_BiCopCDF <- testCheck3(results_BiCopCDF)
if(!all(check_BiCopCDF)){
  print(check_BiCopCDF)
} else {
  rm(results_BiCopCDF)
  gc()
}



# BiCopPDF
results_BiCopPDF <- testRunBiCop("BiCopPDF")
check_BiCopPDF <- testCheck3(results_BiCopPDF)
if(!all(check_BiCopPDF)){
    print(check_BiCopPDF)
} else {
    rm(results_BiCopPDF)
    gc()
}



# BiCopHfunc
results_BiCopHfunc <- testRunBiCop("BiCopHfunc")
check_BiCopHfunc <- testCheck3(results_BiCopHfunc)
if(!all(check_BiCopHfunc)){
    print(check_BiCopHfunc)
} else {
    rm(results_BiCopHfunc)
    gc()
}


# BiCopHinv
results_BiCopHinv <- testRunBiCop("BiCopHinv")
check_BiCopHinv <- testCheck3(results_BiCopHinv)
if(!all(check_BiCopHinv)){
    print(check_BiCopHinv)
} else {
    rm(results_BiCopHinv)
    gc()
}
