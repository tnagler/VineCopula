##' Testsuite - Run
##'
##' Run several tests for the BiCop-functions of the VineCopula-package
##'
##' @author Dr. Ulf Schepsmeier
##' @param FUN function name
##' @return results list of results for each family


## testRun for BiCopPar2Tau, BiCopPar2Beta
## BiCopPar2TailDep geht so leider noch nicht, da lower und upper als return

testRunBiCopPar <- function(FUN){
  ## familyset
  familyset <- c(1:10,13:20,23:30,33:40,104,114,124,134,204,214,224,234)
  #familyset <- c(1:10,13:20,23:30)
  familyset <- familyset[-which(familyset %in% c(15,25,35))]

  if(FUN == "BiCopPar2Beta") familyset <- familyset[-which(familyset == 2)]

  ## parameter sets
  parset3 <- seq(0, 0.999, 0.001)
  parset3a <- seq(1, 1.999, 0.001)
  parset1 <- c(parset3, seq(1, 10, 0.01))
  parset2 <- c(parset3a, seq(2, 10, 0.01))
  parset4 <- seq(0, 50, 1)

  ## return the results in a list
  results <- list()

  k <- 1
  for(fam in familyset){  # run over all families
    ## set the correct parameter set
    if(fam == 1){
      res <- rep(0, length(parset3))
      par <- parset3
    } else if(fam == 2){
      res <- matrix(0, length(parset3), length(parset4))
      par <- parset3
      par2 <- parset4
    } else if(fam %in% c(3, 13, 23, 33)){
      res <- rep(0,length(parset1)-1)
      par <- parset1[-1]
    } else if(fam %in% c(4, 14, 24, 34)){
      res <- rep(0,length(parset2))
      par <- parset2
    } else if(fam == 5){
      res <- rep(0,length(parset1)-1)
      par <- parset1[-1]
    } else if(fam %in% c(6, 16, 26, 36)){
      res <- rep(0,length(parset2)-1)
      par <- parset2[-1]
    } else if(fam %in% c(7, 17, 27, 37, 8, 18, 28, 38)){
      res <- matrix(0, length(parset1)-1, length(parset2))
      par <- parset1[-1]
      par2 <- parset2
    } else if(fam %in% c(9, 19, 29, 39)){
      res <- matrix(0, length(parset2), length(parset1)-1)
      par <- parset2
      par2 <- parset1[-1]
    } else if(fam %in% c(10, 20, 30, 40)){
      res <- matrix(0, length(parset2), length(parset3)-1)
      par <- parset2
      par2 <- parset3[-1]
    } else if(fam > 100){
      res <- matrix(0, length(parset2), length(parset3))
      par <- parset2
      par2 <- parset3
    }

    ## length of results (depending on the parameter set)
    n1 <- ifelse(is.null(dim(res)), length(res), nrow(res))
    n2 <- ifelse(is.null(dim(res)), 0, ncol(res))

    ## for BiCopPar2TailDep double res
    if(FUN == "BiCopPar2TailDep"){
      if(n2 == 0){
        res <- c(res,res)
      } else {
        res <- rbind(res, res)
      }
    }

    ## for rotated copulas switch sign
    if(fam > 20 && fam < 100){
      par <- -par
      par2 <- -par2
    } else if(fam %in% c(124,134,224,234)){
      par <- -par
    }

    for(i in 1:n1){
      if(n2 == 0){
        if(FUN == "BiCopPar2TailDep"){
          tmp <- do.call(what=FUN, args=list(family=fam, par=par[i], par2=0))
          res[i] <- tmp$lower
          res[i+n1] <- tmp$upper
        } else {
          res[i] <- do.call(what=FUN, args=list(family=fam, par=par[i], par2=0))
        }
      } else {
        for(j in n2){
          if(FUN == "BiCopPar2TailDep"){
            tmp <- do.call(what=FUN, args=list(family=fam, par=par[i], par2=par2[j]))
            res[i,j] <- tmp$lower
            res[(i+n1),j] <- tmp$upper
          } else {
            res[i,j] <- do.call(what=FUN, args=list(family=fam, par=par[i], par2=par2[j]))
          }
        }
      }
    }

    ## save the results and give it the name of teh family
    results[[k]] <- res
    names(results)[[k]] <- as.character(fam)

    k <- k+1

  } # end familyset

  return(results)
}



## test for BiCopTau2Par
testRunBiCopTau <- function(FUN){
  ## familyset
  familyset <- c(1:6, 13,14,16,23,24,26,33,34,36)

  ## tau
  tauset <- seq(0.001, 0.999, 0.001)
  ntauset <- -tauset

  ## return the results in a list
  results <- list()

  k <- 1
  for(fam in familyset){  # run over all families
    if(fam %in% c(1,2,5)){
      tau <- c(ntauset[length(ntauset):1], tauset)
    } else if(fam %in% c(3,4,6,13,14,16)){
      tau <- tauset
    } else {
      tau <- ntauset[length(ntauset):1]
    }

    res <- do.call(what=FUN, args=list(family=fam, tau=tau))  # vectorized function

    ## save the results and give it the name of teh family
    results[[k]] <- res
    names(results)[[k]] <- as.character(fam)

    k <- k+1

  } # end familyset

  return(results)
}



testRunBiCop <- function(FUN){
  ## familyset
  #familyset <- c(1:10,13:20,23:30,33:40,104,114,124,134,204,214,224,234)
  familyset <- c(1:10, 104, 204)
  #familyset <- familyset[-which(familyset %in% c(15,25,35))]

  if(FUN == "BiCopCDF") familyset <- familyset[-which(familyset == 2)]

  ## parameter sets
  parset3 <- seq(0, 0.99, 0.01)
  parset3a <- seq(1, 1.99, 0.01)
  parset1 <- c(parset3, seq(1, 10, 0.25))
  parset2 <- c(parset3a, seq(2, 10, 0.25))
  parset4 <- seq(2, 30, 2)

  ## copula data
  u1 <- c(seq(0.001,0.01,0.002), seq(0.01,0.99,0.02), seq(0.99,0.999,0.002))
  u2 <- u1

  ## return the results in a list
  results <- list()

  k <- 1
  for(fam in familyset){  # run over all families
    ## set the correct parameter set
    if(fam == 1){
      res <- array(0, dim=c(length(parset3), length(u1), length(u2)))
      par <- parset3
    } else if(fam == 2){
      res <- array(0, dim=c(length(parset3), length(parset4), length(u1), length(u2)))
      par <- parset3
      par2 <- parset4
    } else if(fam %in% c(3, 13, 23, 33)){
      res <- array(0, dim=c(length(parset1)-1, length(u1), length(u2)))
      par <- parset1[-1]
    } else if(fam %in% c(4, 14, 24, 34)){
      res <- array(0, dim=c(length(parset2), length(u1), length(u2)))
      par <- parset2
    } else if(fam == 5){
      res <- array(0, dim=c(length(parset1)-1, length(u1), length(u2)))
      par <- parset1[-1]
    } else if(fam %in% c(6, 16, 26, 36)){
      res <- array(0, dim=c(length(parset2)-1, length(u1), length(u2)))
      par <- parset2[-1]
    } else if(fam %in% c(7, 17, 27, 37, 8, 18, 28, 38)){
      res <- array(0, dim=c(length(parset1)-1, length(parset2), length(u1), length(u2)))
      par <- parset1[-1]
      par2 <- parset2
    } else if(fam %in% c(9, 19, 29, 39)){
      res <- array(0, dim=c(length(parset2), length(parset1)-1,length(u1), length(u2)))
      par <- parset2
      par2 <- parset1[-1]
    } else if(fam %in% c(10, 20, 30, 40)){
      res <- array(0, dim=c(length(parset2), length(parset3)-1, length(u1), length(u2)))
      par <- parset2
      par2 <- parset3[-1]
    } else if(fam > 100){
      res <- array(0, dim=c(length(parset2), length(parset3), length(u1), length(u2)))
      par <- parset2
      par2 <- parset3
    }

    ## length of results (depending on the parameter set)
    n1 <- dim(res)[1]
    n2 <- ifelse(length(dim(res))==4, dim(res)[2], 0)

    ## for rotated copulas switch sign
    if(fam > 20 && fam < 100){
      par <- -par
      par2 <- -par2
    } else if(fam %in% c(124,134,224,234)){
      par <- -par
    }

    ## for loops are not the best
    iu <- 1
    for(u in u1){
      #iv <- 1
      #for(v in u2){
      uu <- rep(u,length(u2))
        for(i in 1:n1){
            if(n2 == 0){
              if(FUN %in% c("BiCopHfunc", "BiCopHinv")){
                  ## At the moment just test hfunc1
                  res[i,iu,] <- do.call(what=FUN,
                                        args=list(u1=uu, u2=u2, family=fam, par=par[i],
                                                  par2=0, check.pars=FALSE))[[1]]
              } else {
                res[i,iu,] <- do.call(what=FUN, args=list(u1=uu, u2=u2, family=fam, par=par[i], par2=0, check.pars=FALSE))
              }
            } else {
            for(j in n2){
                if(FUN %in% c("BiCopHfunc", "BiCopHinv")){
                    ## At the moment just test hfunc1
                    res[i,j,iu,] <- do.call(what=FUN,
                                          args=list(u1=uu, u2=u2, family=fam, par=par[i],
                                                    par2=0, check.pars=FALSE))[[1]]
                } else {
                    res[i,j,iu,] <- do.call(what=FUN, args=list(u1=uu, u2=u2, family=fam, par=par[i], par2=par2[j], check.pars=FALSE))
                }
            }
          }
        }
        #iv <- iv + 1
      #}
      iu <- iu + 1
    }



    ## save the results and give it the name of teh family
    results[[k]] <- res
    names(results)[[k]] <- as.character(fam)

    k <- k+1

  } # end familyset

  return(results)

}
