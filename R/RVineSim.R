#' Simulation from an R-Vine Copula Model
#'
#' This function simulates from a given R-vine copula model.
#'
#'
#' @param N Number of d-dimensional observations to simulate.
#' @param RVM An \code{\link{RVineMatrix}} object containing the information of
#' the R-vine copula model.
#' @param U If not \code{\link{NULL}}, an (N,d)-matrix of U[0,1] random
#' variates to be transformed to the copula sample.
#' @return An \code{N} x d matrix of data simulated from the given R-vine
#' copula model.
#' @author Jeffrey Dissmann
#' @seealso \code{\link{RVineMatrix}}, \code{\link{BiCopSim}}
#' @references Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka
#' (2013). Selecting and estimating regular vine copulae and application to
#' financial returns. Computational Statistics & Data Analysis, 59 (1), 52-69.
#' @examples
#'
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#'
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#'
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#'
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' # simulate a sample of size 300 from the R-vine copula model
#' set.seed(123)
#' simdata <- RVineSim(300, RVM)
#'
RVineSim <- function(N, RVM, U = NULL) {

    ## sanity checks
    stopifnot(N >= 1)
    if (!is(RVM, "RVineMatrix"))
        stop("'RVM' has to be an RVineMatrix object.")

    ## reorder matrix and U (if provided)
    n <- dim(RVM)
    o <- diag(RVM$Matrix)
    RVM <- normalizeRVineMatrix(RVM)
    takeU <- !is.null(U)
    if (takeU) {
        if (!is.matrix(U))
            U <- rbind(U, deparse.level = 0L)
        if ((d <- ncol(U)) < 2)
            stop("U should be at least bivariate")  # should be an (N, n) matrix
        U <- U[, rev(o)]
    }

    ## create objects for C-call
    matri <- as.vector(RVM$Matrix)
    w1 <- as.vector(RVM$family)
    th <- as.vector(RVM$par)
    th2 <- as.vector(RVM$par2)
    maxmat <- as.vector(RVM$MaxMat)
    conindirect <- as.vector(RVM$CondDistr$indirect)
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0
    tmp <- rep(0, n * N)

    ## simulate R-Vine
    tmp <- .C("SimulateRVine",
              as.integer(N),
              as.integer(n),
              as.integer(w1),
              as.integer(maxmat),
              as.integer(matri),
              as.integer(conindirect),
              as.double(th),
              as.double(th2),
              as.double(tmp),
              as.double(U),
              as.integer(takeU),
              PACKAGE = "VineCopula")[[9]]

    ## store results, bring back to initial order and return
    out <- matrix(tmp, ncol = n, byrow = TRUE)
    if (!is.null(RVM$names)) {
        colnames(out) <- RVM$names
    }
    out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    return(out)
}


#RVineSim <- function(N, RVM, U = NULL) {
RVineSim2 <- function(N, RVM, U = NULL) {

    if (is(RVM, "RVineMatrix")) {
        RVM <- list(RVM)
    }

    ## sanity checks
    stopifnot(N >= 1)
    if (length(RVM) == 0)
        stop("'RVM' is an empty list")
    #if (!is(RVM, "RVineMatrix"))
    if (!all(vapply(RVM, is, F, "RVineMatrix")))
        #stop("'RVM' has to be an RVineMatrix object.")
        stop("'RVM' has to be an RVineMatrix object or a list thereof.")

    . <- dim(RVM[[1]])
    if (any(vapply(RVM, function(RVM) dim(RVM) != ., F)))
        stop("'RVM' is a list of RVineMatrix objects of different dimensions")

    . <- RVM[[1]]$Matrix
    if (!all(vapply(RVM, function(RVM) identical(RVM$Matrix, .), F)))
        stop("'RVM' is a list of RVineMatrix objects with different structures")

    ## reorder matrix and U (if provided)
    #n <- dim(RVM)
    n <- dim(RVM[[1]])
    #o <- diag(RVM$Matrix)
    o <- diag(RVM[[1]]$Matrix)
    #RVM <- normalizeRVineMatrix(RVM)
    RVM[[1]] <- normalizeRVineMatrix(RVM[[1]])  # WARNING creates dependency on normalizeRVineMatrix
    takeU <- !is.null(U)
    if (takeU) {
        if (!is.matrix(U))
            U <- rbind(U, deparse.level = 0L)
        if ((d <- ncol(U)) < 2)
            stop("U should be at least bivariate")  # should be an (N, n) matrix
        U <- U[, rev(o)]
    }

    ## create objects for C-call
    . <- function(name) {
        . <- vapply(RVM, function(RVM) RVM[[name]], RVM[[1]][[name]])  # extract
        array(., dim=c(dim(.)[1:2], N))  # recycle
    }
    #matri <- as.vector(RVM$Matrix)
    matri <- RVM[[1]]$Matrix
    #w1 <- as.vector(RVM$family)
    w1 <- .('family')
    #th <- as.vector(RVM$par)
    th <- .('par')
    #th2 <- as.vector(RVM$par2)
    th2 <- .('par2')
    #maxmat <- as.vector(RVM$MaxMat)
    maxmat <- RVM[[1]]$MaxMat
    #conindirect <- as.vector(RVM$CondDistr$indirect)
    conindirect <- RVM[[1]]$CondDistr$indirect
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0
    #tmp <- rep(0, n * N)

    ## simulate R-Vine
    #tmp <- .C("SimulateRVine",
    #            as.integer(N),
    #            as.integer(n),
    #            as.integer(w1),
    #            as.integer(maxmat),
    #            as.integer(matri),
    #            as.integer(conindirect),
    #            as.double(th),
    #            as.double(th2),
    #            as.double(tmp),
    #            as.double(U),
    #            as.integer(takeU),
    #            PACKAGE = "VineCopula")[[9]]
    tmp <- SimulateRVine(N, n, w1, maxmat, matri, conindirect, th, th2, NULL, U, takeU)

    ## store results, bring back to initial order and return
    out <- matrix(tmp, ncol = n, byrow = TRUE)
    #if (!is.null(RVM$names)) {
    if (!is.null(RVM[[1]]$names)) {
        #colnames(out) <- RVM$names
        colnames(out) <- RVM[[1]]$names
    }
    out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    return(out)
}


#void SimulateRVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out, double* U, int* takeU)
SimulateRVine <- function(T_, d, family, maxmat, matrix_, conindirect, par, par2, out, U, takeU) {
#{
    #int i, j, k, m, one, **fam, **cindirect, **mat, **mmat, **fam2, **cindirect2, **mat2, **mmat2;
    #double **theta, **nu, **theta2, **nu2, ***vdirect, ***vindirect, **U2;

    #one = 1;
    #//Allocate memory
    #theta=create_matrix(*d,*d);
    #nu=create_matrix(*d,*d);
    #fam=create_intmatrix(*d,*d);
    #mmat=create_intmatrix(*d,*d);
    #cindirect=create_intmatrix(*d,*d);
    #mat=create_intmatrix(*d,*d);
    #theta2=create_matrix(*d,*d);
    #nu2=create_matrix(*d,*d);
    #fam2=create_intmatrix(*d,*d);
    #mmat2=create_intmatrix(*d,*d);
    #cindirect2=create_intmatrix(*d,*d);
    #mat2=create_intmatrix(*d,*d);
    #vdirect = create_3darray(*d,*d,1);
    vdirect <- array(dim=c(d, d, T_))
    #vindirect = create_3darray(*d,*d,1);
    vindirect <- array(dim=c(d, d, T_))
    #U2 = create_matrix(*T, *d);

    #//Initialize random number generator:
    #GetRNGstate();

    #//Initialize
    #k=0;
    #for(i=0;i<(*d);i++)
    #{
        #for(j=0;j<(*d);j++)
        #{
            #theta2[i][j]=par[(i+1)+(*d)*j-1] ;
            #nu2[i][j]=par2[(i+1)+(*d)*j-1]    ;
            #mmat2[i][j]=maxmat[(i+1)+(*d)*j-1] ;
            #mat2[i][j]=matrix[(i+1)+(*d)*j-1] ;
            #cindirect2[i][j]=conindirect[(i+1)+(*d)*j-1] ;
            #fam2[i][j]=family[(i+1)+(*d)*j-1] ;
        #}
    #}
    theta2 <- par
    nu2 <- par2
    mmat2 <- maxmat
    mat2 <- matrix_
    cindirect2 <- conindirect
    fam2 <- family
    #if(*takeU == 1)
    if (takeU) {
    #{
        #for(j=0;j<(*d);j++) for(i=0;i<(*T);i++) U2[i][j]=U[(*T)*j+i]; // (T [=N], d)-matrix
        U2 <- U
    #}
    }

    #// Matrizen rotieren f?r den Algo
    #for(i=0;i<(*d);i++)
    #{
        #for(j=0;j<(*d);j++)
        #{
            #theta[(*d-i-1)][(*d-j-1)]=theta2[i][j];
            #nu[(*d-i-1)][(*d-j-1)]=nu2[i][j];
            #mmat[(*d-i-1)][(*d-j-1)]=mmat2[i][j];
            #mat[(*d-i-1)][(*d-j-1)]=mat2[i][j];
            #cindirect[(*d-i-1)][(*d-j-1)]=cindirect2[i][j];
            #fam[(*d-i-1)][(*d-j-1)]=fam2[i][j];
        #}
    #}
    . <- function(.) .[nrow(.):1, ncol(.):1,, drop=F]
    fam <- .(fam2)
    theta <- .(theta2)
    nu <- .(nu2)
    . <- function(.) .[nrow(.):1, ncol(.):1, drop=F]
    mmat <- .(mmat2)
    mat <- .(mat2)
    cindirect <- .(cindirect2)

    #free_matrix(theta2,*d);
    #free_matrix(nu2,*d);
    #free_intmatrix(fam2,*d);
    #free_intmatrix(mmat2,*d);
    #free_intmatrix(cindirect2,*d);
    #free_intmatrix(mat2, *d);

    #// Der eigentliche Algo
    #int nn = 0;
    #for(int n = 0; n < *T; n++) // sample size
    #{
        #if(*takeU == 1) for(i=0;i<*d;i++) vdirect[i][i][0] = U2[n][i]; // j = 'sample size'; i = 'copula dimension'
        #else for(i=0;i<*d;i++) vdirect[i][i][0] = runif(0,1);
    vdirect[1:d + d*(0:(d - 1)) + rep((0:(T_ - 1))*d^2, each=d)] <- if (takeU) t(U2) else runif(d*T_)
        #vindirect[0][0][0] = vdirect[0][0][0];
    vindirect[1, 1,] <- vdirect[1, 1,]

        #for(i=1;i<*d;i++)
    for (i in 2:d) {
        #{
            #for(k=(i-1);k>(-1);k--)
        for (k in (i - 1):1) {
            #{
                #m = mmat[k][i];
            m <- mmat[k, i]
                #if(mat[k][i]==m)
                #{
                    #Hinv1(&fam[k][i],&one,vdirect[k+1][i],vdirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
                #}
                #else
                #{
                    #Hinv1(&fam[k][i],&one,vdirect[k+1][i],vindirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
                #}
            vdirect[k, i,] <- BiCopHinv1((if (mat[k, i] == m) vdirect else vindirect)[k, m,], vdirect[k + 1, i,], fam[k, i,], theta[k, i,], par2=nu[k, i,], check.pars=F)

                #if(i+1<(*d))
            if (i < d) {
                #{
                    #if(cindirect[k+1][i]==1)
                if (cindirect[k + 1, i]) {
                    #{
                        #if(mat[k][i]==m)
                        #{
                            #Hfunc2(&fam[k][i],&one,vdirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
                        #}
                        #else
                        #{
                            #Hfunc2(&fam[k][i],&one,vindirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
                        #}
                    vindirect[k + 1, i,] <- BiCopHfunc2((if (mat[k, i] == m) vdirect else vindirect)[k, m,], vdirect[k, i,], fam[k, i,], theta[k, i,], par2=nu[k, i,], check.pars=F)
                    #}
                }
                #}
            }
            #}
        }
        #}
    }

        #for(i=0;i<(*d);i++)
        #{
            #out[nn]=vdirect[0][i][0];
            #nn++;
        #}
    out <- vdirect[1,,]
    #}


    #//Free memory:
    #free_matrix(theta,*d);
    #free_matrix(nu,*d);
    #free_intmatrix(fam,*d);
    #free_intmatrix(mmat,*d);
    #free_intmatrix(cindirect,*d);
    #free_intmatrix(mat, *d);
    #free_3darray(vdirect,*d,*d);
    #free_3darray(vindirect,*d,*d);
    #free_matrix(U2,*T);
    #PutRNGstate();
    return(out)
}


transform <- function(M) {
    n <- dim(M)[1]

    M.new <- matrix(rep(0, n * n), n, n)
    for (i in 1:n) {
        for (j in 1:i) {
            M.new[(n - i + 1), (n - j + 1)] <- M[i, j]
        }
    }

    return(M.new)
}

