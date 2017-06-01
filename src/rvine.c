/*
 ** rvine.c - C code of the package CDRVine
 **
 ** with contributions from Carlos Almeida, Aleksey Min,
 ** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
 **
 ** A first version was based on code
 ** from Daniel Berg <daniel at danielberg.no>
 ** provided by personal communication.
 **
 */

#include "vine.h"
#include "memoryhandling.h"
#include "likelihood.h"
#include "rvine.h"
#include "hfunc.h"



// Code from Jakob St?ber and Ulf Schepsmeier for R-vine log-likelihood calculation

//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for the pair-copula construction (Rvine)
// (by J.S.)
// Input:
// T        sample size
// d        dimension (>=2)
// family   copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// par      parameter values (at least d*(d-1)/2 parameters
// par2		second set of parameter values (f.e. for student copulas)
// data     data set for which to compute log-likelihood
// matrix   an RVineMatrix in vector form
// condirect, conindirect Matrizes which tell us where we find the right values
// seperate	Control Parameter, do we want to seperate the likelihoods for each data point?
// calcupdate matrix which tells us for which parameters we need to redo the calculations, not newly computed values are taken from ll, vv, vv2
// Output:
// out      Loglikelihood
// ll       array with the contribution to LL (for each copula)
// vv,vv2       array for the transformation operated (Hfunc)
/////////////////////////////////////////////////////////////
void VineLogLikRvine(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
                     double* out, double* ll, double* vv, double* vv2, int* calcupdate, int* seperate)
{
    int i, j, k, t, m, **fam;
    int kk=1;
    double loglik=1.0, sumloglik=0.0, **x, **theta, **nu, ***vdirect, ***vindirect;

    double *sumsplitlog;
    sumsplitlog=(double*) Calloc(*T,double);

    for (t=0;t<*T;t++ )
    {
        sumsplitlog[t] = 0;
    }



    double ***value2;
    value2=create_3darray(*d,*d,*T);
    double **value;
    value=create_matrix(*d,*d);



    //Allocate memory
    x = create_matrix(*d,*T);
    vdirect = create_3darray(*d,*d,*T);
    vindirect = create_3darray(*d,*d,*T);
    theta=create_matrix(*d,*d);
    nu=create_matrix(*d,*d);
    fam=create_intmatrix(*d,*d);


    //Initialize
    k=0;
    for(i=0;i<(*d);i++)
    {
        for (t=0;t<*T;t++ )
        {
            x[i][t] = data[k];
            k++;
        }
    }

    k=0;
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta[i][j]=par[(i+1)+(*d)*j-1] ;
            nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
            fam[i][j]=family[(i+1)+(*d)*j-1] ;
            if(*seperate==0)
            {
                value[i][j]=ll[(i+1)+(*d)*j-1] ;
            }
        }
    }


    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            for(t=0;t<*T;t++ )
            {
                vdirect[i][j][t]=vv[(i+1)+(*d)*j+(*d)*(*d)*t-1];
                vindirect[i][j][t]=vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1];
                if(*seperate==1)
                {
                    value2[i][j][t]=ll[(i+1)+(*d)*j+(*d)*(*d)*t-1];
                }
            }
        }
    }



    for(i=0;i<(*d);i++)
    {
        for(t=0;t<*T;t++ )
        {
            vdirect[*d-1][i][t]=x[*d-1-i][t];
        }
    }


    for(i=*d-2; i>-1; i--)
    {
        for(k=*d-1;k>i;k--)
        {
            if(calcupdate[(k+1)+(*d)*i-1]==1)
            {
                m=maxmat[(k+1)+(*d)*i-1];

                if(m == matrix[(k+1)+(*d)*i-1])
                {
                    if(*seperate==1)
                    {
                        kk = 1;
                        for(t=0;t<*T;t++ )
                        {
                            LL_mod2(&fam[k][i],&kk,&vdirect[k][(*d-m)][t],&vdirect[k][i][t],&theta[k][i],&nu[k][i],&loglik);
                            value2[k][i][t]=loglik;
                        }
                    }
                    else
                    {

                        LL_mod2(&fam[k][i],T,vdirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],&loglik);
                        value[k][i]=loglik;
                    }
                    if(condirect[k+(*d)*i-1]==1)
                    {
                        Hfunc1(&fam[k][i],T,vdirect[k][i],vdirect[k][*d-m],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                    }
                    if(conindirect[k+(*d)*i-1]==1)
                    {
                        Hfunc2(&fam[k][i],T,vdirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                    }

                }
                else
                {
                    if(*seperate==1)
                    {
                        kk = 1;
                        for(t=0;t<*T;t++ )
                        {
                            LL_mod2(&fam[k][i],&kk,&vindirect[k][(*d-m)][t],&vdirect[k][i][t],&theta[k][i],&nu[k][i],&loglik);
                            value2[k][i][t]=loglik;
                        }
                    }
                    else
                    {
                        LL_mod2(&fam[k][i],T,vindirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],&loglik);
                        value[k][i]=loglik;
                    }
                    if(condirect[k+(*d)*i-1]==1)
                    {
                        Hfunc1(&fam[k][i],T,vdirect[k][i],vindirect[k][(*d-m)],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                    }
                    if(conindirect[k+(*d)*i-1]==1)
                    {
                        Hfunc2(&fam[k][i],T,vindirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                    }
                }
            }
            if(*seperate==0)
            {
                sumloglik += value[k][i];
            }
            else
            {
                for(t=0;t<*T;t++ )
                {
                    sumsplitlog[t] += value2[k][i][t];
                }
            }
        }
    }



    if(*seperate==0)
    {
        *out = sumloglik;
    }
    else
    {
        for(t=0;t<*T;t++ )
        {
            out[t] =  sumsplitlog[t];
        }
    }
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            for(t=0;t<*T;t++ )
            {
                vv[(i+1)+(*d)*j+(*d)*(*d)*t-1]=vdirect[i][j][t];
                vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1]=vindirect[i][j][t];
                if(*seperate==1)
                {
                    ll[(i+1)+(*d)*j+(*d)*(*d)*t-1]=value2[i][j][t];
                }
            }
        }
    }

    if(*seperate==0)
    {
        for(i=0;i<(*d);i++)
        {
            for(j=0;j<(*d);j++)
            {
                ll[(i+1)+(*d)*j-1]=value[i][j];
            }
        }
    }


    //Free memory:
    free_matrix(x,*d);
    free_3darray(vdirect,*d,*d);
    free_matrix(theta,*d);
    free_matrix(nu,*d);
    free_intmatrix(fam,*d);
    free_matrix(value,*d);
    free_3darray(value2,*d,*d);
    free_3darray(vindirect,*d,*d);
    Free(sumsplitlog);
}

// memory efficient version, where outer loop is over sample size
void VineLogLikRvine2(int* T, int* d, int* family, int* maxmat, int* matrix,
                      int* condirect, int* conindirect,
                      double* par, double* par2, double* data,
                      double* out)
{
    int i, j, k, t, m, **fam;
    int kk=1;
    double loglik=1.0, **x, **theta, **nu, ***vdirect, ***vindirect;

    double *sumsplitlog;
    sumsplitlog=(double*) Calloc(*T,double);

    for (t=0;t<*T;t++ )
    {
        sumsplitlog[t] = 0;
    }

    double ***value2;
    value2=create_3darray(*d,*d,1);

    //Allocate memory
    x = create_matrix(*d,*T);
    vdirect = create_3darray(*d,*d,1);
    vindirect = create_3darray(*d,*d,1);
    theta=create_matrix(*d,*d);
    nu=create_matrix(*d,*d);
    fam=create_intmatrix(*d,*d);

    //Initialize
    k=0;
    for(i=0;i<(*d);i++)
    {
        for (t=0;t<*T;t++ )
        {
            x[i][t] = data[k];
            k++;
        }
    }

    k=0;
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta[i][j]=par[(i+1)+(*d)*j-1] ;
            nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
            fam[i][j]=family[(i+1)+(*d)*j-1] ;
        }
    }

    for(t=0;t<*T;t++ )
    {
        for(i=0;i<(*d);i++)
        {
            vdirect[*d-1][i][0]=x[*d-1-i][t];
        }


        for(i=*d-2; i>-1; i--)
        {
            for(k=*d-1;k>i;k--)
            {
                m=maxmat[(k+1)+(*d)*i-1];
                kk = 1;
                if(m == matrix[(k+1)+(*d)*i-1])
                {
                    LL_mod2(&fam[k][i],&kk,&vdirect[k][(*d-m)][0],&vdirect[k][i][0],&theta[k][i],&nu[k][i],&loglik);
                    value2[k][i][0]=loglik;

                    if(condirect[k+(*d)*i-1]==1)
                    {
                        Hfunc1(&fam[k][i],&kk,vdirect[k][i],vdirect[k][*d-m],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                    }
                    if(conindirect[k+(*d)*i-1]==1)
                    {
                        Hfunc2(&fam[k][i],&kk,vdirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                    }

                }
                else
                {
                    LL_mod2(&fam[k][i],&kk,&vindirect[k][(*d-m)][0],&vdirect[k][i][0],&theta[k][i],&nu[k][i],&loglik);
                    value2[k][i][0]=loglik;
                    if(condirect[k+(*d)*i-1]==1)
                    {
                        Hfunc1(&fam[k][i],&kk,vdirect[k][i],vindirect[k][(*d-m)],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                    }
                    if(conindirect[k+(*d)*i-1]==1)
                    {
                        Hfunc2(&fam[k][i],&kk,vindirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                    }
                }

                sumsplitlog[t] += value2[k][i][0];
            }
        }
        out[t] =  sumsplitlog[t];
    }


    //Free memory:
    free_matrix(x,*d);
    free_3darray(vdirect,*d,*d);
    free_matrix(theta,*d);
    free_matrix(nu,*d);
    free_intmatrix(fam,*d);
    free_3darray(value2,*d,*d);
    free_3darray(vindirect,*d,*d);
    Free(sumsplitlog);

}


//////////////////////////////////////////////////////////////
// Function to compute likelihood for the pair-copula construction (Rvine)
// (by J.S.)
// Input:
// T        sample size
// d        dimension (>=2)
// family   copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// par      parameter values (at least d*(d-1)/2 parameters
// par2		second set of parameter values (f.e. for student copulas)
// data     data set for which to compute log-likelihood
// matrix   an RVineMatrix in vector form
// condirect, conindirect Matrizes which tell us where we find the right values
// seperate	Control Parameter, do we want to seperate the likelihoods for each data point?
// calcupdate matrix which tells us for which parameters we need to redo the calculations, not newly computed values are taken from ll, vv, vv2
// Output:
// out      Loglikelihood
// ll       array with the contribution to LL (for each copula)
// vv,vv2       array for the transformation operated (Hfunc)
/////////////////////////////////////////////////////////////


/// seperate = 1 not implemented, highly experimental
void VineLikRvine(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
                  double* out, double* ll, double* vv, double* vv2, int* calcupdate, int* seperate)
{
    int i, j, k, t, m, **fam, **cdirect, **cindirect, **mat, **mmat, **calc;
    int kk=1;
    double loglik=1.0, sumloglik=1.0, **x, **theta, **nu, ***vdirect, ***vindirect;

    double *sumsplitlog;
    sumsplitlog=(double*) Calloc(*T,double);

    for (t=0;t<*T;t++ )
    {
        sumsplitlog[t] = 0;
    }



    double ***value2;
    value2=create_3darray(*d,*d,*T);
    double **value;
    value=create_matrix(*d,*d);



    //Allocate memory
    x = create_matrix(*d,*T);
    vdirect = create_3darray(*d,*d,*T);
    vindirect = create_3darray(*d,*d,*T);
    theta=create_matrix(*d,*d);
    nu=create_matrix(*d,*d);
    fam=create_intmatrix(*d,*d);
    mmat=create_intmatrix(*d,*d);
    cdirect=create_intmatrix(*d,*d);
    cindirect=create_intmatrix(*d,*d);
    mat=create_intmatrix(*d,*d);
    calc=create_intmatrix(*d,*d);


    //Initialize
    k=0;
    for(i=0;i<(*d);i++)
    {
        for (t=0;t<*T;t++ )
        {
            x[i][t] = data[k];
            k++;
        }
    }

    k=0;
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta[i][j]=par[(i+1)+(*d)*j-1] ;
            nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
            mmat[i][j]=maxmat[(i+1)+(*d)*j-1] ;
            mat[i][j]=matrix[(i+1)+(*d)*j-1] ;
            cdirect[i][j]=condirect[(i+1)+(*d)*j-1];
            cindirect[i][j]=conindirect[(i+1)+(*d)*j-1] ;
            fam[i][j]=family[(i+1)+(*d)*j-1] ;
            if(*seperate==0){value[i][j]=ll[(i+1)+(*d)*j-1] ;}
            calc[i][j]=calcupdate[(i+1)+(*d)*j-1] ;
        }
    }


    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++){
            for(t=0;t<*T;t++ ) {
                vdirect[i][j][t]=vv[(i+1)+(*d)*j+(*d)*(*d)*t-1];
                vindirect[i][j][t]=vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1];
                if(*seperate==1){ value2[i][j][t]=ll[(i+1)+(*d)*j+(*d)*(*d)*t-1];}
            }}}



    for(i=0;i<(*d);i++)
    {
        for(t=0;t<*T;t++ )
        {
            vdirect[*d-1][i][t]=x[*d-1-i][t];
        }
    }


    for(i=*d-2; i>-1; i--)
    {
        for(k=*d-1;k>i;k--)
        {

            if(calc[k][i]==1){
                m=mmat[k][i];


                if(m == mat[k][i])
                {
                    if(*seperate==1){kk = 1;
                        for(t=0;t<*T;t++ ) {

                            copLik_mod(&fam[k][i],&kk,&vdirect[k][i][t],&vdirect[k][(*d-m)][t],&theta[k][i],&nu[k][i],&loglik);
                            /* sumloglik += loglik; */
                            value2[k][i][t]=loglik;

                        }}else{

                            copLik_mod(&fam[k][i],T,vdirect[k][i],vdirect[k][(*d-m)],&theta[k][i],&nu[k][i],&loglik);
                            /* sumloglik += loglik; */
                            value[k][i]=loglik;}
                        if(cdirect[k-1][i]==1) {
                            Hfunc1(&fam[k][i],T,vdirect[k][i],vdirect[k][*d-m],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                        }
                        if(cindirect[k-1][i]==1)  {
                            Hfunc2(&fam[k][i],T,vdirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                        }



                }else{
                    if(*seperate==1){kk = 1;
                        for(t=0;t<*T;t++ ) {
                            copLik_mod(&fam[k][i],&kk,&vdirect[k][i][t],&vindirect[k][(*d-m)][t],&theta[k][i],&nu[k][i],&loglik);
                            /* sumloglik += loglik; */
                            value2[k][i][t]=loglik;

                        }}else{

                            copLik_mod(&fam[k][i],T,vdirect[k][i],vindirect[k][(*d-m)],&theta[k][i],&nu[k][i],&loglik);
                            /* sumloglik += loglik; */
                            value[k][i]=loglik;}
                        if(cdirect[k-1][i]==1) {
                            Hfunc1(&fam[k][i],T,vdirect[k][i],vindirect[k][(*d-m)],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
                        }
                        if(cindirect[k-1][i]==1)  {
                            Hfunc2(&fam[k][i],T,vindirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
                        }

                }
            }
            if(*seperate==0){sumloglik = sumloglik*value[k][i];}else{ for(t=0;t<*T;t++ ) {  sumsplitlog[t] += value2[k][i][t];}  };
        }
    }





    if(*seperate==0){  *out = sumloglik;}else{for(t=0;t<*T;t++ ) {out[t] =  sumsplitlog[t];}};
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++){
            for(t=0;t<*T;t++ ) {
                vv[(i+1)+(*d)*j+(*d)*(*d)*t-1]=vdirect[i][j][t];
                vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1]=vindirect[i][j][t];
                if(*seperate==1){ll[(i+1)+(*d)*j+(*d)*(*d)*t-1]=value2[i][j][t];}
            }}}
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            if(*seperate==0){ll[(i+1)+(*d)*j-1]=value[i][j];}}}


    //Free memory:
    free_matrix(x,*d); free_3darray(vdirect,*d,*d); free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d); free_matrix(value,*d);free_3darray(value2,*d,*d);
    free_intmatrix(mmat,*d); free_intmatrix(cdirect,*d); free_intmatrix(cindirect,*d);   free_3darray(vindirect,*d,*d); free_intmatrix(calc, *d); free_intmatrix(mat, *d);
    Free(sumsplitlog);
}



////////////////////////////////////////////
// PCC for R-vine
////////////////////////////////////////////

void SimulateRVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out, double* U, int* takeU)
{
    int i, j, k, m, one, **fam, **cindirect, **mat, **mmat, **fam2, **cindirect2, **mat2, **mmat2;
    double **theta, **nu, **theta2, **nu2, ***vdirect, ***vindirect, **U2;

    one = 1;
    //Allocate memory
    theta=create_matrix(*d,*d);
    nu=create_matrix(*d,*d);
    fam=create_intmatrix(*d,*d);
    mmat=create_intmatrix(*d,*d);
    cindirect=create_intmatrix(*d,*d);
    mat=create_intmatrix(*d,*d);
    theta2=create_matrix(*d,*d);
    nu2=create_matrix(*d,*d);
    fam2=create_intmatrix(*d,*d);
    mmat2=create_intmatrix(*d,*d);
    cindirect2=create_intmatrix(*d,*d);
    mat2=create_intmatrix(*d,*d);
    vdirect = create_3darray(*d,*d,1);
    vindirect = create_3darray(*d,*d,1);
    U2 = create_matrix(*T, *d);

    //Initialize random number generator:
    GetRNGstate();

    //Initialize
    k=0;
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta2[i][j]=par[(i+1)+(*d)*j-1] ;
            nu2[i][j]=par2[(i+1)+(*d)*j-1]    ;
            mmat2[i][j]=maxmat[(i+1)+(*d)*j-1] ;
            mat2[i][j]=matrix[(i+1)+(*d)*j-1] ;
            cindirect2[i][j]=conindirect[(i+1)+(*d)*j-1] ;
            fam2[i][j]=family[(i+1)+(*d)*j-1] ;
        }
    }
    if(*takeU == 1)
    {
        for(j=0;j<(*d);j++) for(i=0;i<(*T);i++) U2[i][j]=U[(*T)*j+i]; // (T [=N], d)-matrix
    }

    // Matrizen rotieren f?r den Algo
    for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta[(*d-i-1)][(*d-j-1)]=theta2[i][j];
            nu[(*d-i-1)][(*d-j-1)]=nu2[i][j];
            mmat[(*d-i-1)][(*d-j-1)]=mmat2[i][j];
            mat[(*d-i-1)][(*d-j-1)]=mat2[i][j];
            cindirect[(*d-i-1)][(*d-j-1)]=cindirect2[i][j];
            fam[(*d-i-1)][(*d-j-1)]=fam2[i][j];
        }
    }

    free_matrix(theta2,*d);
    free_matrix(nu2,*d);
    free_intmatrix(fam2,*d);
    free_intmatrix(mmat2,*d);
    free_intmatrix(cindirect2,*d);
    free_intmatrix(mat2, *d);

    // Der eigentliche Algo
    int nn = 0;
    for(int n = 0; n < *T; n++) // sample size
    {
        if(*takeU == 1) for(i=0;i<*d;i++) vdirect[i][i][0] = U2[n][i]; // j = 'sample size'; i = 'copula dimension'
        else for(i=0;i<*d;i++) vdirect[i][i][0] = runif(0,1);
        vindirect[0][0][0] = vdirect[0][0][0];

        for(i=1;i<*d;i++)
        {
            for(k=(i-1);k>(-1);k--)
            {
                m = mmat[k][i];
                if(mat[k][i]==m)
                {
                    Hinv1(&fam[k][i],&one,vdirect[k+1][i],vdirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
                }
                else
                {
                    Hinv1(&fam[k][i],&one,vdirect[k+1][i],vindirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
                }

                if(i+1<(*d))
                {
                    if(cindirect[k+1][i]==1)
                    {
                        if(mat[k][i]==m)
                        {
                            Hfunc2(&fam[k][i],&one,vdirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
                        }
                        else
                        {
                            Hfunc2(&fam[k][i],&one,vindirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
                        }
                    }
                }
            }
        }

        for(i=0;i<(*d);i++)
        {
            out[nn]=vdirect[0][i][0];
            nn++;
        }
    }


    //Free memory:
    free_matrix(theta,*d);
    free_matrix(nu,*d);
    free_intmatrix(fam,*d);
    free_intmatrix(mmat,*d);
    free_intmatrix(cindirect,*d);
    free_intmatrix(mat, *d);
    free_3darray(vdirect,*d,*d);
    free_3darray(vindirect,*d,*d);
    free_matrix(U2,*T);
    PutRNGstate();
}
