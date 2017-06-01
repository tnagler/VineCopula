#include "vine.h"
#include "memoryhandling.h"
#include "likelihood.h"
#include "deriv.h"
#include "tcopuladeriv.h"
#include "tcopuladeriv_new.h"
#include "rvinederiv.h"
#include "rvinederiv2.h"


// Code from Jakob Stoeber and Ulf Schepsmeier for R-vine log-likelihood derivative calculation

//////////////////////////////////////////////////////////////
// Function to compute the derivative of log-likelihood for the pair-copula construction (Rvine) (one-element of the gradient)
// (by J.S.)
// Input:
// T        		sample size
// d        		dimension (>=2)
// family   		copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// kk				row number of the parameter wrt to which we want to calculate the derivative
// ii				column number of the parameter wrt to which we want to calculate the derivative
// par      		parameter values (at least d*(d-1)/2 parameters
// par2				second set of parameter values (f.e. for student copulas)
// data     		data set for which to compute the derivative of the log-likelihood
// matrix   		an RVineMatrix in vector form
// condirect, conindirect Matrices which tell us where we find the right values 
// calcupdate 		matrix which tells which terns we need to consider for the calculation of the derivative
// ll       		array with the contribution to the derivative (for each copula)
// vv,vv2       	array of the h-functions  (given as by-product of the log-likelihood calculation)
// tcop				a special marker for the Student's t-copula (1=first parameter, 2=second parameter)
// margin			derivative wrt to the margins as well? (TRUE/FALSE) (needed by Jakob for some of his calculations)
// 
// Output:
// out      		Log-likelihood derivative
// tilde_value		array of separated derivatives in the R-vine construction
// tilde_vdirect	array of separated derivatives of the h-functions needed
// tilde_vindirect	array of separated derivatives of the h-functions needed
/////////////////////////////////////////////////////////////

// Reference:
// St?ber, J. and U. Schepsmeier (2013)
// Estimating standard errors in regular vine copula models
// Computational Statistics, 28 (6), 2679-2707


void VineLogLikRvineDeriv(int* T, int* d, int* family, int* kk, int* ii, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* out, double* ll, double* vv, double* vv2, int* calcupdate, double* tilde_vdirect, double* tilde_vindirect, double* tilde_value, int* tcop, int* margin)
{
	int i, j, k, t, m, **fam, **calc;
	
	double sumloglik=0.0, **theta, **nu, ***tildevdirect, ***tildevindirect, *zr1, *zr2, *tildezr1, *tildezr2, *cop;
	double *handle1;
	double param[2];
	
	zr1=(double*) Calloc(*T,double);
	zr2=(double*) Calloc(*T,double);
	handle1=(double*) Calloc(*T,double);
	tildezr1=(double*) Calloc(*T,double);
	tildezr2=(double*) Calloc(*T,double);
	cop=(double*) Calloc(*T,double);
	
	double ***tildevalue;
	tildevalue=create_3darray(*d,*d,*T);
	
	
	//Allocate memory
	tildevdirect = create_3darray(*d,*d,*T);
	tildevindirect = create_3darray(*d,*d,*T);
	theta=create_matrix(*d,*d);
	nu=create_matrix(*d,*d);
	fam=create_intmatrix(*d,*d);
	calc=create_intmatrix(*d,*d);

	
	
	//Initialize
    
	k=0;
	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++)
		{
			theta[i][j]=par[(i+1)+(*d)*j-1] ;
			nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
			fam[i][j]=family[(i+1)+(*d)*j-1] ;
			calc[i][j]=calcupdate[(i+1)+(*d)*j-1] ;
		}
	}       
	
	
	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++)
		{	
			for(t=0;t<*T;t++ ) 
			{
				tildevdirect[i][j][t]=0;
				tildevindirect[i][j][t]=0;
				tildevalue[i][j][t]=0;
			}
		}
	}
	
	
	
	m=maxmat[*kk+(*d)*(*ii-1)-1];
	for(t=0;t<*T;t++ ) 
	{
		zr1[t]=vv[*kk+(*d)*(*ii-1)+(*d)*(*d)*t-1];
		cop[t]=exp(ll[*kk+(*d)*(*ii-1)+(*d)*(*d)*t-1]);
	}
	if(m == matrix[*kk+*d*(*ii-1)-1])
	{	
		for(t=0;t<*T;t++ ) 
		{
			zr2[t]=vv[*kk+(*d)*(*d-m)+(*d)*(*d)*t-1];
		}
		
	}
	else 
	{
		for(t=0;t<*T;t++)
		{
			zr2[t]=vv2[*kk+(*d)*(*d-m)+(*d)*(*d)*t-1];
		}	
	}
	

	param[0]=theta[*kk-1][*ii-1];
	param[1]=nu[*kk-1][*ii-1];
	if(*tcop==1)		//For the t-copula (first parameter)
	{
		diffhfunc_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],tildevdirect[*kk-2][*ii-1]);
		diffhfunc_rho_tCopula(zr2,zr1,T,param,&fam[*kk-1][*ii-1],tildevindirect[*kk-2][*ii-1]);
		diffPDF_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],tildevalue[*kk-1][*ii-1]);
		for(t=0;t<*T;t++ ) 
		{
			tildevalue[*kk-1][*ii-1][t]=tildevalue[*kk-1][*ii-1][t]/cop[t];
		}
	}
	else if(*tcop==2)  // for the t-copula (second parameter)
	{
		diffhfunc_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],tildevdirect[*kk-2][*ii-1]);
		diffhfunc_nu_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],tildevindirect[*kk-2][*ii-1]);
		diffPDF_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],tildevalue[*kk-1][*ii-1]);
		for(t=0;t<*T;t++ ) 
		{
			tildevalue[*kk-1][*ii-1][t]=tildevalue[*kk-1][*ii-1][t]/cop[t];
		}
	}
	else
	{
		if( *margin == 0 )		//Das ist unser bisheriger Fall mit stetigen Variablen (ohne t-copula)
		{		
			diffhfunc_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevdirect[*kk-2][*ii-1]);
			diffhfunc_mod2(zr2,zr1,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevindirect[*kk-2][*ii-1]);
			diffPDF_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevalue[*kk-1][*ii-1]);
			for(t=0;t<*T;t++ ) 
			{
				tildevalue[*kk-1][*ii-1][t]=tildevalue[*kk-1][*ii-1][t]/cop[t];
			}
		}
		else if( *margin== 1)	// Ableitung nach dem ersten Argument = margin1
		{
			diffhfunc_v_mod2(zr2,zr1,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevindirect[*kk-2][*ii-1]);
			diffPDF_u_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevalue[*kk-1][*ii-1]); // hier k?nnte difflPDF stehen
			for(t=0;t<*T;t++ ) 
			{
				tildevdirect[*kk-2][*ii-1][t]=cop[t];
				tildevalue[*kk-1][*ii-1][t]=tildevalue[*kk-1][*ii-1][t]/cop[t];
			}
		}
		else					// Ableitung nach dem zweiten Argument = margin2
		{
			diffhfunc_v_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevdirect[*kk-2][*ii-1]);
			diffPDF_v_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],tildevalue[*kk-1][*ii-1]); // hier k?nnte difflPDF stehen
			for(t=0;t<*T;t++ ) 
			{
				tildevindirect[*kk-2][*ii-1][t]=cop[t];
				tildevalue[*kk-1][*ii-1][t]=tildevalue[*kk-1][*ii-1][t]/cop[t];
			}
		
		}
	}
	
		
	// add up for the final derivative
	for(t=0;t<*T;t++ ) 
	{
		sumloglik+=tildevalue[*kk-1][*ii-1][t];
	}
	
	for(i=*ii-1; i>-1; i--)
    {
		for(k=*kk-2;k>i;k--)
        {   
			if(calc[k][i]==1)
			{
				m=maxmat[(k+1)+*d*i-1];
				
				for(t=0;t<*T;t++ ) 
				{
					zr1[t]=vv[(k+1)+(*d)*i+(*d)*(*d)*t-1];
					tildezr1[t]=tildevdirect[k][i][t];
					cop[t]=exp(ll[(k+1)+(*d)*i+(*d)*(*d)*t-1]);
				}
				if(m == matrix[(k+1)+*d*i-1])
				{	
					for(t=0;t<*T;t++ ) 
					{
						zr2[t]=vv[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
						tildezr2[t]=tildevdirect[k][(*d-m)][t];
					}
				}
				else 
				{
					for(t=0;t<*T;t++ ) 
					{
						zr2[t]=vv2[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
						tildezr2[t]=tildevindirect[k][(*d-m)][t];
					}	
				}
				for(t=0;t<*T;t++ ) 
				{
					tildevdirect[k-1][i][t]=0;
					tildevindirect[k-1][i][t]=0;
					tildevalue[k][i][t]=0;
				}
				if(calc[k+1][i]==1)
				{
					param[0]=theta[k][i];
					param[1]=nu[k][i];
					if(fam[k][i]==2)		//For the t-copula
					{
						diffPDF_u_tCopula_new(zr1,zr2,T,param,&fam[k][i],handle1);
					}
					else
					{
						diffPDF_u_mod(zr1,zr2,T,&theta[k][i],&fam[k][i],handle1);
					}
					for(t=0;t<*T;t++ ) 
					{
						tildevalue[k][i][t]+=handle1[t]/cop[t]*tildezr1[t];
					}
					
					if(condirect[k+(*d)*i-1]==1)
					{
						for(t=0;t<*T;t++ ) 
						{
							tildevdirect[k-1][i][t]+=cop[t]*tildezr1[t];
						}
					}
					if(conindirect[k+(*d)*i-1]==1)
					{
						param[0]=theta[k][i];
						param[1]=nu[k][i];
						if(fam[k][i]==2)		//For the t-copula
						{
							diffhfunc_v_tCopula_new(zr2,zr1,T,param,&fam[k][i],handle1);
						}
						else
						{
							diffhfunc_v_mod2(zr2,zr1,T,&theta[k][i],&fam[k][i],handle1);
						}
						for(t=0;t<*T;t++ ) 
						{
							tildevindirect[k-1][i][t]+=handle1[t]*tildezr1[t];
						}
					}
					
				}
				if(calc[k+1][(*d-m)]==1)
				{
					param[0]=theta[k][i];
					param[1]=nu[k][i];
					if(fam[k][i]==2)		//For the t-copula
					{
						diffPDF_u_tCopula_new(zr2,zr1,T,param,&fam[k][i],handle1);
					}
					else
					{
						diffPDF_v_mod(zr1,zr2,T,&theta[k][i],&fam[k][i],handle1);
					}
					for(t=0;t<*T;t++ ) 
					{
						tildevalue[k][i][t]+=handle1[t]/cop[t]*tildezr2[t];
					}
						
					if(condirect[k+(*d)*i-1]==1)
					{
						param[0]=theta[k][i];
						param[1]=nu[k][i];
						if(fam[k][i]==2)		//For the t-copula
						{
							diffhfunc_v_tCopula_new(zr1,zr2,T,param,&fam[k][i],handle1);
						}
						else
						{
							diffhfunc_v_mod(zr1,zr2,T,&theta[k][i],&fam[k][i],handle1);
						}
						for(t=0;t<*T;t++ ) 
						{
							tildevdirect[k-1][i][t]+=handle1[t]*tildezr2[t];
						}
					}
					if(conindirect[k+(*d)*i-1]==1)
					{
						for(t=0;t<*T;t++ ) 
						{
							tildevindirect[k-1][i][t]+=cop[t]*tildezr2[t];
						}
					}
					
				}
			}
			for(t=0;t<*T;t++ ) 
			{
				sumloglik += tildevalue[k][i][t];
			}  
		}
	}//for loops closed

	*out = sumloglik;

	for(i=0;i<(*d);i++)
	{
		for(j=0;j<(*d);j++)
		{	
			for(t=0;t<*T;t++ ) 
			{
				tilde_vdirect[(i+1)+(*d)*j+(*d)*(*d)*t-1]=tildevdirect[i][j][t];
				tilde_vindirect[(i+1)+(*d)*j+(*d)*(*d)*t-1]=tildevindirect[i][j][t];
				tilde_value[(i+1)+(*d)*j+(*d)*(*d)*t-1]=tildevalue[i][j][t];		
			}
		}
	}



	//Free memory:
	free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d); 
	free_intmatrix(calc, *d); 
	free_3darray(tildevindirect,*d,*d); 
	free_3darray(tildevdirect,*d,*d); 
	free_3darray(tildevalue,*d,*d); 
	Free(zr1); Free(zr2); Free(tildezr1); Free(tildezr2); Free(handle1); Free(cop);

}


/////////////////////////////////////////////////
// Calculate the gradient
// (uses the function VineLogLikRvineDeriv)
//
// Input:
// see above
//
// Output:
// out		gradient vector
///////////////////////////////////////////////////


void VineLogLikRvineGradient(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* out, double* ll, double* vv, double* vv2, int* posParams) 
						  
{
	int kk, ii, tt, i, j, tcop=0, dd=1, aa=0, margin=0;
	int *calc;
	calc = Calloc(((*d)*(*d)),int);
	double *tilde_vdirect, *tilde_vindirect, *tilde_value;
	tilde_vdirect = Calloc(((*d)*(*d)*(*T)),double);
	tilde_vindirect = Calloc(((*d)*(*d)*(*T)),double);
	tilde_value = Calloc(((*d)*(*d)*(*T)),double);
	int **pospar, **fam;
	pospar=create_intmatrix(*d,*d);
	fam=create_intmatrix(*d,*d);

	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++)
		{
			pospar[i][j]=posParams[(i+1)+(*d)*j-1] ;
			fam[i][j]=family[(i+1)+(*d)*j-1] ;
			if( j < i && pospar[i][j]==1){
				aa++;
			}
		}
	}

    tt=0;
	for(ii=(*d-1);ii>0;ii--)
	{
		for(kk=(*d);kk>ii;kk--)
		{
			if(pospar[kk-1][ii-1]==1)
			{
				
				calcupdate_func(d, matrix, &kk, &ii, calc);
				
				if(fam[kk-1][ii-1]==2)		// for the t-copula
				{
					tcop=1;		// first parameter
					VineLogLikRvineDeriv(T, d, family, &kk, &ii, maxmat, matrix, condirect, conindirect, par, par2, data, &out[tt], ll, vv, vv2, calc, tilde_vdirect, tilde_vindirect, tilde_value, &tcop, &margin);
					tcop=2;		// second parameter
					VineLogLikRvineDeriv(T, d, family, &kk, &ii, maxmat, matrix, condirect, conindirect, par, par2, data, &out[aa-1+dd], ll, vv, vv2, calc, tilde_vdirect, tilde_vindirect, tilde_value, &tcop, &margin);		// important: position in the gradient out[aa-1+dd]
					dd++;
				}
				else
				{
					tcop=0;
					VineLogLikRvineDeriv(T, d, family, &kk, &ii, maxmat, matrix, condirect, conindirect, par, par2, data, &out[tt], ll, vv, vv2, calc, tilde_vdirect, tilde_vindirect, tilde_value, &tcop, &margin);
				}
				
				
				tt+=1;
			}
		}
	}

Free(calc);free_intmatrix(pospar,*d);free_intmatrix(fam,*d);
Free(tilde_vdirect);Free(tilde_vindirect);Free(tilde_value);
}

