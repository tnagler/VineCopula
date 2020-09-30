//////////////////////////////////////////////////
// PIT - Probability integral transform			//
//												//
// by Ulf Schepsmeier (2012)					//
//////////////////////////////////////////////////

#include "VineCopula/vine.h"
#include "VineCopula/memoryhandling.h"
#include "VineCopula/likelihood.h"
#include "VineCopula/pit.h"
#include "VineCopula/hfunc.h"

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4

//////////////////////////////////////////////////////////////
// Probability integral transform for the C- and D-vine
// Input:
// n         sample size
// d         dimension (>= 2)
// type      vine type (1=Canonical vine, 2=D-vine)
// family    copula family
// par       parameter values (at least d*(d-1)/2 parameters)
////////////////////////////////////////////////////////////////

// The algorithm is based on the pseudo algorithm of Aas et al. (2009)

void pit(int* T, int* d, int* family, int* type, double* par, double* nu, double* data, double* out)
{
	int i, j, in=1, k, **fam, tt;
	double **v, **theta, **z, **ny, **x;

	x = create_matrix(*d+1,*T);
	v = create_matrix(*d+1,2*(*d)-1);
	theta = create_matrix(*d+1,*d+1);
	z = create_matrix(*d+1,*T);
	ny = create_matrix(*d+1,*d+1);
	fam = create_intmatrix(*d+1,*d+1);

	k = 0;
	for(i=0;i<*d;i++)
    {
		for (tt=0;tt<=*T-1;tt++ )
		{
			x[i][tt] = data[k];
			k++;
		}
    }

	//Initialize dependency parameters
	k = 0;
	for(i=0;i<*d-1;i++)
	{
		for(j=0;j<(*d-i-1);j++)
		{
		  fam[i][j] = family[k];
		  ny[i][j] = nu[k];
		  theta[i][j] = par[k];
		  k++;
		}
	}

	// Transform
	if(*type==1) //Canonical vine
	{
		for(j=0;j<*T;j++)
		{
		  z[0][j] = x[0][j];
		  for(i=1;i<*d;i++)
		  {
			z[i][j]=x[i][j];
			for(k=0;k<=(i-1);k++)
			{
			  Hfunc1(&fam[k][i-k-1],&in, &z[i][j],&z[k][j],&theta[k][i-k-1],&ny[k][i-k-1],&z[i][j]);
			}
		  }
		}
	}
	else if(*type==2) //D-vine
	{
		for(j=0;j<*T;j++)
		{
			z[0][j] = x[0][j];
			Hfunc1(&fam[0][0],&in, &x[1][j],&x[0][j],&theta[0][0],&ny[0][0],&z[1][j]);
			v[1][0] = x[1][j];
			Hfunc2(&fam[0][0],&in, &x[0][j],&x[1][j],&theta[0][0],&ny[0][0],&v[1][1]);
			for(i=2;i<*d;i++)
			{
				Hfunc1(&fam[0][i-1],&in, &x[i][j],&x[i-1][j],&theta[0][i-1],&ny[0][i-1],&z[i][j]);
				for(k=1;k<=(i-1);k++)
				{
					Hfunc1(&fam[k][i-k-1],&in, &z[i][j],&v[i-1][2*(k-1)+1],&theta[k][i-k-1],&ny[k][i-k-1],&z[i][j]);
				}
				if(i==(*d-1))
					break;

				v[i][0] = x[i][j];
				Hfunc2(&fam[0][i-1],&in, &v[i-1][0],&v[i][0],&theta[0][i-1],&ny[0][i-1],&v[i][1]);
				Hfunc1(&fam[0][i-1],&in, &v[i][0],&v[i-1][0],&theta[0][i-1],&ny[0][i-1],&v[i][2]);
				if(i>2)
				{
					for(k=0;k<=(i-3);k++)
					{
					  Hfunc2(&fam[k+1][i-k-2],&in, &v[i-1][2*k+1],&v[i][2*k+2],&theta[k+1][i-k-2],&ny[k+1][i-k-2],&v[i][2*k+3]);
					  Hfunc1(&fam[k+1][i-k-2],&in, &v[i][2*k+2],&v[i-1][2*k+1],&theta[k+1][i-k-2],&ny[k+1][i-k-2],&v[i][2*k+4]);
					}
				}
				Hfunc2(&fam[i-1][0],&in, &v[i-1][2*i-3],&v[i][2*i-2],&theta[i-1][0],&ny[i-1][0],&v[i][2*i-1]);
			}
		}
	}

	//Write to output vector:
	k = 0;
	for(i=0;i<*d;i++)
	{
		for(j=0;j<*T;j++)
		{
		  out[k] = z[i][j];
		  k ++;
		}
	}

	//Free memory:
	free_matrix(x,*d+1); free_matrix(v,*d+1); free_matrix(theta,*d+1); free_matrix(ny,*d+1); free_intmatrix(fam,*d+1); free_matrix(z,*d+1);
}



//////////////////////////////////////////////////////////////
// Probability integral transform for the R-vine
//
// Input:
// T, d			dimensions of the data
// family,...	RVM objects
// data			data
// vv, vv2		h-functions derived bei the likelihood function
// calcupdate	which h-functions, inverse h-functions have to be derived
//
// Output:
// out			PIT
//////////////////////////////////////////////////////////////

// Reference: Schepsmeier (2015)

// Ulf Schepsmeier, Efficient information based goodness-of-fit tests for vine copula models with fixed margins: A comprehensive review,
// Journal of Multivariate Analysis, Available online 14 January 2015, ISSN 0047-259X, https://dx.doi.org/10.1016/j.jmva.2015.01.001.


void RvinePIT(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		 double* out, double* vv, double* vv2, int* calcupdate)
{
	int i, j, k, t, m, **fam;
	double **x, **theta, **nu, ***vdirect, ***vindirect, **z;

	//Allocate memory
	x = create_matrix(*d,*T);
	vdirect = create_3darray(*d,*d,*T);
	vindirect = create_3darray(*d,*d,*T);
	theta=create_matrix(*d,*d);
	nu=create_matrix(*d,*d);
	fam=create_intmatrix(*d,*d);
	z = create_matrix(*d,*T);

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

	// From vector to array
	k=0;
	for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta[i][j]=par[(i+1)+(*d)*j-1] ;
            nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
            fam[i][j]=family[(i+1)+(*d)*j-1] ;
			for(t=0;t<*T;t++ )
			{
				vdirect[i][j][t]=vv[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				vindirect[i][j][t]=vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1];
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

	// First column is easy; it's the data
	for(t=0;t<*T;t++)
	{
		z[0][t]=x[0][t];
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
					//if(condirect[k+(*d)*i-1]==1 || (i==0 && k==1)) // oder man ist im letzten Fall (i==0 && k==1), der nicht mehr in CondDist$direkt drin ist
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
					//if(condirect[k+(*d)*i-1]==1 || (i==0 && k==1))
					{
						Hfunc1(&fam[k][i],T,vdirect[k][i],vindirect[k][(*d-m)],&theta[k][i],&nu[k][i],vdirect[k-1][i]);
					}
					if(conindirect[k+(*d)*i-1]==1)
					{
						Hfunc2(&fam[k][i],T,vindirect[k][(*d-m)],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k-1][i]);
					}
				}

				//Rprintf("d-i: %d\t k-1: %d\t i: %d\n",*d-i, k-1, i);
				//Rprintf("vindirect[k-1][i]: %f\t conindirect[k+(*d)*i-1]: %d \n", vindirect[k-1][i][1], conindirect[k+(*d)*i-1]);
				for(t=0;t<*T;t++)
				{
					z[*d-i-1][t]=vdirect[k-1][i][t];
				}

			}

		}
	}

	//Write to output vector:

	k = 0;
	for(i=0;i<*d;i++)
	{
		for(j=0;j<*T;j++)
		{
		  out[k] = z[i][j];
		  k ++;
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
			}
		}
	}

	//Free memory:
	free_matrix(x,*d);
	free_matrix(z,*d);
	free_3darray(vdirect,*d,*d);
	free_matrix(theta,*d);
	free_matrix(nu,*d);
	free_intmatrix(fam,*d);
	free_3darray(vindirect,*d,*d);
}

