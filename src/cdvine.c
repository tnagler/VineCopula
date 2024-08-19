/*
** cdvine.c - C code of the package CDRVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/


// Include all the head files
#include "VineCopula/vine.h"			// general one
#include "VineCopula/memoryhandling.h"	// for creating two and three dimensional arrays
#include "VineCopula/likelihood.h"		// formally main functionality; log-likelihood with help functions; bivariate densities
#include "VineCopula/cdvine.h"			// Header file for this C-file
#include "VineCopula/hfunc.h"			// h-functions, i.e. conditional densities; also inverse h-functions

#define UMAX  1-1e-12

#define UMIN  1e-12

#define XEPS 1e-4


//////////////////////////////////////////////////////////////
// Function to simulate from a C- or D-vine
// Input:
// n         sample size
// d         dimension (>= 2)
// type      vine type (1=Canonical vine, 2=D-vine)
// family    copula family (see help pages which families are now included)
// par       parameter values (at least d*(d-1)/2 parameters)
// nu		 second parameter for t-copula, BB-copulas and Tawn
//
// Output:
// out		 two dimensional array of simulated data
////////////////////////////////////////////////////////////////

void pcc(int* n, int* d, int* family, int* type, double* par, double* nu, double* out)
{
  int i, j, in=1, k, **fam;
  double *w, **v, t, **theta, **x, **ny;

  GetRNGstate();		//Init random number generator
  //Allocate memory:
  w = R_Calloc((*d+1),double);

  v = create_matrix(*d+1,2*(*d)-1);
  theta = create_matrix(*d,*d);
  x = create_matrix(*n+1,*d+1);
  ny = create_matrix(*d,*d);
  fam = create_intmatrix(*d,*d);
  //Initialize dependency parameters
  
  // The function arguments are one-dimensional vectors; for better understanding the transform them back to matrices (see theory)
  // This step may be updated in the future to optimize the algorithms
  k = 0;
  for(i=1;i<=*d-1;i++)
  {
    for(j=1;j<=*d-i;j++)
    {
      fam[i][j] = family[k];
      ny[i][j] = nu[k];
      theta[i][j] = par[k];
      k ++;
    }
  }
  //Simulate: (it follows the theoretical algorithm)
  if(*type==1) //Canonical vine
  {
    for(j=1;j<=*n;j++)		// run over all observations (rows)
    {
      for(i=1;i<=*d;i++) w[i] = runif(0,1);
      x[j][1] = w[1];
      for(i=2;i<=*d;i++)	// run over all dimensions (cols)
      {
        t = w[i];
        for(k=i-1;k>=1;k--)
        {
		  Hinv1(&fam[k][i-k],&in, &t,&w[k],&theta[k][i-k],&ny[k][i-k],&t);
        }
        x[j][i] = t;
      }
    }
  }
  else if(*type==2) //D-vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) { w[i] = runif(0,1);} 
      v[1][1] = w[1];
      v[2][1] = w[2];
      
      Hinv1(&fam[1][1],&in,&w[2],&v[1][1],&theta[1][1],&ny[1][1],&v[2][1]);
      Hfunc2(&fam[1][1],&in, &v[1][1],&v[2][1],&theta[1][1],&ny[1][1],&v[2][2]);
      for(i=3;i<=*d;i++)
      {
        v[i][1] = w[i];
	
        for(k=i-1;k>=2;k--) 
		{ 
			Hinv1(&fam[k][i-k],&in, &v[i][1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][1]);
		}
        Hinv1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][1]);
		
        // Compute conditional cdf's needed in next step:
        if(i<*d)
        {
          Hfunc2(&fam[1][i-1],&in, &v[i-1][1],&v[i][1],&theta[1][i-1],&ny[1][i-1],&v[i][2]);
          Hfunc1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][3]);
          if(i>3)
          {
            for(k=2;k<=(i-2);k++)
            {
              Hfunc2(&fam[k][i-k],&in, &v[i-1][2*k-2],&v[i][2*k-1],&theta[k][i-k],&ny[k][i-k],&v[i][2*k]);
              Hfunc1(&fam[k][i-k],&in, &v[i][2*k-1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][2*k+1]);
            }
          }
          Hfunc2(&fam[i-1][1],&in, &v[i-1][2*i-4],&v[i][2*i-3],&theta[i-1][1],&ny[i-1][1],&v[i][2*i-2]);
        }
      }
      for(i=1;i<=*d;i++) x[j][i] = v[i][1];
    }
  }
  //Write to output vector:
  k = 0;
  for(i=1;i<=*d;i++)
  {
    for(j=1;j<=*n;j++)
    {
      out[k] = x[j][i];
      k ++;
    }
  }
  PutRNGstate();		// Function for the random number generator
  //Free memory:
  R_Free(w); free_matrix(v,*d+1); free_matrix(theta,*d); free_matrix(ny,*d); free_intmatrix(fam,*d); free_matrix(x,*n+1);
}




//////////////////////////////////////////////////////////////
// Function to compute -log-likelihood for C- and D-vine
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula families
// par      parameter values (at least d*(d-1)/2 parameters ) The second parameter is added at the end of par
// data     data set for which to compute log-likelihood
// Output:
// out      Log-likelihood
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////
void VineLogLikm(int* T, int* d, int* type, int* family, double* par, double* data, 
		 double* out, double* ll, double* vv)
{
  int i, j, k, t, kk, **fam;
  double loglik=0.0, sumloglik=0.0, **x, **theta, **nu, ***v;
  
  //Allocate memory:
  x = create_matrix(*d+1,*T);
  // By Ulf Schepsmeier
  if(*type==1) //C-vine
	{
		v = create_3darray(*d-1,*d,*T);
	}
  else //D-vine
	{
		v = create_3darray(*d,2*(*d)-3,*T);
	}
  
  theta = create_matrix(*d,*d);
  nu = create_matrix(*d,*d);
  fam = create_intmatrix(*d+1,*d+1);
  
  //Initialize:
  k = 0;
  for(i=1;i<=*d;i++)
    { 
      for (t=0;t<=*T-1;t++ ) 
	{
	  x[i][t] = data[k];	//transform the data back into a 2-dim array
	  k++;
	}
    }
  k = 0;
  for(i=1;i<=(*d-1);i++)
    {
      for(j=1;j<=(*d-i);j++)
		{
			theta[i][j] = par[k];
			fam[i][j] = family[k];
			nu[i][j] = par[*d*(*d-1)/2+k];		// the second parameter is added at the end of par (not the best solution but was practise at the beginning)
			k++;
		}
    }

  if(*type==1) //C-vine
    {
      // By Ulf Schepsmeier
	  kk=0;
	  //Compute likelihood at level 1:
		for(i=1;i<*d;i++)
		{
			LL_mod2(&fam[1][i],T,x[1],x[i+1],&theta[1][i],&nu[1][i],&loglik);	// call the bivariate log-likelihood function 
			//(with the correct rotation for 90, 180 and 270 degrees)
			sumloglik += loglik;	// sum up
			ll[kk] = loglik; 		// store all bivariate log-likelihoods too
			++kk;
			if(*d>2)
			{
			//Compute variables for next level:
			Hfunc1(&fam[1][i],T,x[i+1],x[1],&theta[1][i],&nu[1][i],v[1][i]);
			}
		}
		//Compute likelihood at next levels:
		if(*d>2)
		{
			for(k=2;k<=(*d-1);k++)
			{
				for(i=1;i<=(*d-k);i++)
				{
					LL_mod2(&fam[k][i],T,v[k-1][1],v[k-1][i+1],&theta[k][i],&nu[k][i],&loglik);
					sumloglik += loglik;
					ll[kk] = loglik; 
					++kk;
				}
				if(k<(*d-1))
				{
					for (i=1;i<=(*d-k);i++)
					{
						Hfunc1(&fam[k][i],T,v[k-1][i+1],v[k-1][1],&theta[k][i],&nu[k][i],v[k][i]);
					}
				}
			}
		}
    }
  else if(*type==2) //D-vine
  {
    kk=0;
    //Compute the likelihood at level 1:
    for(i=1;i<*d;i++)
      {
        LL_mod2(&fam[1][i],T,x[i],x[i+1],&theta[1][i],&nu[1][i],&loglik);
        sumloglik += loglik;
		ll[kk] = loglik; 
		++kk;
      }
    //Compute variables for next level:
	if(*d>2)
	{
		Hfunc2(&fam[1][1],T,x[1],x[2],&theta[1][1],&nu[1][1],v[1][1]);
		for(k=1;k<=(*d-3);k++)
		  {
			Hfunc1(&fam[1][k+1],T,x[k+2],x[k+1],&theta[1][k+1],&nu[1][k+1],v[1][2*k]);
			Hfunc2(&fam[1][k+1],T,x[k+1],x[k+2],&theta[1][k+1],&nu[1][k+1],v[1][2*k+1]);
		  }
		Hfunc1(&fam[1][*d-1],T,x[*d],x[*d-1],&theta[1][*d-1],&nu[1][*d-1],v[1][2*(*d)-4]);
		//Compute likelihood at next levels:
		for(k=2;k<=(*d-1);k++)
		  {
			for(i=1;i<=(*d-k);i++)
			  {
				LL_mod2(&fam[k][i],T,v[k-1][2*i-1],v[k-1][2*i],&theta[k][i],&nu[k][i],&loglik);
				sumloglik += loglik;
				ll[kk] = loglik; ++kk;
			  }
			if(k<(*d-1))
			  {
				Hfunc2(&fam[k][1],T,v[k-1][1],v[k-1][2],&theta[k][1],&nu[k][1],v[k][1]);
				if((*d)>4)
				  {
					for(i=1;i<=(*d-k-2);i++)
					{
						Hfunc1(&fam[k][i+1],T,v[k-1][2*i+2],v[k-1][2*i+1],&theta[k][i+1],&nu[k][i+1],v[k][2*i]);
						Hfunc2(&fam[k][i+1],T,v[k-1][2*i+1],v[k-1][2*i+2],&theta[k][i+1],&nu[k][i+1],v[k][2*i+1]);
					}
				  }
				Hfunc1(&fam[k][*d-k],T,v[k-1][2*(*d)-2*k],v[k-1][2*(*d)-2*k-1],&theta[k][*d-k],&nu[k][*d-k],v[k][2*(*d)-2*k-2]);
			  }
		  }
	}
  }
  //Write to output:
  *out = -sumloglik;
  kk=00;
  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		if(*d>2)
		{
		for(k=1;k<(*d-1);k++)
			for(i=1;i<=(*d-k);i++)
				for(t=0;t<*T;t++)
				{
					vv[kk] = v[k][i][t];		// transformation from a 3-dim array to a vector
  					++kk;
  				}
		}
		//Free memory:
		free_3darray(v,*d-1,*d);
	}
	else //D-Vine
	{
		if(*d>2)
		{
		for(k=1;k<*d;k++)
			for(i=1;i<=2*(*d-k-1);i++)
				for(t=0;t<*T;t++)
  				{
  					vv[kk] = v[k][i][t];
  					++kk;
  				}
		}
		//Free memory:
		free_3darray(v,*d,2*(*d)-3);
	}
  //Free memory:
  free_matrix(x,*d+1); free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d+1);
}

//////////////////////////////////////////////////////////////
// Function to compute an update of the log-likelihood for C- and D-vine
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula families
// par      parameter values (at least d*(d-1)/2 parameters
// mpar     index of modified parameter (related to previous computation)
// data     data set for which to compute log-likelihood
// ll       array with the stored contribution of the likelihood in a previous computation
// vv       3d array  array with the stored transformations in a previous computation
// Output:
// out		log-likelihood (updated)
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////
void VineLogLikmP(int* T, int* d, int* type, int* family, double* par, int* mpar, double* data, 
		  double* out, double* ll, double* vv)
{
  int i, j, ii, jj, k, t, kk,**fam, **ind;
  double sumloglik=0.0,loglik=0.0,  **x, **theta, **nu, ***v;
  //Allocate memory:
  x = create_matrix(*d+1,*T);

  // By Ulf Schepsmeier
  if(*type==1) //C-vine
	{
		v = create_3darray(*d-1,*d,*T);
	}
  else //D-vine
	{
		v = create_3darray(*d,2*(*d)-3,*T);
	}
  
  theta = create_matrix(*d,*d);
  nu = create_matrix(*d,*d);
  fam = create_intmatrix(*d,*d);
  ind = create_intmatrix(*d,*d);
  //Initialize:
  k = 0;
  for(i=1;i<=*d;i++)
    { 
      for (t=0;t<=*T-1;t++ ) 
	{
	  x[i][t] = data[k];
	  k++;
	}
    }
  k = 0;
  jj = *d;
  ii = *d;
  kk=00;

  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		for(i=1;i<=(*d-1);i++)
		{
		  for(j=1;j<=(*d-1);j++)
			{
				ind[i][j] = 0;
				if(j <= *d-i)
				{
					++k;
					if (k == *mpar)
					{
						ii=i;
						jj=j;
					}
					if(ii+jj-i>0)
					{
						if(i>=ii && j==ii+jj-i)
						{
							ind[i][j]=1;
						}
					}
					else if(ii+jj-i<=0)
					{
						ind[i][j]=1;
					}
				}
			}
		}

		for(k=1;k<*d-1;k++)
			for(i=1;i<=*d-k;i++)
				for(t=0;t<*T;t++)
				{
					v[k][i][t] = vv[kk];
  					++kk;
  				}
	}
	else //D-Vine
	{
	  for(i=1;i<=(*d-1);i++)
		{
		  for(j=1;j<=(*d-1);j++)
			{
			ind[i][j] = 0;
			if(j <= *d-i)
			{
				++k;
				if (k == *mpar)
				{
					ii=i;
					jj=j;
				}            
				if(i >= ii &&  j >= ii+jj-i &&  j <= *d-i && j <= jj)
				{
					ind[i][j]=1;
				}
			//	    //printf("%d ", ind[i][j]);
			}
			}
		  //      //printf("\n");
		}

		for(k=1;k<*d;k++)
			for(i=1;i<=2*(*d-k-1);i++)
				for(t=0;t<*T;t++)
				{
					v[k][i][t] = vv[kk];
					++kk;
				}
	}
  k=0;
  for(i=1;i<=(*d-1);i++)
    {
      for(j=1;j<=(*d-i);j++)
	{
	  theta[i][j] = par[k];
	  fam[i][j] = family[k];
	  nu[i][j] = par[*d*(*d-1)/2+k];
	  k ++;
      }
    }
 
  

  if(*type==1) //C-vine
    {
      // By Ulf Schepsmeier
	  kk=0;
	  //Compute likelihood at level 1:
		for(i=1;i<*d;i++)
		{
			if(ind[1][i]==1) 
			{
				LL_mod2(&fam[1][i],T,x[1],x[i+1],&theta[1][i],&nu[1][i],&loglik);
				ll[kk] = loglik;
				if(*d>2)
				{
				//Compute variables for next level:
				Hfunc1(&fam[1][i],T,x[i+1],x[1],&theta[1][i],&nu[1][i],v[1][i]);
				}
			}
			sumloglik += ll[kk];
			++kk;
		}
		//Compute likelihood at next levels:
		if(*d>2)
		{
		for(k=2;k<=(*d-1);k++)
		{
			for(i=1;i<=(*d-k);i++)
			{
				if(ind[k][i]==1)
				{
					LL_mod2(&fam[k][i],T,v[k-1][1],v[k-1][i+1],&theta[k][i],&nu[k][i],&loglik);
					ll[kk] = loglik; 
					if(k<(*d-1))
					{
						Hfunc1(&fam[k][i],T,v[k-1][i+1],v[k-1][1],&theta[k][i],&nu[k][i],v[k][i]);
					}
				}
				sumloglik += ll[kk];
				++kk;
			}
		}
		}
    }
  else if(*type==2) //D-vine
  {
    kk=0;
    //Compute the likelihood at level 1:
    for(i=1;i<*d;i++)
      {
        if(ind[1][i]==1) 
		{
			LL_mod2(&fam[1][i],T,x[i],x[i+1],&theta[1][i],&nu[1][i],&loglik);
			ll[kk] = loglik;
		}
		sumloglik += ll[kk];
		++kk;
      }
    //Compute variables for next level:
	if(*d>2)
	{
    if(ind[1][1]==1) Hfunc2(&fam[1][1],T,x[1],x[2],&theta[1][1],&nu[1][1],v[1][1]);
    for(k=1;k<=(*d-3);k++)
      if(ind[1][k+1]==1)
		{
			Hfunc1(&fam[1][k+1],T,x[k+2],x[k+1],&theta[1][k+1],&nu[1][k+1],v[1][2*k]);
			Hfunc2(&fam[1][k+1],T,x[k+1],x[k+2],&theta[1][k+1],&nu[1][k+1],v[1][2*k+1]);
		}
    if(ind[1][*d-1]) Hfunc1(&fam[1][*d-1],T,x[*d],x[*d-1],&theta[1][*d-1],&nu[1][*d-1],v[1][2*(*d)-4]);
    //Compute likelihood at next levels:
    for(k=2;k<=(*d-1);k++)
      {
        for(i=1;i<=(*d-k);i++)
    	  {
	    if(ind[k][i] == 1)
	      {
		LL_mod2(&fam[k][i],T,v[k-1][2*i-1],v[k-1][2*i],&theta[k][i],&nu[k][i],&loglik);
		ll[kk] = loglik;
	      }
	    sumloglik += ll[kk];
	    ++kk;
    	  }
        if(k<(*d-1))
    	  {
    	    if(ind[k][1] == 1) Hfunc2(&fam[k][1],T,v[k-1][1],v[k-1][2],&theta[k][1],&nu[k][1],v[k][1]);
    	    if((*d)>4)
    	      {
    		for(i=1;i<=(*d-k-2);i++)
    		  {
		    if(ind[k][i+1]==1)
		      {
			Hfunc1(&fam[k][i+1],T,v[k-1][2*i+2],v[k-1][2*i+1],&theta[k][i+1],&nu[k][i+1],v[k][2*i]);
			Hfunc2(&fam[k][i+1],T,v[k-1][2*i+1],v[k-1][2*i+2],&theta[k][i+1],&nu[k][i+1],v[k][2*i+1]);
		      }
		  }
    	      }
    	    if(ind[k][*d-k]==1) Hfunc1(&fam[k][*d-k],T,v[k-1][2*(*d)-2*k],v[k-1][2*(*d)-2*k-1],&theta[k][*d-k],&nu[k][*d-k],v[k][2*(*d)-2*k-2]);
    	  }
      }
	}
  }
  //Write to output:
  *out = -sumloglik;
  kk=00;
  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		if(*d>2)
		{
		for(k=1;k<*d-1;k++)
			for(i=1;i<=*d-k;i++)
				for(t=0;t<*T;t++)
				{
					vv[kk] = v[k][i][t];
  					++kk;
  				}
		}
		//Free memory:
		free_3darray(v,*d-1,*d);
	}
	else //D-Vine
	{
		if(*d>2)
		{
		for(k=1;k<*d;k++)
			for(i=1;i<=2*(*d-k-1);i++)
				for(t=0;t<*T;t++)
  				{
  					vv[kk] = v[k][i][t];
  					++kk;
  				}
		}
		//Free memory:
		free_3darray(v,*d,2*(*d)-3);
	}
  //Free memory:
  free_matrix(x,*d+1); free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d);free_intmatrix(ind,*d);
}
