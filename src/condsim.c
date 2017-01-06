/*
** condsim.c - C code of the package CDRVine  
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
#include "condsim.h"
#include "hfunc.h"

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4

///////////////////////////////////////////////////////////////////

void condsim(int* n, int* d, int* d1, double* u1, int* family, double* par, double* nu, double* out)
{
  int i,j, k;
  double **uf,**ub,**th,**nuu;
  double aux;
  int **fam;
  uf = create_matrix(*d,*d);
  ub = create_matrix(*d,*d);
  th = create_matrix(*d+1,*d+1);
  nuu = create_matrix(*d+1,*d+1);
  fam = create_intmatrix(*d+1,*d+1);
  // param in matrices:
  k = 0;
  for(i=0;i<((*d)-1);i++)
    {
      for(j=0;j<((*d)-i-1);j++)
	{
	  fam[i][j] = family[k];
	  nuu[i][j] = nu[k];
	  th[i][j] = par[k];
	  k++;
	  //printf("%d \t",fam[i][j]);
	}
      //printf("\n");
    }
  // Simulation
  GetRNGstate();

	/*
	Declare variable to hold seconds on clock.
*/
//time_t seconds;
/*
Get value from system clock and
place in seconds variable.
*/
//time(&seconds);
/*
Convert seconds to a unsigned
integer.
*/
//srand((unsigned int) seconds);


  // for i = 0
  uf[0][0] = u1[0];
  ub[0][0] = u1[0];
  // for i = 1,... d1-1
  // compute uf and ub
  for (int i = 1; i < (*d1); ++i)
    {
      uf[i][i] = u1[i];
      ub[i][i] = u1[i];
      for (int j = (i-1); j >= 0; --j)
	{
	  Hfunc(&fam[i-j-1][j],n,&ub[i][j+1], &uf[i-1][j],&th[i-j-1][j],&nuu[i-j-1][j],&ub[i][j]); //backward
  	  //printf("ub: %d,%d : %d, %5.2f : %10.8f   \n",i,j,fam[i-j-1][j], th[i-j-1][j], ub[i][j]);
	} 
      //printf("\n");
      for (int j = 0; j <= i-1; ++j)
	{
	  Hfunc(&fam[i-j-1][j],n, &uf[i-1][j], &ub[i][j+1],&th[i-j-1][j],&nuu[i-j-1][j],&uf[i][j]); //forward
  	  //printf("uf: %d,%d : %d, %5.2f : %10.8f   \n",i,j,fam[i-j-1][j], th[i-j-1][j], uf[i][j]);
	}
      //printf("\n");
    }
  // for  i= d1,.. d-1
  for (int i = (*d1); i < (*d); ++i)
    {
      // (a) Simulate uniform
      //out[i-(*d1)] =  rand()/(RAND_MAX+1.0);
	  out[i-(*d1)]=runif(0,1);
      // (b) inverse transformation:
      for (int j = 0; j < i; ++j)
	{
  	  //printf("inv: %d,%d : %d, %5.2f : %10.8f   \t",i-j-1,j,fam[i-j-1][j], th[i-j-1][j], uf[i-1][j]);
	  Hinv(&fam[i-j-1][j], n, &out[i-*d1] , &uf[i-1][j], &th[i-j-1][j], &nuu[i-j-1][j],&aux );
	  out[i-(*d1)]  = aux;
	  //printf("%10.8f   \n ", aux);
	}
      //printf("\n");
      if (i <((*d)-1))
	{
	  // forward and backward steps:
	  uf[i][i] = out[i-(*d1)];
	  ub[i][i] = out[i-(*d1)];
	  for (int j = i-1; j >= 0; --j)
	    {
	      Hfunc(&fam[i-j-1][j],n,&ub[i][j+1], &uf[i-1][j],&th[i-j-1][j],&nuu[i-j-1][j],&ub[i][j]); //backward
	      //printf("ub: %d,%d : %d, %5.2f : %10.8f   \n",i-j-1,j,fam[i-j-1][j], th[i-j-1][j], ub[i][j]);
	    } 
	  //printf("\n");
	  for (int j = 0; j <= i-1; ++j)
	    {
	      Hfunc(&fam[i-j-1][j],n, &uf[i-1][j], &ub[i][j+1],&th[i-j-1][j],&nuu[i-j-1][j],&uf[i][j]); //forward
	      //printf("uf: %d,%d : %d, %5.2f : %10.8f   \n",i-j-1,j,fam[i-j-1][j], th[i-j-1][j], uf[i][j]);
	    } 
	  //printf("\n");
	}
    }
  // free memory
  free_matrix(th,*d);    
  free_matrix(ub,*d);    
  free_matrix(uf,*d);    
  free_matrix(nuu,*d);    
  free_intmatrix(fam,*d);
  PutRNGstate();
}
