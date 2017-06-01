/*
** cdvine.h - C code of the package CDRVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#if !defined(CDVINE_H)
#define CDVINE_H

//////////////////////////////////////////////////////////////
// Function to simulate from a pair-copula construction (vine)
// Input:
// n         sample size
// d         dimension (>= 2)
// type      vine type (1=Canonical vine, 2=D-vine)
// family    copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB6, 9=BB7, 10=BB8)
// par       parameter values (at least d*(d-1)/2 parameters)
////////////////////////////////////////////////////////////////

void pcc(int* n, int* d, int* family, int* type, double* par, double* nu, double* out);


//////////////////////////////////////////////////////////////
// Function to compute -log-likelihood for the pair-copula construction (vine)
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB6, 9=BB7, 10=BB8)
// par      parameter values (at least d*(d-1)/2 parameters
// data     data set for which to compute log-likelihood
// Output:
// out      Loglikelihood
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////

void VineLogLikm(int* T, int* d, int* type, int* family, double* par, double* data, double* out, double* ll, double* vv);


//////////////////////////////////////////////////////////////
// Function to compute -log-likelihood for the pair-copula construction (vine) 
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB6, 9=BB7, 10=BB8)
// par      parameter values (at least d*(d-1)/2 parameters
// mpar     index of modified parameter (related to previous computation)
// data     data set for which to compute log-likelihood
// ll       array with the stored contribution of the likelihood in a previous computation
// vv       3d array  array with the stored transformations in a previous computation
// Output:
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////

void VineLogLikmP(int* T, int* d, int* type, int* family, double* par, int* mpar, double* data,  double* out, double* ll, double* vv);

#endif
