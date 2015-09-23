#if !defined(RVINEDERIV_H)
#define RVINEDERIV_H

//////////////////////////////////////////////////////////////
// Function to compute the derivative of log-likelihood for the pair-copula construction (Rvine)
// (by J.S. and U.S.)
//
// Input:
// T						sample size
// d						dimension (>=2)
// family					copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// kk						row number of the parameter wrt to which we want to calculate the derivative
// ii						column number of the parameter wrt to which we want to calculate the derivative
// par						parameter values (at least d*(d-1)/2 parameters
// par2						second set of parameter values (f.e. for student copulas)
// data						data set for which to compute the derivative of the log-likelihood
// matrix					an RVineMatrix in vector form
// condirect, conindirect	Matrizes which tell us where we find the right values 
// calcupdate				matrix which tells which terns we need to consider for the calculation of the derivative
// tcop						variable for the t-copula
//
// Output:
// out			Loglikelihood
// ll			array with the contribution to the derivative (for each copula)
// vv,vv2       array for the derivatives of the h-functions  
/////////////////////////////////////////////////////////////

void VineLogLikRvineDeriv(int* T, int* d, int* family, int* kk, int* ii, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* out, double* ll, double* vv, double* vv2, int* calcupdate, double* tilde_vdirect, double* tilde_vindirect, double* tilde_value, int* tcop, int* margin);


//////////////////////////////////////////////////////////////
// Function to compute the gradient for an R-vine
//////////////////////////////////////////////////////////////

void VineLogLikRvineGradient(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* out, double* ll, double* vv, double* vv2, int* posParams);
						  //double* tilde_vdirect, double* tilde_vindirect, double* tilde_value);

void VineLogLikRvineGradient2(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* out, double* ll, double* vv, double* vv2, int* posParams);

#endif
