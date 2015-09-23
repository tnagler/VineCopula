#if !defined(RVINEDERIV2_H)
#define RVINEDERIV2_H

//////////////////////////////////////////////////////////////
// Function to compute the second derivative of log-likelihood for the pair-copula construction (Rvine)
// (by J.S. and U.S.)
//
// Input:
// T				sample size
// d				dimension (>=2)
// family			copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// kk				row number of the parameter wrt to which we want to calculate the derivative
// ii				column number of the parameter wrt to which we want to calculate the derivative
// kkk				row number of the parameter wrt to which we want to calculate the derivative (second parameter)
// iii				column number of the parameter wrt to which we want to calculate the derivative (second parameter)
// maxmat			Maximums-Matrix from RVM object
// matrix			an RVineMatrix in vector form
// condirect, conindirect	Matrizes which tell us where we find the right values 
// par				parameter values (at least d*(d-1)/2 parameters
// par2				second set of parameter values (f.e. for student copulas)
// data				data set for which to compute the derivative of the log-likelihood
// calcupdate, calcupdate2	matrix which tells which terns we need to consider for the calculation of the derivative
// ll_tilde
// vv_tilde, vv2_tilde
// ll_hat
// vv_hat, vv2_hat
// ll
// vv, vv2
// tcop				variable for the t-copula
//
// Output:
// out					second deriv of Loglikelihood
// barvalue_out				array with the contribution to the derivative (for each copula)
// barvdirect_out, barvindirect_out	array for the derivatives of the h-functions  
/////////////////////////////////////////////////////////////

void VineLogLikRvineDeriv2(int* T, int* d, int* family, int* kk, int* ii, int* kkk, int* iii, int* maxmat, int* matrix, int* condirect, 
			int* conindirect, double* par, double* par2, double* data, double* ll_tilde, double* vv_tilde, double* vv2_tilde, 
			double* ll_hat, double* vv_hat, double* vv2_hat, int* calcupdate, int* calupdate2, double* out, double* ll, 
			double* vv, double* vv2, double* barvalue_out, double* barvdirect_out, double* barvindirect_out, int* tcop, int* kk_second);

void hesse_step(int* T, int* d, int* family, int* kk, int* ii, int* kkk, int* iii, int* maxmat, int* matrix, int* condirect, int* conindirect, 
				double* par, double* par2, double* data, int* calcupdate, int* calcupdate2, double* out, double* ll, double* vv, double* vv2,
				double* tilde_value, double* tilde_vdirect, double* tilde_vindirect, double* hat_value, double* hat_vdirect, double* hat_vindirect,
				double* barvalue_out, double* barvdirect_out, double* barvindirect_out, int* kk_second, int* kkk_second);

void hesse(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
			double* out, double* subhess, double* der, double* subder);

void calcupdate_func(int* d, int* matrix, int* i, int* j, int* calc);

#endif
