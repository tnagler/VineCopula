#if !defined(RVINE_H)
#define RVINE_H

////////////////////////////////////////////
// PCC for R-vine
////////////////////////////////////////////

void SimulateRVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out, double* U, int* takeU);


//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for the pair-copula construction (Rvine)
// (by J.S.)
//
// Input:
// T				sample size
// d				dimension (>=2)
// family			copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// par				parameter values (at least d*(d-1)/2 parameters
// par2				second set of parameter values (f.e. for student copulas)
// data				data set for which to compute log-likelihood
// matrix			an RVineMatrix in vector form
// condirect, conindirect	Matrizes which tell us where we find the right values
// seperate			Control Parameter, do we want to seperate the likelihoods for each data point?
// calcupdate			matrix which tells us for which parameters we need to redo the calculations, not newly computed values are taken from ll, vv, vv2
// seperate			Control Parameter, do we want to seperate the likelihoods for each data point?
//
// Output:
// out		Loglikelihood
// ll		array with the contribution to LL (for each copula)
// vv,vv2       array for the transformation operated (Hfunc)
/////////////////////////////////////////////////////////////

void VineLogLikRvine(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		 double* out, double* ll, double* vv, double* vv2, int* calcupdate, int* seperate);


//////////////////////////////////////////////////////////////
// Function to compute likelihood for the pair-copula construction (Rvine)
// (by J.S.)
//
// Input:
// T				sample size
// d				dimension (>=2)
// family			copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// par				parameter values (at least d*(d-1)/2 parameters
// par2				second set of parameter values (f.e. for student copulas)
// data				data set for which to compute log-likelihood
// matrix			an RVineMatrix in vector form
// condirect, conindirect	Matrizes which tell us where we find the right values
// seperate			Control Parameter, do we want to seperate the likelihoods for each data point?
// calcupdate			matrix which tells us for which parameters we need to redo the calculations, not newly computed values are taken from ll, vv, vv2
//
// Output:
// out		Loglikelihood
// ll		array with the contribution to LL (for each copula)
// vv,vv2       array for the transformation operated (Hfunc)
/////////////////////////////////////////////////////////////

void VineLikRvine(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		double* out, double* ll, double* vv, double* vv2, int* calcupdate, int* seperate);

//////////////////////////////////////////////////////////////
// Function to construct the R-vine matrix from the binary matrix
// Input:
// b        the binary matrix
// Output:
// out      an Rvine matrix
/////////////////////////////////////////////////////////////
void getRVM(int* b, int* d, int* RVM);

#endif
