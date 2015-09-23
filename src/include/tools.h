#if !defined(TOOLS_H)
#define TOOLS_H

// File tools.c
void printError(char *text, char filename[200]);
double StableGammaDivision(double x1, double x2);
void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);
void ktau_matrix(double *data, int *d, int *N, double *out);		//New

void ADtest(double* cdf, int* n, double* out);		// Daniel Berg
void CumDist(double* x, int* i_n, int* i_m, double* out);		// Daniel Berg

#endif
