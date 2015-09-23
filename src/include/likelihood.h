#if !defined(LIKELIHOOD_H)
#define LIKELIHOOD_H

// File likelihood.c
//void gen(double* u, int* n, double* param, int* copula, double* out);
//void genInv(double* u, int* n, double* param, int* copula, double* out);
//void genDrv(double* u, int* n, double* param, int* copula, double* out);
//void genDrv2(double* u, int* n, double* param, int* copula, double* out);
//void copCdf(double* u, double* v, int* n, double* param, int* copula, double* out);
//void copPdf(double* u, double* v, int* n, double* param, int* copula, double* out);
void archCDF(double* u, double* v, int* n, double* param, int* copula, double* out);
void LL_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void LL_mod2(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void LL(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void copLik(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik);
void copLik_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik);
void LL_mod_seperate(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void LL_mod_seperate_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void dbb1(double* u, double* v, int* n, double* param, double* out);
void dbb6(double* u, double* v, int* n, double* param, double* out);
void dbb7(double* u, double* v, int* n, double* param, double* out);
void dbb8(double* u, double* v, int* n, double* param, double* out);

#endif
