#if !defined(LOGDERIV_H)
#define LOGDERIV_H

void difflPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void difflPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void difflPDF_rho_tCopula_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void difflPDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2lPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2lPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2lPDF_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);

void difflPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void difflPDF_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void difflPDF(double* u, double* v, int* n, double* param, int* copula, double* out);

void diff2lPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2lPDF(double* u, double* v, int* n, double* param, int* copula, double* out);

#endif
