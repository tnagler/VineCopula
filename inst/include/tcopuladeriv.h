#if !defined(TCOPULA_H)
#define TCOPULA_H

// File tcopuladeriv.c
void diffPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_rho_tCopula_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diffPDF_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_rho_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_rho_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out);

#endif
