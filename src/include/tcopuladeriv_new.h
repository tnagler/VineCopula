#if !defined(TCOPULA_NEW_H)
#define TCOPULA_NEW_H

// File tcopuladeriv_new.c

void diffX_nu_tCopula(double* x, double* param, double* out);
void integral1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
void integral2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
void integral3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
void bigI_nu(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void diff_int_nu(double* x, double* nu, double* out);
void diff_t_nu(double* x, double* nu, double* out);
void diff_dt_x(double* x, double* nu, double* out);
void diff_dt_nu(double* x, double* nu, double* out);
void diff_dt_u(double* x, double* nu, double* out);
void diff_t_nu_nu(double* x, double* nu, double* out);
void diff2_x_nu(double* x, double* nu, double* out);

void diffPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diffPDF_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_rho_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2PDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2PDF_rho_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_u_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2PDF_nu_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_v_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2PDF_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);

void diffhfunc_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diffhfunc_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2hfunc_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_rho_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2hfunc_rho_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_nu_v_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out);
void diff2PDF_rho_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_nu_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out);

#endif
