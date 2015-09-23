#if !defined(EVCOPULA_H)
#define EVCOPULA_H

// Some function for the Tawn copula

//CDF
void ta(double* t, int* n, double* par, double* par2, double* par3, double* out);
void Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out);
void TawnCDF(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);

///////////////////////////////////////////////////
//PDF
void ta2(double* t, int* n, double* par, double* par2, double* par3, double* out);
void d1ta(double* t, int* n, double* par, double* par2, double* par3, double* out);
void d2ta(double* t, int* n, double* par, double* par2, double* par3, double* out);

void Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out);
void d1Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out);
void d2Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out);
void dA_du(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);
void dA_dv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);
void dA_dudv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);

void TawnC(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);
void dC_du(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);
void TawnPDF(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);

// h-function
void dC_dv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out);

#endif
