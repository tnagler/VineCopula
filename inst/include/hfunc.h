#if !defined(HFUNC_H)
#define HFUNC_H

// File hfunc.c
void Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out);
void Hfunc1_vec(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out);
void Hfunc2(int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);
void Hfunc2_vec(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out);
void Hfunc(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv1(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv1_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv2(int* family, int* n, double* v, double* u, double* theta, double* nu, double* out);
void Hinv2_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out);
void pcondbb1(double* u, double* v, int* n, double* param, double* out);
void pcondbb6(double* u, double* v, int* n, double* param, double* out);
void pcondbb7(double* u, double* v, int* n, double* param, double* out);
void pcondbb8(double* u, double* v, int* n, double* param, double* out);
void qcondgum(double* q, double* u, double* de, double* out);
void qcondjoe(double* q, double* u, double* de, double* out);
//void qcondbb1(double* q, double* u, double* de, double* th, double* out);

#endif
