#if !defined(INCOMPLETEBETA_H)
#define INCOMPLETEBETA_H

void incompleBeta_an1_bn1_p(double* x, double p, double q, double* an, double* bn);
void incompleBeta_an1_bn1_q(double* x, double p, double q, double* an, double* bn);
void incompleBeta_an_bn_p(double* x, double p, double q, int n, double* an, double* bn);
void incompleBeta_an_bn_q(double* x, double p, double q, int n, double* an, double* bn);
void inbeder(double* x_in, double* p_in, double* q_in, double* der);

#endif
