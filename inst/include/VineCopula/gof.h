#if !defined(GOF_H)
#define GOF_H

void ADtest(double* cdf, int* n, double* out);					// Daniel Berg
void CumDist(double* x, int* i_n, int* i_m, double* out);		// Daniel Berg
void CvMtest(double* cdf, int* n, double* out); 				// Daniel Berg
void KStest(double* cdf, int* n, double* out); 					// Daniel Berg

void White(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, double* D, double* V);	

void Bj(int *T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
		 double* out, double* vv, double* vv2, int* calcupdate, int* method, int *alpha);
void SimulateBj(double* S, int *T, int* d, int* B, int* method, int *alpha, double* p);
void gofPIT_AD(int *T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
		 double* statistic, double* vv, double* vv2, int* calcupdate, int* method, int *alpha, int* B, int *statisticName);
void gofPIT_AD_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		double* statistic, double* vv, double* vv2, int* calcupdate, int* method, int* alpha, int* B, double* pvalue, int *statisticName);

void MySample(int *k, int *n, int *y);

void gofECP(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* data, double* statistic, int* statisticName);
void gofECP_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* data, double* statistic, int* statisticName, double* pvalue, int* B);
void gofECP2(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
		double* vv, double* vv2, int* calcupdate, double* statistic, int* statisticName);
void gofECP2_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
		double* vv, double* vv2, int* calcupdate, double* statistic, double* pvalue,int* statisticName, int* B);

void ChatZj(double* data, double* u, int* T, int* d, int* m, double* Chat);
void C_ind(double* data, int* n, int* d, double* C);

#endif
