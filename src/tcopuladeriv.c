/*
** tcopuladeriv.c - C code of the package CDRVine
**
** by Ulf Schepsmeier
**
** derivatives of the t-copula (density and h-function)
**
*/

#include "VineCopula/vine.h"
#include "VineCopula/likelihood.h"
#include "VineCopula/hfunc.h"
#include "VineCopula/deriv.h"
#include "VineCopula/tcopuladeriv.h"

#define UMAX  1-1e-12
#define UMIN  1e-12

///////////////////////////////////
// Ableitungen f?r die t-copula
///////////////////////////////////

// Ableitung von c nach rho

void diffPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j, k=1;
	double t1, t2, t3, t4, t5, t6, t7, t10, t11, c;

	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		t1 = qt(u[j],nu,1,0);
		t2 = qt(v[j],nu,1,0);
		t3 = -(nu+2.0)/2.0;
		t10 = nu*(1.0-rho*rho);
		t4 = -2.0*t1*t2/t10;
		t11 = (t1*t1+t2*t2-2.0*rho*t1*t2);
		t5 = 2.0*t11*rho/t10/(1.0-rho*rho);
		t6 = 1.0+(t11/t10);
		t7 = rho/(1.0-rho*rho);
		out[j]=c*(t3*(t4+t5)/t6 + t7);
	}
}

// vectorized version
void diffPDF_rho_tCopula_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffPDF_rho_tCopula(&u[i], &v[i], &nn, ipars, copula, &out[i]);
    };
    free(ipars);
}

// Ableitung von c nach nu mit finiten Differenzen

void diffPDF_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, nu1, nu2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;
	nu1=param[1]-ep;
	nu2=param[1]+ep;

	for(j=0;j<*n;j++)
	{
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &nu1, &out1[j]);
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &nu2, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


/////////////////////////////////

//Ableitung von c nach u

void diffPDF_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, u1, u2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		u1=u[j]-ep;
		u2=u[j]+ep;
		copLik_mod(copula, &k, &u1, &v[j], &param[0], &param[1], &out1[j]);
		copLik_mod(copula, &k, &u2, &v[j], &param[0], &param[1], &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


////////////////////////////


//2. Ableitung nach rho

void diff2PDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j, k=1;
	double c, diff, M, t1, t2, t3, t4, t5, t6, t7;
	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param, copula, &diff);
		t1 = qt(u[j],nu,1,0);
		t2 = qt(v[j],nu,1,0);
		t4 = 1.0-rho*rho;
		M = ( nu*t4 + t1*t1 + t2*t2 - 2.0*rho*t1*t2 );

		t3 = -(nu+1.0)*(1.0+rho*rho)/t4/t4;
		t5 = (nu+2.0)*nu/M;
		t6 = 2.0*(nu+2.0)*pow(nu*rho+t1*t2,2.0)/M/M;
		t7 = diff/c;

		out[j] = c * (t3 + t5 + t6 + t7*t7);
	}
}


// 2. Ableitung von c nach nu

void diff2PDF_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, *out3, nu1=param[1]-ep, nu2=param[1]+ep;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	out3 = (double*) Calloc(*n,double);
	int j, k=1;

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	for(j=0;j<*n;j++)
	{
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &nu1, &out1[j]);
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &nu2, &out2[j]);
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &param[1], &out3[j]);

		out[j] = (out2[j]-2.0*out3[j]+out1[j])/(ep*ep);
	}

	Free(out1);Free(out2);Free(out3);
}




// 2. Ableitung von c nach u

void diff2PDF_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, *out3, u1, u2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	out3 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		u1=u[j]-ep;
		u2=u[j]+ep;
		copLik_mod(copula, &k, &u1, &v[j], &param[0], &param[1], &out1[j]);
		copLik_mod(copula, &k, &u2, &v[j], &param[0], &param[1], &out2[j]);
		copLik_mod(copula, &k, &u[j], &v[j], &param[0], &param[1], &out3[j]);
		out[j] = (out2[j]-2.0*out3[j]+out1[j])/(ep*ep);
	}

	Free(out1);Free(out2);Free(out3);
}


// 2. Ableitung von c nach u und rho

void diff2PDF_rho_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, u1, u2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;


	for(j=0;j<*n;j++)
	{
		u1=u[j]-ep;
		u2=u[j]+ep;
		diffPDF_rho_tCopula(&u1, &v[j], &k, param, copula, &out1[j]);
		diffPDF_rho_tCopula(&u2, &v[j], &k, param, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}

// 2. Ableitung von c nach u und v

void diff2PDF_u_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, u1, u2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		u1=v[j]-ep;
		u2=v[j]+ep;
		diffPDF_u_tCopula(&u[j], &u1, &k, param, copula, &out1[j]);
		diffPDF_u_tCopula(&u[j], &u2, &k, param, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


// 2. Ableitung von c nach u und nu

void diff2PDF_nu_u_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, u1, u2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		u1=u[j]-ep;
		u2=u[j]+ep;
		diffPDF_nu_tCopula(&u1, &v[j], &k, param, copula, &out1[j]);
		diffPDF_nu_tCopula(&u2, &v[j], &k, param, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


// 2. Ableitung von c nach rho und nu mit finiten Differenzen

void diff2PDF_rho_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, param1[2], param2[2];
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;
	param1[0]=param[0];
	param2[0]=param[0];
	param1[1]=param[1]-ep;
	param2[1]=param[1]+ep;

	for(j=0;j<*n;j++)
	{
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param1, copula, &out1[j]);
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param2, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}



// Ableitung von h nach rho

void diffhfunc_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;

	double rho = param[0];
	double nu = param[1];

	for(j=0;j<*n;j++)
	{
		t1 = qt(u[j],nu,1,0);
		t2 = qt(v[j],nu,1,0);
		t3 = t1-rho*t2;
		t4 = nu+t2*t2;
		t5 = 1.0-rho*rho;
		t6 = nu+1.0;
		t7 = t4*t5/t6;
		t8 = sqrt(t7);
		t9 = dt(t3/t8,nu+1.0,0);
		t10 = -t2/t8;
		t11 = t3*t4*rho/t6/t7/t8;
		out[j] = t9*(t10+t11);
	}
}

// Ableitung von h nach nu mit finiten Differenzen

void diffhfunc_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, nu1, nu2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;
	nu1=param[1]-ep;
	nu2=param[1]+ep;

	for(j=0;j<*n;j++)
	{
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &nu1, &out1[j]);
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &nu2, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


// Ableitung von h nach v mit finiten Differenzen

void diffhfunc_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, v1, v2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		v1=v[j]-ep;
		v2=v[j]+ep;
		Hfunc1(copula, &k, &u[j], &v1, &param[0], &param[1], &out1[j]);
		Hfunc1(copula, &k, &u[j], &v2, &param[0], &param[1], &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


// 2. Ableitung von h nach rho

void diff2hfunc_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, rho1, rho2, param1[2], param2[2];
	int j, k=1;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	rho1=param[0]-ep;
	rho2=param[0]+ep;
	param1[0]=rho1;
	param1[1]=param[1];
	param2[0]=rho2;
	param2[1]=param[1];

	for(j=0;j<*n;j++)
	{
		diffhfunc_rho_tCopula(&u[j], &v[j], &k, param1, copula, &out1[j]);
		diffhfunc_rho_tCopula(&u[j], &v[j], &k, param2, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}
	Free(out1);Free(out2);
}


// 2. Ableitung von h nach nu

void diff2hfunc_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, *out3, nu1, nu2;
	int j, k=1;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	out3 = (double*) Calloc(*n,double);
	nu1=param[1]-ep;
	nu2=param[1]+ep;

	for(j=0;j<*n;j++)
	{
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &nu1, &out1[j]);
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &nu2, &out2[j]);
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &param[1], &out3[j]);
		out[j] = (out2[j]-2.0*out3[j]+out1[j])/(ep*ep);
	}
	Free(out1);Free(out2);Free(out3);
}


// 2. Ableitung von h nach v mit finiten Differenzen

void diff2hfunc_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, *out3, v1, v2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	out3 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		v1=v[j]-ep;
		v2=v[j]+ep;
		Hfunc1(copula, &k, &u[j], &v1, &param[0], &param[1], &out1[j]);
		Hfunc1(copula, &k, &u[j], &v2, &param[0], &param[1], &out2[j]);
		Hfunc1(copula, &k, &u[j], &v[j], &param[0], &param[1], &out3[j]);
		out[j] = (out2[j]-2.0*out3[j]+out1[j])/(ep*ep);
	}

	Free(out1);Free(out2);Free(out3);
}


// 2. Ableitung von h nach rho und v mit finiten Differenzen

void diff2hfunc_rho_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, v1, v2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		v1=v[j]-ep;
		v2=v[j]+ep;
		diffhfunc_rho_tCopula(&u[j], &v1, &k, param, copula, &out1[j]);
		diffhfunc_rho_tCopula(&u[j], &v2, &k, param, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}


// 2. Ableitung von h nach nu und v mit finiten Differenzen

void diff2hfunc_nu_v_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, v1, v2;
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;

	for(j=0;j<*n;j++)
	{
		v1=v[j]-ep;
		v2=v[j]+ep;
		diffhfunc_nu_tCopula(&u[j], &v1, &k, param, copula, &out1[j]);
		diffhfunc_nu_tCopula(&u[j], &v2, &k, param, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}

// 2. Ableitung von h nach rho und nu mit finiten Differenzen

void diff2hfunc_rho_nu_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double *out1, *out2, param1[2], param2[2];
	out1 = (double*) Calloc(*n,double);
	out2 = (double*) Calloc(*n,double);
	int j, k=1;
	param1[0]=param[0];
	param2[0]=param[0];
	param1[1]=param[1]-ep;
	param2[1]=param[1]+ep;

	for(j=0;j<*n;j++)
	{
		diffhfunc_rho_tCopula(&u[j], &v[j], &k, param1, copula, &out1[j]);
		diffhfunc_rho_tCopula(&u[j], &v[j], &k, param2, copula, &out2[j]);
		out[j] = (out2[j]-out1[j])/(2.0*ep);
	}

	Free(out1);Free(out2);
}
