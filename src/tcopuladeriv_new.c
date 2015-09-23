// Neue Ableitungen der t-Copula
// Diesmal analytisch und nicht mit finiten differenzen

/*
** tcopuladeriv.c - C code of the package CDRVine  
** 
** by Ulf Schepsmeier
** 
** derivatives of the t-copula (density and h-function)
**
*/

#include "include/vine.h"
#include "include/likelihood.h"
#include "include/hfunc.h"
#include "include/deriv.h"
#include "include/tcopuladeriv.h"
#include "include/tcopuladeriv_new.h"
#include "include/incompleteBeta.h"


void diffX_nu_tCopula(double* x, double* param, double* out)
{
	//double xmin=0, val=0, err=0, val2;
	double xmax=0;
	double *inbeder_out;
	inbeder_out=Calloc(3,double);
	double t1, t2, t3, t4, t5, t6, t7;
	double x_help;

	double nu = param[1];

	if(*x>=0)
	{
		x_help=*x;
	}
	else
	{
		x_help=-*x;
	}

	xmax=nu/(nu+x_help*x_help);


	t1=dt(x_help,nu,0);
	t2=nu/2.0;
	t3=0.5;
	inbeder(&xmax,&t2,&t3,inbeder_out);
	t4=(nu+1.0)/2.0;
	t5=pow(nu,nu/2.0-1.0)*x_help;
	t6=pow(1.0/(x_help*x_help+nu),t4);
	t7=beta(nu/2.0,0.5);
	
	out[0]=1.0/(2.0*t1)*( 0.5*inbeder_out[1]+(t5*t6)/t7 );
	
	if(*x<0)
	{
		out[0]=-out[0];
	}

Free(inbeder_out);
}

void integral1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double nu = *((double *) fdata);
	fval[0] = pow((*x),nu/2.0-1.0)*pow((1.0-(*x)),-0.5)*log((*x));
}

void integral2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double nu = *((double *) fdata);
	fval[0] = pow((*x),nu/2.0-1.0)*pow((1.0-(*x)),-0.5);
}

void integral3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double nu = *((double *) fdata);
	fval[0] = pow((*x),nu/2.0-1.0)*pow((1.0-(*x)),-0.5)*log((*x))*log((*x));
}

void bigI_nu(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	double nu = *((double *) fdata);
	double xmax;

	xmax=nu/(nu+(*x)*(*x));


	fval[0]=pbeta(xmax,nu/2,0.5,1,0);
}


void diff_int_nu(double* x, double* nu, double* out)
{
	double t1, t2, t3, t4, t5, t6, t7, t8;

	t1=(*x)*(*x)+(*nu);
	t2=*nu/t1;
	t3=0.5*(*nu)-1.0;
	t4=pow(t2,t3);
	t5=1.0-t2;
	t6=pow(t5,-0.5);
	t7=log(t2);
	t8=2.0*(*nu)*(*x)/t1/t1;

	out[0]=t4*t6*t7*t8;
}


void diff_t_nu(double* x, double* nu, double* out)
{
	double xmax=0;
	double *inbeder_out;
	inbeder_out=Calloc(3,double);
	double t2, t3, t4, t5, t6, t7;
	double x_help;

	if(*x>=0)
	{
		x_help=*x;
	}
	else
	{
		x_help=-*x;
	}

	xmax=(*nu)/((*nu)+x_help*x_help);

	t2=*nu/2.0;
	t3=0.5;
	inbeder(&xmax,&t2,&t3,inbeder_out);
	t4=(*nu+1.0)/2.0;
	t5=pow(*nu,*nu/2.0-1.0)*x_help;
	t6=pow(1.0/(x_help*x_help+*nu),t4);
	t7=beta(t2,0.5);
	
	out[0]=-0.5*( 0.5*inbeder_out[1]+(t5*t6)/t7 );
		
	if(*x<0)
	{
		out[0]=-out[0];
	}

	Free(inbeder_out);
}



void diff_dt_x(double* x, double* nu, double* out)
{
	double t2, t3, t4, t5, t6, t7;

	t2=((*nu)+1.0)/(*nu);
	t3=sqrt((*nu));
	t4=1.0/t3/beta((*nu)*0.5,0.5);
	t5=1.0+((*x)*(*x))/(*nu);
	t6=((*nu)+3.0)/2.0;
	t7=pow(t5,-t6);
	out[0]=-t4*t2*(*x)*t7;
}


void diff_dt_nu(double* x, double* nu, double* out)
{
	double t1, t2, t3, t4, t6, t10, t11, t13, t14, t15, t16;

	t1=((*nu)+1.0)/2.0;
	t2=digamma(t1);
	t3=beta((*nu)*0.5,0.5);
	t4=sqrt(*nu);
	t6=digamma(0.5*(*nu));
	t10=-0.5/t3/t4*(t6-t2+1.0/(*nu));
	t11=1.0+((*x)*(*x))/(*nu);
	t13=pow(t11,-t1);
	t14=1.0/t3/t4;
	t15=log(t11);
	t16=-t1*(*x)*(*x)/(*nu)/(*nu)/t11;

	out[0]=t10*t13 + t14*(t13*(-0.5*t15-t16));

}

void diff_dt_u(double* x, double* nu, double* out)
{
	*out = -((*x)*(*nu+1.0)/(*nu))/(1.0+((*x)*(*x))/(*nu));
}


void diff_t_nu_nu(double* x, double* nu, double* out)
{
	double xmax=0;
	double *inbeder_out;
	inbeder_out=Calloc(3,double);
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
	double x_help;

	if(*x>=0)
	{
		x_help=*x;
	}
	else
	{
		x_help=-*x;
	}

	xmax=(*nu)/((*nu)+x_help*x_help);

	t1=1.0/(x_help*x_help+*nu);
	t2=*nu/2.0;
	t3=0.5;
	inbeder(&xmax,&t2,&t3,inbeder_out);
	t4=(*nu+1.0)/2.0;
	t5=pow(*nu,*nu/2.0-1.0)*x_help;
	t6=pow(t1,t4);
	t7=beta(t2,0.5);
	t8=t5*t6;
	t9=*nu*t1;
	
	t11=digamma(0.5*(*nu));
	t12=digamma(0.5*(*nu)+0.5);
	t13=t11-t12;
	t14=1.0/t7;
	
	t10=-t1*t4 + (t2-1.0)/(*nu) + 0.5*log(t1) + 0.5*log(*nu);
	
	t1=inbeder_out[2];
	t2=t8*log(t9)/t7;
	t3=t13*t8/t7;
	t4=t8/t7*t10;
	
	out[0]= - 1.0/8.0*inbeder_out[2] + t8*t14*( -0.25 * log(t9) +0.5* t13 - 0.5*t10  );
	
	
	if(*x<0)
	{
		out[0]=-out[0];
	}

	Free(inbeder_out);
}


void diff2_x_nu(double* x, double* nu, double* out)
{
	double t1, t2, t3, t4, t5;
	double *param;
	param=Calloc(2,double);
	param[0]=0;
	param[1]=*nu;

	t1=dt(*x,*nu,0);
	diff_t_nu_nu(x,nu,&t2);
	diff_dt_nu(x,nu,&t3);
	diffX_nu_tCopula(x, param, &t4);
	diff_dt_x(x,nu,&t5);

	out[0]=(-t5*t4*t4-t2-2.0*t3*t4)/t1;
	
	Free(param);
}


/////////////////////////////////////////////////////////


// Ableitung von c nach nu

void diffPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, x1, x2;
	int j=0, k=1;

	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, M, c;

	double rho = param[0];
	double nu = param[1];

	

	t1=digamma((nu+1.0)/2.0);
	t2=digamma(nu/2.0);
	t14=rho*rho;
	t3=0.5*log(1.0-t14);
	t4=(nu-2.0)/(2.0*nu);
	t5=0.5*log(nu);
	t6=-t1+t2+t3-t4-t5;
	t10=(nu+2.0)/2.0;

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		t7=1.0+2.0*x1*out1;
		t8=1.0+2.0*x2*out2;
		t15=x2*x2;
		t16=x1*x1;
		t9=(nu+1.0)/2.0*( t7/(nu+x1*x1) + t8/(nu+t15) );
		M=nu*(1.0-t14) + t16 + t15 - 2.0*rho*x1*x2;
		t11=1.0 - t14 + 2.0*x1*out1 + 2.0*x2*out2 - 2.0*rho*(x1*out2+x2*out1);
		t12=0.5*log((nu+t16)*(nu+t15));
		t13=0.5*log(M);

		out[j]=c*(t6 + t9 + t12 - t10*t11/M - t13);	
	}

}

// vectorized version
void diffPDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffPDF_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


// Ableitung von c nach u

void diffPDF_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double c, out1, t1, t2, t3, t4, t6, t7, t8, x1, x2;

	int j=0, k=1;

	double rho = param[0];
	double nu = param[1];


	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diff_dt_u(&x1, &nu, &out1);
		t1=c/dt(x1,nu,0);
		t2=(nu+2.0)*(x1-rho*x2);
		t6=rho*rho;
		t7=x1*x1;
		t8=x2*x2;
		t3=nu*(1.0-t6)+t7+t8-2*rho*x1*x2;
		t4=t2/t3;
		out[j]=-t1*(t4+out1);
	}
}



// 2. Ableitung von c nach nu

void diff2PDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, out3=0, out4=0, x1, x2, diff_nu=0;
	int j=0, k=1;

	double t1, t2, t3, t4, t5, t6, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, M_nu, M, M_nu_nu, c;

	double rho = param[0];
	double nu = param[1];

	
	t1=(nu+1.0)/2.0;
	t2=nu/2.0;
	t23=nu*nu;
	t3=1.0/t23;
	t4=1.0/(2.0*nu);
	t5=0.5*trigamma(t1);
	t6=(1.0-rho*rho);
	t9=0.5*trigamma(t2);
	t10=-t5+t9-t3-t4;
	

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		M = ( nu*t6 + x1*x1 + x2*x2 - 2.0*rho*x1*x2 );

		t8=(x1*out2+out1*x2);
		M_nu=t6+2.0*x1*out1+2.0*x2*out2-2.0*rho*t8;

		t24=x1*x1;
		t25=x2*x2;

		t11=1.0+2.0*x1*out1;
		t12=nu+t24;
		t13=t11/t12;

		t14=1.0+2.0*x2*out2;
		t15=nu+t25;
		t16=t14/t15;

		diff2_x_nu(&x1,&nu,&out3);
		diff2_x_nu(&x2,&nu,&out4);

		t17=2.0*out1*out1 + 2.0*x1*out3;
		t18=t17/t12;

		t19=2.0*out2*out2 + 2.0*x2*out4;
		t20=t19/t15;

		t21=t13*t13;
		t22=t16*t16;

		M_nu_nu=2.0*out1*out1 + 2.0*x1*out3 + 2.0*out2*out2 + 2.0*x2*out4 - 4.0*rho*out1*out2 - 2.0*rho*(x2*out3 + x1*out4);
		
		diffPDF_nu_tCopula_new(&u[j], &v[j], &k, param, copula, &diff_nu);

		out[j]=c*( t10+0.5*(t13+t16) + t1*(t18-t21+t20-t22) + 0.5*t13 + 0.5*t16 - M_nu/M - (nu/2.0+1.0)*(M_nu_nu/M-M_nu*M_nu/M/M )) + diff_nu*diff_nu/c;
	}
}

// vectorized version
void diff2PDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2PDF_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


// 2. Ableitung von c nach rho und nu

void diff2PDF_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, x1, x2;
	int j=0, k=1;

	double t3, t4, t5, t6, t7, t8, t9, t10, t11, M_rho, M_nu, M, c;

	double rho = param[0];
	double nu = param[1];

	
	t4=1.0-rho*rho;
	t3=rho/t4;
	t5=nu+2.0;

	

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		t10=x1*x1;
		t11=x2*x2;
		M = ( nu*t4 + t10 + t11 - 2.0*rho*x1*x2 );
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param, copula, &t6);
		diffPDF_nu_tCopula_new(&u[j], &v[j], &k, param, copula, &t7);
		M_rho=-2.0*(nu*rho+x1*x2);
		t8=(x1*out2+out1*x2);
		M_nu=t4+2.0*x1*out1+2.0*x2*out2-2.0*rho*t8;
		t9=-t3+t5/M*(rho+t8+0.5*M_nu*M_rho/M)-0.5*M_rho/M;

		out[j]=c*t9+t6*t7/c;	
	}
}

// vectorized version
void diff2PDF_rho_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2PDF_rho_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

// 2. Ableitung von c nach rho und u

void diff2PDF_rho_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double x1, x2;
	int j=0, k=1;

	double t1, t2, t3, t4,  M, c=0, out1=0, out2=0;

	double rho = param[0];
	double nu = param[1];
	
	t1=nu+2.0;
	t4=1.0-rho*rho;

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		M = ( nu*t4 + x1*x1 + x2*x2 - 2.0*rho*x1*x2 );
		t2=nu*rho+x1*x2;
		t3=x1-rho*x2;
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param, copula, &out1);
		diffPDF_u_tCopula_new(&u[j], &v[j], &k, param, copula, &out2);
		
		out[j]=c*(t1/M/dt(x1,nu,0) * (x2 - 2.0*t2/M*t3)) + out1/c*out2;
	}
}

// 2. Ableitung von c nach nu und u

void diff2PDF_nu_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double x1, x2;
	int j=0, k=1;

	double t1, t2, t3, t4, t6, t7, t8, t9, t10, t11, t12, t13, M, c=0, out1=0, out2=0, diffPDF=0, diff_dt=0, diff_dt2=0, diff_dt3=0, M_nu=0;
	double t14, t15, t16, t17;

	double rho = param[0];
	double nu = param[1];
	
	t1=nu+2.0;
	t3=(nu+1.0)/nu;
	t14=rho*rho;
	t4=1.0-t14;
	t17=nu*nu;

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		t15=x1*x1;
		t16=x2*x2;
		M = ( nu*t4 + t15 + t16 - 2.0*rho*x1*x2 );
		t2=dt(x1,nu,0);
		diffPDF_nu_tCopula_new(&u[j], &v[j], &k, param, copula, &diffPDF);
		diff_dt_nu(&x1, &nu, &diff_dt);
		diff_dt_u(&x1, &nu, &diff_dt2);
		diff_dt_x(&x1, &nu, &diff_dt3);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		t8=(x1*out2+out1*x2);
		M_nu=t4+2.0*x1*out1+2.0*x2*out2-2.0*rho*t8;
		
		t6=1.0+t15/nu;
		t7=x1-rho*x2;
		t9=-diffPDF/t2 + c/t2/t2*(diff_dt+diff_dt3*out1);
		t10=t1*t7/M + diff_dt2;
		t11=c/t2;
		t12=t7/M - t1*t7/M/M*M_nu + t1*(out1-rho*out2)/M;
		t13=-out1*t3/t6 + x1/(t17+nu*t15) + x1*t3/t6/t6 * (2.0*x1*out1/nu - t15/t17);
		
		out[j]=t9*t10 - t11*(t12+t13);
	}
	
}

// vectorized version

void diff2PDF_nu_u_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2PDF_nu_u_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


// 2. Ableitung von c nach u

void diff2PDF_u_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double x1, x2;
	int j=0, k=1;

	double t1, t2, t4, t5, t6, t7, t8, t9, t11, t12, t13, M, c=0, diffPDF=0, diff_dt2=0;
	double t14, t15, t16;

	double rho = param[0];
	double nu = param[1];
	
	t1=nu+2.0;
	t14=rho*rho;
	t4=1.0-t14;

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		t15=x1*x1;
		t16=x2*x2;
		M = ( nu*t4 + t15 + t16 - 2.0*rho*x1*x2 );
		t2=dt(x1,nu,0);
		diffPDF_u_tCopula_new(&u[j], &v[j], &k, param, copula, &diffPDF);
		diff_dt_u(&x1, &nu, &diff_dt2);
		
		t7=x1-rho*x2;
		
		t5=-diffPDF/t2 + diff_dt2/t2/t2*c;
		t6=t1*t7/M + diff_dt2;
		
		t11=c/t2;
		t8=1.0/t2;
		t9=t1/M - 2.0*t1*t7*t7/M/M;
		t13=1.0+t15/nu;
		
		t12=t8*( -(nu+1.0)/(nu+t15) + 2.0*t15* (nu+1.0)/nu/nu / t13/t13);
		
		out[j]=t5*t6 - t11*(t8*t9 + t12);
	}
}


// 2. Ableitung von c nach u und v

void diff2PDF_u_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double x1, x2;
	int j=0, k=1;

	double t1, t2, t4, t5, t6, t7, t8, t9, t10, t11, t12, M, c=0, diff_dt1=0, diff_dt2=0;

	double rho = param[0];
	double nu = param[1];
	
	t1=nu+2.0;
	t4=1.0-rho*rho;

	for(j=0;j<*n;j++)
	{
		LL(copula, &k, &u[j], &v[j], &rho, &nu, &c);
		c=exp(c);
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		M = ( nu*t4 + x1*x1 + x2*x2 - 2.0*rho*x1*x2 );
		t2=dt(x1,nu,0);
		t5=dt(x2,nu,0);
		
		t6=x1-rho*x2;
		t7=x2-rho*x1;
		
		diff_dt_u(&x1, &nu, &diff_dt1);
		diff_dt_u(&x2, &nu, &diff_dt2);
		
		t8=c/t2/t5;
		t9=t1*t6/M + diff_dt1;
		t10=t1*t7/M + diff_dt2;
		t11=t1*rho/M;
		t12=2.0*t1*t6*t7/M/M;
		
		out[j]=t8*(t9*t10 + t11 + t12);
	}
}


// Ableitung nach par und v

void diff2PDF_rho_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
       diff2PDF_rho_u_tCopula_new(v, u, n, param, copula, out);
}


// Ableitung nach nu und v

void diff2PDF_nu_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
       diff2PDF_nu_u_tCopula_new(v, u, n, param, copula, out);
}

void diff2PDF_nu_v_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    diff2PDF_nu_u_tCopula_new_vec(v, u, n, par, par2, copula, out);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


// Ableitung von h nach nu

void diffhfunc_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
	double diff_t=0, out1=0, out2=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,t6,0);
		t11=nu+1.0;
		
		diff_t_nu(&t10,&t11,&diff_t);
		
		diffX_nu_tCopula(&t1, param, &out1);
		diffX_nu_tCopula(&t2, param, &out2);
		
		t12=out1-rho*out2;
		t13=t12/t8;
		t14=1.0+2.0*t2*out2;
		t15=t14/t6;
		t16=t4/t6/t6;
		
		out[j]=diff_t + t9*(t13 - 0.5*t10/t7*t5*(t15-t16) );
		
	}
}

// vectorized version

void diffhfunc_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffhfunc_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

// Ableitung von h nach v

void diffhfunc_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;

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
		
		t10=t3/t8;
		t9 = dt(t10,nu+1.0,0);
		t11=nu+1.0;

		t12=1.0/dt(t2,nu,0);
		t13=-rho/t8;
		t14=t10/t7;
		t15=t5/t11*t2;
		
		out[j]=t9*t12*(t13-t14*t15);
		
	}
}



// 2.Ableitung von h nach rho

void diff2hfunc_rho_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
	double diff_t=0;

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
		t14=t4/t6;
		
		t10=t3/t8;
		t9 = dt(t10,nu+1.0,0);
		t11=nu+1.0;
		
		diff_dt_x(&t10,&t11,&diff_t);

		t12 = -t2/t8;
		t13 = t10/t7*rho*t14;
		
		out[j]=diff_t*(t12+t13)*(t12+t13) + t9*((t2/t7/t8 - 1.5*t13/t7) * (-2.0*rho)*t14 + t10/t7*t14);
		
	}
}


// 2.Ableitung von h nach nu

void diff2hfunc_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21;
	double diff_t=0, diff_t2=0, diff_t3=0;
	double out1=0, out2=0, out3=0, out4=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,t6,0);
		t11=nu+1.0;
		
		diff_dt_x(&t10,&t11,&diff_t);
		diff_t_nu_nu(&t10,&t11,&diff_t2);
		diff_dt_nu(&t10,&t11,&diff_t3);
		
		diffX_nu_tCopula(&t1, param, &out1);
		diffX_nu_tCopula(&t2, param, &out2);
		diff2_x_nu(&t1, &nu, &out3);
		diff2_x_nu(&t2, &nu, &out4);
		
		t12=out1-rho*out2;
		t13=t12/t8;
		t14=1.0+2.0*t2*out2;
		t15=t14/t6;
		t16=t4/t6/t6;
		t17=t15-t16;
		
		t18=(t13-0.5*t10/t7*t5*t17);
		t19=-t13/t7 + 0.75*t10/t7/t7*t5*t17;
		t20=(2.0*out2*out2 + 2.0*t2*out4)/t6 - 2.0*t15/t6 + 2.0*t16/t6;
		
		t21=(out3-rho*out4)/t8;
		
		out[j]=diff_t2+2.0*diff_t3*t18 + diff_t*t18*t18 + t9*(t19*t5*t17 + t21 - 0.5*t10/t7*t5*t20);
		
	}
}

// vectorized version
void diff2hfunc_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2hfunc_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


// 2.Ableitung von h nach v

void diff2hfunc_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
	double diff_t=0, diff_t2=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,t6,0);
		t11=nu+1.0;
		
		diff_dt_x(&t10,&t11,&diff_t);
		diff_dt_x(&t2,&nu,&diff_t2);
		t12=dt(t2,nu,0);
		
		t13=-rho/t8 - t10/t7*t5/t6*t2;
		
		t14=0.5*rho/t8/t7 + 1.5*t10/t7/t7*t5/t6*t2;
		t15=t5/t6*2.0*t2/t12;
		t16=t5/t6/t12*(t3-rho*t2)/t7/t8;
		
		out[j]=(diff_t/t12/t12*t13 - t9*diff_t2/t12/t12/t12) * t13 + t9/t12*(t14*t15 - t16);
	}
}



// 2.Ableitung von h nach rho und nu

void diff2hfunc_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23;
	double diff_t=0, diff_t3=0;
	double out1=0, out2=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,nu+1.0,0);
		t11=nu+1.0;
		
		diff_dt_x(&t10,&t11,&diff_t);
		diff_dt_nu(&t10,&t11,&diff_t3);
		
		diffX_nu_tCopula(&t1, param, &out1);
		diffX_nu_tCopula(&t2, param, &out2);
		
		t12=out1-rho*out2;
		t13=t12/t8;
		t14=1.0+2.0*t2*out2;
		t15=t14/t6;
		t16=t4/t6/t6;
		t17=t15-t16;
		
		t18=t13 - 0.5*t10/t7*t5*t17;
		t19=-t2/t8 + t10/t7*rho*t4/t6;
		t20=0.5*t2/t7/t8 - 1.5*t10/t7/t7*rho*t4/t6;
		t21=out2/t8;
		t22=t13/t7*rho*t4/t6;
		t23=t10/t7*rho*(t6*t14-t4)/t6/t6;
		
		out[j]=(diff_t3 + diff_t*t18) * t19 + t9*(t20*t5*t17 - t21 + t22 + t23);
	}
}

void diff2hfunc_rho_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2hfunc_rho_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

// 2.Ableitung von h nach rho und v

void diff2hfunc_rho_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
	double diff_t=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,t6,0);
		t11=nu+1.0;
		
		diff_dt_x(&t10,&t11,&diff_t);
		
		t12=dt(t2,nu,0);
		t13=-t2/t8 + t10/t7*rho*t4/t6;
		t14=-rho/t8 - t10/t7*t2*t5/t6;
		t15=-1.0/t8;
		t16=t2*t2/t8/t7*t5/t6;
		t17=t10/t7*2.0*rho*t2/t6;
		t18=1.5*t10/t7/t7*t5/t6*t2 + 0.5*rho/t8/t7;
		t19=-2.0*rho*t4/t6;
		
		out[j]=diff_t/t12*t13*t14 + t9/t12*(t15+t16+t17+t18*t19);
	}
}


// 2. Ableitung von h nach nu und v

void diff2hfunc_nu_v_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22;
	double diff_t=0, diff_t2=0, diff_t3=0, diff_t4=0;
	double out1=0, out2=0;

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
		
		t10=t3/t8;
		t9 = dt(t10,nu+1.0,0);
		t11=nu+1.0;
		
		diff_dt_nu(&t2,&nu,&diff_t);
		diff_dt_x(&t10,&t11,&diff_t2);
		diff_dt_nu(&t10,&t11,&diff_t3);
		diff_dt_x(&t2,&nu,&diff_t4);
		
		diffX_nu_tCopula(&t1, param, &out1);
		diffX_nu_tCopula(&t2, param, &out2);
		
		t12=dt(t2,nu,0);
		
		t13=out1-rho*out2;
		t14=t13/t8;
		
		t16=(1.0+2.0*t2*out2)/t6 - t4/t6/t6;
		t15=t14 - 0.5*t10/t7*t5*t16;
		
		t17=-rho/t8 - t10/t7*t5/t6*t2;
		t18=0.5*rho/t8/t7*t5*t16;
		t19=t14/t7*t5/t6*t2;
		t20=t10/t7*t5/t6*out2;
		t21=t10/t7*t5/t6/t6*t2;
		t22=1.5*t10/t7/t7*t5*t5/t6*t16*t2;
		
		out[j]=(diff_t3/t12 + diff_t2/t12*t15 - t9*(diff_t/t12/t12 + diff_t4/t12/t12*out2)) * t17 + t9/t12*(t18-t19-t20+t21+t22);
		
	}
}

// vectorized version
void diff2hfunc_nu_v_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2hfunc_nu_v_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

