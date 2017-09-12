#include "VineCopula/vine.h"
#include "VineCopula/likelihood.h"
#include "VineCopula/hfunc.h"
#include "VineCopula/deriv.h"
#include "VineCopula/tcopuladeriv.h"
#include "VineCopula/tcopuladeriv_new.h"
#include "VineCopula/incompleteBeta.h"
#include "VineCopula/logderiv.h"

#define UMAX  1-1e-10
#define UMIN  1e-10


//////////////////////////////////
// we calculated the derivatives of the copula density in deriv.c and deriv2.c
// Further, the derivatives of the Student's t-copula were derived in separate files due to their complexity
// here some derivatives of log(c) are calculated, since sometime it is numerical advantageous to use the log(c) instead of c,
// in particular for the t-copula.
//
// In most cases the calculation is almost the same as for the derivatives of c
//////////////////////////////////



// Ableitung von log(c) nach rho

/////////////////////////////////////////
// Derivative of log(c) wrt to the first parameter rho of the Student's t-copula
//
// Input:
// u, v			copula arguments (vectors)
// n			length of u and v
// param		two-dimensional parameter vector
// copula		copula family (not needed here)
//
// Output:
// out			derivative
///////////////////////////////////////////

void difflPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t10, t11;

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
		t1 = qt(u[j],nu,1,0);
		t2 = qt(v[j],nu,1,0);
		t3 = -(nu+2.0)/2.0;
		t10 = nu*(1.0-rho*rho);
		t4 = -2.0*t1*t2/t10;
		t11 = (t1*t1+t2*t2-2.0*rho*t1*t2);
		t5 = 2.0*t11*rho/t10/(1.0-rho*rho);
		t6 = 1.0+(t11/t10);
		t7 = rho/(1.0-rho*rho);
		out[j]=(t3*(t4+t5)/t6 + t7);
	}
}

// vectorized version
void difflPDF_rho_tCopula_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        difflPDF_rho_tCopula(&u[i], &v[i], &nn, ipars, copula, &out[i]);
    };
}


// Ableitung von log(c) nach nu

void difflPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, x1, x2;
	int j=0;

	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, M;

	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	t1=digamma((nu+1.0)/2.0);
	t2=digamma(nu/2.0);
	t3=0.5*log(1.0-rho*rho);
	t4=(nu-2.0)/(2.0*nu);
	t5=0.5*log(nu);
	t6=-t1+t2+t3-t4-t5;
	t10=(nu+2.0)/2.0;

	for(j=0;j<*n;j++)
	{
		x1=qt(u[j],nu,1,0);

		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		t7=1.0+2.0*x1*out1;
		t8=1.0+2.0*x2*out2;
		t9=(nu+1.0)/2.0*( t7/(nu+x1*x1) + t8/(nu+x2*x2) );
		M=nu*(1.0-rho*rho) + x1*x1 + x2*x2 - 2.0*rho*x1*x2;
		t11=1.0 - rho*rho + 2.0*x1*out1 + 2.0*x2*out2 - 2.0*rho*(x1*out2+x2*out1);
		t12=0.5*log((nu+x1*x1)*(nu+x2*x2));
		t13=0.5*log(M);

		out[j]=(t6 + t9 + t12 - t10*t11/M - t13);
	}

}

// vectorized version
void difflPDF_nu_tCopula_new_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        difflPDF_nu_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}



//2. Ableitung von log(c) nach rho

void diff2lPDF_rho_tCopula(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double M, t1, t2, t3, t4, t5, t6;
	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	t4 = 1.0-rho*rho;
	t3 = -(nu+1.0)*(1.0+rho*rho)/t4/t4;

	for(j=0;j<*n;j++)
	{
		t1 = qt(u[j],nu,1,0);
		t2 = qt(v[j],nu,1,0);

		M = ( nu*t4 + t1*t1 + t2*t2 - 2.0*rho*t1*t2 );

		t5 = (nu+2.0)*nu/M;
		t6 = 2.0*(nu+2.0)*pow(nu*rho+t1*t2,2.0)/M/M;

		out[j] = (t3 + t5 + t6);
	}
}



// 2. Ableitung von log(c) nach nu

void diff2lPDF_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, out3=0, out4=0, x1, x2, diff_nu=0;
	int j=0, k=1;

	double t1, t2, t3, t4, t5, t6, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, M_nu, M, M_nu_nu;

	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}


	t1=(nu+1.0)/2.0;
	t2=nu/2.0;
	t3=1.0/(nu*nu);
	t4=1.0/(2.0*nu);
	t5=0.5*trigamma(t1);
	t6=(1.0-rho*rho);
	t9=0.5*trigamma(t2);
	t10=-t5+t9-t3-t4;


	for(j=0;j<*n;j++)
	{
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		M = ( nu*t6 + x1*x1 + x2*x2 - 2.0*rho*x1*x2 );
		t8=(x1*out2+out1*x2);
		M_nu=t6+2.0*x1*out1+2.0*x2*out2-2.0*rho*t8;

		t11=1.0+2.0*x1*out1;
		t12=nu+x1*x1;
		t13=t11/t12;

		t14=1.0+2.0*x2*out2;
		t15=nu+x2*x2;
		t16=t14/t15;

		diff2_x_nu(&x1,&nu,&out3);
		diff2_x_nu(&x2,&nu,&out4);

		t17=2.0*out1*out1 + 2.0*x1*out3;
		t18=t17/t12;

		t19=2.0*out2*out2 + 2.0*x2*out4;
		t20=t19/t15;

		t21=t11*t11/t12/t12;
		t22=t14*t14/t15/t15;

		M_nu_nu=2.0*out1*out1 + 2.0*x1*out3 + 2.0*out2*out2 + 2.0*x2*out4 - 4.0*rho*out1*out2 - 2.0*rho*(x2*out3 + x1*out4);

		diffPDF_nu_tCopula_new(&u[j], &v[j], &k, param, copula, &diff_nu);

		out[j]=( t10+0.5*(t13+t16) + t1*(t18-t21+t20-t22) + 0.5*t13 + 0.5*t16 - M_nu/M - (nu/2.0+1.0)*(M_nu_nu/M-M_nu*M_nu/M/M ));
	}
}



// 2. Ableitung von log(c) nach rho und nu

void diff2lPDF_rho_nu_tCopula_new(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	double out1=0, out2=0, x1, x2;
	int j=0, k=1;

	double t3, t4, t5, t6, t7, t8, t9, M_rho, M_nu, M;

	double rho = param[0];
	double nu = param[1];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}


	t4=1.0-rho*rho;
	t3=rho/t4;
	t5=nu+2.0;



	for(j=0;j<*n;j++)
	{
		x1=qt(u[j],nu,1,0);
		x2=qt(v[j],nu,1,0);
		diffX_nu_tCopula(&x1, param, &out1);
		diffX_nu_tCopula(&x2, param, &out2);
		M = ( nu*t4 + x1*x1 + x2*x2 - 2.0*rho*x1*x2 );
		diffPDF_rho_tCopula(&u[j], &v[j], &k, param, copula, &t6);
		diffPDF_nu_tCopula_new(&u[j], &v[j], &k, param, copula, &t7);
		M_rho=-2.0*(nu*rho+x1*x2);
		t8=(x1*out2+out1*x2);
		M_nu=t4+2.0*x1*out1+2.0*x2*out2-2.0*rho*t8;
		t9=-t3+t5/M*(rho+t8+0.5*M_nu*M_rho/M)-0.5*M_rho/M;

		out[j]=t9;
	}

}


/////////////////////////////////////////////////////////////
//
// Ableitung der loglik nach dem Parameter
//
/////////////////////////////////////////////////////////////

void difflPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
{
  double* negv = (double *) malloc(*n*sizeof(double));
  double* negu = (double *) malloc(*n*sizeof(double));
  double* nparam = (double *) malloc(2*sizeof(double));

  int ncopula;
  nparam[0]= -param[0];
  nparam[1]= -param[1];
  int i;


  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  difflPDF(u, negv, n, nparam, &ncopula, out);
	  for(i=0;i<*n;i++){out[i]=-out[i];}
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  difflPDF(negu, v, n, nparam, &ncopula, out);
	  for(i=0;i<*n;i++){out[i]=-out[i];}
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i)
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		difflPDF(negu, negv, n, param, &ncopula, out);
	}
  else
	{
		difflPDF(u, v, n, param, copula, out);
	}

  free(negv);
  free(negu);
  free(nparam);
}

//vectorized_version
void difflPDF_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    for (int i = 0; i < (*n); ++i) {
        if (copula[i] == 2) {
            ipars[0] = par[i];
            ipars[1] = par2[i];
            difflPDF_rho_tCopula(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            difflPDF_mod(&u[i], &v[i], &nn, &par[i], &copula[i], &out[i]);
        }
    };
    free(ipars);
}


void difflPDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t14, t15, t16, t17, t18, t19, t20, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t45, t46;
	t3=0;
	t4=0;

	double theta = param[0];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	for(j=0;j<*n;j++)
	{
		if(*copula==0)
		{
			out[j]=0;
		}
		else if(*copula==1)
		{
			t1 = qnorm(u[j],0.0,1.0,1,0);
			t2 = qnorm(v[j],0.0,1.0,1,0);
			t3 = theta*theta;
			t4 = 1.0-t3;
			t5 = t1*t1;
			t6 = t2*t2;
			t7 = t4*t4;
			out[j] = (theta*t4-theta*(t5+t6)+(1.0+t3)*t1*t2)/t7;
		}
		else if(*copula==3)
		{
			t4 = log(u[j]*v[j]);
			t5 = pow(u[j],-1.0*theta);
			t6 = pow(v[j],-1.0*theta);
			t7 = t5+t6-1.0;
			t8 = log(t7);
			t9 = theta*theta;
			t14 = log(u[j]);
			t16 = log(v[j]);
			out[j] = 1/(1.0+theta)-t4+t8/t9+(1/theta+2.0)*(t5*t14+t6*t16)/t7;
		}
		else if(*copula==4)
		{
			t1 = log(u[j]);
			t2 = pow(-t1,1.0*theta);
			t3 = log(v[j]);
			t4 = pow(-t3,1.0*theta);
			t5 = t2+t4;
			t6 = 1/theta;
			t7 = pow(t5,1.0*t6);
			t8 = theta*theta;
			t10 = log(t5);
			t11 = 1/t8*t10;
			t12 = log(-t1);
			t14 = log(-t3);
			t16 = t2*t12+t4*t14;
			t18 = 1/t5;
			t20 = -t11+t6*t16*t18;
			t22 = exp(-t7);
			t23 = -1.0+t6;
			t24 = pow(t5,2.0*t23);
			t25 = t22*t24;
			t27 = t1*t3;
			t28 = theta-1.0;
			t29 = pow(t27,1.0*t28);
			t30 = pow(t5,-1.0*t6);
			t31 = t28*t30;
			t32 = 1.0+t31;
			t34 = 1/u[j];
			t35 = 1/v[j];
			t36 = t34*t35;
			t37 = t29*t32*t36;
			t45 = t25*t29;
			t46 = log(t27);
			out[j] = (-t7*t20*t25*t37+t25*(-2.0*t11+2.0*t23*t16*t18)*t37+t45*t46*t32*t36
					+t45*(t30-t31*t20)*t34*t35)/t22/t24/t29/t32*u[j]*v[j];

		}
		else if(*copula==5)
		{
			t1 = exp(-theta);
			t2 = 1.0-t1;
			t3 = u[j]+v[j];
			t5 = exp(-theta*t3);
			t8 = exp(-theta*u[j]);
			t9 = 1.0-t8;
			t11 = exp(-theta*v[j]);
			t12 = 1.0-t11;
			t14 = 1.0-t1-t9*t12;
			t15 = t14*t14;
			t16 = 1.0/t15;
			t17 = theta*t2;
			out[j] = (t2*t5*t16+theta*t1*t5*t16-t17*t3*t5*t16-2.0*t17*t5/t15/t14*(t1-u[j]*
						t8*t12-t9*v[j]*t11))/theta/t2/t5*t15;

		}
		else if(*copula==6)
		{
			t1 = 1.0-u[j];
			  t2 = pow(t1,1.0*theta);
			  t3 = 1.0-v[j];
			  t4 = pow(t3,1.0*theta);
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/theta-2.0;
			  t9 = pow(t6,1.0*t8);
			  t10 = theta*theta;
			  t12 = log(t6);
			  t14 = log(t1);
			  t15 = t2*t14;
			  t16 = log(t3);
			  t17 = t4*t16;
			  t18 = t15*t4;
			  t19 = t5*t16;
			  t26 = theta-1.0;
			  t27 = pow(t1,1.0*t26);
			  t28 = pow(t3,1.0*t26);
			  t30 = theta-1.0+t2+t4-t5;
			  t33 = t9*t27;
			  out[j] = (t9*(-1/t10*t12+t8*(t15+t17-t18-t19)/t6)*t27*t28*t30+t33*t14*t28*
					t30+t33*t28*t16*t30+t33*t28*(1.0+t15+t17-t18-t19))/t9/t27/t28/t30;

		}
	}

}



////////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach dem Parameter
//
////////////////////////////////////////////////////////////////////

void diff2lPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
{
  double* negv;
  double* negu;
  double* nparam;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  nparam = (double *) malloc(2*sizeof(double));
  int ncopula, i;
  nparam[0]=-param[0];
  nparam[1]=-param[1];



  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2lPDF(u, negv, n, nparam, &ncopula, out);
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2lPDF(negu, v, n, nparam, &ncopula, out);
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i)
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2lPDF(negu, negv, n, param, &ncopula, out);
	}
  else
	{
		diff2lPDF(u, v, n, param, copula, out);
	}

  free(negv);
  free(negu);
  free(nparam);
}


void diff2lPDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
	double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67, t68;
	double t70, t73, t74, t75, t76, t78, t80, t84, t85, t87, t88, t89, t90, t92, t93, t95, t96, t98;
	double t108, t109, t112, t116, t117, t132, t133, t136, t138, t139, t153, t154, t165;

	double theta = param[0];

	for(int i=0;i<*n;i++)
	{
	    if(u[i]<UMIN) u[i]=UMIN;
	    else if(u[i]>UMAX) u[i]=UMAX;
	    if(v[i]<UMIN) v[i]=UMIN;
	    else if(v[i]>UMAX) v[i]=UMAX;
	}

	for(j=0;j<*n;j++)
	{
		if(*copula==0)
		{
			out[j]=0;
		}
		else if(*copula==1)
		{
			t5 = qnorm(u[j],0.0,1.0,1,0);
			t6 = qnorm(v[j],0.0,1.0,1,0);
			t1 = theta*theta;
			t3 = t5*t5;
			t4 = t6*t6;
			t9 = 1.0-t1;
			t10 = t9*t9;
			out[j] = (1.0-3.0*t1-t3-t4+2.0*theta*t5*t6)/t10+4.0*(theta*t9-theta*(t3+t4)+
					(1.0+t1)*t5*t6)/t10/t9*theta;
		}
		else if(*copula==3)
		{
			t1 = u[j]*v[j];
			t2 = -theta-1.0;
			t3 = pow(t1,1.0*t2);
			t4 = log(t1);
			t6 = pow(u[j],-1.0*theta);
			t7 = pow(v[j],-1.0*theta);
			t8 = t6+t7-1.0;
			t10 = -2.0-1/theta;
			t11 = pow(t8,1.0*t10);
			t14 = t3*t11;
			t15 = theta*theta;
			t16 = 1/t15;
			t17 = log(t8);
			t19 = log(u[j]);
			t21 = log(v[j]);
			t23 = -t6*t19-t7*t21;
			t25 = 1/t8;
			t27 = t16*t17+t10*t23*t25;
			t30 = -t2*t3;
			t31 = t4*t4;
			t34 = t4*t11;
			t38 = t27*t27;
			t48 = t19*t19;
			t50 = t21*t21;
			t55 = t23*t23;
			t57 = t8*t8;
			t64 = -1/t2;
			t68 = 1/t3/t11;
			t73 = t14-t30*t34+t30*t11*t27;
			t74 = t2*t2;
			t78 = t73*t64;
			out[j] = (-2.0*t3*t4*t11+2.0*t14*t27+t30*t31*t11-2.0*t30*t34*t27+t30*t11*t38
						+t30*t11*(-2.0/t15/theta*t17+2.0*t16*t23*t25+t10*(t6*t48+t7*t50)*t25-t10*t55/
						t57))*t64*t68-t73/t74*t68+t78*t68*t4-t78*t68*t27;
		}
		else if(*copula==4)
		{
			t3 = log(u[j]);
			  t4 = pow(-t3,1.0*theta);
			  t5 = log(v[j]);
			  t6 = pow(-t5,1.0*theta);
			  t7 = t4+t6;
			  t8 = 1/theta;
			  t9 = pow(t7,1.0*t8);
			  t10 = theta*theta;
			  t11 = 1/t10;
			  t12 = log(t7);
			  t13 = t11*t12;
			  t14 = log(-t3);
			  t16 = log(-t5);
			  t18 = t4*t14+t6*t16;
			  t20 = 1/t7;
			  t22 = -t13+t8*t18*t20;
			  t23 = t22*t22;
			  t25 = exp(-t9);
			  t27 = t25/u[j];
			  t29 = 1/v[j];
			  t30 = -1.0+t8;
			  t31 = pow(t7,2.0*t30);
			  t32 = t29*t31;
			  t33 = t3*t5;
			  t34 = theta-1.0;
			  t35 = pow(t33,1.0*t34);
			  t36 = pow(t7,-1.0*t8);
			  t37 = t34*t36;
			  t38 = 1.0+t37;
			  t39 = t35*t38;
			  t40 = t32*t39;
			  t44 = 1/t10/theta*t12;
			  t47 = t11*t18*t20;
			  t49 = t14*t14;
			  t51 = t16*t16;
			  t53 = t4*t49+t6*t51;
			  t56 = t18*t18;
			  t58 = t7*t7;
			  t59 = 1/t58;
			  t61 = 2.0*t44-2.0*t47+t8*t53*t20-t8*t56*t59;
			  t65 = t9*t9;
			  t70 = t9*t22*t27;
			  t74 = -2.0*t13+2.0*t30*t18*t20;
			  t75 = t74*t35;
			  t80 = log(t33);
			  t87 = t36-t37*t22;
			  t88 = t35*t87;
			  t92 = t27*t29;
			  t93 = t74*t74;
			  t108 = t80*t38;
			  t112 = t31*t74;
			  t116 = t31*t35;
			  t117 = t80*t80;
			  t132 = -t9*t23*t27*t40-t9*t61*t27*t40+t65*t23*t27*t40-2.0*t70*t32*t75*t38
						-2.0*t70*t32*t35*t80*t38-2.0*t70*t32*t88+t92*t31*t93*t39+t92*t31*(4.0*t44-4.0*
						t47+2.0*t30*t53*t20-2.0*t30*t56*t59)*t39+2.0*t27*t32*t75*t108+2.0*t92*t112*t88+
						t92*t116*t117*t38+2.0*t92*t116*t80*t87+t92*t116*(-2.0*t36*t22+t37*t23-t37*t61);
			  t133 = 1/t25;
			  t136 = 1/t31;
			  t138 = 1/t35;
			  t139 = 1/t38;
			  t153 = (-t70*t40+t92*t112*t39+t92*t116*t108+t92*t116*t87)*t133*u[j]*v[j];
			  t154 = t136*t138;
			  t165 = t38*t38;
			  out[j] = t132*t133*u[j]*v[j]*t136*t138*t139+t153*t154*t139*t9*t22-t153*t154*t139*
						t74-t153*t154*t139*t80-t153*t154/t165*t87;
		}
		else if(*copula==5)
		{
			t1 = exp(theta);
			  t2 = theta*v[j];
			  t3 = theta*u[j];
			  t5 = exp(t2+t3+theta);
			  t8 = exp(t2+t3);
			  t10 = exp(t2+theta);
			  t12 = exp(t3+theta);
			  t13 = t8-t10-t12+t1;
			  t14 = t13*t13;
			  t15 = 1/t14;
			  t18 = t1-1.0;
			  t19 = v[j]+u[j]+1.0;
			  t21 = t5*t15;
			  t24 = t18*t5;
			  t26 = 1/t14/t13;
			  t27 = v[j]+u[j];
			  t29 = v[j]+1.0;
			  t31 = u[j]+1.0;
			  t33 = t27*t8-t29*t10-t31*t12+t1;
			  t37 = theta*t1;
			  t38 = t37*t21;
			  t40 = t19*t5*t15;
			  t43 = t5*t26;
			  t44 = t43*t33;
			  t47 = theta*t18;
			  t48 = t19*t19;
			  t55 = t14*t14;
			  t58 = t33*t33;
			  t62 = t27*t27;
			  t64 = t29*t29;
			  t66 = t31*t31;
			  t73 = 1/theta;
			  t75 = 1/t18;
			  t76 = 1/t5;
			  t78 = t75*t76*t14;
			  t84 = t24*t15+t38+t47*t40-2.0*t47*t44;
			  t85 = theta*theta;
			  t89 = t84*t73;
			  t90 = t18*t18;
			  t93 = t76*t14;
			  t96 = t89*t75;
			  out[j] = (2.0*t1*t5*t15+2.0*t18*t19*t21-4.0*t24*t26*t33+t38+2.0*t37*t40-4.0
						*t37*t44+t47*t48*t5*t15-4.0*t47*t19*t44+6.0*t47*t5/t55*t58-2.0*t47*t43*(t62*t8-
						t64*t10-t66*t12+t1))*t73*t78-t84/t85*t78-t89/t90*t93*t1-t96*t93*t19+2.0*t96*t76
						*t13*t33;
		}
		else if(*copula==6)
		{
			t1 = 1.0-u[j];
			  t2 = pow(t1,1.0*theta);
			  t3 = 1.0-v[j];
			  t4 = pow(t3,1.0*theta);
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/theta-2.0;
			  t9 = pow(t6,1.0*t8);
			  t10 = theta*theta;
			  t11 = 1/t10;
			  t12 = log(t6);
			  t14 = log(t1);
			  t15 = t2*t14;
			  t16 = log(t3);
			  t17 = t4*t16;
			  t18 = t15*t4;
			  t19 = t5*t16;
			  t20 = t15+t17-t18-t19;
			  t22 = 1/t6;
			  t24 = -t11*t12+t8*t20*t22;
			  t25 = t24*t24;
			  t27 = theta-1.0;
			  t28 = pow(t1,1.0*t27);
			  t29 = pow(t3,1.0*t27);
			  t30 = t28*t29;
			  t31 = theta-1.0+t2+t4-t5;
			  t32 = t30*t31;
			  t41 = t14*t14;
			  t42 = t2*t41;
			  t43 = t16*t16;
			  t49 = t42+t4*t43-t42*t4-2.0*t15*t17-t5*t43;
			  t52 = t20*t20;
			  t54 = t6*t6;
			  t60 = t9*t24;
			  t61 = t60*t28;
			  t62 = t14*t29;
			  t63 = t62*t31;
			  t66 = t29*t16;
			  t67 = t66*t31;
			  t70 = 1.0+t15+t17-t18-t19;
			  t74 = t9*t28;
			  t92 = t9*t25*t32+t9*(2.0/t10/theta*t12-2.0*t11*t20*t22+t8*t49*t22-t8*t52/
						t54)*t32+2.0*t61*t63+2.0*t61*t67+2.0*t60*t30*t70+t74*t41*t29*t31+2.0*t74*t14*
						t67+2.0*t74*t62*t70+t74*t29*t43*t31+2.0*t74*t66*t70+t74*t29*t49;
			  t93 = 1/t9;
			  t95 = 1/t28;
			  t96 = 1/t29;
			  t98 = 1/t31;
			  t108 = (t60*t32+t74*t63+t74*t67+t74*t29*t70)*t93*t95;
			  t109 = t96*t98;
			  t116 = t31*t31;
			  out[j] = t92*t93*t95*t96*t98-t108*t109*t24-t108*t109*t14-t108*t109*t16-t108
						*t96/t116*t70;
		}
	}
}

