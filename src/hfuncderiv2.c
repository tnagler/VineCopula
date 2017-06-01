/*
** hfuncderiv2.c - C code of the package CDRVine  
** 
** by Ulf Schepsmeier
** 
** Second derivatives of the h-function with respect to different parameters
**
*/

#include "vine.h"
#include "deriv.h"
#include "deriv2.h"
#include "tcopuladeriv.h"
#include "tcopuladeriv_new.h"

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4


////////////////////////////////////////////////////////////////////
//
// 2. Ableitung der h-function nach dem Parameter
//
////////////////////////////////////////////////////////////////////

void diff2hfunc_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
{
  double* negv;
  double* negu;
  double* nparam;
  double* out2;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  nparam = (double *) malloc(2*sizeof(double));
  out2 = (double *) malloc(*n*sizeof(double));
  int ncopula, i;
  nparam[0]=-param[0];
  nparam[1]=-param[1];

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc(u, v, n, nparam, &ncopula, out);
		diffhfunc(u, v, n, nparam, &ncopula, out2);
		for (i = 0; i < *n; i++) {out[i]=out[i]*4/pow(1-param[0],4)+out2[i]*4/pow(1-param[0],3);}
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc(u, negv, n, nparam, &ncopula, out);
			diffhfunc(u, negv, n, nparam, &ncopula, out2);
			for (i = 0; i < *n; i++) {out[i]=out[i]*4/pow(1+param[0],4)+out2[i]*4/pow(1+param[0],3);}
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc(u, v, n, nparam, &ncopula, out);
		diffhfunc(u, v, n, nparam, &ncopula, out2);
		for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1-param[0],4)+out2[i]*2/pow(1-param[0],3);}
		}else{
			nparam[0]=1/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc(u, negv, n, nparam, &ncopula, out);
			diffhfunc(u, negv, n, nparam, &ncopula, out2);
			for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1+param[0],4)+out2[i]*2/pow(1+param[0],3);}
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc(u, negv, n, nparam, &ncopula, out);
	  //for (i = 0; i < *n; i++) {out[i]=-out[i];}
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc(negu, v, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=-out[i];}
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc(negu, negv, n, param, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=-out[i];};
	}
  else 
	{
		diff2hfunc(u, v, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
  free(out2);
}

// vectorized version
void diff2hfunc_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    for (int i = 0; i < (*n); ++i) {
        if (copula[i] == 2) {
            ipars[0] = par[i];
            ipars[1] = par2[i];
            diff2hfunc_rho_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diff2hfunc_mod(&u[i], &v[i], &nn, &par[i], &copula[i], &out[i]);
        }
    };
    free(ipars);
}


void diff2hfunc_mod2(double* v, double* u, int* n, double* param, int* copula, double* out)  // Achtung u und v vertauscht; Notaion aus Hfunc1 und Hfunc2
{
  double* negv;
  double* negu;
  double* nparam;
  double* out2;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  nparam = (double *) malloc(2*sizeof(double));
  out2 = (double *) malloc(*n*sizeof(double));
  int ncopula, i;
  nparam[0]=-param[0];
  nparam[1]=-param[1];

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc(v, u, n, nparam, &ncopula, out);
		diffhfunc(v, u, n, nparam, &ncopula, out2);
		for (i = 0; i < *n; i++) {out[i]=out[i]*4/pow(1-param[0],4)+out2[i]*4/pow(1-param[0],3);}
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc(negv, u, n, nparam, &ncopula, out);
	  diffhfunc(negv, u, n, nparam, &ncopula, out2);
	  for (i = 0; i < *n; i++) {out[i]=-out[i]*4/pow(1+param[0],4)-out2[i]*4/pow(1+param[0],3);};
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc(v, u, n, nparam, &ncopula, out);
		diffhfunc(v, u, n, nparam, &ncopula, out2);
		for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1-param[0],4)+out2[i]*2/pow(1-param[0],3);}
		}else{
			nparam[0]=1/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc(negv, u, n, nparam, &ncopula, out);
			diffhfunc(negv, u, n, nparam, &ncopula, out2);
			for (i = 0; i < *n; i++) {out[i]=-out[i]/pow(1+param[0],4)-out2[i]*2/pow(1+param[0],3);};
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc(negv, u, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=-out[i];};
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc(v, negu, n, nparam, &ncopula, out);
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc(negv, negu, n, param, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=-out[i];};
	}
  else 
	{
		diff2hfunc(v, u, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
  free(out2);
}



void diff2hfunc(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t18, t20, t21, t23, t25, t26, t28, t29, t30, t31, t32;
	double t36, t39, t40, t42, t44, t47, t49, t50, t56, t61, t62, t65, t70, t73;

	double theta = param[0];
	//double delta = param[1];

	for(j=0;j<*n;j++)
	{
		if(*copula==0)
		{
			out[j]=0;
		}
		else if(*copula==1)
		{
			t1=qnorm(u[j],0.0,1.0,1,0);
			t2=qnorm(v[j],0.0,1.0,1,0);
			t4=dnorm(t2,0.0,1.0,0);
			t5=1/t4;
			t6=t1-theta*t2;
			t7=1-pow(theta,2.0);
			t8=sqrt(t7);
			t9=t6/t8;
			t10=dnorm(t9,0.0,1.0,0);
			t11=t6*t2/t7;
			t12=pow(t6,2.0)*theta/pow(t7,2.0);
			t13=t10*(t11-t12);
			t14=-t2*t8+(t1-theta*t2)*theta/t8;
			t15=t14/t7;
			t16=2.0*t1*theta*theta - 3.0*theta*t2 + t1;
			t18=pow(t7,2.5);
			t20=t16/t18;
			out[j]=t13*t15 + t10*t20;
		}
		else if(*copula==3)
		{
			t2 = pow(v[j],-1.0*theta-1.0);
			t3 = log(v[j]);
			t4 = t3*t3;
			t6 = pow(u[j],-1.0*theta);
			t7 = pow(v[j],-1.0*theta);
			t8 = t6+t7-1.0;
			t10 = -1.0-1/theta;
			t11 = pow(t8,1.0*t10);
			t14 = theta*theta;
			t15 = 1/t14;
			t16 = log(t8);
			t18 = log(u[j]);
			t21 = -t6*t18-t7*t3;
			t23 = 1/t8;
			t25 = t15*t16+t10*t21*t23;
			t29 = t2*t11;
			t32 = t25*t25;
			t12 = t18*t18;
			t9 = t21*t21;
			t5 = t8*t8;
			out[j] = t2*t4*t11-2.0*t2*t3*t11*t25+t29*t32+t29*(-2.0/t14/theta*t16+2.0*t15
					*t21*t23+t10*(t6*t12+t7*t4)*t23-t10*t9/t5);
		}
		else if(*copula==4)
		{
			  t1 = log(v[j]);
			  t2 = pow(-t1,1.0*theta);
			  t3 = log(u[j]);
			  t4 = pow(-t3,1.0*theta);
			  t5 = t2+t4;
			  t6 = 1/theta;
			  t7 = pow(t5,1.0*t6);
			  t8 = theta*theta;
			  t9 = 1/t8;
			  t10 = log(t5);
			  t11 = t9*t10;
			  t12 = log(-t1);
			  t13 = t2*t12;
			  t14 = log(-t3);
			  t16 = t13+t4*t14;
			  t18 = 1/t5;
			  t20 = -t11+t6*t16*t18;
			  t21 = t20*t20;
			  t23 = exp(-t7);
			  t25 = t6-1.0;
			  t26 = pow(t5,1.0*t25);
			  t28 = 1/v[j];
			  t29 = 1/t1;
			  t30 = t28*t29;
			  t31 = t26*t2*t30;
			  t36 = 2.0/t8/theta*t10;
			  t39 = 2.0*t9*t16*t18;
			  t40 = t12*t12;
			  t42 = t14*t14;
			  t44 = t2*t40+t4*t42;
			  t47 = t16*t16;
			  t49 = t5*t5;
			  t50 = 1/t49;
			  t56 = t7*t7;
			  t61 = t23*t26;
			  t62 = t7*t20*t61;
			  t65 = -t11+t25*t16*t18;
			  t70 = t13*t30;
			  t73 = t65*t65;
			  t15 = t2*t28*t29;
			  out[j] = t7*t21*t23*t31+t7*(t36-t39+t6*t44*t18-t6*t47*t50)*t23*t31-t56*t21*
					t23*t31+2.0*t62*t65*t2*t30+2.0*t62*t70-t61*t73*t15-t61*(t36-t39+t25*t44*t18-t25
					*t47*t50)*t15-2.0*t61*t65*t70-t61*t2*t40*t28*t29;
		}
		else if(*copula==5)
		{
			t1 = exp(theta);
			t2 = theta*u[j];
			t3 = exp(t2);
			t5 = t1*(t3-1.0);
			t6 = theta*v[j];
			t8 = exp(t6+t2);
			t10 = exp(t6+theta);
			t12 = exp(t2+theta);
			t13 = t8-t10-t12+t1;
			t14 = 1/t13;
			t16 = t1*u[j];
			t18 = t3*t14;
			t20 = t13*t13;
			t21 = 1/t20;
			t23 = v[j]+u[j];
			t25 = v[j]+1.0;
			t26 = u[j]+1.0;
			t28 = t23*t8-t25*t10-t26*t12+t1;
			t32 = u[j]*u[j];
			t42 = t28*t28;
			t44 = t23*t23;
			t47 = t25*t25;
			t49 = t26*t26;
			out[j] = -t5*t14-2.0*t16*t18+2.0*t5*t21*t28-t1*t32*t18+2.0*t16*t3*t21*t28
					-2.0*t5/t20/t13*t42+t5*t21*(t44*t8-t47*t10-t49*t12+t1);
		}
		else if(*copula==6)
		{
			  t1 = 1.0-u[j];
			  t2 = pow(t1,1.0*theta);
			  t3 = 1.0-v[j];
			  t4 = pow(t3,1.0*theta);
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/theta-1.0;
			  t9 = pow(t6,1.0*t8);
			  t10 = theta*theta;
			  t11 = 1/t10;
			  t12 = log(t6);
			  t14 = log(t1);
			  t15 = t2*t14;
			  t16 = log(t3);
			  t18 = t4*t16;
			  t20 = t15+t18-t15*t4-t5*t16;
			  t21 = 1/t6;
			  t23 = -t11*t12+t8*t20*t21;
			  t25 = t23*t23;
			  t28 = pow(t3,1.0*theta-1.0);
			  t29 = 1.0-t2;
			  t30 = t28*t29;
			  t39 = t14*t14;
			  t40 = t2*t39;
			  t42 = t16*t16;
			  t50 = t20*t20;
			  t56 = t6*t6;
			  t13 = t9*t23;
			  t7 = t9*t28;
			  out[j] = t9*t25*t30+t9*(2.0/t10/theta*t12-2.0*t11*t20*t21+t8*(t40+t4*t42-t40
					*t4-2.0*t15*t18-t5*t42)*t21-t8*t50/t56)*t30+2.0*t13*t28*t16*t29-2.0*t13*t28*t2*
					t14+t7*t42*t29-2.0*t7*t16*t2*t14-t7*t40;
		}
	}
}

////////////////////////////////////////////////////////////////////
//
// 2. Ableitung der h-function nach v
//
////////////////////////////////////////////////////////////////////

void diff2hfunc_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc_v(u, v, n, nparam, &ncopula, out);
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_v(u, negv, n, nparam, &ncopula, out);
			//for(i=0;i<*n;i++){out[i]=-out[i];}
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc_v(u, v, n, nparam, &ncopula, out);
		}else{
			nparam[0]=1/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_v(u, negv, n, nparam, &ncopula, out);
			//for(i=0;i<*n;i++){out[i]=-out[i];}
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_v(u, negv, n, nparam, &ncopula, out);
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc_v(negu, v, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=-out[i];};
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc_v(negu, negv, n, param, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=1-out[i];};
	}
  else 
	{
		diff2hfunc_v(u, v, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
}

// vectorized version
void diff2hfunc_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2hfunc_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


void diff2hfunc_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out)
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

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc_v(v, u, n, nparam, &ncopula, out);
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_v(negv, u, n, nparam, &ncopula, out);
			//for(i=0;i<*n;i++){out[i]=-out[i];}
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc_v(v, u, n, nparam, &ncopula, out);
		}else{
			nparam[0]=1/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_v(negv, u, n, nparam, &ncopula, out);
			//for(i=0;i<*n;i++){out[i]=-out[i];}
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_v(negv, u, n, nparam, &ncopula, out);
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc_v(v, negu, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=-out[i];};
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc_v(negv, negu, n, param, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=1-out[i];};
	}
  else 
	{
		diff2hfunc_v(v, u, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
}


void diff2hfunc_v(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j, k=1;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t18, t20, t21, t22, t23, t24, t26, t27, t29, t30, t31, t32, t34, t36;
	double t37, t38, t40, t41, t48, t49, t50, t51, t55, t59, t62, t65, t69, t70, t79, t86;

	double theta = param[0];

	for(j=0;j<*n;j++)
	{
		if(*copula==0)
		{
			out[j]=0;
		}
		else if(*copula==1)
		{
			t1=qnorm(u[j],0.0,1.0,1,0);
			t2=qnorm(v[j],0.0,1.0,1,0);
			t4=dnorm(t2,0.0,1.0,0);
			t5=1-theta*theta;
			t6=sqrt(t5);
			t7=(t1-theta*t2)/t6;
			t8=dnorm(t7,0.0,1.0,0);
			t9=pow(t4,2.0);
			out[j]=-theta*t8/t6/t9*(t7/t6*theta + t2);
		}
		else if(*copula==2)
		{
			diff2hfunc_v_tCopula_new(&u[j], &v[j], &k, param, copula, &out[j]);
		}
		else if(*copula==3)
		{
			t1 = -theta-1.0;
			t2 = pow(v[j],1.0*t1);
			t3 = t1*t1;
			t5 = v[j]*v[j];
			t6 = 1/t5;
			t7 = pow(u[j],-1.0*theta);
			t8 = pow(v[j],-1.0*theta);
			t9 = t7+t8-1.0;
			t11 = -1.0-1/theta;
			t12 = pow(t9,1.0*t11);
			t13 = t6*t12;
			t16 = t2*t1*t13;
			t18 = 1/t9;
			t23 = t2*t12;
			t24 = t11*t11;
			t26 = t8*t8;
			t27 = theta*theta;
			t29 = t9*t9;
			t32 = t26*t27*t6/t29;
			t34 = t23*t11;
			t36 = t6*t18;
			out[j] = t2*t3*t13-t16-2.0*t16*t11*t8*theta*t18+t23*t24*t32+t34*t8*t27*t36+
					t34*t8*theta*t36-t34*t32;
		}
		else if(*copula==4)
		{
			  t1 = log(v[j]);
			  t2 = pow(-t1,1.0*theta);
			  t3 = log(u[j]);
			  t4 = pow(-t3,1.0*theta);
			  t5 = t2+t4;
			  t6 = 1/theta;
			  t7 = pow(t5,1.0*t6);
			  t8 = exp(-t7);
			  t9 = t6-1.0;
			  t10 = pow(t5,1.0*t9);
			  t11 = t8*t10;
			  t12 = v[j]*v[j];
			  t14 = 1/t12/v[j];
			  t15 = t2*t14;
			  t20 = t1*t1;
			  t21 = 1/t20;
			  t26 = 1/t20/t1;
			  t30 = t7*t7;
			  t31 = t2*t2;
			  t32 = t31*t2;
			  t34 = t5*t5;
			  t36 = 1/t34;
			  t37 = t26*t36;
			  t38 = t37*t11;
			  t40 = t11*t2;
			  t41 = theta*t14;
			  t48 = t7*t31;
			  t49 = t48*t14;
			  t50 = 1/t5;
			  t51 = t21*t50;
			  t55 = t26*t50;
			  t59 = t7*t32;
			  t62 = t14*t26;
			  t65 = t10*theta;
			  t69 = t59*t62;
			  t70 = t36*t8;
			  t79 = t11*t9*t31;
			  t86 = t9*t9;
			  t18 = theta*theta;
			  t16 = t18*t14;
			  t13 = t16*t37;
			  out[j] = -2.0*t11*t15/t1-3.0*t11*t15*t21-2.0*t11*t15*t26-t30*t32*t14*t38+
					3.0*t40*t41*t21+3.0*t40*t41*t26-3.0*t49*t51*t11-3.0*t49*t55*t11+t59*t14*t38+3.0
					*t48*t62*t50*t8*t65-t69*t70*t65+2.0*t69*t70*t10*t9*theta+3.0*t79*t41*t51+3.0*
					t79*t41*t55-t11*t86*t32*t13-3.0*t79*t16*t55+t11*t9*t32*t13-t40*t16*t26;
		}
		else if(*copula==5)
		{
			t1 = exp(theta);
			t2 = theta*u[j];
			t3 = exp(t2);
			t5 = t1*(t3-1.0);
			t6 = theta*v[j];
			t8 = exp(t6+t2);
			t10 = exp(t6+theta);
			t12 = exp(t2+theta);
			t13 = t8-t10-t12+t1;
			t14 = t13*t13;
			t20 = pow(theta*t8-theta*t10,2.0);
			t24 = theta*theta;
			out[j] = -2.0*t5/t14/t13*t20+t5/t14*(t24*t8-t24*t10);
		}
		else if(*copula==6)
		{
			t2 = pow(1.0-u[j],1.0*theta);
			t3 = 1.0-v[j];
			t4 = pow(t3,1.0*theta);
			t5 = t2*t4;
			t6 = t2+t4-t5;
			t8 = 1/theta-1.0;
			t9 = pow(t6,1.0*t8);
			t10 = t8*t8;
			t12 = t4*theta;
			t13 = 1/t3;
			t14 = -t12*t13+t5*theta*t13;
			t18 = t14*t14;
			t20 = t6*t6;
			t22 = theta-1.0;
			t23 = pow(t3,1.0*t22);
			t24 = 1.0-t2;
			t26 = 1/t20*t23*t24;
			t27 = t9*t8;
			t29 = theta*theta;
			t31 = t3*t3;
			t32 = 1/t31;
			t41 = 1/t6;
			t51 = t9*t23;
			t55 = t22*t22;
			out[j] = t9*t10*t18*t26+t27*(t4*t29*t32-t12*t32-t5*t29*t32+t5*theta*t32)*t41
					*t23*t24-t27*t18*t26-2.0*t27*t14*t41*t23*t22*t13*t24+t51*t55*t32*t24-t51*t22*
					t32*t24;
		}
	}
}


////////////////////////////////////////////////////////////////////
//
// Ableitung der h-function nach dem Parameter und v
//
////////////////////////////////////////////////////////////////////

void diff2hfunc_par_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc_par_v(u, v, n, nparam, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=out[i]*2/pow(1-param[0],2);}
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_par_v(u, negv, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=-out[i]*2/pow(1+param[0],2);};
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc_par_v(u, v, n, nparam, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1-param[0],2);}
		}else{
			nparam[0]=1/(1+param[0]);
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_par_v(u, negv, n, nparam, &ncopula, out);
	  for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1+param[0],2);};
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_par_v(u, negv, n, nparam, &ncopula, out);
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc_par_v(negu, v, n, nparam, &ncopula, out);
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc_par_v(negu, negv, n, param, &ncopula, out);
	}
  else 
	{
		diff2hfunc_par_v(u, v, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
}

// vectorized version
void diff2hfunc_par_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    
    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        if (copula[i] == 2) {
            diff2hfunc_rho_v_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diff2hfunc_par_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        }
    };
    free(ipars);
}


void diff2hfunc_par_v(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t28, t29, t30, t31, t32, t33, t34, t35, t36;
	double t37, t38, t39, t41, t42, t46, t47, t50, t51, t53, t57, t60, t64, t67, t69, t75, t77, t83, t84, t86, t91, t92;
	double t104, t106, t109, t117, t122;

	double theta = param[0];

	for(j=0;j<*n;j++)
	{
		if(*copula==0)
		{
			out[j]=0;
		}
		else if(*copula==1)
		{
			t1=qnorm(u[j],0.0,1.0,1,0);
			t2=qnorm(v[j],0.0,1.0,1,0);
			t4=dnorm(t2,0.0,1.0,0);
			t5=1/t4;
			t6=t1-theta*t2;
			t7=1-pow(theta,2.0);
			t8=sqrt(t7);
			t9=t6/t8;
			t10=dnorm(t9,0.0,1.0,0);
			t11=t6*t2/t7;
			t12=pow(t6,2.0)*theta/pow(t7,2.0);
			t13=-theta/t8;
			t14=-1.0/t8;
			t15=pow(theta,2)/pow(t7,1.5);
			out[j]=t5*(t10*(t11-t12)*t13+t10*(t14-t15));
		}
		else if(*copula==3)
		{
			  t1 = -theta-1.0;
			  t2 = pow(v[j],1.0*t1);
			  t3 = t2*t1;
			  t4 = 1/v[j];
			  t5 = log(v[j]);
			  t6 = t4*t5;
			  t7 = pow(u[j],-1.0*theta);
			  t8 = pow(v[j],-1.0*theta);
			  t9 = t7+t8-1.0;
			  t10 = 1/theta;
			  t11 = -1.0-t10;
			  t12 = pow(t9,1.0*t11);
			  t20 = t8*theta;
			  t21 = 1/t9;
			  t22 = t4*t21;
			  t26 = theta*theta;
			  t28 = log(t9);
			  t30 = log(u[j]);
			  t34 = t11*(-t7*t30-t8*t5);
			  t36 = 1/t26*t28+t34*t21;
			  t39 = t2*t12;
			  t53 = t9*t9;
			  out[j] = -t3*t6*t12-t2*t4*t12+t2*t5*t12*t11*t20*t22+t3*t4*t12*t36-t39*t11*t8
					*theta*t4*t21*t36+t39*(-t10*t8*t22+t11*(t20*t6-t8*t4)*t21+t34/t53*t20*t4);
		}
		else if(*copula==4)
		{
			  t1 = log(v[j]);
			  t2 = pow(-t1,1.0*theta);
			  t3 = log(u[j]);
			  t4 = pow(-t3,1.0*theta);
			  t5 = t2+t4;
			  t6 = 1/theta;
			  t7 = pow(t5,1.0*t6);
			  t8 = t2*t2;
			  t10 = v[j]*v[j];
			  t11 = 1/t10;
			  t12 = t1*t1;
			  t13 = 1/t12;
			  t14 = t11*t13;
			  t15 = t7*t8*t14;
			  t16 = 1/t5;
			  t17 = theta*theta;
			  t19 = log(t5);
			  t20 = 1/t17*t19;
			  t21 = log(-t1);
			  t23 = log(-t3);
			  t25 = t2*t21+t4*t23;
			  t28 = -t20+t6*t25*t16;
			  t30 = exp(-t7);
			  t31 = t6-1.0;
			  t32 = pow(t5,1.0*t31);
			  t33 = t30*t32;
			  t37 = 1/v[j];
			  t38 = 1/t1;
			  t39 = t37*t38;
			  t41 = t6*t2*t39*t16;
			  t42 = t2*theta;
			  t46 = t2*t37*t38;
			  t47 = t42*t39*t21+t46;
			  t50 = t5*t5;
			  t51 = 1/t50;
			  t57 = t32*t2;
			  t60 = t7*t7;
			  t64 = t13*t16;
			  t67 = t7*t28;
			  t75 = t42*t14;
			  t77 = t67*t30;
			  t83 = t16*t30;
			  t84 = t31*t25;
			  t86 = -t20+t84*t16;
			  t91 = t33*t31*t8;
			  t92 = theta*t11;
			  t104 = t33*t86;
			  t106 = t2*t11;
			  t109 = t106*t13;
			  t117 = t33*t2;
			  t122 = t21*t11;
			  out[j] = t15*t16*t28*t33+t7*(-t41+t6*t47*t16-t25*t51*t46)*t30*t57*t39-t60*
					t28*t8*t11*t64*t33+t67*t33*t31*t8*theta*t14*t16+t67*t33*t75-t77*t57*t11*t38-t77
					*t57*t14+t15*t83*t32*t86-t91*t92*t64*t86-t33*(-t41+t31*t47*t16-t84*t51*t42*t39)
					*t46-t104*t75+t104*t106*t38+t104*t109+t15*t83*t32*t21-t91*t92*t64*t21-t117*t92*
					t13*t21-t33*t109+t117*t122*t38+t117*t122*t13;
		}
		else if(*copula==5)
		{
			  t1 = exp(theta);
			  t2 = theta*u[j];
			  t3 = exp(t2);
			  t5 = t1*(t3-1.0);
			  t6 = theta*v[j];
			  t8 = exp(t6+t2);
			  t10 = exp(t6+theta);
			  t12 = exp(t2+theta);
			  t13 = t8-t10-t12+t1;
			  t14 = t13*t13;
			  t15 = 1/t14;
			  t18 = theta*t8-theta*t10;
			  t28 = v[j]+u[j];
			  t29 = v[j]+1.0;
			  out[j] = t5*t15*t18+t1*u[j]*t3*t15*t18-2.0*t5/t14/t13*(t28*t8-t29*t10-(u[j]+1.0)*
					t12+t1)*t18+t5*t15*(t8+t28*theta*t8-t10-t29*theta*t10);
		}
		else if(*copula==6)
		{
			  t1 = 1.0-u[j];
			  t2 = pow(t1,1.0*theta);
			  t3 = 1.0-v[j];
			  t4 = pow(t3,1.0*theta);
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/theta-1.0;
			  t9 = pow(t6,1.0*t8);
			  t11 = t4*theta;
			  t12 = 1/t3;
			  t13 = t11*t12;
			  t14 = theta*t12;
			  t16 = -t13+t5*t14;
			  t17 = t9*t8*t16;
			  t18 = 1/t6;
			  t19 = theta*theta;
			  t20 = 1/t19;
			  t21 = log(t6);
			  t23 = log(t1);
			  t24 = t2*t23;
			  t25 = log(t3);
			  t30 = t8*(t24+t4*t25-t24*t4-t5*t25);
			  t32 = -t20*t21+t30*t18;
			  t34 = theta-1.0;
			  t35 = pow(t3,1.0*t34);
			  t36 = 1.0-t2;
			  t37 = t35*t36;
			  t42 = t12*t25;
			  t57 = t6*t6;
			  t64 = t18*t35;
			  t67 = t9*t35;
			  t69 = t67*t34;
			  out[j] = t17*t18*t32*t37+t9*(-t20*t16*t18+t8*(-t11*t42-t4*t12+t24*t13+t5*t14
					*t25+t5*t12)*t18-t30/t57*t16)*t37-t9*t32*t35*t34*t12*t36+t17*t64*t25*t36-t69*
					t42*t36-t67*t12*t36-t17*t64*t24+t69*t12*t2*t23;
		}
	}

}


void diff2hfunc_par_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out)
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

if((*copula)==43)
	{
		ncopula=3;
		if(param[0] > 0){
			  nparam[0]=2*(param[0])/(1-param[0]);
		diff2hfunc_par_v(v, u, n, nparam, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=out[i]*2/pow(1-param[0],2);}
		}else{
			nparam[0]=-2*(param[0])/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_par_v(negv, u, n, nparam, &ncopula, out);
			for (i = 0; i < *n; i++) {out[i]=out[i]*2/pow(1+param[0],2);}
		}
	}else if((*copula)==44)
	{
		ncopula=4;
		if(param[0] > 0){
			nparam[0]=1/(1-param[0]);
		diff2hfunc_par_v(v, u, n, nparam, &ncopula, out);
		for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1-param[0],2);}
		}else{
			nparam[0]=1/(1+param[0]);
			for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
			diff2hfunc_par_v(negv, u, n, nparam, &ncopula, out);
			for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1+param[0],2);}
		}
	}else{
  if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
	  ncopula = (*copula)-20;
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
	  diff2hfunc_par_v(negv, u, n, nparam, &ncopula, out);
    }
  else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
	  ncopula = (*copula)-30;
      for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
	  diff2hfunc_par_v(v, negu, n, nparam, &ncopula, out);
    }
  else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
	{
		ncopula = (*copula)-10;
		for (i = 0; i < *n; ++i) 
		{
			negv[i] = 1 - v[i];
			negu[i] = 1 - u[i];
		}
		diff2hfunc_par_v(negv, negu, n, param, &ncopula, out);
	}
  else 
	{
		diff2hfunc_par_v(v, u, n, param, copula, out);
	}
	}
  free(negv);
  free(negu);
  free(nparam);
}


