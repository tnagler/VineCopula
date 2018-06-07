/*
 ** hfuncderiv.c - C code of the package CDRVine
 **
 ** by Ulf Schepsmeier
 **
 **
 **
 */

#include "VineCopula/vine.h"
#include "VineCopula/deriv.h"
#include "VineCopula/tcopuladeriv.h"
#include "VineCopula/tcopuladeriv_new.h"

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4


////////////////////////////////////////////////////////////////////
//
// Ableitung der h-function nach dem Parameter
//
////////////////////////////////////////////////////////////////////

void diffhfunc_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diffhfunc(negu, v, n, nparam, &ncopula, out);
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diffhfunc(negv, u, n, nparam, &ncopula, out);
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
        diffhfunc(negu, negv, n, param, &ncopula, out);
        for (i = 0; i < *n; i++) {out[i]=-out[i];};
    }
    else
    {
        diffhfunc(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diffhfunc_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    for (int i = 0; i < (*n); ++i) {
        if (copula[i] == 2) {
            ipars[0] = par[i];
            ipars[1] = par2[i];
            diffhfunc_rho_tCopula(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diffhfunc_mod(&u[i], &v[i], &nn, &par[i], &copula[i], &out[i]);
        }
    };
    free(ipars);
}


void diffhfunc_mod2(double* v, double* u, int* n, double* param, int* copula, double* out)  // Achtung u und v vertauscht; Notaion aus Hfunc1 und Hfunc2
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
            diffhfunc(v, u, n, nparam, &ncopula, out);
            for (i = 0; i < *n; i++) {out[i]=out[i]*2/pow(1-param[0],2);}
        }else{
            nparam[0]=-2*(param[0])/(1+param[0]);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc(negv, u, n, nparam, &ncopula, out);
            for (i = 0; i < *n; i++) {out[i]=out[i]*2/pow(1+param[0],2);};
        }
    }else if((*copula)==44)
    {
        ncopula=4;
        if(param[0] > 0){
            nparam[0]=1/(1-param[0]);
            diffhfunc(v, u, n, nparam, &ncopula, out);
            for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1-param[0],2);}
        }else{
            nparam[0]=1/(1+param[0]);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc(negv, u, n, nparam, &ncopula, out);
            for (i = 0; i < *n; i++) {out[i]=out[i]/pow(1+param[0],2);};
        }
    }else{
        if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
        {
            ncopula = (*copula)-20;
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc(negv, u, n, nparam, &ncopula, out);
        }
        else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
        {
            ncopula = (*copula)-30;
            for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            diffhfunc(v, negu, n, nparam, &ncopula, out);
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
            diffhfunc(negv, negu, n, param, &ncopula, out);
            for (i = 0; i < *n; i++) {out[i]=-out[i];};
        }
        else
        {
            diffhfunc(v, u, n, param, copula, out);
        }
    }
    free(negv);
    free(negu);
    free(nparam);
}


void diffhfunc(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t14, t15, t16, t18, t22, t24, t25, t27, t28, t32;

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
            t3=t1-theta*t2;
            t4=1.0-pow(theta,2);
            t5=sqrt(t4);
            t6=t3/t5;
            t7=dnorm(t6,0.0,1.0,0);
            t8=-1.0*t2*t5+1.0*t3*theta/t5;
            t9=t8/t4;
            out[j]=t7*t9;
        }
        else if(*copula==3)
        {
            t1 = pow(v[j],-1.0*theta-1.0);
            t2 = log(v[j]);
            t3 = pow(u[j],-1.0*theta);
            t4 = pow(v[j],-1.0*theta);
            t5 = t3+t4-1.0;
            t6 = -1.0-1/theta;
            t7 = pow(t5,1.0*t6);
            t8 = theta*theta;
            t9 = log(t5);
            t10 = log(u[j]);
            out[j] = -t1*t2*t7+t1*t7*(1/t8*t9+t6*(-t3*t10-t4*t2)/t5);
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
            t9 = log(t5);
            t10 = 1/t8*t9;
            t11 = log(-t1);
            t14 = log(-t3);
            t16 = t2*t11+t4*t14;
            t18 = 1/t5;
            t22 = exp(-t7);
            t24 = t6-1.0;
            t25 = pow(t5,1.0*t24);
            t27 = 1/v[j];
            t28 = 1/t1;
            t32 = t22*t25;
            out[j] = t7*(-t10+t6*t16*t18)*t22*t25*t2*t27*t28-t32*(-t10+t24*t16*t18)*t2*t27*t28-t32*t2*t11*t27*t28;
        }
        else if(*copula==5)
        {
            t1 = exp(theta);
            t2 = theta*u[j];
            t3 = exp(t2);
            t5 = t1*(t3-1.0);
            t6 = theta*v[j];
            t8 = exp(t6+t2);
            t9 = exp(t6+theta);
            t10 = exp(t2+theta);
            t11 = t8-t9-t10+t1;
            t14 = 1/t11;
            t18 = t11*t11;
            out[j] = -t5*t14-t1*u[j]*t3*t14+t5/t18*((v[j]+u[j])*t8-(v[j]+1.0)*t9-(u[j]+1.0)*t10+t1);
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
            t12 = log(t6);
            t14 = log(t1);
            t15 = t2*t14;
            t16 = log(t3);
            t27 = pow(t3,1.0*theta-1.0);
            t7 = 1.0-t2;
            t11 = t9*t27;
            out[j] = t9*(-1.0/t10*t12+t8*(t15+t4*t16-t15*t4-t5*t16)/t6)*t27*t7+t11*t16*t7-t11*t15;
        }
    }

}


////////////////////////////////////////////////////////////////////
//
// Ableitung der h-function nach v
//
////////////////////////////////////////////////////////////////////

void diffhfunc_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    double* negv;
    double* negu;
    double* nparam;
    negv = (double *) malloc(*n*sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    nparam = (double *) malloc(2*sizeof(double));
    int ncopula;
    nparam[0]=-param[0];
    nparam[1]=-param[1];

    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
        ncopula = (*copula)-20;
        for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diffhfunc_v(negu, v, n, nparam, &ncopula, out);
        for (int i = 0; i < *n; i++) {out[i]=-out[i];};
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diffhfunc_v(u, negv, n, nparam, &ncopula, out);
        for (int i = 0; i < *n; i++) {out[i]=-out[i];};
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (int i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diffhfunc_v(negu, negv, n, param, &ncopula, out);
    }
    else
    {
        diffhfunc_v(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diffhfunc_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffhfunc_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


void diffhfunc_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out)	// Achtung u und v vertauscht
{
    double* negv;
    double* negu;
    double* nparam;
    negv = (double *) malloc(*n*sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    nparam = (double *) malloc(2*sizeof(double));
    int ncopula;
    nparam[0]=-param[0];
    nparam[1]=-param[1];

    if((*copula)==43)
    {
        ncopula=3;
        if(param[0] > 0){
            nparam[0]=2*(param[0])/(1-param[0]);
            diffhfunc_v(v, u, n, nparam, &ncopula, out);
        }else{
            nparam[0]=-2*(param[0])/(1+param[0]);
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc_v(negv, u, n, nparam, &ncopula, out);
            for (int i = 0; i < *n; i++) {out[i]=-out[i];};
        }
    }else if((*copula)==44)
    {
        ncopula=4;
        if(param[0] > 0){
            nparam[0]=1/(1-param[0]);
            diffhfunc_v(v, u, n, nparam, &ncopula, out);
        }else{
            nparam[0]=1/(1+param[0]);
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc_v(negv, u, n, nparam, &ncopula, out);
            for (int i = 0; i < *n; i++) {out[i]=-out[i];};
        }
    }else{
        if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
        {
            ncopula = (*copula)-20;
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diffhfunc_v(negv, u, n, nparam, &ncopula, out);
            for (int i = 0; i < *n; i++) {out[i]=-out[i];};
        }
        else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
        {
            ncopula = (*copula)-30;
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            diffhfunc_v(v, negu, n, nparam, &ncopula, out);
            for (int i = 0; i < *n; i++) {out[i]=-out[i];};
        }
        else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
        {
            ncopula = (*copula)-10;
            for (int i = 0; i < *n; ++i)
            {
                negv[i] = 1 - v[i];
                negu[i] = 1 - u[i];
            }
            diffhfunc_v(negv, negu, n, param, &ncopula, out);
        }
        else
        {
            diffhfunc_v(v, u, n, param, copula, out);
        }
    }
    free(negv);
    free(negu);
    free(nparam);
}


void diffhfunc_v(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j, k=1;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t12, t13, t15, t16, t18, t19, t20, t21, t22, t27, t33;

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
            t3=t1-theta*t2;
            t4=1.0-pow(theta,2);
            t5=sqrt(t4);
            t6=t3/t5;
            t7=dnorm(t6,0.0,1.0,0);
            t8=sqrt(2.0*pi);
            t9=pow(t2,2);
            t10=exp(-t9/2.0);
            out[j]=t7*t8*(-theta)/t5/t10;
        }
        else if(*copula==2)
        {
            diffhfunc_v_tCopula_new(&u[j], &v[j], &k, param, copula, &out[j]);
        }
        else if(*copula==3)
        {
            t1 = -theta-1.0;
            t2 = pow(v[j],1.0*t1);
            t4 = 1/v[j];
            t5 = pow(u[j],-1.0*theta);
            t6 = pow(v[j],-1.0*theta);
            t7 = t5+t6-1.0;
            t9 = -1.0-1/theta;
            t10 = pow(t7,1.0*t9);
            out[j] = t10*t4*t1*t2-1/t7*t4*theta*t6*t9*t10*t2;
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
            t10 = t6*t6;
            t12 = v[j]*v[j];
            t13 = 1/t12;
            t15 = t5*t5;
            t16 = 1/t15;
            t18 = t16/t7;
            t19 = exp(-t9);
            t20 = t8-1.0;
            t21 = pow(t7,1.0*t20);
            t22 = t19*t21;
            t27 = theta*t13;
            t33 = t6*t13;
            out[j] = t9*t10*t13*t18*t22-t22*t20*t10*t27*t18-t22*t6*t27*t16+t22*t33/t5+t22*t33*t16;
        }
        else if(*copula==5)
        {
            t1 = exp(theta);
            t2 = theta*u[j];
            t3 = exp(t2);
            t6 = theta*v[j];
            t8 = exp(t6+t2);
            t10 = exp(t6+theta);
            t12 = exp(t2+theta);
            t13 = pow(t8-t10-t12+t1,2.0);
            out[j] = t1*(t3-1.0)/t13*(theta*t8-theta*t10);
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
            t12 = 1/t3;
            t19 = theta-1.0;
            t20 = pow(t3,1.0*t19);
            t22 = 1.0-t2;
            out[j] = t9*t8*(-t4*theta*t12+t5*theta*t12)/t6*t20*t22-t9*t20*t19*t12*t22;
        }
    }

}

