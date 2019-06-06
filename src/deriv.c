/*
 ** deriv.c - C code of the package CDRVine
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



/////////////////////////////////////////////////////////////
//
// Ableitung der Copula nach dem Parameter
// Derivative of bivariate copulas with respect to the (first) parameter
//
// Input:
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family
//
// Output:
// out			derivative
/////////////////////////////////////////////////////////////

void diffPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    // for the rotated copulas we need some help variables
    double* negv;
    double* negu;
    double* nparam;
    negv = (double *) malloc(*n*sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    nparam = (double *) malloc(2*sizeof(double));
    int ncopula;
    nparam[0]=-param[0];
    nparam[1]=-param[1];
    int i;
    // for the rotation see the master thesis of Jakob Stoeber
    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90 rotated copulas
    {
        ncopula = (*copula)-20;
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diffPDF(negu, v, n, nparam, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270 rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diffPDF(u, negv, n, nparam, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180 rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diffPDF(negu, negv, n, param, &ncopula, out);
    }
    else
    {
        diffPDF(u, v, n, param, copula, out);		// eigentliche Ableitungsfunktion
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diffPDF_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    for (int i = 0; i < (*n); ++i) {
        if (copula[i] == 2) {
            ipars[0] = par[i];
            ipars[1] = par2[i];
            diffPDF_rho_tCopula(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diffPDF_mod(&u[i], &v[i], &nn, &par[i], &copula[i], &out[i]);
        }
    };
    free(ipars);
}


//////////////////////////////////////////////////
// Derivative of bivariate copulas with respect to the parameter (standard copula form without rotations, see above)
//
// Input:
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family (1,3,4,5,6)
//
// Output:
// out			derivative
//
// Reference: Schepsmeier and Stoeber (2012, 2013)
/////////////////////////////////////////////////////////////

// the stepwise calculation is due to performance and numerical stability reasons (t1,t2,...)
// for Gauss some of the step can be found in the reference
// for the archimedean copulas one gets this optimization with Maple

void diffPDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t24, t26, t27, t28, t29, t30, t32, t33, t34, t35, t36;
    t3=0;
    t4=0;

    double theta = param[0];
    //double delta = param[1];

    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }

    for(j=0;j<*n;j++)
    {
        if(*copula==0)	// independence copulas
        {
            out[j]=0;
        }
        else if(*copula==1)		// gauss, formula see reference
        {
            t1 = qnorm(u[j],0.0,1.0,1,0);
            t2 = qnorm(v[j],0.0,1.0,1,0);
            t4 = t1*t1;
            t5 = t2*t2;
            t3 = t4+t5;
            t7 = theta*theta;
            t8 = 1.0-t7;
            t9 = 1/t8/2.0;
            t15 = t7*t3-2.0*theta*t1*t2;
            t22 = exp(-t15*t9);
            t24 = sqrt(t8);
            out[j] = (-2.0*(theta*t3-t1*t2)*t9-t15/(t8*t8)*theta)*t22/t24+t22/t24/t8*theta;
        }
        // t-copula is separate; very complicated
        else if(*copula==3)	// the archimedean copula derivatives are derived by Maple
        {
            t1 = u[j]*v[j];
            t2 = -theta-1.0;
            t3 = pow(t1,1.0*t2);
            t4 = pow(u[j],-1.0*theta);
            t5 = pow(v[j],-1.0*theta);
            t6 = t4+t5-1.0;
            t7 = -2.0-1/theta;
            t8 = pow(t6,1.0*t7);
            t9 = -t2*t3;
            t10 = log(t1);
            t11 = theta*theta;
            t12 = log(t6);
            t13 = log(u[j]);
            t14 = log(v[j]);
            out[j] = t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);
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
            t12 = log(t7);
            t13 = 1/t10*t12;
            t14 = log(-t3);
            t16 = log(-t5);
            t18 = t4*t14+t6*t16;
            t20 = 1/t7;
            t22 = -t13+t8*t18*t20;
            t24 = exp(-t9);
            t26 = t24/u[j];
            t28 = 1/v[j];
            t29 = -1.0+t8;
            t30 = pow(t7,2.0*t29);
            t32 = t3*t5;
            t33 = theta-1.0;
            t34 = pow(t32,1.0*t33);
            t35 = pow(t7,-1.0*t8);
            t36 = t33*t35;
            t17 = 1.0+t36;
            t15 = t34*t17;
            t11 = t26*t28;
            t2 = t30*t34;
            t1 = log(t32);
            out[j] = -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*
                t2*t1*t17+t11*t2*(t35-t36*t22);

        }
        else if(*copula==5)
        {
            t2 = exp(theta);
            t3 = t2-1.0;
            t4 = theta*v[j];
            t5 = theta*u[j];
            t7 = exp(t4+t5+theta);
            t10 = exp(t4+t5);
            t12 = exp(t4+theta);
            t14 = exp(t5+theta);
            t15 = t10-t12-t14+t2;
            t16 = t15*t15;
            t17 = 1/t16;
            t21 = theta*t3;
            out[j] = t3*t7*t17+theta*t2*t7*t17+t21*(v[j]+u[j]+1.0)*t7*t17-2.0*t21*t7/t15/t16*((v[j]+u[j])*t10-(v[j]+1.0)*t12-(u[j]+1.0)*t14+t2);
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
            t11 = log(t6);
            t12 = log(t1);
            t13 = t2*t12;
            t14 = log(t3);
            t15 = t4*t14;
            t16 = t13*t4;
            t19 = t5*t14;
            t21 = theta-1.0;
            t27 = pow(t1,1.0*t21);
            t28 = pow(t3,1.0*t21);
            t30 = theta-1.0+t2+t4-t5;
            t33 = t9*t27;
            out[j] = t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30
                +t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
        }

    }

}

////////////////////////////////////////////////////////////////////
//
// 1. Ableitung von c nach u
// First derivative of the bivariate copula density with respect to u (first argument)
// Input:
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family
//
// Output:
// out			derivative
//
////////////////////////////////////////////////////////////////////

void diffPDF_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
        diffPDF_u(negu, v, n, nparam, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diffPDF_u(u, negv, n, nparam, &ncopula, out);
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diffPDF_u(negu, negv, n, param, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else
    {
        diffPDF_u(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diffPDF_u_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out) {
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffPDF_u_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}


////////////////////////////////////////////////////////////////////
//
// 1. Ableitung von c nach u (eigentliche Funktion ohne die Rotationen)
// First derivative of the bivariate copula density with respect to u (first argument)
// Input:
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family (1,2,3,4,5,6)
//
// Output:
// out			derivative
//
////////////////////////////////////////////////////////////////////

void diffPDF_u(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j, k=1;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t18, t19, t20, t21, t22, t23, t24, t25, t27, t28, t29;
    double t30, t33, t36;

    double theta = param[0];
    //double delta = param[1];

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
            t5 = sqrt(t4);
            t6 = pow(t1,2.0);
            t7 = pow(t2,2.0);
            t8 = t3*(t6+t7)-(2.0*theta*t1*t2);
            t9 = exp(-t8/t4/2.0);
            t10 = t9/t5;
            t11 = sqrt(2.0*pi);
            t12 = exp(-t6/2.0);
            t13 = t11/t12;
            out[j] = -t10*(theta*t13/t4)*(theta*t1-t2);
        }
        else if(*copula==2)
        {
            diffPDF_u_tCopula_new(&u[j], &v[j], &k, param, copula, &out[j]);	// special function for t-copula
        }
        else if(*copula==3)
        {
            t1 = 1.0+theta;
            t3 = pow(u[j]*v[j],-1.0*t1);
            t4 = t1*t3;
            t5 = 1/u[j];
            t7 = pow(u[j],-1.0*theta);
            t8 = pow(v[j],-1.0*theta);
            t9 = t7+t8-1.0;
            t11 = -2.0-1/theta;
            t12 = pow(t9,1.0*t11);
            out[j] = -t4*t1*t5*t12-t4*t12*t11*t7*theta*t5/t9;
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
            t11 = u[j]*u[j];
            t12 = 1/t11;
            t13 = 1/t3;
            t15 = 1/t7;
            t18 = exp(-t9);
            t19 = 1/v[j];
            t21 = -1.0+t8;
            t22 = pow(t7,2.0*t21);
            t24 = theta-1.0;
            t25 = pow(t3*t5,1.0*t24);
            t27 = pow(t7,-1.0*t8);
            t28 = t24*t27;
            t29 = 1.0+t28;
            t30 = t22*t25*t29;
            t33 = t18*t12;
            t36 = t19*t22;
            out[j] = -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*theta
                *t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15;
        }
        else if(*copula==5)
        {
            t1 = theta*theta;
            t2 = exp(theta);
            t3 = t2-1.0;
            t5 = theta*v[j];
            t6 = theta*u[j];
            t8 = exp(t5+t6+theta);
            t10 = exp(t5+t6);
            t12 = exp(t5+theta);
            t14 = exp(t6+theta);
            t15 = t10-t12-t14+t2;
            t16 = t15*t15;
            out[j] = t1*t3*t8/t16-2.0*theta*t3*t8/t16/t15*(theta*t10-theta*t14);
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
            t11 = t2*theta;
            t12 = 1/t1;
            t16 = -t11*t12+t11*t12*t4;
            t19 = theta-1.0;
            t20 = pow(t1,1.0*t19);
            t22 = pow(t3,1.0*t19);
            t23 = theta-1.0+t2+t4-t5;
            t27 = t9*t20;
            out[j] = t9*t8*t16/t6*t20*t22*t23-t27*t19*t12*t22*t23+t27*t22*t16;
        }
    }

}


////////////////////////////////////////////////////////////////////
//
// 1. Ableitung von c nach v
// First derivative of the bivariate copula density with respect to v (second argument)
// Input:
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family (1,2,3,4,5,6)
//
// Output:
// out			derivative
//
////////////////////////////////////////////////////////////////////

void diffPDF_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
        diffPDF_u(v, negu, n, nparam, &ncopula, out);		// we can use again the function for the derivative of c wrt u but change the arguments
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diffPDF_u(negv, u, n, nparam, &ncopula, out);
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
        diffPDF_u(negv, negu, n, param, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else
    {
        diffPDF_u(v, u, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diffPDF_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out) {
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diffPDF_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

