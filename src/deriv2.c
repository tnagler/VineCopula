/*
 ** deriv2.c - C code of the package CDRVine
 **
 ** by Ulf Schepsmeier
 **
 ** Second derivatives of c (density)
 **
 */

#include "VineCopula/vine.h"
#include "VineCopula/likelihood.h"
#include "VineCopula/deriv.h"
#include "VineCopula/deriv2.h"
#include "VineCopula/tcopuladeriv.h"
#include "VineCopula/tcopuladeriv_new.h"

#define UMAX  1-1e-12

#define UMIN  1e-12

#define XEPS 1e-4


////////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach dem Parameter
// Second derivative of the bivariate copula density with respect to the parameter
//
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family (1,2,3,4,5,6)
//
// Output:
// out			derivative
//
////////////////////////////////////////////////////////////////////

void diff2PDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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

    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
        ncopula = (*copula)-20;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diff2PDF(u, negv, n, nparam, &ncopula, out);
        //for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diff2PDF(negu, v, n, nparam, &ncopula, out);
        //for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diff2PDF(negu, negv, n, param, &ncopula, out);
    }
    else
    {
        diff2PDF(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
    free(out2);
}

// vectorized version
void diff2PDF_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));
    for (int i = 0; i < (*n); ++i) {
        if (copula[i] == 2) {
            ipars[0] = par[i];
            ipars[1] = par2[i];
            diff2PDF_rho_tCopula(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diff2PDF_mod(&u[i], &v[i], &nn, &par[i], &copula[i], &out[i]);
        }
    };
    free(ipars);
}

////////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach dem Parameter
// Second derivative of the bivariate copula density with respect to the parameter (main function)
//
// u,v			copula arguments (data vectors)
// n			length of u,v
// param		parameter vector (par,par2)
// copula		copula family (1,2,3,4,5,6)
//
// Output:
// out			derivative
//
// Reference: Schepsmeier and Stoeber (2012, 2013)
////////////////////////////////////////////////////////////////////

// the structure is the same as for the first derivative
// see also the comments for the first derivatives

void diff2PDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t47, t48, t49, t51, t53, t56, t58, t59, t60, t61, t62, t65, t66;
    double t67, t70, t74, t75, t80, t87, t88;

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
            t6 = qnorm(u[j],0.0,1.0,1,0);
            t7 = qnorm(v[j],0.0,1.0,1,0);
            t1 = t6*t7;
            t2 = theta*theta;
            t3 = 1.0-t2;
            t4 = 4.0*t3*t3;
            t5 = 1/t4;
            t12 = t6*t6;
            t13 = t7*t7;
            t14 = 2.0*theta*t6*t7-t12-t13;
            t21 = t14*t5;
            t26 = 1/t3/2.0;
            t29 = exp(t12/2.0+t13/2.0+t14*t26);
            t31 = sqrt(t3);
            t32 = 1/t31;
            t38 = 2.0*t1*t26+4.0*t21*theta;
            t39 = t38*t38;
            t44 = 1/t31/t3;
            t48 = t3*t3;
            out[j] = (16.0*t1*t5*theta+16.0*t14/t4/t3*t2+4.0*t21)*t29*t32+t39*t29*t32+2.0*
                t38*t29*t44*theta+3.0*t29/t31/t48*t2+t29*t44;
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
            t15 = theta*theta;
            t16 = 1/t15;
            t17 = log(t8);
            t19 = log(u[j]);
            t21 = log(v[j]);
            t24 = -t6*t19-t7*t21;
            t26 = 1/t8;
            t27 = t16*t17+t10*t24*t26;
            t30 = -t2*t3;
            t32 = t4*t4;
            t14 = t27*t27;
            t13 = t19*t19;
            t12 = t21*t21;
            t9 = t24*t24;
            t5 = t8*t8;
            out[j] = -2.0*t3*t4*t11+2.0*t3*t11*t27+t30*t32*t11-2.0*t30*t4*t11*t27+t30*
                t11*t14+t30*t11*(-2.0/t15/theta*t17+2.0*t16*t24*t26+t10*(t6*t13+t7*t12)*t26-t10*t9/t5);
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
            t17 = t27*t29;
            t15 = t74*t74;
            t2 = t31*t35;
            t1 = t80*t80;
            out[j] = -t9*t23*t27*t40-t9*t61*t27*t40+t65*t23*t27*t40-2.0*t70*t32*t75*t38
                -2.0*t70*t32*t35*t80*t38-2.0*t70*t32*t88+t17*t31*t15*t39+t17*t31*(4.0*t44-4.0*
                    t47+2.0*t30*t53*t20-2.0*t30*t56*t59)*t39+2.0*t27*t32*t75*t80*t38+2.0*t17*t31*
                    t74*t88+t17*t2*t1*t38+2.0*t17*t2*t80*t87+t17*t2*(-2.0*t36*t22+t37*t23-
                    t37*t61);

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
            t26 = 1/t14/t13;
            t27 = v[j]+u[j];
            t29 = v[j]+1.0;
            t31 = u[j]+1.0;
            t33 = t27*t8-t29*t10-t31*t12+t1;
            t37 = theta*t1;
            t43 = t5*t26;
            t44 = t43*t33;
            t47 = theta*t18;
            t48 = t19*t19;
            t11 = t14*t14;
            t9 = t33*t33;
            t7 = t27*t27;
            t6 = t29*t29;
            t4 = t31*t31;
            out[j] = 2.0*t1*t5*t15+2.0*t18*t19*t21-4.0*t18*t5*t26*t33+t37*t21+2.0*t37*
                t19*t5*t15-4.0*t37*t44+t47*t48*t5*t15-4.0*t47*t19*t44+6.0*t47*t5/t11*t9-2.0*
                t47*t43*(t7*t8-t6*t10-t4*t12+t1);
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
            t51 = t20*t20;
            t53 = t6*t6;
            t60 = t9*t24;
            t61 = t60*t28;
            t62 = t14*t29;
            t66 = t29*t16;
            t67 = t66*t31;
            t70 = 1.0+t15+t17-t18-t19;
            t74 = t9*t28;
            out[j] = t9*t25*t32+t9*(2.0/t10/theta*t12-2.0*t11*t20*t22+t8*t49*t22-t8*t51/
                t53)*t32+2.0*t61*t62*t31+2.0*t61*t67+2.0*t60*t30*t70+t74*t41*t29*t31+2.0*t74*
                    t14*t67+2.0*t74*t62*t70+t74*t29*t43*t31+2.0*t74*t66*t70+t74*t29*t49;
        }
    }

}


///////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach u (2mal)
// second derivative with respect to u (two times)
//
////////////////////////////////////////////////////////////////////

void diff2PDF_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
    int i;

    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
        ncopula = (*copula)-20;
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diff2PDF_u(negu, v, n, nparam, &ncopula, out);
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diff2PDF_u(u, negv, n, nparam, &ncopula, out);
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diff2PDF_u(negu, negv, n, param, &ncopula, out);
    }
    else
    {
        diff2PDF_u(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diff2PDF_u_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2PDF_u_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

//////////////////////////////////////////
// Ableitung von c nach v (2 mal)
// Second derivative with respect to v (two times)
//////////////////////////////////////////

void diff2PDF_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
    int i;

    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
        ncopula = (*copula)-20;
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diff2PDF_u(v, negu, n, nparam, &ncopula, out);
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diff2PDF_u(negv, u, n, nparam, &ncopula, out);
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diff2PDF_u(negv, negu, n, param, &ncopula, out);
    }
    else
    {
        diff2PDF_u(v, u, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diff2PDF_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        diff2PDF_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
    };
    free(ipars);
}

////////////////////////////////
// Main function to calculate the derivative with respect to u
////////////////////////////////


void diff2PDF_u(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t33, t34, t35, t36, t37, t39, t40, t41, t42, t43, t45, t46, t47, t52, t53, t54, t55, t56, t57, t66, t68, t71, t79, t83, t87, t91;
    double t92, t95, t100, t106, t113, t118, t125, t128, t129 ;

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
            double nu=0, c=0, diffc=0;
            int k=1;

            LL(copula, &k, &u[j], &v[j], &theta, &nu, &c);	// one needs the density
            c=exp(c);
            diffPDF_u_mod(&u[j],&v[j],&k,param,copula,&diffc); // and also the derivative with respect to u (first derivative)

            t1=qnorm(u[j],0.0,1.0,1,0);
            t2=qnorm(v[j],0.0,1.0,1,0);
            t3=dnorm(t1,0.0,1.0,0);
            //t4=dnorm(t2,0.0,1.0,0);
            t5=1.0-theta*theta;
            t6=theta*t1-t2;
            t7=theta+t6*t1;
            t8=t7/t3/t3;
            out[j]=-theta/t5 * (diffc*t6/t3 + c*t8);
        }
        else if(*copula==2)
        {
            int k=1;
            diff2PDF_u_tCopula_new(&u[j], &v[j], &k, param, copula, &out[j]);  // special function for t-copula
        }
        else if(*copula==3)
        {
            t1 = 1.0+theta;
            t3 = pow(u[j]*v[j],-1.0*t1);
            t4 = t1*t3;
            t5 = t1*t1;
            t6 = u[j]*u[j];
            t7 = 1/t6;
            t9 = pow(u[j],-1.0*theta);
            t10 = pow(v[j],-1.0*theta);
            t11 = t9+t10-1.0;
            t13 = -2.0-1/theta;
            t14 = pow(t11,1.0*t13);
            t17 = -t1*t7;
            t21 = t14*t13;
            t22 = t9*theta;
            t23 = 1/t11;
            t28 = t13*t13;
            t31 = t9*t9;
            t32 = theta*theta;
            t34 = t11*t11;
            t37 = t31*t32*t7/t34;
            t39 = t4*t21;
            t41 = t7*t23;
            out[j] = t4*t5*t7*t14-t4*t17*t14-2.0*t4*t17*t21*t22*t23+t4*t14*t28*t37+t39*
                t9*t32*t41+t39*t22*t41-t39*t37;
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
            t10 = exp(-t9);
            t11 = u[j]*u[j];
            t13 = 1/t11/u[j];
            t14 = t10*t13;
            t15 = 1/v[j];
            t16 = t14*t15;
            t17 = -1.0+t8;
            t18 = pow(t7,2.0*t17);
            t20 = theta-1.0;
            t21 = pow(t3*t5,1.0*t20);
            t23 = pow(t7,-1.0*t8);
            t24 = t20*t23;
            t25 = t24+1.0;
            t26 = t18*t21*t25;
            t29 = t9*t4;
            t30 = 1/t3;
            t32 = 1/t7;
            t35 = t10*t15;
            t36 = t35*t26;
            t39 = t3*t3;
            t40 = 1/t39;
            t41 = t13*t40;
            t43 = t29*t41*t32;
            t45 = t15*t18;
            t46 = t14*t45;
            t47 = t21*t20;
            t52 = t4*t4;
            t53 = t9*t52;
            t54 = t7*t7;
            t55 = 1/t54;
            t56 = t41*t55;
            t57 = t53*t56;
            t66 = t35*t18;
            t68 = t21*t25*theta;
            t71 = t9*t9;
            t79 = 2.0*t45*t17;
            t83 = t47*t25;
            t87 = t47*t23;
            t91 = t14*t79;
            t92 = t4*theta;
            t95 = t32*t21*t25;
            t100 = t14*t45*t21;
            t106 = 2.0*t16*t26+3.0*t29*t13*t30*t32*t36+t43*t36-3.0*t46*t47*t30*t25-
                t57*t36-t29*theta*t13*t40*t32*t10*t15*t26+t57*t66*t68+t71*t52*t56*t36-2.0*t53*
                t13*t40*t55*t10*t79*t68-2.0*t43*t66*t83+2.0*t57*t66*t87-3.0*t91*t92*t30*t95+3.0
            *t100*t24*t4*t30*t32;
            t113 = theta*theta;
            t118 = t52*t113*t40*t55*t21*t25;
            t125 = 2.0*t18*t17;
            t128 = theta*t40;
            t129 = t128*t32;
            t19 = t128*t55;
            t12 = t40*t25;
            t2 = t20*t20;
            t1 = -t91*t92*t40*t95+4.0*t14*t45*t17*t17*t118+t91*t4*t113*t40*t95-t91*
                t118+2.0*t16*t125*t4*t129*t83-2.0*t16*t125*t52*t19*t87-t46*t47*t12+t46*t21*
                t2*t12-2.0*t100*t2*t40*t23*t4*t32+t100*t24*t4*t40*t32+t100*t24*t52*t40*t55
                -t100*t24*t4*t129+t100*t24*t52*t19;
            out[j] = t106+t1;
        }
        else if(*copula==5)
        {
            t1 = theta*theta;
            t3 = exp(theta);
            t4 = t3-1.0;
            t6 = theta*v[j];
            t7 = theta*u[j];
            t9 = exp(t6+t7+theta);
            t11 = exp(t6+t7);
            t13 = exp(t6+theta);
            t15 = exp(t7+theta);
            t16 = t11-t13-t15+t3;
            t17 = t16*t16;
            t24 = t9/t17/t16;
            t27 = theta*t11-theta*t15;
            t31 = theta*t4;
            t32 = t17*t17;
            t35 = t27*t27;
            out[j] = t1*theta*t4*t9/t17-4.0*t1*t4*t24*t27+6.0*t31*t9/t32*t35-2.0*t31*t24
                *(t1*t11-t1*t15);
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
            t10 = t8*t8;
            t12 = t2*theta;
            t13 = 1/t1;
            t17 = -t12*t13+t12*t13*t4;
            t18 = t17*t17;
            t20 = t6*t6;
            t22 = theta-1.0;
            t23 = pow(t1,1.0*t22);
            t25 = pow(t3,1.0*t22);
            t26 = theta-1.0+t2+t4-t5;
            t27 = t25*t26;
            t28 = 1/t20*t23*t27;
            t30 = t9*t8;
            t31 = theta*theta;
            t32 = t2*t31;
            t33 = t1*t1;
            t34 = 1/t33;
            t37 = t34*t4;
            t40 = t32*t34-t12*t34-t32*t37+t12*t37;
            t42 = 1/t6;
            t43 = t42*t23;
            t46 = t30*t18;
            t16 = t13*t25;
            t15 = t9*t23;
            t14 = t22*t22;
            t11 = t34*t25*t26;
            t7 = t15*t22;
            out[j] = t9*t10*t18*t28+t30*t40*t43*t27-t46*t28-2.0*t30*t17*t42*t23*t22*t16*
                t26+2.0*t46*t43*t25+t15*t14*t11-t7*t11-2.0*t7*t16*t17+t15*t25*t40;
        }
    }
}



///////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach u und v
// Second derivative with respect to u and v (first and second argument)
//
////////////////////////////////////////////////////////////////////

void diff2PDF_u_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
    int i;

    if((*copula)==43)
    {
        ncopula=3;
        if(param[0] > 0){
            nparam[0]=2*(param[0])/(1-param[0]);
            diff2PDF_u_v(u, v, n, nparam, &ncopula, out);
        }else{
            nparam[0]=-2*(param[0])/(1+param[0]);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diff2PDF_u_v(u, negv, n, nparam, &ncopula, out);
            for(i=0;i<*n;i++){out[i]=-out[i];}
        }
    }else if((*copula)==44)
    {
        ncopula=4;
        if(param[0] > 0){
            nparam[0]=1/(1-param[0]);
            diff2PDF_u_v(u, v, n, nparam, &ncopula, out);
        }else{
            nparam[0]=1/(1+param[0]);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diff2PDF_u_v(u, negv, n, nparam, &ncopula, out);
            for(i=0;i<*n;i++){out[i]=-out[i];}
        }
    }else{
        if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
        {
            ncopula = (*copula)-20;
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            diff2PDF_u_v(u, negv, n, nparam, &ncopula, out);
            for(i=0;i<*n;i++){out[i]=-out[i];}
        }
        else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
        {
            ncopula = (*copula)-30;
            for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            diff2PDF_u_v(negu, v, n, nparam, &ncopula, out);
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
            diff2PDF_u_v(negu, negv, n, param, &ncopula, out);
        }
        else
        {
            diff2PDF_u_v(u, v, n, param, copula, out);
        }
    }
    free(negv);
    free(negu);
    free(nparam);
}


void diff2PDF_u_v(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t34, t35, t36, t37, t38, t39, t40, t41, t42, t45, t48, t49, t51, t52, t53, t55, t59, t64, t70, t71, t72, t84, t85;
    double t94, t98, t101, t102, t104, t106, t107, t110, t114, t118, t128, t129 ;

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
            double c=0, diffc=0, nu=0;
            int k=1;

            LL(copula, &k, &u[j], &v[j], &theta, &nu, &c);
            c=exp(c);
            diffPDF_v_mod(&u[j],&v[j],&k,param,copula,&diffc);

            t1=qnorm(u[j],0.0,1.0,1,0);
            t2=qnorm(v[j],0.0,1.0,1,0);
            t3=dnorm(t1,0.0,1.0,0);
            t4=dnorm(t2,0.0,1.0,0);
            t5=1.0-theta*theta;
            //t6=t1-theta*t2;
            t7=theta*t1-t2;
            out[j]=-theta/t3/t5 * (diffc*t7 - c/t4);
        }
        else if(*copula==2)
        {
            int k=1;
            diff2PDF_u_v_tCopula_new(&u[j], &v[j], &k, param, copula, &out[j]);
        }
        else if(*copula==3)
        {
            t1 = 1.0+theta;
            t3 = pow(u[j]*v[j],-1.0*t1);
            t4 = t1*t3;
            t5 = t1*t1;
            t7 = 1/v[j];
            t8 = 1/u[j];
            t10 = pow(u[j],-1.0*theta);
            t11 = pow(v[j],-1.0*theta);
            t12 = t10+t11-1.0;
            t14 = -2.0-1/theta;
            t15 = pow(t12,1.0*t14);
            t23 = 1/t12;
            t35 = t14*t14;
            t39 = theta*theta;
            t41 = t12*t12;
            t42 = 1/t41;
            out[j] = t4*t5*t7*t8*t15+t4*t1*t8*t15*t14*t11*theta*t7*t23+t4*t1*t7*t15*t14*
                t10*theta*t8*t23+t4*t15*t35*t11*t39*t7*t42*t10*t8-t4*t15*t14*t10*t39*t8*t42*t11*t7;
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
            t10 = exp(-t9);
            t11 = u[j]*u[j];
            t12 = 1/t11;
            t13 = t10*t12;
            t14 = v[j]*v[j];
            t15 = 1/t14;
            t16 = t13*t15;
            t17 = -1.0+t8;
            t18 = pow(t7,2.0*t17);
            t20 = theta-1.0;
            t21 = pow(t3*t5,1.0*t20);
            t22 = t18*t21;
            t23 = pow(t7,-1.0*t8);
            t24 = t20*t23;
            t25 = 1.0+t24;
            t26 = t22*t25;
            t28 = t9*t4;
            t29 = 1/t3;
            t30 = t12*t29;
            t31 = 1/t7;
            t34 = t10*t15;
            t37 = t9*t6;
            t38 = t37*t15;
            t39 = 1/t5;
            t40 = t7*t7;
            t41 = 1/t40;
            t48 = t28*t12;
            t49 = t29*t41;
            t51 = t48*t49*t10;
            t52 = t15*t18;
            t53 = t52*t21;
            t55 = theta*t39;
            t59 = t9*t9;
            t64 = t15*t39;
            t70 = 2.0*t18*t17;
            t71 = t70*t6;
            t72 = t21*t25;
            t84 = t6*t39;
            t85 = t24*t84;
            t94 = 2.0*t13*t52*t17;
            t98 = t31*t21*t25;
            t101 = t13*t52;
            t102 = t21*t20;
            t104 = t102*t39*t25;
            t106 = t13*t53;
            t107 = t84*t31;
            t110 = t4*theta;
            t114 = t16*t26+t28*t30*t31*t34*t26-t38*t39*t41*t4*t30*t10*t26+t51*t53*t25
                *t6*t55+t59*t4*t12*t49*t6*t64*t10*t26-2.0*t48*t49*t34*t71*t55*t72-t48*t29*t31*
                    t10*t53*t20*t39*t25+2.0*t51*t53*t85+t37*t64*t31*t13*t26-t94*t6*theta*t39*t98-
                    t101*t104+t106*t24*t107-t94*t110*t29*t98;
            t118 = theta*theta;
            t128 = t4*t29;
            t129 = t16*t70*t4;
            t32 = theta*t29*t31*t104;
            t27 = t20*t20;
            t19 = t128*t31;
            t2 = t16*t22*t20;
            t1 = 4.0*t16*t18*t17*t17*t6*t118*t39*t41*t128*t72-t129*t118*t29*t41*t72
                *t84+t129*t32-2.0*t16*t70*t110*t49*t21*t85-t101*t102*t29*t25-t38*t39*t31*t10*
                    t12*t18*t21*t20*t29*t25+t16*t71*t32+t101*t21*t27*t39*t29*t25-t106*t27*t29*
                    t23*t107+t106*t24*t19-t106*t27*t39*t23*t19+t2*t23*t6*t39*t41*t4*t29+t2*
                    t23*t4*t29*t41*t6*t55;
            out[j] = t114+t1;
        }
        else if(*copula==5)
        {
            t1 = theta*theta;
            t3 = exp(theta);
            t4 = t3-1.0;
            t5 = t1*theta*t4;
            t6 = theta*v[j];
            t7 = theta*u[j];
            t9 = exp(t6+t7+theta);
            t11 = exp(t6+t7);
            t13 = exp(t6+theta);
            t15 = exp(t7+theta);
            t16 = t11-t13-t15+t3;
            t17 = t16*t16;
            t21 = t1*t4;
            t24 = t9/t17/t16;
            t25 = theta*t11;
            t27 = t25-theta*t13;
            t32 = t25-theta*t15;
            t38 = t17*t17;
            out[j] = t5*t9/t17-2.0*t21*t24*t27-2.0*t21*t24*t32+6.0*theta*t4*t9/t38*t32*
                t27-2.0*t5*t24*t11;
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
            t10 = t8*t8;
            t13 = 1/t3;
            t17 = -t4*theta*t13+t5*theta*t13;
            t18 = t6*t6;
            t19 = 1/t18;
            t22 = t2*theta;
            t23 = 1/t1;
            t27 = -t22*t23+t22*t23*t4;
            t28 = theta-1.0;
            t29 = pow(t1,1.0*t28);
            t31 = pow(t3,1.0*t28);
            t32 = theta-1.0+t2+t4-t5;
            t36 = t9*t8;
            t37 = theta*theta;
            t41 = t4*t13;
            t42 = 1/t6;
            t45 = t29*t31;
            t55 = t28*t13;
            t71 = t23*t31;
            t72 = t9*t29;
            t7 = t28*t28;
            out[j] = t9*t10*t17*t19*t27*t29*t31*t32-t36*t2*t37*t23*t41*t42*t45*t32-t36*
                t27*t19*t45*t32*t17-t36*t27*t42*t45*t55*t32+2.0*t36*t27*t42*t29*t31*t17-t36*t17
                *t42*t29*t28*t71*t32+t72*t7*t71*t13*t32-t72*t28*t71*t17-t72*t31*t55*t27-t72*
                    t31*t2*t37*t23*t41;
        }
    }
}


///////////////////////////////////////////////////////////////////
//
// 2. Ableitung von c nach par und u
// Second derivative with respect to the parameter and the first argument
//
////////////////////////////////////////////////////////////////////

void diff2PDF_par_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
    int i;

    if(((*copula==23) | (*copula==24) | (*copula==26) | (*copula==27) | (*copula==28) | (*copula==29) | (*copula==30)))	// 90? rotated copulas
    {
        ncopula = (*copula)-20;
        for (i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
        diff2PDF_par_u(negu, v, n, nparam, &ncopula, out);
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diff2PDF_par_u(u, negv, n, nparam, &ncopula, out);
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
        diff2PDF_par_u(negu, negv, n, param, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else
    {
        diff2PDF_par_u(u, v, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diff2PDF_par_u_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        if (copula[i] == 2) {
            diff2PDF_rho_u_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diff2PDF_par_u_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        }
    };
    free(ipars);
}


void diff2PDF_par_u(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t46, t47, t48, t49, t50, t51, t52, t55, t56, t59, t60, t61, t62, t63, t64, t65, t67, t68, t69, t71, t80, t81, t87;
    double t94, t97, t100, t103, t104, t106, t107,t109, t111, t112, t113, t121, t123, t126, t137, t140, t144, t146 ;

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
            double c=0, diffc=0, nu=0;
            int k=1;

            LL(copula, &k, &u[j], &v[j], &theta, &nu, &c);
            c=exp(c);
            diffPDF(&u[j],&v[j],&k,param,copula,&diffc);

            t1=qnorm(u[j],0.0,1.0,1,0);
            t2=qnorm(v[j],0.0,1.0,1,0);
            t3=1/dnorm(t1,0.0,1.0,0);
            t4=theta*(theta*t1-t2);
            t5=1.0-theta*theta;
            t6=-t4/t5;
            t7=-2.0*theta*t1+t2+t2*theta*theta;
            t8=pow(-1+theta*theta,2);
            t9=t7/t8;
            out[j]=diffc*t6*t3 + c*t3*t9;
        }
        else if(*copula==3)
        {
            t1 = u[j]*v[j];
            t2 = -theta-1.0;
            t3 = pow(t1,1.0*t2);
            t5 = 1/u[j];
            t6 = pow(u[j],-1.0*theta);
            t7 = pow(v[j],-1.0*theta);
            t8 = t6+t7-1.0;
            t9 = 1/theta;
            t10 = -2.0-t9;
            t11 = pow(t8,1.0*t10);
            t12 = t5*t11;
            t16 = t6*theta;
            t17 = 1/t8;
            t18 = t5*t17;
            t21 = -t3*t2;
            t22 = t21*t2;
            t23 = log(t1);
            t35 = theta*theta;
            t37 = log(t8);
            t39 = log(u[j]);
            t41 = log(v[j]);
            t44 = t10*(-t6*t39-t7*t41);
            t46 = 1/t35*t37+t44*t17;
            t62 = t8*t8;
            out[j] = t3*t2*t12-t3*t11*t10*t16*t18-t22*t5*t23*t11-t21*t12+t21*t23*t11*t10
                *t6*theta*t5*t17+t22*t12*t46-t21*t11*t10*t16*t18*t46+t21*t11*(-t9*t6*t18+t10*(
                        t16*t5*t39-t6*t5)*t17+t44/t62*t16*t5);
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
            t17 = t9*t4*t12*t13*t15;
            t18 = theta*theta;
            t20 = log(t7);
            t21 = 1/t18*t20;
            t22 = log(-t3);
            t24 = log(-t5);
            t26 = t4*t22+t6*t24;
            t29 = -t21+t8*t26*t15;
            t30 = exp(-t9);
            t32 = 1/v[j];
            t34 = -1.0+t8;
            t35 = pow(t7,2.0*t34);
            t36 = t3*t5;
            t37 = theta-1.0;
            t38 = pow(t36,1.0*t37);
            t39 = t35*t38;
            t40 = pow(t7,-1.0*t8);
            t41 = t37*t40;
            t42 = t41+1.0;
            t43 = t39*t42;
            t47 = 1/u[j];
            t48 = t47*t13;
            t49 = t48*t15;
            t50 = t8*t4*t49;
            t51 = t4*theta;
            t55 = t4*t47*t13;
            t56 = t51*t48*t22+t55;
            t59 = t7*t7;
            t60 = 1/t59;
            t63 = -t50+t8*t56*t15-t26*t60*t55;
            t65 = t30*t47;
            t67 = t32*t35;
            t68 = t38*t42;
            t69 = t67*t68;
            t71 = t9*t9;
            t80 = t9*t29;
            t81 = t30*t12;
            t87 = t80*t30*t12*t32*t35;
            t94 = t81*t32;
            t97 = t37*t13*t42;
            t100 = t38*t37;
            t103 = t4*t13*t15;
            t104 = t100*t40*t103;
            t106 = t30*t32;
            t107 = t106*t35;
            t109 = 2.0*t34*t26;
            t111 = -2.0*t21+t109*t15;
            t112 = t111*t38;
            t113 = t112*t42;
            t121 = 2.0*t94*t35*t34*t4;
            t123 = theta*t13*t15;
            t126 = t65*t32;
            t137 = t81*t67;
            t140 = -t17*t29*t30*t32*t43-t9*t63*t65*t69+t71*t29*t4*t12*t13*t15*t30*t32
                *t43+t80*t81*t69-2.0*t87*t34*t4*theta*t13*t15*t68-t80*t94*t39*t97+t87*t104-t17*
                    t107*t113-t94*t35*t111*t68+t121*t123*t113+t126*t35*(-2.0*t50+2.0*t34*t56*t15-
                    t109*t60*t51*t48)*t68+t137*t112*t97;
            t144 = log(t36);
            t146 = t38*t144*t42;
            t10 = t40-t41*t29;
            t2 = t39*t10;
            t1 = -t81*t67*t111*t104-t17*t107*t146-t94*t39*t144*t42+t121*t123*t146+
                t137*t100*t13*t144*t42+t94*t39*t13*t42-t81*t67*t38*t144*t37*t40*t103-t17*t106*
                t2-t94*t2+2.0*t81*t67*t34*t51*t13*t15*t38*t10+t137*t100*t13*t10+t126*t39*
                (-t40*t4*t49+t41*t4*t48*t15*t29-t41*t63);
            out[j] = t140+t1;
        }
        else if(*copula==5)
        {
            t1 = exp(theta);
            t2 = t1-1.0;
            t3 = theta*t2;
            t4 = theta*v[j];
            t5 = theta*u[j];
            t7 = exp(t4+t5+theta);
            t9 = exp(t4+t5);
            t11 = exp(t4+theta);
            t13 = exp(t5+theta);
            t14 = t9-t11-t13+t1;
            t15 = t14*t14;
            t16 = 1/t15;
            t17 = t7*t16;
            t22 = 1/t15/t14;
            t25 = theta*t9-theta*t13;
            t29 = theta*theta;
            t33 = t7*t22;
            t34 = t33*t25;
            t37 = t29*t2;
            t38 = v[j]+u[j]+1.0;
            t46 = v[j]+u[j];
            t49 = u[j]+1.0;
            t51 = t46*t9-(v[j]+1.0)*t11-t49*t13+t1;
            t56 = t15*t15;
            out[j] = 2.0*t3*t17-2.0*t2*t7*t22*t25+t29*t1*t17-2.0*theta*t1*t34+t37*t38*t7
                *t16-2.0*t3*t38*t34-2.0*t37*t33*t51+6.0*t3*t7/t56*t51*t25-2.0*t3*t33*(t9+t46*
                    theta*t9-t13-t49*theta*t13);
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
            t10 = t9*t8;
            t11 = t2*theta;
            t12 = 1/t1;
            t14 = t12*t4;
            t16 = -t11*t12+t11*t14;
            t17 = 1/t6;
            t19 = t10*t16*t17;
            t20 = theta*theta;
            t21 = 1/t20;
            t22 = log(t6);
            t24 = log(t1);
            t25 = t2*t24;
            t26 = log(t3);
            t27 = t4*t26;
            t28 = t25*t4;
            t29 = t5*t26;
            t31 = t8*(t25+t27-t28-t29);
            t33 = -t21*t22+t31*t17;
            t34 = theta-1.0;
            t35 = pow(t1,1.0*t34);
            t37 = pow(t3,1.0*t34);
            t38 = theta-1.0+t2+t4-t5;
            t39 = t37*t38;
            t44 = t12*t24;
            t46 = t2*t12;
            t52 = -t11*t44-t46+t11*t44*t4+t46*t4+t11*t14*t26;
            t55 = t6*t6;
            t61 = t35*t37;
            t64 = t9*t33;
            t71 = t9*t35;
            t80 = t71*t34;
            t81 = t12*t37;
            t87 = t26*t38;
            t94 = 1.0+t25+t27-t28-t29;
            out[j] = t19*t33*t35*t39+t9*(-t21*t16*t17+t8*t52*t17-t31/t55*t16)*t61*t38-
                t64*t35*t34*t12*t39+t64*t61*t16+t19*t35*t24*t39-t80*t44*t39-t71*t81*t38+t71*t24
                *t37*t16+t19*t61*t87-t80*t81*t87+t71*t37*t26*t16+t10*t16*t17*t35*t37*t94-t80*
                    t81*t94+t71*t37*t52;
        }
    }
}



// The same with respect to the parameter and v (second argument)

void diff2PDF_par_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out)
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
        diff2PDF_par_u(v, negu, n, nparam, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else if(((*copula==33) | (*copula==34) | (*copula==36) | (*copula==37) | (*copula==38) | (*copula==39) | (*copula==40)))	// 270? rotated copulas
    {
        ncopula = (*copula)-30;
        for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
        diff2PDF_par_u(negv, u, n, nparam, &ncopula, out);
    }
    else if(((*copula==13) | (*copula==14) | (*copula==16) | (*copula==17) | (*copula==18) | (*copula==19) | (*copula==20)))	// 180? rotated copulas
    {
        ncopula = (*copula)-10;
        for (i = 0; i < *n; ++i)
        {
            negv[i] = 1 - v[i];
            negu[i] = 1 - u[i];
        }
        diff2PDF_par_u(negv, negu, n, param, &ncopula, out);
        for(i=0;i<*n;i++){out[i]=-out[i];}
    }
    else
    {
        diff2PDF_par_u(v, u, n, param, copula, out);
    }
    free(negv);
    free(negu);
    free(nparam);
}

// vectorized version
void diff2PDF_par_v_mod_vec(double* u, double* v, int* n, double* par, double* par2, int* copula, double* out)
{
    int nn = 1;
    double* ipars = (double *) malloc(2*sizeof(double));

    for (int i = 0; i < (*n); ++i) {
        ipars[0] = par[i];
        ipars[1] = par2[i];
        if (copula[i] == 2) {
            diff2PDF_rho_v_tCopula_new(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        } else {
            diff2PDF_par_v_mod(&u[i], &v[i], &nn, ipars, &copula[i], &out[i]);
        }
    };
    free(ipars);
}

