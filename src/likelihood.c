/*
 ** likelihood.c - C code of the package CDRVine
 **
 ** with contributions from Carlos Almeida, Aleksey Min,
 ** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
 **
 ** A first version was based on code
 ** from Daniel Berg <daniel at danielberg.no>
 ** provided by personal communication.
 **
 */

#include "include/vine.h"
#include "include/memoryhandling.h"
#include "include/tools.h"
#include "include/likelihood.h"
#include "include/evCopula.h"
#include <math.h>

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4

#define XINFMAX DBL_MAX


///////////////////////////////////////////////////////
// New

double log1mexp(double a)
{
    double result;
    if (a<log(2)) {
        result=log(-expm1(-a));
    }else{
        result=log1p(-exp(-a));
    }
    return result;
}



void archCDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
    int j;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;

    for(j=0;j<*n;j++)
    {
        if(u[j]>UMAX && v[j]>UMAX){ out[j]=1;}
        else if(u[j]>UMAX){ out[j]=v[j];}
        else if(v[j]>UMAX){ out[j]=u[j];}
        else if(u[j]<UMIN || v[j]<UMIN){ out[j]=0;}
        else
        {
            if(*copula==3)	//Clayton
            {
                t1 = pow(u[j],-param[0]);
                t2 = pow(v[j],-param[0]);
                t3 = t1+t2-1;
                out[j] = pow(t3,-1/param[0]);
            }
            else if(*copula==4)	//Gumbel
            {
                t1 = -log(u[j]);
                t2 = -log(v[j]);
                t3 = pow(t1,param[0]);
                t4 = pow(t2,param[0]);
                t5 = t3+t4;
                t6 = -pow(t5,1/param[0]);
                out[j] = exp(t6);
            }
            else if(*copula==5)	//Frank
            {
                if (param[0]>0) {
                    t1=-log1p(exp(-param[0]) * expm1(param[0]-u[j]*param[0])/expm1(-param[0]));
                    t2=-log1p(exp(-param[0]) * expm1(param[0]-v[j]*param[0])/expm1(-param[0]));
                    out[j] = -log1mexp(t1+t2-log1mexp(param[0]))/param[0];
                } else {
                    out[j] =-1/param[0] * log(1 + exp(-(-log((exp(-param[0] * u[j]) - 1)/(exp(-param[0]) - 1)) + -log((exp(-param[0] * v[j]) - 1)/(exp(-param[0]) - 1)))) * (exp(-param[0]) - 1));
                }
            }
            else if(*copula==6)	//Joe
            {
                t1 = 1-u[j];
                t2 = 1-v[j];
                t3 = pow(t1,param[0]);
                t4 = pow(t2,param[0]);
                t5 = t3*t4;
                out[j] = 1-pow(t3+t4-t5,1/param[0]);
            }
            else if(*copula==7)	//BB1
            {
                t1 = pow(u[j],-param[0]);
                t2 = pow(v[j],-param[0]);
                t3 = t1-1;
                t4 = t2-1;
                t5 = pow(t3,param[1]);
                t6 = pow(t4,param[1]);
                t7 = t5+t6;
                t8 = pow(t7,1/param[1]);
                out[j] = pow(1+t8,-1/param[0]);
            }
            else if(*copula==8)	//BB6
            {
                t1 = 1-u[j];
                t2 = 1-v[j];
                t3 = pow(t1,param[0]);
                t4 = pow(t2,param[0]);
                t5 = 1-t3;
                t6 = 1-t4;
                t7 = -log(t5);
                t8 = -log(t6);
                t9 = pow(t7,param[1]);
                t10 = pow(t8,param[1]);
                t11 = t9+t10;
                t12 = pow(t11,1/param[1]);
                t13 = exp(-t12);
                t14 = 1-t13;
                out[j] = 1-pow(t14,1/param[0]);
            }
            else if(*copula==9)	//BB7
            {
                t1 = 1-u[j];
                t2 = 1-v[j];
                t3 = pow(t1,param[0]);
                t4 = pow(t2,param[0]);
                t5 = 1-t3;
                t6 = 1-t4;
                t7 = pow(t5,-param[1]);
                t8 = pow(t6,-param[1]);
                t9 = t7+t8-1;
                t10 = pow(t9,-1/param[1]);
                t11 = 1-t10;
                t12 = pow(t11,1/param[0]);
                out[j] = 1-t12;
            }
            else if(*copula==10)    //BB8
            {
                double nu;
                t1 = param[1]*u[j];
                t2 = param[1]*v[j];
                t3 = 1-t1;
                t4 = 1-t2;
                t5 = pow(t3,param[0]);
                t6 = pow(t4,param[0]);
                t7 = 1-t5;
                t8 = 1-t6;
                nu = 1-param[1];
                nu = pow(nu,param[0]);
                nu = 1-nu;
                nu = 1/nu;
                t9 = 1-nu*t7*t8;
                t10 = pow(t9,1/param[0]);
                out[j] = 1/param[1]*(1-t10);
            }
            else if(*copula==41)
            {
                t1=qgamma(1.0-u[j],param[0],1,1,0);
                t2=qgamma(1.0-v[j],param[0],1,1,0);
                t3=pow(pow(t1,param[0])+pow(t2,param[0]),(1.0/param[0]));
                out[j]=1.0-pgamma(t3,param[0],1,1,0);
            }
        }
    }

}





void dbb1(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t16, t17, t38, t39, t4, t5, t6, t7, t9, t10, t12, t13, t20, t24, t25, t27, t29, t32, t33, t34, t36, t43, t59;

    th = param[0];
    de = param[1];

    for(i=0;i<*n;i++)
    {
        t1 = pow(u[i],(-th));
        t2 = t1-1.0;
        t3 = pow(t2,de);
        t16 = 1./u[i];
        t17 = 1./t2;
        t38 = t1*t16;
        t39 = t38*t17;
        t4 = pow(v[i],(-th));
        t5 = t4-1.0;
        t6 = pow(t5,de);
        t7 = t3+t6;
        t9 = pow(t7,(1./de));
        t10 = 1.0+t9;
        t12 = pow(t10,(-1./th));
        t13 = t12*t9;
        t20 = 1./t10;
        t24 = t9*t9;
        t25 = t12*t24;
        t27 = 1./v[i];
        t29 = 1./t5;
        t32 = t7*t7;
        t33 = 1./t32;
        t34 = t10*t10;
        t36 = t33/t34;
        t43 = t4*th;
        t59 = t43*t27*t29;

        out[i] = t25*t6*t27*t4*t29*t36*t3*t39-t13*t6*t43*t27*t29*t33*t3*t38*t17*t20+
            t13*t3*t38*t17*t33*t20*t6*de*t59+t25*t3*t39*t36*t6*t59;
    }

}


void dbb6(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t4, t5, t12, t16, t32, t38, t39, t40, t47, t50, t61, t90, t6, t7, t8, t9, t10, t11, t13, t14, t35, t36, t37, t42, t48, t53, t56, t57, t59, t78, t80, t87, t93;

    th = param[0];
    de = param[1];

    for(i=0;i<*n;i++)
    {
        t1 = 1.0-u[i];
        t2 = pow(t1,th);
        t3 = 1.0-t2;
        t4 = log(t3);
        t5 = pow(-t4,de);
        t12 = 1/de;
        t16 = 1/th;
        t32 = de-1.0;
        t38 = 2.0*de;
        t39 = -1.0+t38;
        t40 = pow(-t4,t39);
        t47 = 3.0*de-1.0;
        t50 = pow(-t4,t32);
        t61 = pow(-t4,t47);
        t90 = pow(-t4,t38);
        t6 = 1.0-v[i];
        t7 = pow(t6,th);
        t8 = 1.0-t7;
        t9 = log(t8);
        t10 = pow(-t9,de);
        t11 = t5+t10;
        t13 = pow(t11,t12);
        t14 = exp(-t13);
        t35 = pow(t11,-2.0*t32*t12);
        t36 = t35*th;
        t37 = exp(t13);
        t42 = pow(-t9,t39);
        t48 = pow(-t9,t47);
        t53 = t13*de;
        t56 = pow(-t9,t32);
        t57 = t37*t50*t56;
        t59 = t13*th;
        t78 = t37-1.0;
        t80 = pow(t78*t14,t16);
        t87 = t78*t78;
        t93 = pow(-t9,t38);

        out[i] = (2.0*t36*t37*t40*t42+t36*t37*t48*t50+t53*th*t57-t59*t57+
            t36*t37*t61*t56-2.0*t35*t40*t42-t35*t61*t56-t53*th*t50*t56+t59*t50*t56-
            t35*t48*t50) *t80*t7*t2/t3/t8/t87/(t90+2.0*t5*t10+t93)/t1/t6;
    }

}


void dbb7(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t11, t12, t14, t15, t16, t18, t20, t23, t24, t25, t27, t30, t31, t32, t35, t37, t42, t54;

    th = param[0];
    de = param[1];

    for(i=0;i<*n;i++)
    {
        t1 = 1.0-u[i];
        t2 = pow(t1,th);
        t3 = 1.0-t2;
        t4 = pow(t3,-de);
        t5 = 1.0-v[i];
        t6 = pow(t5,th);
        t7 = 1.0-t6;
        t8 = pow(t7,-de);
        t9 = t4+t8-1.0;
        t11 = pow(t9,-1.0/de);
        t12 = 1.0-t11;
        t14 = pow(t12,1.0/th);
        t15 = t11*t11;
        t16 = t14*t15;
        t18 = 1./t5;
        t20 = 1./t7;
        t23 = t9*t9;
        t24 = 1./t23;
        t25 = t12*t12;
        t27 = t24/t25;
        t30 = t2/t1;
        t31 = 1./t3;
        t32 = t30*t31;
        t35 = t14*t11;
        t37 = t6*th;
        t42 = 1./t12;
        t54 = t37*t18*t20;

        out[i] = -t16*t8*t6*t18*t20*t27*t4*t32 + t35*t8*t37*t18*t20*t24*t4*t30*t31*t42+
            t35*t4*t30*t31*t24*t42*t8*de*t54+t16*t4*t32*t27*t8*t54;
    }

}


void dbb8(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t2, t3, t12, t16, t6, t7, t10, t11, t33, t38, t39, t49, t59, t69, t25, t26, t29, t44, t45, t50, t54, t62, t67;

    th = param[0];
    de = param[1];

    for(i=0;i<*n;i++)
    {
        t2 = 1.0-de*u[i];
        t3 = pow(t2,th);
        t10 = 1.0-de;
        t11 = pow(t10,th);
        t12 = 1.0-t11;
        t16 = 1/th;
        t33 = th*t3;
        t38 = 2.0*th;
        t39 = pow(t10,t38);
        t49 = pow(t2,t38);
        t59 = pow(t10,3.0*th);
        t69 = t12*t12;
        t6 = 1.0-de*v[i];
        t7 = pow(t6,th);
        t25 = t3*t7;
        t26 = t11-t7-t3+t25;
        t29 = pow(-t26/t12,t16);
        t44 = pow(t6,t38);
        t45 = t3*t44;
        t50 = t49*t7;
        t54 = t49*t44;
        t62 = -2.0*t25*t11+t25-t33*t7+3.0*t33*t7*t11-3.0*t33*t7*t39+t25*t39+
            2.0* t45*t11-t45*t39+2.0*t50*t11-t50*t39-2.0*t54*t11+t54*t39+t54-
            t50-t45+t33*t7*t59;
        t67 = t26*t26;
        out[i] = -de*t29*t62/t6/t2/t67/t69;
    }

}



void LL_mod2(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    double* negv;
    double* negu;
    negv = (double *) malloc(*n*sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    double ntheta, nnu;
    int nfamily;
    ntheta = -*theta;
    nnu = -*nu;

    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }
    if((*family)==43)
    {
        nfamily=3;
        if(*theta > 0){
            ntheta=2*(*theta)/(1-*theta);
            LL(&nfamily, n, u,  v, &ntheta, nu, loglik);
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
        }
    }else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            LL(&nfamily, n, u,  v, &ntheta, nu, loglik);
        }else{
            ntheta=1/(1+*theta);
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
        }
    }else{

        if((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61))	// 90? rotated copulas
        {
            nfamily = (*family)-20;
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
        }
        else if((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71))	// 270? rotated copulas
        {
            nfamily = (*family)-30;
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            LL(&nfamily, n, u,  negv, &ntheta, &nnu, loglik);
        }
        else if((*family==124) | (*family==224))
        {
            nfamily = (*family)-20;
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            LL(&nfamily, n, v, negu, &ntheta, nu, loglik);
        }
        else if((*family==134) | (*family==234))
        {
            nfamily = (*family)-30;
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            LL(&nfamily, n, negv, u, &ntheta, nu, loglik);
        }
        else {
            LL(family, n, u,  v, theta, nu, loglik);
        }
    }
    free(negv);
    free(negu);
}


//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// loglik    output
//////////////////////////////////////////////////////////////
void LL(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int j;
    double *dat, rho, ll=0.0, t1=0.0, t2=0.0, f;
    //Allocate memory:
    dat = Calloc(2,double);

    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }

    //Compute log-likelihood:
    if(*family==0) //independent
        ll = 0;
    else if(*family==1) //Gaussian
    {
        rho=*theta;
        for(j=0;j<*n;j++)
        {
            dat[0]=u[j]; dat[1]=v[j];
            t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
            f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==2) //Student
    {
        rho=*theta;
        for(j=0;j<*n;j++)
        {
            dat[0] = u[j]; dat[1] = v[j];
            t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
            f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==3) //Clayton
    {
        if(*theta == 0) ll = 0;
        else if(*theta < 1e-10) ll = 0;
        else
        {
            for(j=0;j<*n;j++)
            {
                dat[0] = u[j]; dat[1] = v[j];
                f=log1p(*theta)-(1.0+*theta)*log(dat[0]*dat[1])-(2.0+1.0/(*theta))*log(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0);
                if(f>XINFMAX) ll += log(XINFMAX);
                else if(f<log(DBL_MIN)) ll += log(DBL_MIN);
                else ll += f;
            }
        }
    }
    else if(*family==4) //Gumbel
    {
        for(j=0;j<*n;j++)
        {
            dat[0] = u[j]; dat[1] = v[j];
            t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
            f= -pow(t1,1.0/(*theta))+(2.0/(*theta)-2.0)*log(t1)+(*theta-1.0)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log1p((*theta-1.0)*pow(t1,-1.0/(*theta)));

            if(f>XINFMAX) ll += log(XINFMAX);
            else if(f<log(DBL_MIN))ll += log(DBL_MIN);
            else ll += f;
        }
    }
    else if(*family==5) // Frank
    {
        for(j=0;j<*n;j++)
        {
            if (fabs(*theta) < 1e-10) {
                ll = 0;
            } else {
                dat[0] = u[j]; dat[1] = v[j];
                f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
                if(log(f)>XINFMAX) ll += log(XINFMAX);
                else if(f < DBL_MIN) ll += log(DBL_MIN);
                else ll += log(f);
            }
        }
    }
    else if(*family==6)	//Joe
    {
        for(j=0;j<*n;j++)
        {
            f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==7)	//BB1
    {
        if(*theta == 0){
            for(j=0;j<*n;j++)
            {
                dat[0] = u[j]; dat[1] = v[j];
                t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
                f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
                if(f>XINFMAX) ll += log(XINFMAX);
                else if(f<log(DBL_MIN))ll += log(DBL_MIN);
                else ll += f;
            }
        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            dbb1(u, v, n, param, fuc);
            for(j=0;j<*n;j++)
            {
                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
                else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
                else ll += log(fuc[j]);
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==8)	//BB6
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        dbb6(u, v, n, param, fuc);
        for(j=0;j<*n;j++)
        {
            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
            else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
            else ll += log(fuc[j]);
        }
        Free(fuc); Free(param);
    }
    else if(*family==9)	//BB7
    {
        if(*nu==0)
        {
            for(j=0;j<*n;j++)
            {
                f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
                if(log(f)>XINFMAX) ll += log(XINFMAX);
                else if(f < DBL_MIN) ll += log(DBL_MIN);
                else ll += log(f);
            }
        }
        else
        {
            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            dbb7(u, v, n, param, fuc);
            for(j=0;j<*n;j++)
            {
                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
                else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
                else ll += log(fuc[j]);
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==10) //BB8
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        dbb8(u, v, n, param, fuc);
        for(j=0;j<*n;j++)
        {
            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
            else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
            else ll += log(fuc[j]);
        }
        Free(fuc); Free(param);

    }
    else if(*family==13) //rotated Clayton (180?)
    {
        if(*theta == 0) ll = 0;
        else if(*theta < XEPS) ll = 0;
        else
        {
            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];
                f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
                f = MAX(f,0);
                if(log(f)>XINFMAX) ll += log(XINFMAX);
                else if(f < DBL_MIN) ll += log(DBL_MIN);
                else ll += log(f);
            }
        }
    }
    else if(*family==14) //rotated Gumbel (180?)
    {
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
            t2 = exp(-pow(t1,1.0/(*theta)));
            f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==16) //rotated Joe (180?)
    {
        for(j=0;j<*n;j++)
        {
            u[j]=1-u[j]; v[j]=1-v[j];
            f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
            u[j]=1-u[j]; v[j]=1-v[j];
        }
    }
    else if(*family==17) //rotated BB1
    {
        if(*theta == 0){
            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];
                t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
                f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
                if(f>XINFMAX) ll += log(XINFMAX);
                else if(f<log(DBL_MIN))ll += log(DBL_MIN);
                else ll += f;
            }
        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            int k=1;
            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];

                dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);

                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
                else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
                else ll += log(fuc[j]);
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==18) //rotated BB6
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        int k=1;
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];

            dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);

            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
            else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
            else ll += log(fuc[j]);
        }
        Free(fuc); Free(param);
    }
    else if(*family==19) //rotated BB7
    {
        if(*nu==0){
            for(j=0;j<*n;j++)
            {
                f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
                if(log(f)>XINFMAX) ll += log(XINFMAX);
                else if(f < DBL_MIN) ll += log(DBL_MIN);
                else ll += log(f);
            }
        }else{
            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            int k=1;

            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];
                dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);

                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
                else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
                else ll += log(fuc[j]);
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==20) //rotated BB8
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        int k=1;

        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            dbb8(&dat[0], &dat[1], &k, param, &fuc[j]);

            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
            else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
            else ll += log(fuc[j]);
        }
        Free(fuc); Free(param);
    }
    else if(*family==41)		// New: 1-parametric asymmetric copula (from Harry Joe)
    {
        double tem1, tem2, con, sm, tem;

        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            tem1=qgamma(1.0-dat[0],*theta,1,1,0);
            tem2=qgamma(1.0-dat[1],*theta,1,1,0);
            con=gammafn(1.0+(*theta))/(*theta);
            sm=pow(tem1,*theta)+pow(tem2,*theta);
            tem=pow(sm,(1.0/(*theta)));
            f=con*tem*exp(-tem+tem1+tem2)/sm;

            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==51)		// New: rotated 1-parametric asymmetric copula (from Harry Joe)
    {
        double tem1, tem2, con, sm, tem;

        for(j=0;j<*n;j++)
        {
            tem1=qgamma(1.0-u[j],*theta,1,1,0);
            tem2=qgamma(1.0-v[j],*theta,1,1,0);
            con=gammafn(1.0+(*theta))/(*theta);
            sm=pow(tem1,*theta)+pow(tem2,*theta);
            tem=pow(sm,(1.0/(*theta)));
            f=con*tem*exp(-tem+tem1+tem2)/sm;

            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==104)		//New: Tawn
    {
        int T=1; //length of sample, different from T in TawnPDF
        double par3=1.0;
        for(j=0;j<*n;j++)
        {
            TawnPDF(&u[j], &v[j], &T, theta, nu, &par3, &f);
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==114)		//New: rotated Tawn
    {
        int T=1; // length of sample, different from T in TawnPDF
        double par3=1.0;
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            TawnPDF(&dat[0], &dat[1], &T, theta, nu, &par3, &f);
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==204)		//New: Tawn2
    {
        int T=1; // length of sample, different from T in TawnPDF
        double par2=1.0;
        for(j=0;j<*n;j++)
        {
            TawnPDF(&u[j], &v[j], &T, theta, &par2, nu, &f);
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else if(*family==214)		//New: rotated Tawn2
    {
        int T=1; // length of sample, different from T in TawnPDF
        double par2=1.0;
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            TawnPDF(&dat[0], &dat[1], &T, theta, &par2, nu, &f);
            if(log(f)>XINFMAX) ll += log(XINFMAX);
            else if(f < DBL_MIN) ll += log(DBL_MIN);
            else ll += log(f);
        }
    }
    else
    {
        Rprintf("%d\n",*family);
        printError("Error in LL\t:","Unknown copula family");
    }
    //Free memory:
    Free(dat);
    //Write to output vector:
    *loglik = ll;
}

//////////////////////////////////////////////////////////////
// Function to compute likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// coplik    output
//////////////////////////////////////////////////////////////
void copLik_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
    double* negv;
    double* negu;
    negv = (double *) malloc(*n*sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    double ntheta, nnu;
    int nfamily, i;
    ntheta = -*theta;
    nnu = -*nu;

    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }

    if((*family)==43)
    {
        nfamily=3;
        if(*theta > 0){
            ntheta=2*(*theta)/(1-*theta);
            copLik(&nfamily, n, u,  v, &ntheta, nu, coplik);
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
        }
    }else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            copLik(&nfamily, n, u,  v, &ntheta, nu, coplik);
        }else{
            ntheta=1/(1+*theta);
            for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
        }
    }else{

        if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)) )	// 90? rotated copulas
        {
            nfamily = (*family)-20;
            for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
            copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
        }
        else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)) )	// 270? rotated copulas
        {
            nfamily = (*family)-30;
            for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
            copLik(&nfamily, n, negu,  v, &ntheta, &nnu, coplik);
        }
        else {
            copLik(family, n, u,  v, theta, nu, coplik);
        }
    }
    free(negv);
    free(negu);
}


void copLik(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
    int j;
    double *dat, rho, lik=1.0, t1=0.0, t2=0.0, f;
    //Allocate memory:
    dat = Calloc(2,double);

    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }

    //Compute likelihood:
    if(*family==0) //independent
        lik = 1.0;
    else if(*family==1) //Gaussian
    {
        rho=*theta;
        for(j=0;j<*n;j++)
        {
            dat[0]=u[j]; dat[1]=v[j];
            t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
            f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
            lik *= f;
        }
    }
    else if(*family==2) //Student
    {
        rho=*theta;
        for(j=0;j<*n;j++)
        {
            dat[0] = u[j]; dat[1] = v[j];
            t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
            f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
            lik *= f;
        }
    }
    else if(*family==3) //Clayton
    {
        if(*theta == 0) lik = 1.0;
        if(*theta < XEPS) lik = 1.0;
        else
        {
            for(j=0;j<*n;j++)
            {
                dat[0] = u[j]; dat[1] = v[j];
                f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
                f = MAX(f,0);
                lik *= f;
            }
        }
    }
    else if(*family==4) //Gumbel
    {
        for(j=0;j<*n;j++)
        {
            dat[0] = u[j]; dat[1] = v[j];
            t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
            t2 = exp(-pow(t1,1.0/(*theta)));
            f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
            lik *= f;
        }
    }
    else if(*family==5) // Frank
    {
        for(j=0;j<*n;j++)
        {
            dat[0] = u[j]; dat[1] = v[j];
            f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
            lik *= f;
        }
    }
    else if(*family==6)	//Joe
    {
        for(j=0;j<*n;j++)
        {
            f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
            lik *= f;
        }
    }
    else if(*family==7)	//BB1
    {
        if(*theta == 0){

            for(j=0;j<*n;j++)
            {
                dat[0] = u[j]; dat[1] = v[j];
                t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
                t2 = exp(-pow(t1,1.0/(*nu)));
                f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
                lik *= f;
            }

        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            dbb1(u, v, n, param, fuc);
            for(j=0;j<*n;j++)
            {
                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                lik *= fuc[j];
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==8)	//BB6
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        dbb6(u, v, n, param, fuc);
        for(j=0;j<*n;j++)
        {
            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            lik *= fuc[j];
        }
        Free(fuc); Free(param);
    }
    else if(*family==9)	//BB7
    {
        if(*nu == 0){
            for(j=0;j<*n;j++)
            {
                f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
                lik *= f;
            }
        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            dbb7(u, v, n, param, fuc);
            for(j=0;j<*n;j++)
            {
                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                lik *= fuc[j];
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==10)	//BB8
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        dbb8(u, v, n, param, fuc);
        for(j=0;j<*n;j++)
        {
            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            lik *= fuc[j];
        }
        Free(fuc); Free(param);
    }
    else if(*family==13) //rotated Clayton (180?)
    {
        if(*theta == 0) lik = 1.0;
        else
        {
            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];
                f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
                lik *= f;
            }
        }
    }
    else if(*family==14) //rotated Gumbel (180?)
    {
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];
            t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
            t2 = exp(-pow(t1,1.0/(*theta)));
            f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
            lik *= f;
        }
    }
    else if(*family==16) //rotated Joe (180?)
    {
        for(j=0;j<*n;j++)
        {
            dat[0]=1-u[j]; dat[1]=1-v[j];
            f = pow(pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta),1/(*theta)-2)*pow(1-dat[0],*theta-1)*pow(1-dat[1],*theta-1)*(*theta-1+pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta));
            lik *= f;

        }
    }
    else if(*family==17)	//rotated BB1
    {
        if(*theta == 0){

            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];
                t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
                t2 = exp(-pow(t1,1.0/(*nu)));
                f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
                lik *= f;
            }

        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            int k=1;
            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];

                dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);

                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                lik *= fuc[j];
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==18)	//rotated BB6
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        int k=1;
        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];

            dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);

            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            lik *= fuc[j];
        }
        Free(fuc); Free(param);
    }
    else if(*family==19)	//rotated BB7
    {
        if(*nu == 0){
            for(j=0;j<*n;j++)
            {
                f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
                lik *= f;
            }
        }else{

            double *param, *fuc;
            param=Calloc(2,double);
            param[0]=*theta;
            param[1]=*nu;
            fuc = Calloc(*n,double);
            int k=1;

            for(j=0;j<*n;j++)
            {
                dat[0] = 1-u[j]; dat[1] = 1-v[j];

                dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);

                if(!isfinite(fuc[j]) || isnan(fuc[j]))
                {
                    fuc[j]=1;
                }

                lik *= fuc[j];
            }
            Free(fuc); Free(param);
        }
    }
    else if(*family==20)	//rotated BB8
    {
        double *param, *fuc;
        param=Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        fuc = Calloc(*n,double);
        int k=1;

        for(j=0;j<*n;j++)
        {
            dat[0] = 1-u[j]; dat[1] = 1-v[j];

            dbb8(&dat[0], &dat[1], &k, param, &fuc[j]);

            if(!isfinite(fuc[j]) || isnan(fuc[j]))
            {
                fuc[j]=1;
            }

            lik *= fuc[j];
        }
        Free(fuc); Free(param);
    }

    else printError("Error in copLik\t:","Unknown copula family");
    //Free memory:
    Free(dat);
    //Write to output vector:
    *coplik = lik;
}









//// log likelihood for each observation --------
// unvectorized version
void LL_mod_seperate(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int kk=1;
    for(int i=0; i<(*n); i++){
        LL_mod2(family,&kk,&u[i],&v[i],theta,nu,&loglik[i]);
    };
}

// vectorized version
void LL_mod_seperate_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        LL_mod2(&family[i],&nn,&u[i],&v[i],&theta[i],&nu[i],&loglik[i]);
    };
}

//// density for each observation --------
// unvectorized version
void PDF_seperate(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int kk=1;
    for(int i=0; i<(*n); i++){
        LL_mod2(family,&kk,&u[i],&v[i],theta,nu,&loglik[i]);
        loglik[i] = exp(loglik[i]);
    };
}

// vectorized version
void PDF_seperate_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        LL_mod2(&family[i],&nn,&u[i],&v[i],&theta[i],&nu[i],&loglik[i]);
        loglik[i] = exp(loglik[i]);
    };
}
