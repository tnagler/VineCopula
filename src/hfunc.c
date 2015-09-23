/*
 ** hfunc.c - C code of the package CDRVine  
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
#include "include/hfunc.h"
#include "include/evCopula.h"

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4


// h-func for BB1

void pcondbb1(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t16, t17, t4, t5, t6, t7, t9, t10, t12, t13, t20;
    
    th = param[0];
    de = param[1];
    for(i=0;i<*n;i++)
    {
        t1 = pow(u[i],-th);
        t2 = t1-1.;
        t3 = pow(t2,de);
        t16 = 1./u[i];
        t17 = 1./t2;
        t4 = pow(v[i],-th);
        t5 = t4-1.;
        t6 = pow(t5,de);
        t7 = t3+t6;
        t9 = pow(t7,1/de);
        t10 = 1.0+t9;
        t12 = pow(t10,-1/th);
        t13 = t12*t9;
        t20 = 1./t10;
        out[i] = t13*t3*t1*t16*t17/t7*t20;
    }
    
}



void pcondbb6(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t4, t5, t12, t16, t6, t7, t8, t9, t10, t11, t13, t14, t15, t17;
    
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
        t6 = 1.0-v[i];
        t7 = pow(t6,th);
        t8 = 1.0-t7;
        t9 = log(t8);
        t10 = pow(-t9,de);
        t11 = t5+t10;
        t13 = pow(t11,t12);
        t14 = exp(-t13);
        t15 = 1.0-t14;
        t17 = pow(t15,t16);
        
        out[i] = -t17*t13*t5*t2/t1/t3/t4/t11*t14/t15;
    }
    
}


void pcondbb7(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t1, t2, t3, t4, t6, t8, t9, t11, t12, t14;
    
    th = param[0];
    de = param[1];
    
    for(i=0;i<*n;i++)
    {
        t1 = 1.0-u[i];
        t2 = pow(t1,1.0*th);
        t3 = 1.0-t2;
        t4 = pow(t3,-1.0*de);
        t6 = pow(1.0-v[i],1.0*th);
        t8 = pow(1.0-t6,-1.0*de);
        t9 = t4+t8-1.0;
        t11 = pow(t9,-1.0/de);
        t12 = 1.0-t11;
        t14 = pow(t12,1.0/th);
        
        out[i] = t14*t11*t4*t2/t1/t3/t9/t12;
    }
    
}


void pcondbb8(double* u, double* v, int* n, double* param, double* out)
{
    int i;
    double th, de;
    double t2, t3, t12, t16, t6, t7, t8, t10, t11, t13, t15, t17;
    
    th = param[0];
    de = param[1];
    
    for(i=0;i<*n;i++)
    {
        t2 = 1.0-de*u[i];
        t3 = pow(t2,th);
        t10 = 1.0-de;
        t11 = pow(t10,th);
        t12 = 1.0-t11;
        t13 = 1/t12;
        t16 = 1/th;
        t6 = 1.0-de*v[i];
        t7 = pow(t6,th);
        t8 = 1.0-t7;
        t15 = 1.0-(1.0-t3)*t8*t13;
        t17 = pow(t15,t16);
        
        out[i] = t17*t3/t2*t8*t13/t15;
    }
    
}


// Since the h function is not symmetric in case of double Gumbel and double Clayton we have two implement both separately,
// i.e. Hfunc1 and Hfunc2
void  Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
    double *negv, *negu;
    negv = (double *) malloc(*n* sizeof(double));
    negu = (double *) malloc(*n*sizeof(double));
    double ntheta, nnu;
    int nfamily, j, T=1;
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
            Hfunc(&nfamily, n, u, v, &ntheta, &nnu, out);  
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc(&nfamily, n, u, negv, &ntheta, &nnu, out);
        }
    }else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            Hfunc (&nfamily, n, u, v, &ntheta, &nnu, out);  
        }else{
            ntheta=1/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc (&nfamily, n, u, negv, &ntheta, &nnu, out);
        }
    }else{
        
        if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61) ))
        {
            nfamily=(*family)-20;
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc (&nfamily, n, u, negv, &ntheta, &nnu, out);
        }
        else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71) ))
        {
            nfamily=(*family)-30;
            for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
            Hfunc(&nfamily, n, negu, v, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
        // u and v enter in wrong order from BiCopHfunc and have to be treated accordingly 
        else if(*family==104)
        {
            double par3=1;
            dC_du(v,u,n,theta,nu,&par3,out);		
        }
        else if(*family==114)
        {
            double par3=1;
            for(j=0;j<*n;j++)
            {
                negv[j]= 1-v[j];
                negu[j]= 1-u[j];
                dC_du(&negv[j],&negu[j],&T,theta,nu,&par3,&out[j]);
                out[j]= 1-out[j];
            }
        }
        else if(*family==124)
        {
            double par3=1;
            for(j=0;j<*n;j++)
            {
                negv[j]= 1-v[j];
                dC_du(&negv[j],&u[j],&T,&ntheta,nu,&par3,&out[j]);
            }
        }
        else if(*family==134)
        {
            double par3=1;
            for(j=0;j<*n;j++)
            {
                negu[j]= 1-u[j];
                dC_du(&v[j],&negu[j],&T,&ntheta,nu,&par3,&out[j]);
                out[j]=1-out[j];
                
            }
        }
        else if(*family==204)
        {
            double par3=*nu;
            double par2=1;
            dC_du(v,u,n,theta,&par2,&par3,out);
        }
        else if(*family==214)
        {
            double par3=*nu;
            double par2=1;
            for(j=0;j<*n;j++)
            {
                negv[j]= 1-v[j];
                negu[j]= 1-u[j];
                dC_du(&negv[j],&negu[j],&T,theta,&par2,&par3,&out[j]);
                out[j]= 1-out[j];
            }
        }
        else if(*family==224)
        {
            double par3=*nu;
            double par2=1;
            for(j=0;j<*n;j++)
            {
                negv[j]= 1-v[j];
                dC_du(&negv[j],&u[j],&T,&ntheta,&par2,&par3,&out[j]);
            }
        }
        else if(*family==234)
        {
            double par3=*nu;
            double par2=1;
            for(j=0;j<*n;j++)
            {
                negu[j]= 1-u[j];
                dC_du(&v[j],&negu[j],&T,&ntheta,&par2,&par3,&out[j]);
                out[j]=1-out[j];
            }
        }
        else {
            Hfunc (family, n, u, v, theta, nu, out); 
        }
    }
    // ensure that results are in [0,1]
    for(int j=0; j <* n; ++j){out[j] = MIN(MAX(out[j], 0), 1);}
    free(negv);
    free(negu);
}

void  Hfunc2(int* family,int* n,double* v,double* u,double* theta,double* nu,double* out)
{
    double *negv, *negu;
    negv = (double *) malloc(*n * sizeof(double));
    negu = (double *) malloc(*n * sizeof(double));
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
            Hfunc (&nfamily, n, v, u, &ntheta, &nnu, out);  
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
    }else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            Hfunc (&nfamily, n, v, u, &ntheta, &nnu, out);  
        }else{
            ntheta=1/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};		}
    }else{
        
        if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61) ))
        {
            nfamily=(*family)-20;
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
        else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71) ))
        {
            nfamily=(*family)-30;
            for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
            Hfunc(&nfamily, n, v, negu, &ntheta, &nnu, out);
        }
        // else if(*family==104 | *family==204 | *family==114 | *family==214)
        // {
        // u und v vertauschen (Unsauber, aber so sollte es funktionieren in unserer bisherigen Notation)
        // Hfunc(family,n,u,v,theta,nu,out);
        // }
        else if((*family==124) | (*family==224))
        {
            nfamily=(*family)-20;
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hfunc(&nfamily, n, negv, u, &ntheta, nu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
        else if((*family==134) | (*family==234))
        {
            nfamily=(*family)-30;
            for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
            Hfunc(&nfamily, n, v, negu, &ntheta, nu, out);
        }
        else
        { 
            Hfunc(family, n, v, u, theta, nu, out);
        }
    }
    // ensure that results are in [0,1]
    for(int i=0; i < *n; ++i) {out[i] = MIN(MAX(out[i], 0), 1);}
    free(negv);
    free(negu);
}

// vectorized versions
void Hfunc1_vec(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        Hfunc1(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
    };
}

void Hfunc2_vec(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        Hfunc2(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
    };
}


//////////////////////////////////////////////////////////////
// Function to compute h-function for vine simulation and estimation
// Input:
// family   copula family (0=independent,  1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
// out      output
//////////////////////////////////////////////////////////////
void Hfunc(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
    int j;
    double *h;
    h = Calloc(*n,double);
    double x;
    

    for(j=0;j<*n;j++)
    {
        if((v[j]==0) | ( u[j]==0)) h[j] = 0;
        else if (v[j]==1) h[j] = u[j];
        else
        {
            if(*family==0) //independent
            {
                h[j] = u[j]; 
            }
            else if(*family==1) //gaussian
            {
                x = (qnorm(u[j],0.0,1.0,1,0) - *theta*qnorm(v[j],0.0,1.0,1,0))/sqrt(1.0-pow(*theta,2.0));
                if (isfinite(x))
                    h[j] = pnorm(x,0.0,1.0,1,0);
                else if ((qnorm(u[j],0.0,1.0,1,0) - *theta*qnorm(v[j],0.0,1.0,1,0)) < 0)
                    h[j] = 0;
                else 
                    h[j] = 1;
            }
            else if(*family==2) //student
            {
                double t1, t2, mu, sigma2;
                t1 = qt(u[j],*nu,1,0); t2 = qt(v[j],*nu,1,0); mu = *theta*t2; sigma2 = ((*nu+t2*t2)*(1.0-*theta*(*theta)))/(*nu+1.0);
                h[j] = pt((t1-mu)/sqrt(sigma2),*nu+1.0,1,0);
            }
            else if(*family==3) //clayton
            { 
                if(*theta == 0) h[j] = u[j] ;
                if(*theta < XEPS) h[j] = u[j] ;
                else
                { 
                    x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
                    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta));
                    if(*theta < 0)
                    {
                        if(x < 0) h[j] = 0; 
                    }
                }
            }
            else if(*family==4) //gumbel
            {
                if(*theta == 1) h[j] = u[j] ; 
                else
                {
                    h[j] = -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*log(v[j]));
                }
            }
            else if(*family==5) //frank
            {
                if(*theta==0) h[j]=u[j];
                else
                {
                    h[j] = -(exp(*theta)*(exp(*theta*u[j])-1.0))/(exp(*theta*v[j]+*theta*u[j])-exp(*theta*v[j]+*theta)-exp(*theta*u[j]+*theta)+exp(*theta));
                }
            }
            else if(*family==6) //joe
            {
                if(*theta==1) h[j]=u[j];
                else
                {
                    h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
                }
            }
            else if(*family==7)	//BB1
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*nu==1) 
                {
                    if(*theta==0) h[j]=u[j];
                    else h[j]=pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1,-1/(*theta)-1)*pow(v[j],-*theta-1);
                }
                else if(*theta==0) 
                {
                    h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
                }
                else
                {
                    pcondbb1(&v[j],&u[j],&T,param,&h[j]);
                }
                Free(param);
            }
            else if(*family==8) //BB6
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*theta==1) 
                {
                    if(*nu==1) h[j]=u[j];
                    else h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
                }     
                else if(*nu==1) 
                {
                    h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
                }
                else
                {
                    pcondbb6(&v[j],&u[j],&T,param,&h[j]);
                }
                Free(param);			
            }
            else if(*family==9)	//BB7
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*theta==1)
                {
                    if(*nu==0) h[j]=u[j];
                    else h[j]=pow(pow(u[j],-*nu)+pow(v[j],-*nu)-1,-1/(*nu)-1)*pow(v[j],-*nu-1);
                }
                else if(*nu==0)
                {
                    h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
                }
                else
                {
                    pcondbb7(&v[j],&u[j],&T,param,&h[j]);
                }
                Free(param);
            }
            else if(*family==10) //BB8
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*nu==0)
                {
                    h[j]=u[j];
                }			
                else if(*nu==1)
                {
                    if(*theta==1) h[j]=u[j];
                    else h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
                }
                else
                {
                    pcondbb8(&v[j],&u[j],&T,param,&h[j]);
                }
                Free(param);		  
            }
            else if(*family==13) //rotated clayton (180?)
            {
                if(*theta == 0) h[j] = u[j] ;
                if(*theta < XEPS) h[j] = u[j] ;
                else
                {
                    u[j]=1-u[j]; 
                    v[j]=1-v[j];
                    x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
                    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta)); // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
                    h[j]= 1-h[j];
                    u[j]=1-u[j];
                    v[j]=1-v[j];
                }
            }
            else if(*family==14) //rotated gumbel (180?)
            {
                v[j]= 1-v[j];
                u[j]= 1-u[j];
                h[j]= -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*	log(v[j]));
                h[j]= 1-h[j];
                u[j]=1-u[j];
                v[j]=1-v[j];
            }
            else if(*family==16)
            {
                v[j]= 1-v[j];
                u[j]= 1-u[j];
                h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
                h[j]= 1-h[j];
                u[j]=1-u[j];
                v[j]=1-v[j];
            }
            else if(*family==17) //rotated BB1
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*nu==1) 
                {
                    if(*theta==0) h[j]=u[j];
                    else
                    { 
                        h[j]=pow(pow(1-u[j],-*theta)+pow(1-v[j],-*theta)-1,-1/(*theta)-1)*pow(1-v[j],-*theta-1);
                        h[j]= 1-h[j];
                    }
                }
                else if(*theta==0) 
                {
                    h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
                    h[j]= 1-h[j];
                }
                else
                {
                    v[j]= 1-v[j];
                    u[j]= 1-u[j];
                    pcondbb1(&v[j],&u[j],&T,param,&h[j]);
                    u[j]=1-u[j];
                    v[j]=1-v[j];
                    h[j]= 1-h[j];
                }
                Free(param);
            }
            else if(*family==18) //rotated BB6
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*theta==1) 
                {
                    if(*nu==1) h[j]=u[j];
                    else
                    {
                        h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
                        h[j]= 1-h[j];
                    }
                }     
                else if(*nu==1) 
                {
                    h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
                    h[j]= 1-h[j];
                }
                else
                {
                    v[j]= 1-v[j];
                    u[j]= 1-u[j];
                    pcondbb6(&v[j],&u[j],&T,param,&h[j]);
                    u[j]=1-u[j];
                    v[j]=1-v[j];
                    h[j]= 1-h[j];	
                }
                Free(param);	  		  
            }
            else if(*family==19) //rotated BB7
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*theta==1)
                {
                    if(*nu==0) h[j]=u[j];
                    else{
                        h[j]=pow(pow(1-u[j],-*nu)+pow(1-v[j],-*nu)-1,-1/(*nu)-1)*pow(1-v[j],-*nu-1);
                        h[j]= 1-h[j];
                    }
                }
                else if(*nu==0)
                {
                    h[j] = pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
                    h[j]= 1-h[j];
                }				
                else
                {
                    v[j]= 1-v[j];
                    u[j]= 1-u[j];
                    pcondbb7(&v[j],&u[j],&T,param,&h[j]);
                    u[j]=1-u[j];
                    v[j]=1-v[j];
                    h[j]= 1-h[j];
                }
                Free(param);
            }
            else if(*family==20) //rotated BB8
            {
                double* param;
                param = Calloc(2,double);
                param[0]=*theta;
                param[1]=*nu;
                int T=1;
                if(*nu==0)
                {
                    h[j]=u[j];
                }			
                else if(*nu==1)
                {
                    if(*theta==1) h[j]=u[j];
                    else{
                        h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
                        h[j]= 1-h[j];
                    }
                }
                else
                {
                    v[j]= 1-v[j];
                    u[j]= 1-u[j];
                    pcondbb8(&v[j],&u[j],&T,param,&h[j]);
                    u[j]=1-u[j];
                    v[j]=1-v[j];
                    h[j]= 1-h[j];
                }
                Free(param);		  
            }
            else if(*family==41)
            {
                double t1,t2,t3;
                t1=qgamma(1.0-u[j],*theta,1,1,0);
                t2=qgamma(1.0-v[j],*theta,1,1,0);
                t3=pow(pow(t1,*theta)+pow(t2,*theta),(1.0/(*theta)));
                h[j]=exp(-t3+t1);
            }
            else if(*family==104)
            {
                int T=1;
                double par3=1;
                dC_dv(&u[j],&v[j],&T,theta,nu,&par3,&h[j]);
            }
            else if(*family==114)
            {
                int T=1;
                double par3=1;
                v[j]= 1-v[j];
                u[j]= 1-u[j];
                dC_dv(&u[j],&v[j],&T,theta,nu,&par3,&h[j]);
                u[j]=1-u[j];
                v[j]=1-v[j];
                h[j]= 1-h[j];
            }
            else if(*family==204)
            {
                int T=1;
                double par3=*nu, par2=1;
                dC_dv(&u[j],&v[j],&T,theta,&par2,&par3,&h[j]);	
            }
            else if(*family==214)
            {
                int T=1;
                double par3=*nu, par2=1;
                v[j]= 1-v[j];
                u[j]= 1-u[j];
                dC_dv(&u[j],&v[j],&T,theta,&par2,&par3,&h[j]);
                u[j]=1-u[j];
                v[j]=1-v[j];
                h[j]= 1-h[j];
            }
        }	
        out[j] = MAX(MIN(h[j],UMAX),UMIN);
    }	
    Free(h);
}

///////////////////////////////////////////////////////////////
void qcondgum(double* q, double* u, double* de, double* out)
{
    double a,p,g,gp,z1,z2,con,de1,dif;
    double mxdif;
    int iter;
    
    p = 1-*q;
    z1 = -log(*u);
    con=log(1.-p)-z1+(1.-*de)*log(z1); de1=*de-1.;
    a=pow(2.*pow(z1,*de),1./(*de));
    mxdif=1; iter=0;
    dif=.1;  // needed in case first step leads to NaN
    while(mxdif>1.e-6 && iter<20)
    { g=a+de1*log(a)+con;
        gp=1.+de1/a;
        if(isnan(g) || isnan(gp) || isnan(g/gp) ) { dif/=-2.; }  // added for de>50
        else dif=g/gp;
        a-=dif; iter++;
        while(a<=z1) { dif/=2.; a+=dif; }
        mxdif=fabs(dif);
    }
    z2=pow(pow(a,*de)-pow(z1,*de),1./(*de));
    *out = exp(-z2);
}

void qcondjoe(double* q, double* u, double* de, double* out)
{ double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t19,t23,t28,t31;
    double c21,pdf;
    int iter;
    double diff,v,de1,dtem,de1inv,tem;
    
    t1 = 1.0-*u;
    t2 = pow(t1,1.0*(*de));
    t7 = 1./(*de);
    t10 = t2*(*de);
    t11 = 1./t1;
    t19 = (*de)*(*de);
    de1=*de-1;  // may need better modification for large delta
    dtem=-de1/(1.+de1); de1inv=-1./de1;
    
    // v = 0.5 * (q+u); // starting guess
    
    // Use a better starting point based on reflected B4 copula
    // A good starting point is crucial when delta is large because
    //    C_{2|1} will be steep
    // C_{R,2|1}(v|u)=1-C_{2|1}(1-v|1-u),
    // C_{R,2|1}^{-1}(q|u)=1-C_{2|1}^{-1}(1-q|1-u)
    tem=pow(1.-*q,dtem)-1.;
    tem=tem*pow(1.-*u,-de1)+1.;
    v=pow(tem,de1inv); v=1.-v;
    diff=1; iter=0;
    while(fabs(diff)>1.e-6 && iter<20)
    { t3 = 1.-v;
        t4 = pow(t3,*de);
        t5 = t2*t4;
        t6 = t2+t4-t5;
        t8 = pow(t6,t7);
        t9 = t7*t8;
        t13 = t11*t4;
        t15 = -t10*t11+t10*t13;
        t16 = 1./t6;
        t23 = 1./t3;
        t28 = t6*t6;
        t31 = (-t4*(*de)*t23+t5*(*de)*t23)/t28*t15;
        c21 = -t9*t15*t16;
        pdf = -t8/t19*t31+t8*(*de)*t2*t13*t23*t16+t9*t31;
        iter++;
        if(isnan(pdf) || isnan(c21) ) { diff/=-2.; }  // added for de>=30
        else diff=(c21-*q)/pdf;
        v-=diff;
        while(v<=0 || v>=1 || fabs(diff)>0.25 ) { diff/=2.; v+=diff; }
    }
    *out = v;
}




///////////////////////////////////////////////////////////////
// Function to compute inversion of H numerically through Bisection
// Input:
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// out      output
//////////////////////////////////////////////////////////////
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out)
{
    
    int br=0, in=1;
    double ans=0.0, tol=UMIN, x0=UMIN, x1=UMAX, fl, fh, val;
    //Rprintf("family in HNumInv: %d\n", *family);
    Hfunc1(family,&in,&x0,v,theta,nu,&fl); fl -= *u; 
    Hfunc1(family,&in,&x1,v,theta,nu,&fh); fh -= *u;
    if(fabs(fl)<=tol) { ans=x0; br=1; }
    if(fabs(fh)<=tol) { ans=x1; br=1; }
    
    while(!br){
        
        ans = (x0+x1)/2.0;
        Hfunc1(family,&in,&ans,v,theta,nu,&val);
        val -= *u;
        if(fabs(val)<=tol) br=1;
        if(fabs(x0-x1)<=1e-10) br=1; //stop if values become too close (avoid infinite loop)
        
        if(val > 0.0) {x1 = ans; fh = val;}
        else {x0 = ans; fl = val;}
        
    }
    *out = ans;
}  

/////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
/////////////////////////////////////////////
void Hinv1(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
    double *negv, *negu;
    negv = (double*) Calloc(*n,double);
    negu = (double*) Calloc(*n,double);
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
            Hinv(&nfamily, n, u, v, &ntheta, &nnu, out);
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hinv(&nfamily, n, u, negv, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
    }
    else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            Hinv(&nfamily, n, u, v, &ntheta, &nnu, out);
        }else{
            ntheta=1/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hinv(&nfamily, n, u, negv, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
    }
    else if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30) | (*family==61) ))
    {
        nfamily=(*family)-20;
        for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
        Hinv(&nfamily,  n,  u,  negv,  &ntheta,  &nnu,  out);
    }
    else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40) | (*family==71) ))
    {
        nfamily=(*family)-30;
        for (int i = 0; i < *n; i++) {negu[i]=1 - u[i];};
        Hinv(&nfamily,  n,  negu,  v,  &ntheta,  &nnu,  out);
        for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    else if((*family==124) | (*family==224))
    {
        nfamily=(*family)-20;
        for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
        Hinv(&nfamily,  n,  u,  negv,  &ntheta,  nu,  out);
    }
    else if((*family==134) | (*family==234))
    {
        nfamily=(*family)-30;
        for (int i = 0; i < *n; i++) {negu[i]=1 - u[i];};
        Hinv(&nfamily,  n,  negu,  v,  &ntheta,  nu,  out);
        for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    else {
        Hinv( family,  n,  u,  v,  theta,  nu,  out);
    }
    Free(negv);
    Free(negu);
}

void Hinv2(int* family, int* n, double* v, double* u, double* theta, double* nu, double* out)
{
    double *negv, *negu;
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
            Hinv(&nfamily, n, v, u, &ntheta, &nnu, out);
        }else{
            ntheta=-2*(*theta)/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hinv(&nfamily, n, negv, u, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
    }
    else if((*family)==44)
    {
        nfamily=4;
        if(*theta > 0){
            ntheta=1/(1-*theta);
            Hinv(&nfamily, n, v, u, &ntheta, &nnu, out);
        }else{
            ntheta=1/(1+*theta);
            for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
            Hinv(&nfamily, n, negv, u, &ntheta, &nnu, out);
            for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
        }
    }
    else if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30) | (*family==61) ))
    {
        nfamily = (*family)-20;
        for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
        Hinv(&nfamily,  n,  negv, u,  &ntheta,  &nnu,  out);
        for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40) | (*family==71) ))
    {
        nfamily=(*family)-30;
        for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
        Hinv(&nfamily,  n,  v,  negu,  &ntheta,  &nnu,  out);
    }
    else if((*family==124) | (*family==224))
    {
        nfamily = (*family)-20;
        if ((*family)/100 == 1) {
            nfamily = (nfamily) + 100;
        } else if ((*family)/100 == 2) {
            nfamily = nfamily - 100;
        }
        for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
        Hinv(&nfamily,  n,  negv, u,  &ntheta,  nu,  out);
        for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    else if((*family==134) | (*family==234))
    {
        nfamily = (*family)-30;
        if ((*family)/100 == 1) {
            nfamily = (nfamily) + 100;
        } else if ((*family)/100 == 2) {
            nfamily = nfamily - 100;
        }
        for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
        Hinv(&nfamily,  n,  v,  negu,  &ntheta,  nu,  out);
    }
    else {
        nfamily = (*family);
        if ((*family)/100 == 1) {
            nfamily = nfamily + 100;
        } else if ((*family)/100 == 2) {
            nfamily = nfamily - 100;
        }  
        Hinv(&nfamily,  n,  v,  u,  theta,  nu,  out);
    }
    free(negv);
    free(negu);
}

// vectorized versions
void Hinv1_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        Hinv1(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
    };
}

void Hinv2_vec(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
    int nn=1;
    for(int i=0; i<(*n); i++){
        Hinv2(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
    };
}




//////////////////////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
// Input:
// family   copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
//////////////////////////////////////////////////////////////
void Hinv(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
    int j;
    double *hinv;
    hinv = Calloc(*n,double);
    
    for(int i=0;i<*n;i++)
    {
        if(u[i]<UMIN) u[i]=UMIN;
        else if(u[i]>UMAX) u[i]=UMAX;
        if(v[i]<UMIN) v[i]=UMIN;
        else if(v[i]>UMAX) v[i]=UMAX;
    }
    
    for(j=0;j<*n;j++)
    {
        if(*family==0)
        {
            hinv[j]=u[j];
        }
        else if(*family==1) //gaussian
        {
            hinv[j] = pnorm(qnorm(u[j],0.0,1.0,1,0)*sqrt(1.0-pow(*theta,2.0))+*theta*qnorm(v[j],0.0,1.0,1,0),0.0,1.0,1,0);
        }
        else if(*family==2) //student
        {
            double temp1, temp2, mu, var;
            temp1 = qt(u[j],*nu+1.0,1,0); temp2 = qt(v[j],*nu,1,0); mu = *theta*temp2; var=((*nu+(temp2*temp2))*(1.0-(*theta*(*theta))))/(*nu+1.0);
            hinv[j] = pt((sqrt(var)*temp1)+mu,*nu,1,0);
        }
        else if(*family==3) //clayton
        {
            if(*theta < XEPS) hinv[j]=u[j];
            else
                hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
        }
        else if(*family==4) //gumbel - must turn to numerical inversion
        {
            qcondgum(&u[j],&v[j],theta,&hinv[j]);
        }
        else if(*family==5) //frank - numerical inversion
        {
            hinv[j] = -1/(*theta)*log(1-(1-exp(-*theta)) / ((1/u[j]-1)*exp(-*theta*v[j])+1));
        }
        else if(*family==6) //joe - numerical inversion
        {
            if(*theta<40)
            {
                qcondjoe(&u[j],&v[j],theta,&hinv[j]);
            }
            else
            {
                double nu=0.0;
                HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
            }
        }
        else if(*family==7) //BB1
        {
            HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
        }
        else if(*family==8) //BB6
        {
            HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
        }
        else if(*family==9) //BB7
        {
            HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
        }
        else if(*family==10) //BB8
        {
            HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
        }
        else if(*family==13)
        {
            u[j]=1-u[j];
            v[j]=1-v[j];
            hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==14) //rotated gumbel (180?) - must turn to numerical inversion
        {
            u[j]=1-u[j];
            v[j]=1-v[j];
            qcondgum(&u[j],&v[j],theta,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==16) //rotated joe (180?) - must turn to numerical inversion
        {
            u[j]=1-u[j];
            v[j]=1-v[j];			
            if(*theta<40)
            {
                qcondjoe(&u[j],&v[j],theta,&hinv[j]);
            }
            else
            {
                int jj=6;
                double nu=0.0;
                HNumInv(&jj,&u[j],&v[j],theta,&nu,&hinv[j]);
            }
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==17) //rotated BB1 (180?) - must turn to numerical inversion
        {
            int jj=7;
            u[j]=1-u[j];
            v[j]=1-v[j];
            HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==18) //rotated BB6 (180?) - must turn to numerical inversion
        {
            int jj=8;
            u[j]=1-u[j];
            v[j]=1-v[j];
            HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==19) //rotated BB7 (180?) - must turn to numerical inversion
        {
            int jj=9;
            u[j]=1-u[j];
            v[j]=1-v[j];
            HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==20) //rotated BB8 (180?) - must turn to numerical inversion
        {
            int jj=10;
            u[j]=1-u[j];
            v[j]=1-v[j];
            HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==41)	// 1-parametric asymmetric copula (Harry Joe)
        {
            double de=*theta;
            double tem1,y;
            tem1=qgamma(1.0-v[j],de,1,1,0);
            y=pow(tem1-log(u[j]),de) - pow(tem1,de);
            y=pow(y,(1.0/de));
            hinv[j]=1.0-pgamma(y,de,1,1,0);
        }
        else if(*family==51)	// rotated (180) 1-parametric asymmetric copula (Harry Joe)
        {
            double de=*theta;
            double tem1,y;
            u[j]=1-u[j];
            v[j]=1-v[j];
            tem1=qgamma(1.0-v[j],de,1,1,0);
            y=pow(tem1-log(u[j]),de) - pow(tem1,de);
            y=pow(y,(1.0/de));
            hinv[j]=1.0-pgamma(y,de,1,1,0);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        else if(*family==104 || *family==204) //Tawn
        {
            HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
        }
        else if(*family==114 || *family==214) //Tawn
        {
            int jj=*family-10;
            u[j]=1-u[j];
            v[j]=1-v[j];
            HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
            hinv[j]=1-hinv[j];
            u[j]=1-u[j];
            v[j]=1-v[j];
        }
        
        out[j] = MAX(MIN(hinv[j],UMAX),UMIN); 
    }
    Free(hinv);
}

