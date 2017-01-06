// Implementierung der incomplete beta Funktion und ihrer Ableitungen nach p

// The incomplete beta function is needed for the derivative of the Student's t-copula
// Also Its derivative with respect to the parameter p is needed
// For reference see:  Boik and Robinson-Cox (1998).
 
//  Boik, R. J. and J. F. Robinson-Cox (1998).
// Derivatives of the incomplete beta function.
// Journal of Statistical Software 3(1).

// The implementation follows directly their algorithm and is closely related to their published code

#include "vine.h"
#include "incompleteBeta.h"


void incompleBeta_an1_bn1_p(double* x, double p, double q, double* an, double* bn)
{ 
      double t2 = 1.0/(1.0-(*x));
      double t3 = (*x)*t2;
      double t4 = q-1.0;
      double t5 = p+1.0;
      double t9 = t5*t5;
      double t19 = q*(*x)*t2;
      double t20 = 2.0*t19;
      double t21 = 4.0*q;
      double t27 = p*q;
      double t28 = p-2.0-t19;
      double t31 = 1/q;
      double t32 = (t20+t21+2.0*(t19+2.0*q)*(p-1.0)+t27*t28)*t31;
      double t33 = 1.0/p;
      double t34 = p+2.0;
      double t35 = 1/t34;
      double t36 = t33*t35;
      double t40 = (t20+t21+q*t28+t27)*t31;
      double t42 = p*p;
      double t43 = 1/t42;
      double t44 = t43*t35;
      double t46 = t34*t34;
      double t47 = 1/t46;
      double t48 = t33*t47;
      an[0] = t3*t4/t5;
      an[1] = -t3*t4/t9;
      an[2] = 2.0*t3*t4/t9/t5;
      bn[0] = t32*t36;
      bn[1] = t40*t36-t32*t44-t32*t48;
      bn[2] = 2.0*t36-2.0*t40*t44-2.0*t40*t48+2.0*t32/t42/p*t35+2.0*t32*t43*
		  t47+2.0*t32*t33/t46/t34;

}


void incompleBeta_an1_bn1_q(double* x, double p, double q, double* an, double* bn)
{ 
      double t2 = 1.0/(1.0-*x);
      double t3 = *x*t2;
      double t6 = 1.0/(p+1.0);
      double t11 = q*(*x)*t2;
      double t16 = p-1.0;
      double t19 = p*q;
      double t20 = p-2.0-t11;
      double t22 = 2.0*t11+4.0*q+2.0*(t11+2.0*q)*t16+t19*t20;
      double t23 = 1.0/q;
      double t27 = 1.0/(p+2.0);
      double t28 = 1.0/p*t27;
      double t36 = 2.0*t3+4.0+2.0*(t3+2.0)*t16+p*t20-t19*t3;
      double t39 = q*q;
      double t40 = 1.0/t39;
      an[0] = t3*(q-1.0)*t6;
      an[1] = t3*t6;
      an[2] = 0.0;
      bn[0] = t22*t23*t28;
      bn[1] = t36*t23*t28-t22*t40*t28;
      bn[2] = -2.0*t3*t23*t27-2.0*t36*t40*t28+2.0*t22/t39/q*t28;

}

void incompleBeta_an_bn_p(double* x, double p, double q, int n, double* an, double* bn)
{ 
      double t1 = *x*(*x);
      double t2 = 1.0-*x;
      double t3 = t2*t2;
      double t5 = t1/t3;
      double t6 = n-1.0;
      double t9 = t5*t6*(p+q+n-2.0);
      double t10 = p+n-1.0;
      double t11 = q-n;
      double t12 = t10*t11;
      double t13 = 2.0*n;
      double t14 = p+t13-3.0;
      double t15 = 1/t14;
      double t16 = p+t13-2.0;
      double t17 = t16*t16;
      double t18 = 1/t17;
      double t19 = t15*t18;
      double t20 = p+t13-1.0;
      double t21 = 1/t20;
      double t26 = t5*t6*t10;
      double t27 = t11*t15;
      double t28 = t18*t21;
      double t29 = t27*t28;
      double t32 = t14*t14;
      double t33 = 1/t32;
      double t34 = t33*t18;
      double t39 = 1/t17/t16;
      double t40 = t15*t39;
      double t45 = t20*t20;
      double t46 = 1/t45;
      double t55 = t11*t33*t28;
      double t59 = t27*t39*t21;
      double t63 = t27*t18*t46;
      double t88 = t17*t17;
      double t105 = 2.0*t5*t6*t29-2.0*t26*t55-4.0*t26*t59-2.0*t26*t63-2.0*t9*t55-4.0*
t9*t59-2.0*t9*t63+2.0*t9*t12/t32/t14*t18*t21+4.0*t9*t12*t33*t39*t21+2.0*t9*t12*
t34*t46+6.0*t9*t12*t15/t88*t21+4.0*t9*t12*t40*t46+2.0*t9*t12*t19/t45/t20;
      double t108 = q*(*x)/t2;
      double t110 = t108+2.0*q;
      double t111 = n*n;
      double t118 = p*q;
      double t119 = p-2.0-t108;
      double t122 = 1/q;
      double t123 = (2.0*t110*t111+2.0*t110*(p-1.0)*n+t118*t119)*t122;
      double t124 = 1/t16;
      double t125 = p+t13;
      double t126 = 1/t125;
      double t127 = t124*t126;
      double t133 = (2.0*t110*n+q*t119+t118)*t122;
      double t135 = t18*t126;
      double t137 = t125*t125;
      double t138 = 1/t137;
      double t139 = t124*t138;
      an[0] = t9*t12*t19*t21;
      an[1] = t26*t29+t9*t29-t9*t12*t34*t21-2.0*t9*t12*t40*t21-t9*t12*t19*t46;
      an[2] = t105;
      bn[0] = t123*t127;
      bn[1] = t133*t127-t123*t135-t123*t139;
      bn[2] = 2.0*t127-2.0*t133*t135-2.0*t133*t139+2.0*t123*t39*t126+2.0*
		  t123*t18*t138+2.0*t123*t124/t137/t125;

}


void incompleBeta_an_bn_q(double* x, double p, double q, int n, double* an, double* bn)
{
	double t1 = *x*(*x);
      double t2 = 1.0-*x;
      double t3 = t2*t2;
      double t5 = t1/t3;
      double t6 = n-1.0;
      double t9 = t5*t6*(p+q+n-2.0);
      double t10 = p+n-1.0;
      double t11 = q-n;
      double t13 = 2.0*n;
      double t15 = 1/(p+t13-3.0);
      double t16 = p+t13-2.0;
      double t17 = t16*t16;
      double t18 = 1/t17;
      double t21 = 1/(p+t13-1.0);
      double t28 = t18*t21;
      double t32 = t10*t15*t28;
      double t39 = 1/t2;
      double t40 = q*(*x)*t39;
      double t42 = t40+2.0*q;
      double t43 = n*n;
      double t46 = p-1.0;
      double t50 = p*q;
      double t51 = p-2.0-t40;
      double t53 = 2.0*t42*t43+2.0*t42*t46*n+t50*t51;
      double t54 = 1/q;
      double t56 = 1/t16;
      double t58 = 1/(p+t13);
      double t59 = t56*t58;
      double t61 = *x*t39;
      double t62 = t61+2.0;
      double t70 = 2.0*t62*t43+2.0*t62*t46*n+p*t51-t50*t61;
      double t73 = q*q;
      double t74 = 1/t73;
      an[0] = t9*t10*t11*t15*t18*t21;
      an[1] = t5*t6*t10*t11*t15*t28+t9*t32;
      an[2] = 2.0*t5*t6*t32;
      bn[0] = t53*t54*t59;
      bn[1] = t70*t54*t59-t53*t74*t59;
      bn[2] = -2.0*p*(*x)*t39*t54*t56*t58-2.0*t70*t74*t59+2.0*t53/t73/q*t59;

}


///////////////////////////////////////////////////////////////
// incomplete Beta function and the first and second derivative wrt p
// 
// Input:
// x_in 	value at which the beta function is integrated
// p,q 		beta shape parameters
//
// Output:
// der		I (incompete beta function), Ip, Ipp
/////////////////////////////////////////////////////////////////


void inbeder(double* x_in, double* p_in, double* q_in, double* der)
{

  double lbet, pa, pa1, pb, pb1, pab, pab1, err=1e-12;
  double p, q, x;
  int minappx=3, maxappx=200, n=0;

  // falls x>p/(p+q)
  if (*x_in>*p_in/(*p_in+*q_in))
  {
	  x=1-*x_in;
	  p=*q_in;
	  q=*p_in;
  }
  else
  {
	  x=*x_in;
	  p=*p_in;
	  q=*q_in;
  }
  
  // Compute Log Beta, digamma, and trigamma functions
  
  lbet=lbeta(p,q);
  pa=digamma(p);
  pa1=trigamma(p);
  pb=digamma(q);
  pb1=trigamma(q);
  pab=digamma(p+q);
  pab1=trigamma(p+q);


  double omx=1-x;
  double logx=log(x);
  double logomx=log(omx);

  // Compute derivatives of K(x,p,q)=x^p(1-x)^(q-1)/[p beta(p,q)

  double *c;
  double c0, d;
  c=Calloc(3,double);
  c[0]=p*logx+(q-1)*logomx-lbet-log(p);
  c0=exp(c[0]);
  if (*x_in>*p_in/(*p_in+*q_in))
  {
	 c[1]=logomx-pb+pab;
  	c[2]=c[1]*c[1]-pb1+pab1;
  }
  else
  {
  	c[1]=logx-1/p-pa+pab;
 	 c[2]=c[1]*c[1]+1/p/p-pa1+pab1;
  }
  

  int del=1, i=0;
  double *an, *bn, *an1, *an2, *bn1, *bn2, *dr;
  an=Calloc(3,double);
  bn=Calloc(3,double);
  an1=Calloc(3,double);
  bn1=Calloc(3,double);
  an2=Calloc(3,double);
  bn2=Calloc(3,double);
  dr=Calloc(3,double);
  double *dan, *dbn, *der_old, *d1;
  dan=Calloc(3,double);
  dbn=Calloc(3,double);
  der_old=Calloc(3,double);
  d1=Calloc(3,double);

  double Rn=0, pr=0;

  an1[0]=1;
  an2[0]=1;
  bn1[0]=1;
  bn2[0]=0;
  der_old[0]=0;
  for(i=1;i<3;i++)
  {
	  an1[i]=0;
	  an2[i]=0;
	  bn1[i]=0;
	  bn2[i]=0;
	  der_old[i]=0;
  }

	
  while(del==1)
  {
	  n++;
	  if(n==1)
	  {
		  if (*x_in>*p_in/(*p_in+*q_in))
		  {
			  incompleBeta_an1_bn1_q(&x, p, q, an, bn);
		  }
		  else
		  {
			  incompleBeta_an1_bn1_p(&x, p, q, an, bn);

		  }
	  }
	  else
	  {
		  if (*x_in>*p_in/(*p_in+*q_in))
		  {
			  incompleBeta_an_bn_q(&x, p, q, n, an, bn);
		  }
		  else
		  {
			  incompleBeta_an_bn_p(&x, p, q, n, an, bn);
		  }
	  }
	  

	  // Use forward recurrance relations to compute An, Bn, and their derivatives
	  
	  dan[0]=an[0]*an2[0]+bn[0]*an1[0];
	  dbn[0]=an[0]*bn2[0]+bn[0]*bn1[0];
	  dan[1]=an[1]*an2[0]+an[0]*an2[1]+bn[1]*an1[0]+bn[0]*an1[1];
	  dbn[1]=an[1]*bn2[0]+an[0]*bn2[1]+bn[1]*bn1[0]+bn[0]*bn1[1];
	  dan[2]=an[2]*an2[0]+2*an[1]*an2[1]+an[0]*an2[2]+bn[2]*an1[0]+2*bn[1]*an1[1]+bn[0]*an1[2];
	  dbn[2]=an[2]*bn2[0]+2*an[1]*bn2[1]+an[0]*bn2[2]+bn[2]*bn1[0]+2*bn[1]*bn1[1]+bn[0]*bn1[2];
	  
	  
	  // Scale derivatives to prevent overflow
	  
	  Rn=dan[0];
	  if(fabs(dbn[0])>fabs(dan[0]))
	  {
	    Rn=dbn[0];
	  }
	  for(i=0;i<3;i++)
	  {
	      an1[i]=an1[i]/Rn;
	      bn1[i]=bn1[i]/Rn;
	  }
	  dan[1]=dan[1]/Rn;
	  dan[2]=dan[2]/Rn;
	  dbn[1]=dbn[1]/Rn;
	  dbn[2]=dbn[2]/Rn;
	  if(fabs(dbn[0])>fabs(dan[0]))
	  {
	    dan[0]=dan[0]/dbn[0];
	    dbn[0]=1;
	  }
	  else
	  {
	    dbn[0]=dbn[0]/dan[0];
	    dan[0]=1;
	  }
	  
	  // Compute components of derivatives of the nth approximant
	  
	  dr[0]=dan[0]/dbn[0];
	  Rn=dr[0];
	  dr[1]=(dan[1]-Rn*dbn[1])/dbn[0];
	  dr[2]=(-2*dan[1]*dbn[1]+2*Rn*dbn[1]*dbn[1])/dbn[0]/dbn[0]+(dan[2]-Rn*dbn[2])/dbn[0];
	  
	  // Save terms corresponding to approximants n-1 and n-2
	  
	  for(i=0;i<3;i++)
	  {
	    an2[i]=an1[i];
	    an1[i]=dan[i];
	    bn2[i]=bn1[i];
	    bn1[i]=dbn[i];
	  }
	  
	  //  Compute nth approximants
	  pr=0;
	  if(dr[0]>0)
	  {
	    pr=exp(c[0]+log(dr[0]));
	  }
	  der[0]=pr;
	  der[1]=pr*c[1]+c0*dr[1];
	  der[2]=pr*c[2]+2*c0*c[1]*dr[1]+c0*dr[2];
	  
	  
	  // Check for convergence, check for maximum and minimum iterations.
	  
	  for(i=0;i<3;i++)
	  {
	    d1[i]=MAX(err,fabs(der[i]));
	    d1[i]=fabs(der_old[i]-der[i])/d1[i];
	    der_old[i]=der[i];
	  }
	  d=MAX(MAX(d1[0],d1[1]),d1[2]);
	  
	  if(n< minappx)
	  {
	    d=1;
	  }
	  if(n>= maxappx)
	  {
	    d=0;
	  }
	  del=0;
	  if(d> err)
	  {
	    del=1;
	  }
	  
	  
  }
 
  	// Adjust results if I(x,p,q) = 1- I(1-x,q,p) was used
	  
	  if (*x_in>*p_in/(*p_in+*q_in))
	  {
		der[0]=1-der[0];
		der[1]=-der[1];
		der[2]=-der[2];
	  }

  Free(c);
  Free(an);
  Free(bn);
  Free(dan);
  Free(dbn);
  Free(dr);
  Free(an1);
  Free(an2);
  Free(bn1);
  Free(bn2);
  Free(d1);
  Free(der_old);
  
}


