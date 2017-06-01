#include "vine.h"
#include "memoryhandling.h"
#include "gof.h"
#include "rvinederiv2.h"
#include "pit.h"
#include "rvine.h"

/////////////////////////////////////////////////////////////////////
// Code form Daniel Berg, R-package copulaGOF
//     AD: Anderson-Darling  GOF test
//     (Cumulative distribution function test)
//     INPUT:
//          cdf         CDF for which to compute the test
//          n             Length of cdf
/////////////////////////////////////////////////////////////////////
void ADtest(double* cdf, int* n, double* out)
{
  int j;
  double sum=0.0;
  for(j=0;j<*n;j++) sum += (2.0*((double)j+1.0)-1.0)*(log(cdf[j])+log(1.0-cdf[*n-1-j]));
  *out = -(double)*n-(1.0/(double)*n)*sum;
}


///////////////////////////////////////////////////////////////////////////////
// Code form Daniel Berg, R-package copulaGOF
// Function to compute cumulative distribution function of a uniform vector x ($\hat F(x)$)
///////////////////////////////////////////////////////////////////////////////
void CumDist(double* x, int* i_n, int* i_m, double* out)
{
  int i,j,n,m;
  double *y;
  n=*i_n; m=*i_m;
  y = malloc(m*sizeof(double));
  for(i=0;i<m;i++)
  {
    y[i]=0.0;
    for(j=0;j<n;j++)
    {
      if(x[j]<=((double)i+1.0)/((double)m+1.0)) y[i] += 1.0/((double)n+1.0);
    }
    if(y[i]==0.0) y[i] = 1.0/((double)n+1.0);
    out[i] = y[i];
  }
  free(y);
}


/////////////////////////////////////////////////////////////////////
// Code form Daniel Berg, R-package copulaGOF
//     CvM: Cramer-von Mises GOF test
//     (Cumulative distribution function test)
//     INPUT:
//          cdf         CDF for which to compute the test
//          n             Length of cdf
/////////////////////////////////////////////////////////////////////
void CvMtest(double* cdf, int* n, double* out)
{
  int i;
  double sum1=0.0,sum2=0.0;
  for(i=0;i<*n;i++)
  {
    sum1 += pow(cdf[i],2.0);
    sum2 += cdf[i]*(2.0*((double)i+1.0)+1.0);
  }
  *out = (double)*n/3.0 + (double)*n/((double)*n + 1.0)*sum1 - (double)*n/(pow((double)*n+1.0,2.0))*sum2;
}


/////////////////////////////////////////////////////////////////////
// Code form Daniel Berg, R-package copulaGOF
//     KS: Kolmogorov-Smirnof GOF test
//     (Cumulative distribution function test)
//     INPUT:
//          cdf         CDF for which to compute the test
//          n             Length of cdf
/////////////////////////////////////////////////////////////////////
void KStest(double* cdf, int* n, double* out)
{
  int j;
  double tmp, maxdist=0.0;
  for(j=0;j<*n;j++)
  {
    tmp = MAX( fabs(cdf[j]-((double)j+1.0)/((double)*n+1.0)), fabs(cdf[j]-((double)j+2.0)/((double)*n+1.0)) );
    if(tmp>maxdist) maxdist = tmp;
  }
  *out = sqrt((double)*n)*maxdist;
}



////////////////////////////////////////////////////////
// Goodness-of-fit test based on White's information equality
// by U. Schepsmeier
///////////////////////////////////////////////////////////

void White(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, double* D, double* V)
{
	int i=0, dd=0, tt=0, k=1, j=0, kk=0, t=0, mm=0, dd2=0;
	double *Dprime, *hess, *subhess, *der, *subder, *dat, *hess_red, *der_red;

	for(i=0; i<(*d*(*d));i++)
	{
		if(family[i]!=0) dd++;
		if(family[i]==2) tt++;
	}
	mm=(dd+tt)*(dd+tt+1)/2;
	dd2=*d*(*d-1)/2;

	//Allocate memory
	Dprime = malloc((dd+tt)*(dd+tt+1)/2*sizeof(double));
	hess = malloc((dd2+tt)*(dd2+tt)*sizeof(double));
	subhess = malloc((dd2+tt)*(dd2+tt)*sizeof(double));
	der = malloc((dd2+tt)*(dd2+tt)*sizeof(double));
	subder = malloc((dd2+tt)*(dd2+tt)*sizeof(double));
	hess_red = malloc((dd+tt)*(dd+tt)*sizeof(double));
	der_red = malloc((dd+tt)*(dd+tt)*sizeof(double));
	dat = malloc(*d*sizeof(double));

	// initialisieren
	for(i=0;i<mm;i++)
	{
		Dprime[i]=0;
	}


	for(t=0;t<*T;t++)
	{
		for(i=0; i<*d;i++)
		{
			dat[i]=data[(t+1)+(i*(*T))-1];
		}
		for(i=0;i<((dd2+tt)*(dd2+tt));i++)
		{
			hess[i]=0;
			subhess[i]=0;
			der[i]=0;
			subder[i]=0;
		}
		hesse(&k, d, family, maxmat, matrix, condirect, conindirect, par, par2, dat, hess, subhess, der, subder);

		// independence aus der Hesse herausnehmen
		kk=0;
		for(i=0;i<dd2+tt;i++)
		{
			for(j=0;j<dd2+tt;j++)
			{
				if(hess[i+1+((dd2+tt)*j)-1]!=0)
				{
					hess_red[kk]=hess[i+1+((dd2+tt)*j)-1];
					der_red[kk]=der[i+1+((dd2+tt)*j)-1];
					kk++;
				}
			}
		}

		kk=0;
		for(i=0;i<(dd+tt);i++)
		{
			for(j=i;j<(dd+tt);j++)
			{
				Dprime[kk] = hess[(j+1)+((dd+tt)*i)-1] + der[(j+1)+((dd+tt)*i)-1];	// D = H+C
				D[kk] = D[kk] + (Dprime[kk]/(double)(*T));		// Schaetzer, deswegen durch die Anzahl der Beobachtungen teilen
				kk++;
			}
		}

		for(i=0;i<mm;i++)
		{
			for(j=0;j<mm;j++)
			{
				V[(i+1)+mm*j-1]+=(Dprime[i]*Dprime[j]/(double)(*T));
			}
		}
	}

	// Nicht fertig, da hier das Problem D%*%solve(V)%*%t(D) zu loesen ist

	// Free memory
	free(Dprime);
	free(hess);
	free(subhess);
	free(der);
	free(subder);
	free(dat);
	free(hess_red);
	free(der_red);
}


///////////////////////////////////////////////////
// Functions for PIT based GOF
// by U. Schepsmeier
///////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Function to compute probability of observing rank i, given rank 1,...,i-1
// Help function by Daniel Berg
// Input: q (vector to be transformed, must be sorted ascending),
// d(length of vector q), out(output)
///////////////////////////////////////////////////////////////////////////////

void ZStar(double* q, int* d, double* out)
{
  int i;
  double *qprev;
  qprev = malloc(*d*sizeof(double));
  for(i=0;i<*d;i++)
  {
    if(i==0) qprev[i]=0.0;
    else qprev[i] = q[i-1];
    out[i] = 1.0 - pow((1.0-q[i])/(1.0-qprev[i]),*d-i);
    if(out[i]==1.0) { out[i] = 0.9999999999; }
    else if(out[i]==0.0) { out[i] = 1.0-0.9999999999; }
  }
  free(qprev);
}

///////////////////////////////////////////////////////////////////////////////
// Function to compare two numbers
// Help function by Daniel Berg
///////////////////////////////////////////////////////////////////////////////
int comp_nums(const void *num1, const void *num2)
{
  if (*(double*)num1 <  *(double*)num2) return -1;
  else if (*(double*)num1 == *(double*)num2) return  0;
  else return  1;
}

////////////////////////////////////////
// Calculation of the B_j in Berg and Bakken
// Input:
// T, d 		dimensions
// family,...	RVM objects
// data			observed data
// vv vv2		side products of log-likelihood calculation (h-function)
// calcupdate	which h-functions have to be calculated
// method		numeric value (1=Breymann, 2=Berg, 3=Berg2; see my paper)
// alpha		power for the Berg function
//
// Output:
// out			sum of transformed data (PIT) (one step for the Breymann, Berg and Berg2 GOF)
//////////////////////////////////////////

void Bj(int *T, int* d, int* family, int* maxmat, int* matrix, int* condirect,
        int* conindirect, double* par, double* par2, double* data, double* out,
        double* vv, double* vv2, int* calcupdate, int* method, int *alpha)
{
	int i=0, t=0, ii=0, j=0;
	double *udata, **tmp, **u;
	udata = malloc(*d*(*T)*sizeof(double));
	tmp = create_matrix(*T,*d);
	u = create_matrix(*T,*d);

	RvinePIT(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, data, udata, vv, vv2, calcupdate);

	ii=0;
	for(j=0;j<*T;j++)
	{
		if(*method==2 || *method==3)
		{
			for(i=0;i<*d;i++)
			{
				u[j][i] = udata[ii];
				ii += 1;
			}
			qsort(u[j], *d, sizeof(double), comp_nums);
			ZStar(u[j],d,tmp[j]);		//Transformation von Berg and Bakken (2007); ordered PIT
		}
		else	// Im Fall von Breymann ist es besser keine Transformation zu machen
		{
			for(i=0;i<*d;i++)
			{
				tmp[j][i] = udata[ii];
				ii += 1;
			}
		}
	}

	for(t=0;t<*T;t++)
	{
		for(i=0;i<*d;i++)
		{
			if(*method==1)	//Breymann has the standard norma quantile function squared as transformation function
				tmp[t][i]=pow(qnorm(tmp[t][i],0.0,1.0,1,0),2);
			else if(*method==2)  // Berg: absolute value of u-0.5
				tmp[t][i]=fabs(tmp[t][i]-0.5);
			else if(*method==3)  // Berg2: (u-0.5)^alpha
				tmp[t][i]=pow(tmp[t][i]-0.5,*alpha);

			out[t]+=tmp[t][i];
		}
	}

	free(udata);
	free_matrix(tmp,*T);
	free_matrix(u,*T);
}


/////////////////////////
// For bootstrap simulate B_j
//
// Input:
// S			test statistic
// B			number of bootstrap samples
// T,d			dimensions
// method		numeric value (1=Breymann, 2=Berg, 3=Berg2; see my paper)
// alpha		power for the Berg function
//
// p			p-value

void SimulateBj(double* S, int *T, int* d, int* B, int* method, int *alpha, double* p)
{
	int i=0, t=0, m=0;
	double *tmp, Sb=0, *ustar;
	tmp = malloc(*d*sizeof(double));
	ustar = malloc(*d*sizeof(double));

	GetRNGstate();		// random number generator

	for(t=0;t<*T;t++)
	{
		p[t]=0;
	}

	for(m=0;m<*B;m++)
	{
		for(i=1;i<=*d;i++) { tmp[i] = runif(0,1);}
		qsort(tmp, *d, sizeof(double), comp_nums);
		ZStar(tmp,d,ustar);		//Transformation von Berg and Bakken (2007)

		for(i=0;i<*d;i++)
		{
			if(*method==1)
				tmp[i]=pow(qnorm(ustar[i],0.0,1.0,1,0),2);
			else if(*method==2)
				tmp[i]=fabs(ustar[i]-0.5);
			else if(*method==3)
				tmp[i]=pow(ustar[i]-0.5,*alpha);

			Sb+=tmp[i];
		}
		for(t=0;t<*T;t++)
		{
			if(Sb<=S[t]) p[t]+=1.0/(1.0+(double)*B);
		}
		Sb=0;
	}
	for(t=0;t<*T;t++)
	{
		if(p[t] == 0) p[t] = 1.0/(1.0+(double)*B);
	}

	PutRNGstate();

	free(tmp);
	free(ustar);
}


////////////////////////////////
// PIT based GOF tests
// Breymann, Berg, Berg2
//
// Input:
// data with its dimensions
// some help variables like vv and vv2
// method		numeric value (1=Breymann, 2=Berg, 3=Berg2; see my paper)
// alpha		power for the Berg function
// statisticName	numeric value (1=Anderson-Darling, 2=Kolmogorov-Smirnov, 3=Cramer-von Mises)
//
// Output:
// statistic		test statistic
/////////////////////////////////

void gofPIT_AD(int *T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		 double* statistic, double* vv, double* vv2, int* calcupdate, int* method, int *alpha, int* B, int *statisticName)
{
	int t=0;
	double *S, *helpvar, *Bhat;
	S = malloc(*T*sizeof(double));
	helpvar = malloc(*T*sizeof(double));
	Bhat = malloc(*T*sizeof(double));

	for(t=0;t<*T;t++)
	{
		S[t]=0;
		helpvar[t]=0;
		Bhat[t]=0;
	}

	Bj(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, data, S, vv, vv2, calcupdate, method, alpha);

	// Statistic berechnen
	if(*B==0)  // if an asymptotic based test statistic should be returned
	{
		if(*method==1)
		{
			for(t=0;t<*T;t++)
			{
				Bhat[t]=pchisq(S[t],*d,1.0,0.0);  // for Breymann the asymptotic is known (although it is shown that it is not that correct)
			}
		}
		else
			CumDist(S, T, T, Bhat);  // for the other two we need the empirical distribution function

		if(*statisticName==1)		//Anderson-Darling
			ADtest(Bhat, T, statistic);
		else if(*statisticName==2)	//Kolmogorov-Smirnov
			KStest(Bhat, T, statistic);
		else if(*statisticName==3)	//Cramer-von Mises
		{
			CvMtest(Bhat, T, statistic);
		}

	}
	else		//bootstrap
	{
		SimulateBj(S, T, d, B, method, alpha, helpvar);
		CumDist(helpvar, T, T, Bhat);
		if(*statisticName==1)		//Anderson-Darling
			ADtest(Bhat, T, statistic);
		else if(*statisticName==2)	//Kolmogorov-Smirnov
			KStest(Bhat, T, statistic);
		else if(*statisticName==3)	//Cramer-von Mises
		{
			CvMtest(Bhat, T, statistic);
		}
	}

	free(S);
	free(helpvar);
	free(Bhat);
}


///////////////////////////
// p-value estimation for the PIT based GOF tests
//
// Input:	see above
//
// Output:
// pvalue	bootstrapped p-value

void gofPIT_AD_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		double* statistic, double* vv, double* vv2, int* calcupdate, int* method, int* alpha, int* B, double* pvalue, int *statisticName)
{
	int i=0, j=0, m=0, t=0, *f, B2=1000;
	double *bdata, bstat=0;
	double *bvv, *bvv2;

	f = malloc(*T*sizeof(int));
	bdata = malloc(*d*(*T)*sizeof(double));
	bvv = malloc(*d*(*d)*(*T)*sizeof(double));
	bvv2 = malloc(*d*(*d)*(*T)*sizeof(double));

	for(m=0;m<*B;m++)	// B bootstrap steps
	{
		MySample(T, T, f);  // get a sample from my data
		for(t=0;t<*T;t++)
		{
			for(i=0;i<*d;i++)
			{
				bdata[(t+1)+(*T*i)-1]=data[(f[t])+(*T*i)-1];
				// Forget to change vv and vv2, too
				for(j=0; j<*d; j++)
				{
					bvv[(i+1)+(*d)*j+(*d)*(*d)*t-1] = vv[(i+1)+(*d)*j+(*d)*(*d)*(f[t]-1)-1];	// f[t]-1 because C starts to count at 0
					bvv2[(i+1)+(*d)*j+(*d)*(*d)*t-1] = vv2[(i+1)+(*d)*j+(*d)*(*d)*(f[t]-1)-1];
				}
			}
		}
		bstat=0;
		gofPIT_AD(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, bdata,
				&bstat, bvv, bvv2, calcupdate, method, alpha, &B2, statisticName);

		if(bstat>=*statistic)
			*pvalue+=1.0/(*B);
	}

	free(f);
	free(bdata);
	free(bvv);
	free(bvv2);
}


/////////////////////////////////////////////
/* Equal probability sampling; with-replacement case */
// Input:
// k		how many samples
// n		max value of sample
//
// output:
// y		vector of length k returning the samples
///////////////////////////////////////

void MySample(int *k, int *n, int *y)
{
    int i;

	GetRNGstate();
    for (i = 0; i < *k; i++)
	{
		y[i] = (int) *n * unif_rand() + 1;
	}
	PutRNGstate();
}


////////////////////////////////////////////////////////////////
// GOF test based on empirical copula process
//
// Input:
// data				data
// t,d				dimensions
// family,...		RVM object
// statisticName	numerical value (2=Kolmogorov-Smirnov, 3=Cramer-von Mises)
//
// Output:
// statistic		test statistic
////////////////////////////////////////////

void gofECP(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* data, double* statistic, int* statisticName)
{
	double *znull, *Chat1, *Chat2, U=0;
	int T2=1000, i=0, t=0, takeU=0;
	znull = malloc(*d*1000*sizeof(double));
	Chat1 = malloc(*T*sizeof(double));
	Chat2 = malloc(*T*sizeof(double));

	for(t=0;t<T2;t++)
	{
		for(i=0;i<*d;i++)
		{
			znull[t+1+(T2*i)-1]=0;
		}
	}

	SimulateRVine(&T2, d, family, maxmat, matrix, conindirect, par, par2, znull, &U, &takeU);


	ChatZj(data, data, T, d, T, Chat1);		// empirical copula distribution
	ChatZj(znull, data, T, d, &T2, Chat2);

	*statistic=0;
	if(*statisticName==3)	//Cramer-von Mises test statistic
	{
		for(i=0;i<*T;i++)
		{
			*statistic+=pow(Chat1[i]-Chat2[i],2);
		}
	}
	else if(*statisticName==2)	// KS
	{
		for(i=0;i<*T;i++)
		{
			*statistic=MAX(fabs(Chat1[i]-Chat2[i]),*statistic);
		}
		*statistic=*statistic*sqrt(*T);
	}

	free(znull);
	free(Chat1);
	free(Chat2);
}


///////////////////////////
// estimate p-value for ECP based GOF tests
//////////////////////////////

void gofECP_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* data, double* statistic, int* statisticName, double* pvalue, int* B)
{
	int i=0, m=0, t=0, *f;
	double *bdata, bstat=0;

	f = malloc(*T*sizeof(int));
	bdata = malloc(*d*(*T)*sizeof(double));

	for(m=0;m<*B;m++)
	{
		MySample(T, T, f);
		for(t=0;t<*T;t++)
		{
			for(i=0;i<*d;i++)
			{
				bdata[(t+1)+(*T*i)-1]=data[(f[t])+(*T*i)-1];
			}
		}
		bstat=0;
		gofECP(T, d, family, maxmat, matrix, conindirect, par, par2, bdata, &bstat, statisticName);

		if(bstat>=*statistic)
			*pvalue+=1.0/(*B);
	}

	free(f);
	free(bdata);
}

////////////////////////////
// Empirical copula distribution
//
// Input:
// n = dim(u)[1]
// m = dim(data)[1]
//
// Output:
// Chat		empirical copula distribution

void ChatZj(double* data, double* u, int* n, int* d, int* m, double* Chat)
{
	int i,j,k;
	double *helpvar;
	helpvar=malloc(*m*sizeof(double));

	for(j=0;j<*n;j++)
	{
		Chat[j]=0;
		for(k=0;k<*m;k++)
		{
			helpvar[k]=0;
			for(i=0;i<*d;i++)
			{
				if(data[k+1+(*m*i)-1]<=u[j+1+(*n*i)-1])
					helpvar[k]++;
			}
			if(helpvar[k]==*d)
				Chat[j]++;
		}
		Chat[j]/=(*m+1);
	}

	free(helpvar);
}


///////////////////////////
// Copula distribution of the independence copula
////////////////////////////

void C_ind(double* data, int* n, int* d, double* C)
{
	int t=0, i=0;

	for(t=0;t<*n;t++)
	{
		for(i=0;i<*d;i++)
		{
			if(i==0)
				C[t]=data[t+1+(*n*i)-1];
			else
				C[t]=C[t] * data[t+1+(*n*i)-1];
		}

	}
}



////////////////////////////
// GOF test based on ECP and PIT (ECP2)
// rest see above
//////////////////////

void gofECP2(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		double* vv, double* vv2, int* calcupdate, double* statistic, int* statisticName)
{
	double *udata, *Chat1, *Chat2;
	int i=0, t=0;
	udata = malloc(*d*(*T)*sizeof(double));
	Chat1 = malloc(*T*sizeof(double));
	Chat2 = malloc(*T*sizeof(double));

	for(t=0;t<*T;t++)
	{
		for(i=0;i<*d;i++)
		{
			udata[t+1+(*T*i)-1]=0;
		}
	}
	for(t=0;t<*T;t++)
	{
		Chat1[t]=0;
		Chat2[t]=1;
	}

	RvinePIT(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, data, udata, vv, vv2, calcupdate);
	ChatZj(udata, udata, T, d, T, Chat1);

	C_ind(udata,T,d,Chat2);

	*statistic=0;
	if(*statisticName==3)	//Cramer-von Mises test statistic
	{
		for(i=0;i<*T;i++)
		{
			*statistic+=pow(Chat1[i]-Chat2[i],2);
		}
	}
	else if(*statisticName==2)	// KS
	{
		for(i=0;i<*T;i++)
		{
			*statistic=MAX(fabs(Chat1[i]-Chat2[i]),*statistic);
		}
		*statistic=*statistic*sqrt(*T);
	}

	free(udata);
	free(Chat1);
	free(Chat2);
}


///////////////////////////
// p-value for ECP2
///////////////////////////

void gofECP2_pvalue(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data,
		double* vv, double* vv2, int* calcupdate, double* statistic, double* pvalue, int* statisticName, int* B)
{
	int i=0, j=0, m=0, t=0, *f;
	double *bdata, bstat=0;
	double *bvv, *bvv2;

	f = malloc(*T*sizeof(int));
	bdata = malloc(*d*(*T)*sizeof(double));
	bvv = malloc(*d*(*d)*(*T)*sizeof(double));
	bvv2 = malloc(*d*(*d)*(*T)*sizeof(double));

	for(m=0;m<*B;m++)
	{
		MySample(T, T, f);
		for(t=0;t<*T;t++)
		{
			for(i=0;i<*d;i++)
			{
				bdata[(t+1)+(*T*i)-1]=data[(f[t])+(*T*i)-1];
				// Forget to change vv annd vv2, too
				for(j=0; j<*d; j++)
				{
					bvv[(i+1)+(*d)*j+(*d)*(*d)*t-1] = vv[(i+1)+(*d)*j+(*d)*(*d)*(f[t]-1)-1];	// f[t]-1 because C starts to count at 0
					bvv2[(i+1)+(*d)*j+(*d)*(*d)*t-1] = vv2[(i+1)+(*d)*j+(*d)*(*d)*(f[t]-1)-1];
				}
			}
		}
		bstat=0;
		gofECP2(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, bdata, bvv, bvv2, calcupdate, &bstat, statisticName);

		if(bstat>=*statistic)
			*pvalue+=1.0/(*B);
	}

	free(f);
	free(bdata);
	free(bvv);
	free(bvv2);
}
