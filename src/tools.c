/*
** tools.c - C code of the package CDRVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#include "vine.h"
#include "memoryhandling.h"
#include "tools.h"

//////////////////////////////////////////////////////////////////////////////
// Print error text and return to R
//////////////////////////////////////////////////////////////////////////////
void printError(char *text, char filename[200])
{
  Rprintf(text);
  Rprintf(": %s ", filename);
  Rprintf(" !!!\n");
}

//////////////////////////////////////////////////////////////////////////////
// Compute gamma division in a more stable way than gamma(x1)/gamma(x2)
// Input:
// x1 - Divisor
// x2 - Denominator
//-------------------
// a1 - modulus of a
// a2 - integer part of a
// b1 - modulus of b
// b2 - integer part of b
//-------------------
// We are after gamma(x1)/gamma(x2) and this computation will hopefully make it more numerically 
// stable. gamma(x1)/gamma(x2) will sometimes be INF/INF.
//////////////////////////////////////////////////////////////////////////////
double StableGammaDivision(double x1, double x2)
{
  int i;
  double a1, a2, b1, b2, sum=1.0;
  a1 = fmod(MAX(x1,x2),1.0);
  a2 = MAX(x1,x2)-a1;
  b1 = fmod(MIN(x1,x2),1.0);
  b2 = MIN(x1,x2)-b1;
  if(a1==0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
  }
  else if(a1>0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1);
  }
  else if(a1==0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum /= gammafn(b1);
  }
  else if(a1>0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1)/gammafn(b1);
  }
  if(x2 > x1) sum = 1.0/sum;
  return sum;
}


////////////////////////////////////////////////////////////////////////
// Fast computation of Kendall's tau
//
// Given by Shing (Eric) Fu, Feng Zhu, Guang (Jack) Yang, and Harry Joe
// Based on work of the method by Knight (1966)
/////////////////////////////////////////////////////////////////////////


void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V)
{
	// Defining variables
	int K, L, I, J, Iend, Jend;
	int i, j, m;
	double *Y2 = Calloc(*N, double);
	double *X2 = Calloc(*N, double);
	double *xptr,*yptr; // HJ addition for swapping
	boolean Iflag, Jflag, Xflag;
	*S = 0.; *D = 0.; *T = 0; *U = 0; *V = 0;

	/* 1.1 Sort X and Y in X order */
	/* Break ties in X according to Y */
	K=1;
	do
	{
		L=0;
		do
		{
			I = L;
			J = (I+K)<(*N)?(I+K):(*N);
			Iend = J;
			Jend = (J+K)<(*N)?(J+K):(*N);
			do
			{
				Iflag = (I < Iend);
				Jflag = (J < Jend);
				if (Iflag & Jflag) 
				{
				 	Xflag = ((X[I] > X[J]) | ((X[I] == X[J]) & (Y[I] > Y[J])));
				} 
				else
				{
					Xflag = FALSE;
				}
				if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
				{
					X2[L] = X[I];
					Y2[L] = Y[I];
					I++;
					L++;
				};
				if((!Iflag && Jflag) | (Iflag && Jflag && Xflag))
				{
					X2[L] = X[J];
					Y2[L] = Y[J];
					J++;
					L++;
				};
			} 
			while(Iflag | Jflag);
		} 
		while(L < *N);
		
		// Swap lists
		xptr=X; X=X2; X2=xptr;
		yptr=Y; Y=Y2; Y2=yptr;
		#ifdef OLD
		for(i = 0; i < *N; i++)
		{ 
			Xtem = X[i]; Ytem = Y[i];
			X[i] = X2[i]; Y[i] = Y2[i];
			X2[i] = Xtem; Y2[i] = Ytem;
		};
		#endif
		K *= 2;
	} 
	while (K < *N);

	/* 1.2 Count pairs of tied X, T */
	j = 1;
	m = 1;
	for(i = 1; i < *N; i++)
    if(X[i] == X[i-1])
    {
		j++;
		if(Y[i] == Y[i-1])
		m++;
    }
    else if(j > 1)
    {
      *T += j * (j - 1) / 2;
      if(m > 1)
		*V += m * (m - 1) / 2;
      j = 1;
      m = 1;
    };
	*T += j * (j - 1) / 2;
	*V += m * (m - 1) / 2;

	/* 2.1 Sort Y again and count exchanges, S */
	/* Keep original relative order if tied */

	K=1;
	do
	{
		L=0;
		do
		{
			I = L;
			J = (I+K)<(*N)?(I+K):(*N);
			Iend = J;
			Jend = (J+K)<(*N)?(J+K):(*N);
			do
			{
				Iflag = (I < Iend);
				Jflag = (J < Jend);
				if (Iflag & Jflag) 
				{
				 	Xflag = (Y[I] > Y[J]);
				} 
				else
				{
					Xflag = FALSE;
				}
				if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
				{
					X2[L] = X[I];
					Y2[L] = Y[I];
					I++;
					L++;
				};
				if((!Iflag && Jflag) | (Iflag && Jflag && Xflag))
				{
					X2[L] = X[J];
					Y2[L] = Y[J];
					*S += Iend - I;
					J++;
					L++;
				};
			} 
			while((Iflag | Jflag));
		} 
		while(L < *N);
    
		// Swap lists
		xptr=X; X=X2; X2=xptr;
		yptr=Y; Y=Y2; Y2=yptr;
		#ifdef OLD
		for(i = 0; i < *N; i++)
		{ 
			Xtem = X[i]; Ytem = Y[i];
			X[i] = X2[i]; Y[i] = Y2[i];
			X2[i] = Xtem; Y2[i] = Ytem;
		};
		#endif
		K *= 2;
	} 
	while (K < *N);

	/* 2.2 Count pairs of tied Y, U */
	j=1;
	for(i = 1; i < *N; i++)
		if(Y[i] == Y[i-1])
			j++;
		else if(j > 1)
		{
			*U += j * (j - 1) / 2;
			j = 1;
		};
	*U += j * (j - 1) / 2;


	/* 3. Calc. Kendall's Score and Denominator */
	*D = 0.5 * (*N) * (*N - 1);
	*S = *D - (2. * (*S) + *T + *U - *V);
	//if(*T > 0 | *U > 0) // adjust for ties
    *D = sqrt((*D - *T) * (*D - *U));
	*tau = (*S) / (*D);


  Free(Y2);
  Free(X2);
}


//////////////////////////////////////////////////
// ktau_matrix
//
// Input:
// data			data vector
// d			data dimension 1
// N			data dimension 2
//
// Output:
// out			Kendall's tau Matrix (as vector)

void ktau_matrix(double *data, int *d, int *N, double *out)
{
	double **x, S=0.0, D=0.0, *X, *Y;
	int k=0, i, j, t, T=0, U=0, V=0;
	x = create_matrix(*d,*N);
	X = (double*) Calloc(*N,double);
	Y = (double*) Calloc(*N,double);

	for(i=0;i<*d;i++)
    { 
		for (t=0;t<*N;t++ ) 
		{
			x[i][t] = data[k];
			k++;
		}
    }
	
	k=0;
	for(i=0;i<((*d)-1);i++)
	{
		for(j=(i+1);j<(*d);j++)
		{
			for(t=0;t<*N;t++)
			{
				X[t]=x[i][t];
				Y[t]=x[j][t];
			}
			ktau(X, Y, N, &out[k], &S, &D, &T, &U, &V);
			k++;
		}
	}

	Free(X);Free(Y);free_matrix(x, *d);
}


