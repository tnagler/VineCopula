/*
** memoryhandling.c - C code of the package CDRVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#include "VineCopula/vine.h"
#include "VineCopula/memoryhandling.h"

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a double matrix.
//  Input: Dimension of the matrix to be created//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
double **create_matrix(int rows, int columns)
{
  double **a;
  int i=0;
  a = (double**) Calloc(rows, double*);
  for(i=0;i<rows;i++) a[i] = (double*) Calloc(columns,double);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a double matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_matrix(double **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates an int matrix.
//  Input: Dimension of the matrix to be created.
//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
int **create_intmatrix(int rows, int columns)
{
  int **a;
  int i=0;
  a = (int**) Calloc(rows,int*);
  for(i=0;i<rows;i++) a[i] = (int*) Calloc(columns,int);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that an int matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_intmatrix(int **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a 3-d double array.
//  Input: Dimensions of the array to be created.
//  Output: Pointer to the created array.
///////////////////////////////////////////////////////////////////////////////
double ***create_3darray(int d1, int d2, int d3)
{
  double ***a;
  int i=0,j=0;
  a = (double ***) Calloc(d1,double*);
  for(i=0;i<d1;i++)
  {  
    a[i] = (double**) Calloc(d2, double*);
    for(j=0;j<d2;j++)
    {
      a[i][j] = (double*) Calloc(d3,double);
    }
  }
  return a;
}



///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a 3-d double array has been allocated.
//  Input: Dimensions of the array and a pointer to the array.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_3darray(double ***a, int d1, int d2)
{
  int i=0,j=0;
  for(i=0;i<d1;i++)
  {  
    for(j=0;j<d2;j++)
    {
      Free(a[i][j]);
    }
    Free(a[i]);
  }
  Free(a);
}
