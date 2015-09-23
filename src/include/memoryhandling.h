#if !defined(MEMORY_H)
#define MEMORY_H

// File memoryhandling.c
double **create_matrix(int rows, int columns);
void free_matrix(double **a, int rows);
int **create_intmatrix(int rows, int columns);
void free_intmatrix(int **a, int rows);
double ***create_3darray(int d1, int d2, int d3);
void free_3darray(double ***a, int d1, int d2);

#endif
