/*
 ** rvinesample.c - C code of the package VineCopula
 **
 */

#include "VineCopula/vine.h"
#include "VineCopula/memoryhandling.h"
#include "VineCopula/rvine.h"


int largest(int *array, int size)
{
    int i;
    int value = array[0];
    int index = 0;

    for (i = 0; i < size; i++) {
        if (array[i] > value) {
            value = array[i];
            index = i;
        }
    }
    return value;
}

int find_index(int *array, int size, int value)
{
    int i;
    int index = 0;

    for (i = 0; i < size; i++) {
        //printf("%d - %d - %d \n", i, array[i], value);
        if ( array[i] == value ) {
            index = i;
        }
    }
    return index;
}

void remove_element(int *array, int index, int size)
{
    int i;

    for(i = index; i < size; i++)
    {
        array[i] = array[i + 1];
    }
}

//////////////////////////////////////////////////////////////
// Function to construct the R-vine matrix from the binary matrix
// Input:
// b        the binary matrix
// Output:
// out      an Rvine matrix
/////////////////////////////////////////////////////////////
void getRVM(int* b, int* d, int* RVM)
{
    int i, j, ac, size, index, nn, n, **b2, **RVM2;

    b2 = create_intmatrix(*d,*d);
    RVM2 = create_intmatrix(*d,*d);
    //n = (*d)*((*d)-1)/2-((*d)-1)-1;

    //Initialize
    for (i=0;i<(*d);i++)
    {
        for (j=0;j<(*d);j++ )
        {
            b2[i][j]=b[(i+1)+(*d)*j-1] ;
            if (i == j || i == j-1)
            {
                RVM2[i][j] = i+1;
            }
            else
            {
                RVM2[i][j] = 0;
            }
        }
    }
    RVM2[0][2] = 1;

    n = 0;
    nn = 0;
    for (j=3;j<(*d);j++)
    {
        int *toAssign;
        size = j-1;
        toAssign=(int*) R_Calloc(size,double);
        for (i=0;i<size;i++ )
        {
            toAssign[i] = i+1;
        }

        ac = j-2;
        for (i=j-2;i>-1;i--)
        {
            //printf("before: %d - %d - %d - %d \n", i, j, ac, size);
            if (b2[i][j] == 1)
            {
                //printf("b1 \n");
                RVM2[i][j] = ac + 1;
                index = find_index(toAssign,size,ac+1);
                //printf("%d - %d \n", index, ac+1);
                if (size > 1)
                {
                    remove_element(toAssign, index, size);
                    size = size - 1;
                    ac = largest(toAssign, size) - 1;
                    //printf("%d \n", ac+1);
                }
            }
            else
            {
                //printf("b0 \n");
                RVM2[i][j] = RVM2[i-1][ac];
                index = find_index(toAssign,size,RVM2[i-1][ac]);
                remove_element(toAssign, index, size);
                size = size - 1;
            }
            //printf("after: %d - %d - %d - %d \n", i, j, RVM2[i][j], ac);
            //RVM[nn] = RVM2[i][j];
            RVM[n+i+1] = RVM2[i][j];
            nn = nn + 1;
        }
        n = nn;
        R_Free(toAssign);
    }
    RVM[0] = 1;

    free_intmatrix(b2,*d);
    free_intmatrix(RVM2,*d);
}
