#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if 0
double getCellValue(double *arr, int rowLen, int row, int col)
/*
returns value of cell in a 2 dimensional array.
*/
{
    return arr[row * rowLen + col];
}

void setCellValue(double *arr, int rowLen, int row, int col, double val)
/*
sets value of cell in a 2 dimensional array.
*/
{
    arr[row * rowLen + col] = val;
}

int wamC(double *vectors, double *w, int n, int d)
/*
vectors (double *) - one dimentional array, representing a matrix of dim n*d of all data points that were observed.
w - one dimentional array, representing a matrix of dim n*n.
n - number of vectors.
d - dimention of each vector.
this function calculates w, the weighed adjacency matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    int i, j, l;
    double euclid_dist, w_ij, diff;
    /*
    loop over all the vectors,
    for each vector v_i it's enough to calculate the weight w_ij with all the vectors v_j with index j larger then i,
    because w is a symmetric matrix.
    w's diagonal remains all zero.
    */
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            euclid_dist = 0;
            w_ij = 0;
            for (l = 0; l < d; l++)
            {
                diff = getCellValue(vectors, d, i, l) - getCellValue(vectors, d, j, l);

                #if 0
                printf("coord from i: %f\n", getCellValue(vectors, d, i, l));
                printf("coord from j: %f\n", getCellValue(vectors, d, j, l));
                printf("diff: %f\n", diff);
                #endif

                euclid_dist += (diff * diff);
            }
            w_ij = exp((-1.0 / 2) * euclid_dist);
            setCellValue(w, n, i, j, w_ij);
            setCellValue(w, n, j, i, w_ij);

            #if 0
            printf("euclid_dist: %f\n", euclid_dist);
            printf("(-1.0 / 2) * euclid_dist: %f\n", (-1.0 / 2) * euclid_dist);
            printf("w_ij: %f\n", w_ij);
            #endif
        }
        
    }
    return 0;
}
#endif

int wamC(int n, int d, double (*vectors)[d], double (*w)[n])
/*
n - number of vectors.
d - dimention of each vector.
vectors - array representing a matrix of dim n*d of all data points that were observed.
w - array representing a matrix of dim n*n.
this function calculates w, the weighed adjacency matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    int i, j, l;
    double euclid_dist, w_ij, diff;
    /*
    loop over all the vectors,
    for each vector v_i it's enough to calculate the weight w_ij with all the vectors v_j with index j larger then i,
    because w is a symmetric matrix.
    w's diagonal remains all zero.
    */
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            euclid_dist = 0;
            w_ij = 0;
            for (l = 0; l < d; l++)
            {
                diff = vectors[i][l] - vectors[j][l];
                euclid_dist += (diff * diff);
            }
            w_ij = exp((-1.0 / 2) * euclid_dist);
            w[i][j] = w_ij;
            w[j][i] = w_ij;
        }
    }
    return 0;
}

int ddgC(int n, int d, double (*vectors)[d], double (*dg)[n])
/*
n - number of vectors.
d - dimention of each vector.
vectors - array representing a matrix of dim n*d of all data points that were observed.
dg - array representing a matrix of dim n*n.
this function calculates dg, the diagonal degree matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    return 0;
}

int glC(int n, int d, double (*vectors)[d], double (*l)[n])
/*
n - number of vectors.
d - dimention of each vector.
vectors - array representing a matrix of dim n*d of all data points that were observed.
l - array representing a matrix of dim n*n.
this function calculates l, the graph laplacian matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    return 0;
}

int jacobiC(int n, int d, double (*a)[d], double (*jac)[n])
/*
the Jacobi eiganvalue algorithm is an iterative method
for the calculation of eiganvalues and eiganvectors of a real symmetric matrix.
n - number of vectors.
d - dimention of each vector.
a - array representing a symmetric matrix of dim n*n,
for which to calculate eigenvalues and eigenvectors.
jac - 
calculates:
(1) diag - diagonalized matrix where the diagonal values are a's eigenvalues.
(2) vecs - matrix in which the columns are a's eigenvectors.
returns 0 if operation is successful.
*/
{
    return 0;
}


int main(int argc, char *argv[])
{
}
