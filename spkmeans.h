#ifndef __SPKMEANS_H__
#define __SPKMEANS_H__
#define ERROR_OUT_OF_MEMORY 1
#define MAX_ITER 300
#define KMEANS_EPS 0
#define EPS 1.0 * pow(10, -5)
#define MAX_ROTATIONS 100
#define MATRIX_IS_DIAGONAL -1
#define GOAL_IDX 0
#define FILENAME_IDX 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int k;
int n;
int d;

typedef struct vector
{
    double *coords;
    struct vector *next;
    struct vector *prev;
} vector;

typedef struct matrix
/*
one dimentional array that will be used as a two dimentional matrix.
*/
{
    int rowsNum; /* number of rows, or number of vectors. */
    int colsNum; /* number od columns, or the dimention of each vector. */
    double *array;
} matrix;

void printMatrix(matrix *a);
matrix *mallocMatrix(int rowsNum, int colsNum);
double getCellValue(matrix *arr, int row, int col);
void setCellValue(matrix *arr, int row, int col, double val);
int kmeansC(int n, int d, int k, matrix *vectors, matrix *centroids);
int wamC(int n, int d, matrix *vectors, matrix *w);
int ddgC(int n, int d, matrix *vectors, matrix *dg);
int glC(int n, int d, matrix *vectors, matrix *l);
int jacobiC(int n, matrix *a, matrix *eigenvals, matrix *eigenvecs);

#endif /* __SPKMEANS_H__ */