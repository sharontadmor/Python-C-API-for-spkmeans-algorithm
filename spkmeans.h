#ifndef __SPKMEANS_H__
#define __SPKMEANS_H__
#define ERROR_OUT_OF_MEMORY 1
#define GENERAL_ERROR 1
#define OPERATION_SUCCESSFL 0
#define GENERAL_ERROR_MESSAGE "An Error Has Occurred\n"
#define MAX_ITER 300
#define KMEANS_EPS 0
#define EPS 1.0 * pow(10, -5)
#define MAX_ROTATIONS 100
#define ARGS_NUM 3
#define GOAL_IDX 1
#define FILENAME_IDX 2
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct vector
/*
vector node in a linked list.
*/
{
    double *coords;
    struct vector *next;
    struct vector *prev;
} vector;

typedef struct matrix
/*
one dimentional array that will be used as a matrix.
*/
{
    int rowsNum; /* number of rows, or number of vectors. */
    int colsNum; /* number od columns, or the dimention of each vector. */
    double *array;
} matrix;

/* Utils */
void vectorCleanup(vector **partition, vector **clusters, int n);
void matrixCleanup(matrix *a);
vector *mallocVector(matrix *vectors, int vecIdx);
matrix *mallocMatrix(int rowsNum, int colsNum);
double getCellValue(matrix *mat, int row, int col);
void setCellValue(matrix *mat, int row, int col, double val);
/* main */
void printMatrix(matrix *a);
void printJacobi(matrix *eigenvals, matrix *eigenvecs);
/* spkmeans */
int kmeansC(int n, int k, matrix *vectors, matrix *centroids);
int wamC(int n, int d, matrix *vectors, matrix *w);
int ddgC(int n, int d, matrix *vectors, matrix *dg);
int glC(int n, int d, matrix *vectors, matrix *l);
int jacobiC(int n, matrix *a, matrix *eigenvals, matrix *eigenvecs);

#endif /* __SPKMEANS_H__ */