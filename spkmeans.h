#ifndef __SPKMEANS_H__
#define __SPKMEANS_H__
#define ERROR_OUT_OF_MEMORY 1
#define EPS 1.0 * pow(10, -5)
#define MAX_ROTATIONS 100
#define MATRIX_IS_DIAGONAL -1
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
int k;
int iter;
double eps;
int n;
int d;
*/

typedef struct vector
{
    double *coords;
    struct vector *next;
    struct vector *prev;
} Vector;

int wamC(int n, int d, double (*vectors)[d], double (*w)[n]);
int ddgC(int n, int d, double (*vectors)[d], double (*dg)[n]);
int glC(int n, int d, double (*vectors)[d], double (*l)[n]);
int jacobiC(int n, double (*a)[n], double (*eigenvalues)[n], double (*eigenvectors)[n]);

PyMODINIT_FUNC PyInit_mykmeanssp(void);

#endif /* __SPKMEANS_H__ */