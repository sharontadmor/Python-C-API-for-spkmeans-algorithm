#ifndef KMEANS_H_
#define KMEANS_H_
#define ERROR_OUT_OF_MEMORY 1
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

double getCellValue(double *arr, int rowLen, int row, int col);
void setCellValue(double *arr, int rowLen, int row, int col, double val);
int wamC(int n, int d, double (*vectors)[d], double (*w)[n]);
int ddgC(int n, int d, double (*vectors)[d], double (*dg)[n]);
int glC(int n, int d, double (*vectors)[d], double (*l)[n]);
int jacobiC(int n, int d, double (*a)[d], double (*j)[n]);

static PyObject *wam(PyObject *self, PyObject *args);
static PyObject *getList(double *arr, int k, int d);
PyMODINIT_FUNC PyInit_mykmeanssp(void);

#endif /* KMEANS */