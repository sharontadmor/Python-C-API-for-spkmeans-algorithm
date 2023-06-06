#include "spkmeans.h"

static void setVectorData(vector *vec, double *coords, vector *next, vector *prev);
static void setMatrixData(matrix *mat, double *arr, int rowsNum, int colsNum);

void vectorCleanup(vector **partition, vector **clusters, int n)
/*
loop over partition, free nodes of vectors.
free partition.
free clusters.
*/
{
    int i;
    if (partition != NULL)
    {
        for (i = 0; i < n; i++)
        {
            if (partition[i] != NULL)
            {
                free(partition[i]);
            }
        }
        free(partition);
    }
    if (clusters != NULL)
    {
        free(clusters);
    }
}

void matrixCleanup(matrix *a)
/*
free a->array.
free a.
*/
{
    if (a != NULL)
    {
        free(a->array);
        free(a);
    }
}

vector *mallocVector(matrix *vectors, int vecIdx)
/*
returns a vector if vector creation is successful,
or NULL otherwise.
*/
{
    double *coords;
    vector *vec;
    vec = (vector *)malloc(sizeof(vector));
    if (vec == NULL)
    {
        return NULL;
    }
    coords = &vectors->array[vecIdx * vectors->colsNum];
    setVectorData(vec, coords, NULL, NULL);
    return vec;
}

static void setVectorData(vector *vec, double *coords, vector *next, vector *prev)
/*
sets data of vector node.
*/
{
    vec->coords = coords;
    vec->next = next;
    vec->prev = prev;
}

matrix *mallocMatrix(int rowsNum, int colsNum)
/*
returns a matrix if matrix creation is successful,
or NULL otherwise.
*/
{
    matrix *mat;
    double *arr;
    mat = (matrix *)malloc(sizeof(matrix));
    arr = (double *)malloc(sizeof(double) * rowsNum * colsNum);
    if (mat == NULL || arr == NULL)
    {
        return NULL;
    }
    setMatrixData(mat, arr, rowsNum, colsNum);
    return mat;
}

static void setMatrixData(matrix *mat, double *arr, int rowsNum, int colsNum)
/*
sets data of matrix.
*/
{
    mat->rowsNum = rowsNum;
    mat->colsNum = colsNum;
    mat->array = arr;
}

double getCellValue(matrix *mat, int row, int col)
/*
returns value of cell in a 2 dimensional array.
*/
{
    return mat->array[mat->colsNum * row + col];
}

void setCellValue(matrix *mat, int row, int col, double val)
/*
sets value of cell in a 2 dimensional array.
*/
{
    mat->array[mat->colsNum * row + col] = val;
}

