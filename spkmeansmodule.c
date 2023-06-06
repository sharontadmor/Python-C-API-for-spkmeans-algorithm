#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject *spk(PyObject *self, PyObject *args);
static PyObject *wam(PyObject *self, PyObject *args);
static PyObject *ddg(PyObject *self, PyObject *args);
static PyObject *gl(PyObject *self, PyObject *args);
static PyObject *jacobi(PyObject *self, PyObject *args);
static PyObject *buildMatrix(matrix *mat);

static PyObject *spk(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyCentroids, *pyVec, *finalCentroids;
    int n, k, d, i, j, res;
    double coord;
    matrix *vectors, *centroids;
    if (!PyArg_ParseTuple(args, "OOi", &pyVectors, &pyCentroids, &k))
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    n = PyObject_Length(pyVectors);
    if (n <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (d <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    vectors = mallocMatrix(n, d);
    centroids = mallocMatrix(k, d);
    if (vectors == NULL || centroids == NULL)
    {
        matrixCleanup(vectors);
        matrixCleanup(centroids);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(vectors, i, j, coord);
        }
    }
    for (i = 0; i < k; i++)
    {
        pyVec = PyList_GetItem(pyCentroids, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(centroids, i, j, coord);
        }
    }
    /* update centroids array */
    res = kmeansC(n, k, vectors, centroids);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(vectors);
        matrixCleanup(centroids);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    finalCentroids = buildMatrix(centroids);
    matrixCleanup(vectors);
    matrixCleanup(centroids);
    return Py_BuildValue("O", finalCentroids);
}

static PyObject *wam(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyW;
    size_t i, j;
    int n, d, res;
    double coord;
    matrix *vectors, *w;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    n = PyObject_Length(pyVectors);
    if (n <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (d <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    vectors = mallocMatrix(n, d);
    w = mallocMatrix(n, n); /* weighed adjacency matrix */
    if (vectors == NULL || w == NULL)
    {
        matrixCleanup(vectors);
        matrixCleanup(w);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(vectors, i, j, coord);
        }
        for (j = 0; j < n; j++)
        {
            setCellValue(w, i, j, 0);
        }
    }
    res = wamC(n, d, vectors, w);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(vectors);
        matrixCleanup(w);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    pyW = buildMatrix(w);
    matrixCleanup(vectors);
    matrixCleanup(w);
    return Py_BuildValue("O", pyW);
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyDg;
    size_t i, j;
    int n, d, res;
    double coord;
    matrix *vectors, *dg;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    n = PyObject_Length(pyVectors);
    if (n <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (d <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    vectors = mallocMatrix(n, d);
    dg = mallocMatrix(n, n); /* diagonal degree matrix */
    if (vectors == NULL || dg == NULL)
    {
        matrixCleanup(vectors);
        matrixCleanup(dg);
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(vectors, i, j, coord);
        }
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                setCellValue(dg, i, j, 1);
            }
            else
            {
                setCellValue(dg, i, j, 0);
            }
        }
    }
    res = ddgC(n, d, vectors, dg);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(vectors);
        matrixCleanup(dg);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    pyDg = buildMatrix(dg);
    matrixCleanup(vectors);
    matrixCleanup(dg);
    return Py_BuildValue("O", pyDg);
}

static PyObject *gl(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyL;
    size_t i, j;
    int n, d, res;
    double coord;
    matrix *vectors, *l;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    n = PyObject_Length(pyVectors);
    if (n <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (d <= 0)
    {
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    vectors = mallocMatrix(n, d);
    l = mallocMatrix(n, n); /* graph laplacian matrix */
    if (vectors == NULL || l == NULL)
    {
        matrixCleanup(vectors);
        matrixCleanup(l);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(vectors, i, j, coord);
        }
        for (j = 0; j < n; j++)
        {
            setCellValue(l, i, j, 0);
        }
    }
    res = glC(n, d, vectors, l);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(vectors);
        matrixCleanup(l);
        return Py_BuildValue("i", EXIT_FAILURE);
    }
    pyL = buildMatrix(l);
    matrixCleanup(vectors);
    matrixCleanup(l);
    return Py_BuildValue("O", pyL);
}

static PyObject *jacobi(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyVals, *pyVecs;
    size_t i, j;
    int n, res;
    double coord;
    matrix *eigenvals, *eigenvecs, *mat;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("ii", EXIT_FAILURE, EXIT_FAILURE);
    }
    n = PyObject_Length(pyVectors);
    if (n <= 0)
    {
        return Py_BuildValue("ii", EXIT_FAILURE, EXIT_FAILURE);
    }
    mat = mallocMatrix(n, n);
    eigenvals = mallocMatrix(n, n);
    eigenvecs = mallocMatrix(n, n);
    if (mat == NULL || eigenvals == NULL || eigenvecs == NULL)
    {
        matrixCleanup(mat);
        matrixCleanup(eigenvals);
        matrixCleanup(eigenvecs);
        return Py_BuildValue("ii", EXIT_FAILURE, EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < n; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            setCellValue(mat, i, j, coord);
            setCellValue(eigenvals, i, j, coord);
            if (i == j)
            {
                setCellValue(eigenvecs, i, j, 1);
            }
            else
            {
                setCellValue(eigenvecs, i, j, 0);
            }
        }
    }
    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);
    if (pvalue != NULL)
    {
        matrixCleanup(mat);
        matrixCleanup(eigenvals);
        matrixCleanup(eigenvecs);
        return Py_BuildValue("ii", EXIT_FAILURE, EXIT_FAILURE); // ERROR_INVALID_ARGUMENT
    }
    res = jacobiC(n, mat, eigenvals, eigenvecs);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(mat);
        matrixCleanup(eigenvals);
        matrixCleanup(eigenvecs);
        return Py_BuildValue("ii", EXIT_FAILURE, EXIT_FAILURE);
    }
    pyVals = buildMatrix(eigenvals);
    pyVecs = buildMatrix(eigenvecs);
    matrixCleanup(mat);
    matrixCleanup(eigenvals);
    matrixCleanup(eigenvecs);
    return Py_BuildValue("OO", pyVals, pyVecs);
}

static PyObject *buildMatrix(matrix *mat)
{
    int i, j;
    PyObject *pyLst, *item, *num;
#if 0
    pyLst = PyList_New(k);
    for (i = 0; i < k; ++i)
    {
        item = PyList_New(d);
        for (j = 0; j < d; j++)
        {
            num = Py_BuildValue("d", arr[i * d + j]);
            PyList_SetItem(item, j, num);
        }
        PyList_SetItem(pyLst, i, item);
    }
#endif
    pyLst = PyList_New(mat->rowsNum);
    for (i = 0; i < mat->rowsNum; ++i)
    {
        item = PyList_New(mat->colsNum);
        for (j = 0; j < mat->colsNum; j++)
        {
            num = Py_BuildValue("d", mat->array[mat->colsNum * i + j]);
            PyList_SetItem(item, j, num);
        }
        PyList_SetItem(pyLst, i, item);
    }
    return pyLst;
}

static PyMethodDef kmeansMethods[] = {
    {"spk",
     (PyCFunction)spk,
     METH_VARARGS,
     PyDoc_STR(
         "Implementation for K-means clustering algorithm. \
         arguments: \
         (1) vectors - an array of all data points that were observed. \
         (2) centroids - an array of data points chosen as initial centroids. \
         (3) clustersNum - number of required clusters. \
         (4) vectorsNum - number of given vectors. \
         return : 0 if run was successful, 1 otherwise.")},
    {"wam",
     (PyCFunction)wam,
     METH_VARARGS,
     PyDoc_STR(
         "calculates weighed adjacency matrix. \
         arguments: \
         (1) n - number of vectors. \
         (2) d - dimention of each vector. \
         (3) vectors - array representing a matrix of dim n*d of \
         all data points that were observed. \
         (4) w - array representing a matrix of dim n*n. \
         return : 0 if operation is successful.")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR(
         "calculates diagonal degree matrix. \
         arguments: \
         (1) n - number of vectors. \
         (2) d - dimention of each vector. \
         (3) vectors - array representing a matrix of dim n*d of \
         all data points that were observed. \
         (4) dg - array representing a matrix of dim n*n. \
         return : 0 if run was successful, 1 otherwise.")},
    {"gl",
     (PyCFunction)gl,
     METH_VARARGS,
     PyDoc_STR(
         "calculates graph laplacian matrix. \
         arguments: \
         (1) n - number of vectors. \
         (2) d - dimention of each vector. \
         (3) vectors - array representing a matrix of dim n*d of \
         all data points that were observed. \
         (4) l - array representing a matrix of dim n*n. \
         return : 0 if run was successful, 1 otherwise.")},
    {"jacobi",
     (PyCFunction)jacobi,
     METH_VARARGS,
     PyDoc_STR(
         "Implementation for the Jacobi eigenvalue algorithm. \
         arguments: \
         (1) n - size of matrix a. \
         (2) a - array representing a symmetric matrix of dim n*n, \
         for which to calculate eigenvalues and eigenvectors. \
         (3) eigenvals - array representing a matrix of dim n*n, \
         in which the diagonal values are a's eigenvalues. \
         (4) eigenvecs - array representing a matrix of dim n*n, \
         in which the columns are a's eigenvectors. \
         return : 0 if run was successful, 1 otherwise.")},
    {"build_matrix",
     (PyCFunction)buildMatrix,
     METH_VARARGS,
     PyDoc_STR(
         "builds a python list of lists from C array, returns this python list. \
         argument mat : C array.")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "kmeans Module",
    -1,
    kmeansMethods};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}
