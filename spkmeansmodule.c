#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"

static PyObject *spk(PyObject *self, PyObject *args);
static PyObject *wam(PyObject *self, PyObject *args);
static PyObject *ddg(PyObject *self, PyObject *args);
static PyObject *gl(PyObject *self, PyObject *args);
static PyObject *jacobi(PyObject *self, PyObject *args);
static PyObject *getList(double *arr, int k, int d);

static PyObject *spk(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyCentroids, *pyVec, *finalCentroids;
    int n, k, d, i, j, res;
    double coord;
    if (!PyArg_ParseTuple(args, "OOi", &pyVectors, &pyCentroids, &k))
    {
        return NULL;
    }
    n = PyObject_Length(pyVectors);
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (n < 0 || d < 0)
    {
        return NULL;
    }
    double (*vectors)[d], (*centroids)[d];
    vectors = (double (*)[d])malloc(sizeof(double) * n * d);
    centroids = (double (*)[d])malloc(sizeof(double) * k * d);
    if (vectors == NULL || centroids == NULL)
    {
        /* cleanup */
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            vectors[i][j] = coord;
        }
    }
    for (i = 0; i < k; i++)
    {
        pyVec = PyList_GetItem(pyCentroids, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            centroids[i][j] = coord;
        }
    }
    /* update centroids array */
    res = kmeansC(n, d, k, vectors, centroids);
    if (res == ERROR_OUT_OF_MEMORY)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    finalCentroids = getList(centroids[0], k, d);
    /* TODO : cleanup */
    return Py_BuildValue("O", finalCentroids);
}

static PyObject *wam(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyW;
    size_t i, j;
    int n, d, res;
    double coord;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    n = PyObject_Length(pyVectors);
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (n < 0 || d < 0)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    double (*vectors)[d], (*w)[n];
    vectors = (double (*)[d])malloc(sizeof(double) * n * d);
    w = (double (*)[n])malloc(sizeof(double) * n * n);  /* weighed adjacency matrix */
    if (vectors == NULL || w == NULL)
    {
        /* cleanup */
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            vectors[i][j] = coord;
        }
        for (j = 0; j < n; j++)
        {
            w[i][j] = 0;
        }
    }
    res = wamC(n, d, vectors, w);
    if (res == ERROR_OUT_OF_MEMORY)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    pyW = getList(w[0], n, n);
    /* cleanup */
    return Py_BuildValue("O", pyW);
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyDg;
    size_t i, j;
    int n, d, res;
    double coord;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    n = PyObject_Length(pyVectors);
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (n < 0 || d < 0)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    double (*vectors)[d], (*dg)[n];
    vectors = (double (*)[d])malloc(sizeof(double) * n * d);
    dg = (double (*)[n])malloc(sizeof(double) * n * n);  /* diagonal degree matrix */
    if (vectors == NULL || dg == NULL)
    {
        /* cleanup */
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            vectors[i][j] = coord;
        }
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                dg[i][j] = 1;
            }
            else
            {
                dg[i][j] = 0;
            }
        }
    }
    res = ddgC(n, d, vectors, dg);
    if (res == ERROR_OUT_OF_MEMORY)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    pyDg = getList(dg[0], n, n);
    /* cleanup */
    return Py_BuildValue("O", pyDg);
}

static PyObject *gl(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyL;
    size_t i, j;
    int n, d, res;
    double coord;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    n = PyObject_Length(pyVectors);
    d = PyObject_Length(PyList_GetItem(pyVectors, 0));
    if (n < 0 || d < 0)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    double (*vectors)[d], (*l)[n];
    vectors = (double (*)[d])malloc(sizeof(double) * n * d);
    l = (double (*)[n])malloc(sizeof(double) * n * n);  /* graph laplacian matrix */
    if (vectors == NULL || l == NULL)
    {
        /* cleanup */
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < d; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            vectors[i][j] = coord;
        }
        for (j = 0; j < n; j++)
        {
            l[i][j] = 0;
        }
    }
    res = glC(n, d, vectors, l);
    if (res == ERROR_OUT_OF_MEMORY)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    pyL = getList(l[0], n, n);
    /* cleanup */
    return Py_BuildValue("O", pyL);
}

static PyObject *jacobi(PyObject *self, PyObject *args)
{
    PyObject *pyVectors, *pyVec, *pyVals, *pyVecs;
    size_t i, j;
    int n, res;
    double coord;
    if (!PyArg_ParseTuple(args, "O", &pyVectors))
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    n = PyObject_Length(pyVectors);
    if (n < 0)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    double (*mat)[n], (*eigenvalues)[n], (*eigenvectors)[n];
    mat = (double (*)[n])malloc(sizeof(double) * n * n);
    eigenvalues = (double(*)[n])malloc(sizeof(double) * n * n);
    eigenvectors = (double(*)[n])malloc(sizeof(double) * n * n);
    if (mat == NULL || eigenvalues == NULL || eigenvectors == NULL)
    {
        /* cleanup */
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    for (i = 0; i < n; i++)
    {
        pyVec = PyList_GetItem(pyVectors, i);
        for (j = 0; j < n; j++)
        {
            coord = PyFloat_AsDouble(PyList_GetItem(pyVec, j));
            mat[i][j] = coord;
            eigenvalues[i][j] = coord;
            eigenvectors[i][j] = 0;
        }
    }
    res = jacobiC(n, mat, eigenvalues, eigenvectors);
    if (res == ERROR_OUT_OF_MEMORY)
    {
        return Py_BuildValue("i", ERROR_OUT_OF_MEMORY);
    }
    pyVals = getList(eigenvalues[0], n, n);
    pyVecs = getList(eigenvectors[0], n, n);
    /* cleanup */
    return Py_BuildValue("OO", pyVals, pyVecs);
}

static PyObject *getList(double *arr, int k, int d)
{
    int i, j;
    PyObject *pyLst, *item, *num;
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
    return pyLst;
}

static PyMethodDef kmeansMethods[] = {
    {"spk",
     (PyCFunction)spk,
     METH_VARARGS,
     PyDoc_STR(
         "description. \
         arguments description: \
         (1) vectors - an array of all data points that were observed. \
         (2) ... \
         return : .")},
    {"wam",
     (PyCFunction)wam,
     METH_VARARGS,
     PyDoc_STR(
         "description. \
         arguments description: \
         (1) vectors - an array of all data points that were observed. \
         (2) ... \
         return : .")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR(
         "description. \
         arguments description: \
         (1) vectors - an array of all data points that were observed. \
         (2) ... \
         return : .")},
    {"gl",
     (PyCFunction)gl,
     METH_VARARGS,
     PyDoc_STR(
         "description. \
         arguments description: \
         (1) vectors - an array of all data points that were observed. \
         (2) ... \
         return : .")},
    {"jacobi",
     (PyCFunction)jacobi,
     METH_VARARGS,
     PyDoc_STR(
         "description. \
         arguments description: \
         (1) vectors - an array of all data points that were observed. \
         (2) ... \
         return : .")},
    {"get_list",
     (PyCFunction)getList,
     METH_VARARGS,
     PyDoc_STR(
         "builds a python list of lists from C array, returns this python list. \
        arguments description: \
        (1) arr - C array. \
        (2) k (int) - number of lists in list. \
        (3) d (int) - length of lists in list.")},
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
