#include "spkmeans.h"

static double getSumSquares(int n, double (*a)[n]);
static int isBreakCondition(size_t count, double offA, double offVecs);
static int getPivot(int n, double (*a)[n], size_t *idx);
static int calculateP(int n, double (*a)[n], double (*p)[n], size_t i, size_t j);
static double getT(int n, double (*a)[n], size_t i, size_t j);
static double getC(double t);
static double getS(double c, double t);
static int rotate(int n, double (*a)[n], double (*b)[n], double (*p)[n], size_t i, size_t j);
static int matrixMultiplication(int n, double (*a)[n], double (*b)[n], double (*prod)[n]);

void printMatrix(int size, int dim, double (*a)[dim])
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < dim; j++)
        {
            printf("%f, ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int kmeansC(int n, int d, int k, double (*vectors)[d], double (*centroids)[d])
/*
n - number of vectors.
d - dimention of each vector.
k - number of clusters.
vectors - array representing a matrix of dim n*d of all data points that were observed.
centroids - array representing a matrix of dim k*d of the initial centroids.
this function calculates the final centroids of clustering algorithm kmeans++.
returns 0 if operation is successful.
*/
{
    return 0;
}

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
#if 0
                if (i == 0 && j == 2)
                {
                    printf("coord 1: %f\n", vectors[i][l]);
                    printf("coord 2: %f\n", vectors[j][l]);
                    printf("diff: %f\n", diff);
                }
#endif
            }
            w_ij = exp((-1.0 / 2) * euclid_dist);
            w[i][j] = w_ij;
            w[j][i] = w_ij;
#if 0
            if (i == 0 && j == 2)
            {
                printf("---> euclid_dist: %f\n", euclid_dist);
                printf("---> w_ij: %f\n", w_ij);
            }
#endif
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
    size_t i, j;
    double(*w)[n], degree;
    w = (double(*)[n])malloc(sizeof(double) * n * n); /* weighed adjacency matrix */
    if (w == NULL)
    {
        /* cleanup */
        return ERROR_OUT_OF_MEMORY;
    }
    wamC(n, d, vectors, w);
    for (i = 0; i < n; i++)
    {
        degree = 0;
        for (j = 0; j < n; j++)
        {
            degree += w[i][j];
        }
        dg[i][i] = degree;
    }
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
    size_t i, j;
    double(*w)[n], (*dg)[n];
    w = (double(*)[n])malloc(sizeof(double) * n * n);  /* weighed adjacency matrix */
    dg = (double(*)[n])malloc(sizeof(double) * n * n); /* diagonal degree matrix */
    if (w == NULL || dg == NULL)
    {
        /* cleanup */
        return ERROR_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            w[i][j] = 0;
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
    wamC(n, d, vectors, w);
    ddgC(n, d, vectors, dg);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            l[i][j] = dg[i][j] - w[i][j];
        }
    }
    return 0;
}

int jacobiC(int n, double (*a)[n], double (*eigenvalues)[n], double (*eigenvectors)[n])
/*
the Jacobi eigenvalue algorithm is an iterative method
for the calculation of eigenvalues and eigenvectors of a real symmetric matrix.
n - size of matrix a.
a - array representing a symmetric matrix of dim n*n,
for which to calculate eigenvalues and eigenvectors.
eigenvalues -
eigenvectors -
calculates:
(1) diag - diagonalized matrix where the diagonal values are a's eigenvalues.
(2) vecs - matrix in which the columns are a's eigenvectors.
returns 0 if operation is successful.
*/
{
    size_t i, j, count, idx[2], res;
    /* int count; */
    double offA, offVecs;
    double(*p)[n], (*temp)[n], (*prod)[n];
    idx[0] = 0;
    idx[1] = 0;
    offA = 0;
    offVecs = 0;
    count = 0;
    p = (double(*)[n])malloc(sizeof(double) * n * n);
    if (p == NULL)
    {
        /* cleanup */
        return ERROR_OUT_OF_MEMORY;
    }
    offA = getSumSquares(n, a);
    while (!isBreakCondition(count, offA, offVecs))
    {
        count += 1;

        /*
        printf("---> iteration: %d\n", count);
        printf("p 1 - start of iter:\n");
        printMatrix(n, n, p);
        printf("eigenvalues 1 - start of iter:\n");
        printMatrix(n, n, eigenvalues);
        */

        res = getPivot(n, eigenvalues, idx);
        if (res == MATRIX_IS_DIAGONAL)
        {
            if (count == 1)
            {
                /* eigenvectors = p; */
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        eigenvectors[i][j] = p[i][j];
                    }
                }
            }
            else
            {
                prod = (double(*)[n])malloc(sizeof(double) * n * n);
                if (prod == NULL)
                {
                    /* cleanup */
                    return ERROR_OUT_OF_MEMORY;
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        prod[i][j] = 0;
                    }
                }
                matrixMultiplication(n, eigenvectors, p, prod);
                free(eigenvectors);
                eigenvectors = prod;
            }
            return 0;
        }
        temp = (double(*)[n])malloc(sizeof(double) * n * n);
        if (temp == NULL)
        {
            /* cleanup */
            return ERROR_OUT_OF_MEMORY;
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                temp[i][j] = eigenvalues[i][j];
                if (i == j)
                {
                    p[i][j] = 1;
                }
                else
                {
                    p[i][j] = 0;
                }
            }
        }

        /*
        printf("p 2 - after p init:\n");
        printMatrix(n, n, p);
        printf("eigenvalues 2 - after p init:\n");
        printMatrix(n, n, eigenvalues);
        */

        calculateP(n, eigenvalues, p, idx[0], idx[1]);
        rotate(n, eigenvalues, temp, p, idx[0], idx[1]);
        free(eigenvalues);
        eigenvalues = temp;

        /*
        printf("p 3 - after calculations:\n");
        printMatrix(n, n, p);
        printf("eigenvalues 3 - after calculations:\n");
        printMatrix(n, n, eigenvalues);
        */

        if (count == 1)
        {
            /* eigenvectors = p; */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    eigenvectors[i][j] = p[i][j];
                }
            }
        }
        else
        {
            prod = (double(*)[n])malloc(sizeof(double) * n * n);
            if (prod == NULL)
            {
                /* cleanup */
                return ERROR_OUT_OF_MEMORY;
            }
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    prod[i][j] = 0;
                }
            }
            matrixMultiplication(n, eigenvectors, p, prod);
            free(eigenvectors);
            eigenvectors = prod;
        }
        offVecs = getSumSquares(n, eigenvalues);
        /*
        printf("p 4 - end of iter:\n");
        printMatrix(n, n, p);
        */
    }

    /*
    printf("original a:\n");
    printMatrix(n, n, a);
    printf("final eigenvectors:\n");
    printMatrix(n, n, eigenvectors);
    printf("final eigenvalues:\n");
    printMatrix(n, n, eigenvalues);
    */

    return 0;
}

static double getSumSquares(int n, double (*a)[n])
/*
n - size of matrix a.
a - array representing symmetric matrix of dim n*n.
returns the sum of squares of all off-diagonal elements of matrix a.
*/
{
    size_t i, j;
    double sum;
    sum = 0;
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            sum += a[i][j] * a[i][j];
        }
    }
    return 2.0 * sum;
}

static int isBreakCondition(size_t count, double offA, double offVecs)
{
    int convergence;
    convergence = offA - offVecs <= EPS;
    if (count == MAX_ROTATIONS || convergence)
    {
        return 1;
    }
    return 0;
}

static int getPivot(int n, double (*a)[n], size_t *idx)
/*
n - size of matrix a.
a - array representing symmetric matrix of dim n*n.
idx - array of size 2.
this function finds pivot a_ij,
and puts its row and column in array idx, such that the 0 index is the row i of a_ij, and the index 1 is the column j of a_ij.
pivot is the off-diagonal element of matrix a with the largest absolute value.
*/
{
    size_t i, j;
    double pivot;
    pivot = 0;
    /* printMatrix(n, n, a); */
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(a[i][j]) > pivot)
            {
                pivot = fabs(a[i][j]);
                idx[0] = i;
                idx[1] = j;

#if 0
                printf("---> i: %d, j: %d\n", i, j);
                printf("---> pivot: %f\n", pivot);
#endif
            }
        }
    }
    if (pivot == 0) /* all off-diagonal elements are 0 */
    {
        return MATRIX_IS_DIAGONAL;
    }
    return 0;
}

static int calculateP(int n, double (*a)[n], double (*p)[n], size_t i, size_t j)
/*
calculates rotation matrix p.
*/
{
    double t, c, s;
    t = getT(n, a, i, j);
    c = getC(t);
    s = getS(c, t);
    p[i][i] = c;
    p[j][j] = c;
    p[i][j] = s;
    p[j][i] = (-1) * s;

    /*
    printf("---> t: %f\n", t);
    printf("---> c: %f\n", c);
    printf("---> s: %f\n", s);
    */

    return 0;
}

static double getT(int n, double (*a)[n], size_t i, size_t j)
/*
calculates value t out of matrix a.
*/
{
    double theta, denominator, t;
    theta = 0;
    theta = (a[j][j] - a[i][i]) / (2 * a[i][j]);
    /*
    printf("i: %d, j: %d\n", i, j);
    printf("a:\n");
    printMatrix(n, n, a);
    printf("---> theta: %f\n", theta);
    */
    denominator = fabs(theta) + sqrt(theta * theta + 1);
    if (theta >= 0)
    {
        t = 1 / denominator;
    }
    else
    {
        t = (-1) / denominator;
    }
    return t;
}

static double getC(double t)
/*
calculates value c out of matrix a.
*/
{
    double c;
    c = 1 / sqrt(t * t + 1);
    return c;
}

static double getS(double c, double t)
/*
calculates value s out of matrix a.
*/
{
    double s;
    s = c * t;
    return s;
}

static int rotate(int n, double (*a)[n], double (*b)[n], double (*p)[n], size_t i, size_t j)
/*
rotates matrix b to make b_ij = 0.
*/
{
    size_t r;
    double c, s;
    c = p[i][i];
    s = p[i][j];
    b[i][i] = (c * c * a[i][i]) + (s * s * a[j][j]) - (2 * s * c * a[i][j]);
    b[j][j] = (s * s * a[i][i]) + (c * c * a[j][j]) + (2 * s * c * a[i][j]);
    b[i][j] = 0; /* b[i][j] = (c * c - s * s) * a[i][j] + (s * c)(a[i][i] - a[j][j]) = 0 */
    b[j][i] = 0;
    for (r = 0; r < n; r++)
    {
        if ((r != i) && (r != j))
        {
            b[r][i] = c * a[r][i] - s * a[r][j];
            b[r][j] = c * a[r][j] + s * a[r][i];
            b[i][r] = c * a[r][i] - s * a[r][j];
            b[j][r] = c * a[r][j] + s * a[r][i];
        }
    }

    return 0;
}

static int matrixMultiplication(int n, double (*a)[n], double (*b)[n], double (*prod)[n])
/*
calculates matrix multiplication between a and b, square matrices of size n.
stores product into matrix prod of size n.
*/
{
    size_t i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            prod[i][j] = 0;
            for (k = 0; k < n; k++)
            {
                prod[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return 0;
}