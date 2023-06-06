#include "spkmeans.h"

/* helper functions for K-means*/
void cleanup(vector **partition, vector **clusters);
static vector *mallocVector(matrix *vectors, int vecIdx);
static void setVectorData(vector *vec, double *coords, vector *next, vector *prev);
static void setMatrixData(matrix *mat, double *arr, int rowsNum, int colsNul);
static void addNode(vector **clusters, int clIdx, vector *vec);
static void removeNode(vector **clusters, vector **partition, int vecIdx);
static void removeHeadVector(vector **clusters, vector *vec);
static int isBreakConditionKmeans(int count, int convergences);
double distance(matrix *vectors1, matrix *vectors2, int vecIdx1, int vecIdx2);
static void assignVectorToCluster(vector **clusters, int clIdx, vector **partition, vector *vec, int vecIdx);
int updateCenrtoids(matrix *centroids, vector **clusters);
/* helper functions for Jacobi */
static double getSumSquares(matrix *a);
static int isBreakConditionJacobi(int count, double offA, double offVecs);
static int getPivot(matrix *a, int *idx);
static int calculateP(matrix *a, matrix *p, int i, int j);
static double getT(matrix *a, int i, int j);
static double getC(double t);
static double getS(double c, double t);
static int rotate(matrix *a, matrix *b, matrix *p, int i, int j);
static int matrixMultiplication(matrix *a, matrix *b, matrix *prod);

void printMatrix(matrix *a)
{
    int i, j;
    for (i = 0; i < a->rowsNum; i++)
    {
        for (j = 0; j < a->colsNum; j++)
        {
            printf("%.4f, ", getCellValue(a, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void cleanup(vector **partition, vector **clusters)
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

static vector *mallocVector(matrix *vectors, int vecIdx)
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

static void addNode(vector **clusters, int clIdx, vector *vec)
/*
adds a node of type Vector to linked list.
*/
{
    vector *curr;
    if (clusters[clIdx] == NULL)
    {
        clusters[clIdx] = vec;
    }
    else
    {
        curr = clusters[clIdx];
        while (curr->next != NULL)
        {
            curr = curr->next;
        }
        curr->next = vec;
        vec->prev = curr;
    }
}

static void removeNode(vector **clusters, vector **partition, int vecIdx)
/*
removes a given node from linked list.
*/
{
    vector *vec = partition[vecIdx];
    if (vec->prev == NULL) /* head vector */
    {
        removeHeadVector(clusters, vec);
    }
    if (vec->prev != NULL)
    {
        vec->prev->next = vec->next;
    }
    if (vec->next != NULL)
    {
        vec->next->prev = vec->prev;
    }
    vec->prev = NULL;
    vec->next = NULL;
}

static void removeHeadVector(vector **clusters, vector *vec)
/*
removes a given node from linked list,
pre: node is the head of linked list, vec->prev == NULL
*/
{
    int i;
    /* if cluster size is 1 */
    if (vec->prev == NULL && vec->next == NULL)
    {
        for (i = 0; i < k; i++)
        {
            if (clusters[i] == vec)
            {
                clusters[i] = NULL;
                return;
            }
        }
    }
    for (i = 0; i < k; i++)
    {
        if (clusters[i] == vec)
        {
            clusters[i] = vec->next;
        }
    }
}

static int isBreakConditionKmeans(int count, int convergences)
{
    if (count == MAX_ITER || convergences == k)
    {
        return 1;
    }
    return 0;
}

double distance(matrix *vectors1, matrix *vectors2, int vecIdx1, int vecIdx2)
/*
calculates eclidean distance between vec1 and vec2.
pre: len(vec1) == len(vec2).
*/
{
    int i;
    double sum, val1, val2;
    sum = 0;
    for (i = 0; i < vectors1->colsNum; i++)
    {
        val1 = getCellValue(vectors1, vecIdx1, i);
        val2 = getCellValue(vectors2, vecIdx2, i);
        sum += (val1 - val2) * (val1 - val2);
    }
    return sqrt(sum);
}

static void assignVectorToCluster(vector **clusters, int clIdx, vector **partition, vector *vec, int vecIdx)
/*
removes vector vec from current cluster,
inserts vec to the cluster in index clIdx in clusters array.
*/
{
    if (partition[vecIdx] == NULL) /* vector is not yet assinged to a cluster */
    {
        addNode(clusters, clIdx, vec);
        partition[vecIdx] = vec;
    }
    else
    {
        removeNode(clusters, partition, vecIdx);
        addNode(clusters, clIdx, partition[vecIdx]);
    }
}

int updateCenrtoids(matrix *centroids, vector **clusters)
/*
iterates over all the clusters and calculates their new centroids.
*/
{
    int i, j, size, convergences;
    double newCoord;
    matrix *oldCnt;
    vector *curr;
    newCoord = 0;
    convergences = 0;
    /* define a temp centroid */
    oldCnt = mallocMatrix(1, centroids->colsNum);
    if (oldCnt == NULL)
    {
        /* cleanup */
        return -1;
    }
    for (i = 0; i < centroids->rowsNum; i++)
    {
        size = 0;
        /* save current centroid, then zero it */
        for (j = 0; j < centroids->colsNum; j++)
        {
            oldCnt->array[j] = getCellValue(centroids, i, j);
            /* setCellValue(oldCnt, 0, j, getCellValue(centroids, i, j)); */
            setCellValue(centroids, i, j, 0);
        }
        /* get sum vector of all the vectors in cluster */
        curr = clusters[i];
        while (curr != NULL)
        {
            for (j = 0; j < centroids->colsNum; j++)
            {
                newCoord = getCellValue(centroids, i, j) + curr->coords[j];
                setCellValue(centroids, i, j, newCoord);
            }
            size++;
            curr = curr->next;
        }
        /* devide sum vector by the size of cluster */
        if (size > 0)
        {
            for (j = 0; j < centroids->colsNum; j++)
            {
                newCoord = getCellValue(centroids, i, j) / size;
                setCellValue(centroids, i, j, newCoord);
            }
        }
        /* for each cluster calculate dist between prev and curr centroids and count convergances */
        if (distance(centroids, oldCnt, i, 0) < KMEANS_EPS)
        {
            convergences++;
        }
    }
    free(oldCnt->array);
    return convergences;
}

int kmeansC(int vectorsNum, int dim, int clustersNum, matrix *vectors, matrix *centroids)
/*
Implementation for K-means clustering algorithm.
arguments:
(1) vectors - an array of all data points that were observed.
(2) centroids - an array of data points chosen as initial centroids.
(3) clustersNum - number of required clusters.
(4) vectorsNum - number of given vectors.
(5) dim - dimention of a vector.
variables:
(1) partition - an array of pointers to vectors assigned to clusters.
each index in partition is appropriate to the index of same vector in vectors.
(2) clusters - an array of vectors, each index holds a linked list of vectors of the same cluster.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, iterations, convergences, assignedCluster;
    double dist, minDist;
    vector **partition, **clusters, *vec;

    iterations = 0;
    convergences = 0;
    k = clustersNum;
    n = vectorsNum;
    d = dim;
    partition = (vector **)malloc(sizeof(vector *) * n);
    clusters = (vector **)malloc(sizeof(vector *) * k);
    if (partition == NULL || clusters == NULL)
    {
        /* cleanup(partition, clusters); */
        return ERROR_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        partition[i] = NULL;
    }
    for (i = 0; i < k; i++)
    {
        clusters[i] = NULL;
    }

    while (!isBreakConditionKmeans(iterations, convergences))
    {
        iterations++;
        /* iterate over all the vectors, sort vectors to clusters */
        for (i = 0; i < n; i++)
        {
            assignedCluster = 0;
            minDist = distance(centroids, vectors, 0, i);
            /* iterate over all the clusters */
            for (j = 0; j < k; j++)
            {
                dist = distance(centroids, vectors, j, i);
                if (dist < minDist)
                {
                    minDist = dist;
                    assignedCluster = j;
                }
            }
            if (partition[i] == NULL) /* vector is not yet assinged to a cluster */
            {
                vec = mallocVector(vectors, i);
                if (vec == NULL)
                {
                    /* cleanup(partition, clusters); */
                    return ERROR_OUT_OF_MEMORY;
                }
                assignVectorToCluster(clusters, assignedCluster, partition, vec, i);
            }
            else
            {
                assignVectorToCluster(clusters, assignedCluster, partition, NULL, i);
            }
        }
        convergences = updateCenrtoids(centroids, clusters);
    }
    /* cleanup(partition, clusters); */
    return 0;
}

int wamC(int n, int d, matrix *vectors, matrix *w)
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
    init w diagonal:
    */
    for (i = 0; i < n; i++)
    {
        setCellValue(w, i, i, 0);
    }
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
                diff = getCellValue(vectors, i, l) - getCellValue(vectors, j, l);
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
            setCellValue(w, i, j, w_ij);
            setCellValue(w, j, i, w_ij);
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

int ddgC(int n, int d, matrix *vectors, matrix *dg)
/*
n - number of vectors.
d - dimention of each vector.
vectors - array representing a matrix of dim n*d of all data points that were observed.
dg - array representing a matrix of dim n*n.
this function calculates dg, the diagonal degree matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    int i, j;
    double degree;
    matrix *w;
    w = mallocMatrix(n, n); /* weighed adjacency matrix */
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
            degree += getCellValue(w, i, j);
        }
        setCellValue(dg, i, i, degree);
    }
    return 0;
}

int glC(int n, int d, matrix *vectors, matrix *l)
/*
n - number of vectors.
d - dimention of each vector.
vectors - array representing a matrix of dim n*d of all data points that were observed.
l - array representing a matrix of dim n*n.
this function calculates l, the graph laplacian matrix of the graph created from data points in vectors.
returns 0 if operation is successful.
*/
{
    int i, j;
    double diff;
    matrix *w, *dg;
    w = mallocMatrix(n, n); /* weighed adjacency matrix */
    dg = mallocMatrix(n, n); /* diagonal degree matrix */
    if (w == NULL || dg == NULL)
    {
        /* cleanup */
        return ERROR_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            setCellValue(w, i, j, 0);
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
    wamC(n, d, vectors, w);
    ddgC(n, d, vectors, dg);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            diff = getCellValue(dg, i, j) - getCellValue(w, i, j);
            setCellValue(l, i, j, diff);
        }
    }
    return 0;
}

int jacobiC(int n, matrix *a, matrix *eigenvals, matrix *eigenvecs)
/*
the Jacobi eigenvalue algorithm is an iterative method
for the calculation of eigenvalues and eigenvectors of a real symmetric matrix.
n - size of matrix a.
a - array representing a symmetric matrix of dim n*n,
for which to calculate eigenvalues and eigenvectors.
eigenvals - array representing a matrix of dim n*n,
in which the diagonal values are a's eigenvalues.
eigenvecs - array representing a matrix of dim n*n,
in which the columns are a's eigenvectors.
returns 0 if operation is successful.
*/
{
    int i, j, count, idx[2], res;
    /* int count; */
    double offA, offVecs;
    matrix *p, *temp, *prod;
    idx[0] = 0;
    idx[1] = 0;
    offA = 0;
    offVecs = 0;
    count = 0;
    p = mallocMatrix(n, n);
    if (p == NULL)
    {
        /* cleanup */
        return ERROR_OUT_OF_MEMORY;
    }
    offA = getSumSquares(a);

    /*
    printf("a:\n");
    printMatrix(n, n, a);
    */
    while (!isBreakConditionJacobi(count, offA, offVecs))
    {
        count += 1;

        /*
        printf("---> iteration: %d\n", count);
        printf("p 1 - start of iter:\n");
        printMatrix(n, n, p);
        printf("eigenvalues 1 - start of iter:\n");
        printMatrix(n, n, eigenvalues);
        */
        /*
        printf("---> iteration: %d\n", count);
        printf("p 1 - start of iter:\n");
        printMatrix(n, n, p);
        printf("eigenvectors 1 - start of iter:\n");
        printMatrix(n, n, eigenvectors);
        */

        res = getPivot(eigenvals, idx);
        if (res == MATRIX_IS_DIAGONAL)
        {
            if (count == 1)
            {
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        setCellValue(eigenvecs, i, j, getCellValue(p, i, j));
                    }
                }
            }
            else
            {
                prod = mallocMatrix(n, n);
                if (prod == NULL)
                {
                    /* cleanup */
                    return ERROR_OUT_OF_MEMORY;
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        setCellValue(prod, i, j, 0);
                    }
                }
                matrixMultiplication(eigenvecs, p, prod);
                free(eigenvecs);
                eigenvecs = prod;
            }
            return 0;
        }
        temp = mallocMatrix(n, n);
        if (temp == NULL)
        {
            /* cleanup */
            return ERROR_OUT_OF_MEMORY;
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                setCellValue(temp, i, j, getCellValue(eigenvals, i, j));
                if (i == j)
                {
                    setCellValue(p, i, j, 1);
                }
                else
                {
                    setCellValue(p, i, j, 0);
                }
            }
        }

        /*
        printf("p 2 - after p init:\n");
        printMatrix(n, n, p);
        printf("eigenvalues 2 - after p init:\n");
        printMatrix(n, n, eigenvalues);
        */

        calculateP(eigenvals, p, idx[0], idx[1]);
        rotate(eigenvals, temp, p, idx[0], idx[1]);
        free(eigenvals);
        eigenvals = temp;

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
                    setCellValue(eigenvecs, i, j, getCellValue(p, i, j));
                }
            }
        }
        else
        {
            prod = mallocMatrix(n, n);
            if (prod == NULL)
            {
                /* cleanup */
                return ERROR_OUT_OF_MEMORY;
            }
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    setCellValue(prod, i, j, 0);
                }
            }
            matrixMultiplication(eigenvecs, p, prod);
            free(eigenvecs);
            eigenvecs = prod;
        }
        offVecs = getSumSquares(eigenvals);
        /*
        printf("p 4 - end of iter:\n");
        printMatrix(n, n, p);
        printf("eigenvectors 4 - end of iter:\n");
        printMatrix(n, n, eigenvectors);
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

static double getSumSquares(matrix *a)
/*
n - size of matrix a.
a - array representing symmetric matrix of dim n*n.
returns the sum of squares of all off-diagonal elements of matrix a.
*/
{
    int i, j, n;
    double sum;
    n = a->rowsNum;
    sum = 0;
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            sum += getCellValue(a, i, j) * getCellValue(a, i, j);
        }
    }
    return 2.0 * sum;
}

static int isBreakConditionJacobi(int count, double offA, double offVecs)
{
    int convergence;
    convergence = offA - offVecs <= EPS;
    /*
    printf("offA: %.100f\n", offA);
    printf("offVecs: %.100f\n", offVecs);
    printf("diff: %.100f\n", offA - offVecs);
    printf("EPS: %.100f\n", EPS);
    */
    if (count == MAX_ROTATIONS || convergence)
    {
        return 1;
    }
    return 0;
}

static int getPivot(matrix *a, int *idx)
/*
n - size of matrix a.
a - array representing symmetric matrix of dim n*n.
idx - array of size 2.
this function finds pivot a_ij,
and puts its row and column in array idx, such that the 0 index is the row i of a_ij, and the index 1 is the column j of a_ij.
pivot is the off-diagonal element of matrix a with the largest absolute value.
*/
{
    int i, j, n;
    double pivot;
    n = a->rowsNum;
    pivot = 0;
    /* printMatrix(n, n, a); */
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(getCellValue(a, i, j)) > pivot)
            {
                pivot = fabs(getCellValue(a, i, j));
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

static int calculateP(matrix *a, matrix *p, int i, int j)
/*
calculates rotation matrix p.
*/
{
    double t, c, s;
    t = getT(a, i, j);
    c = getC(t);
    s = getS(c, t);
    setCellValue(p, i, i, c);
    setCellValue(p, j, j, c);
    setCellValue(p, i, j, s);
    setCellValue(p, j, i, (-1) * s);

    /*
    printf("---> t: %f\n", t);
    printf("---> c: %f\n", c);
    printf("---> s: %f\n", s);
    */

    return 0;
}

static double getT(matrix *a, int i, int j)
/*
calculates value t out of matrix a.
*/
{
    double theta, denominator, t;
    theta = 0;
    theta = (getCellValue(a, j, j) - getCellValue(a, i, i)) / (2 * getCellValue(a, i, j));
    /* theta = (a[j][j] - a[i][i]) / (2 * a[i][j]); */
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

static int rotate(matrix *a, matrix *b, matrix *p, int i, int j)
/*
rotates matrix b to make b_ij = 0.
*/
{
    int r, n;
    double c, s, val;
    n = a->rowsNum;
    c = getCellValue(p, i, i);
    s = getCellValue(p, i, j);
    val = (c * c * getCellValue(a, i, i)) + (s * s * getCellValue(a, j, j)) - (2 * s * c * getCellValue(a, i, j));
    setCellValue(b, i, i, val);
    val = (s * s * getCellValue(a, i, i)) + (c * c * getCellValue(a, j, j)) + (2 * s * c * getCellValue(a, i, j));
    setCellValue(b, i, j, 0); /* b[i][j] = (c * c - s * s) * a[i][j] + (s * c)(a[i][i] - a[j][j]) = 0 */
    setCellValue(b, j, i, 0);
    for (r = 0; r < n; r++)
    {
        if ((r != i) && (r != j))
        {
            val = c * getCellValue(a, r, i) - s * getCellValue(a, r, j);
            setCellValue(b, r, i, val);
            setCellValue(b, i, r, val);
            val = c * getCellValue(a, r, j) + s * getCellValue(a, r, i);
            setCellValue(b, r, j, val);
            setCellValue(b, j, r, val);
        }
    }
    /*
    c = p[i][i];
    s = p[i][j];
    b[i][i] = (c * c * a[i][i]) + (s * s * a[j][j]) - (2 * s * c * a[i][j]);
    b[j][j] = (s * s * a[i][i]) + (c * c * a[j][j]) + (2 * s * c * a[i][j]);
    b[i][j] = 0;
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
    */

    return 0;
}

static int matrixMultiplication(matrix *a, matrix *b, matrix *prod)
/*
calculates matrix multiplication between a and b, square matrices of size n.
stores product into matrix prod of size n.
*/
{
    int i, j, k, n;
    double res;
    n = a->rowsNum;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            setCellValue(prod, i, j, 0);
            for (k = 0; k < n; k++)
            {
                res = getCellValue(prod, i, j) + getCellValue(a, i, k) * getCellValue(b, k, j);
                setCellValue(prod, i, j, res);
            }
        }
    }
    return 0;
}