/* #include <float.h> */
#include "spkmeans.h"

/* helper functions for K-means*/
static void addNode(vector **clusters, int clIdx, vector *vec);
static void removeNode(vector **clusters, vector **partition, int vecIdx, int clustersNum);
static void removeHeadVector(vector **clusters, vector *vec, int clustersNum);
static void assignVectorToCluster(vector **clusters, int clIdx, vector **partition, vector *vec, int vecIdx, int clustersNum);
static int isBreakConditionKmeans(int count, int convergences, int clustersNum);
static double distance(matrix *vectors1, matrix *vectors2, int vecIdx1, int vecIdx2);
static int updateCenrtoids(matrix *centroids, vector **clusters);
/* helper functions for Jacobi */
static double getSumSquares(matrix *a);
static int isBreakConditionJacobi(int count, int isDiagonal);
static int getPivot(matrix *a, int *idx);
static void buildRotationMatrix(matrix *p, double s, double c, int idx1, int idx2);
static double getT(matrix *a, int i, int j);
static double getC(double t);
static double getS(double c, double t);
static int calculateEigenvectors(matrix *eigenvecs, matrix *p);
static int matrixMultiplication(matrix *a, matrix *b, matrix *prod);
static int calculateEigenvalues(matrix *eigenvals, double s, double c, int idx1, int idx2);
static int rotate(matrix *a, matrix *b, double s, double c, int i, int j);
static int correctMinusZero(matrix *eigenvecs, matrix *eigenvals);

static void addNode(vector **clusters, int clIdx, vector *vec)
/*
adds a node of type Vector to linked list.
adds to the end of the list.
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
        vec->next = NULL;
    }
}

static void removeNode(vector **clusters, vector **partition, int vecIdx, int clustersNum)
/*
removes a given node from linked list.
*/
{
    vector *vec = partition[vecIdx];
    if (vec->prev == NULL) /* head vector */
    {
        removeHeadVector(clusters, vec, clustersNum);
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

static void removeHeadVector(vector **clusters, vector *vec, int clustersNum)
/*
removes a given node from linked list,
pre: node is the head of linked list, vec->prev == NULL
*/
{
    int i;
    for (i = 0; i < clustersNum; i++)
    {
        if (clusters[i] == vec)
        {
            clusters[i] = vec->next;
        }
    }
}

static void assignVectorToCluster(vector **clusters, int clIdx, vector **partition, vector *vec, int vecIdx, int clustersNum)
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
        removeNode(clusters, partition, vecIdx, clustersNum);
        addNode(clusters, clIdx, partition[vecIdx]);
    }
}

static int isBreakConditionKmeans(int count, int convergences, int clustersNum)
{
    if (count == MAX_ITER || convergences == clustersNum)
    {
        return 1;
    }
    return 0;
}

static double distance(matrix *vectors1, matrix *vectors2, int vecIdx1, int vecIdx2)
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

static int updateCenrtoids(matrix *centroids, vector **clusters)
/*
iterates over all the clusters and calculates their new centroids.
return : convergences - number of centroids whose new centroid has a distance of less than epsilon from previus centroid.
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
        return -1;
    }
    for (i = 0; i < centroids->rowsNum; i++)
    {
        size = 0;
        /* save current centroid, then zero it */
        for (j = 0; j < centroids->colsNum; j++)
        {
            oldCnt->array[j] = getCellValue(centroids, i, j);
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
        /* divide sum vector by the size of cluster */
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
    matrixCleanup(oldCnt);
    return convergences;
}

int kmeansC(int vectorsNum, int clustersNum, matrix *vectors, matrix *centroids)
/*
Implementation for K-means clustering algorithm.
arguments:
(1) vectors - an array of all data points that were observed.
(2) centroids - an array of data points chosen as initial centroids.
(3) clustersNum - number of required clusters.
(4) vectorsNum - number of given vectors.
variables:
(1) partition - an array of pointers to vectors assigned to clusters.
each index in partition is appropriate to the index of same vector in vectors.
(2) clusters - an array of vectors, each index holds a linked list of vectors of the same cluster.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, k, n, iterations, convergences, assignedCluster;
    double dist, minDist;
    vector **partition, **clusters, *vec;

    iterations = 0;
    convergences = 0;
    k = clustersNum;
    n = vectorsNum;
    partition = (vector **)malloc(sizeof(vector *) * n);
    if (partition == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < n; i++)
    {
        partition[i] = NULL;
    }
    clusters = (vector **)malloc(sizeof(vector *) * k);
    if (clusters == NULL)
    {
        free(partition);
        return EXIT_FAILURE;
    }
    for (i = 0; i < k; i++)
    {
        clusters[i] = NULL;
    }

    while (!isBreakConditionKmeans(iterations, convergences, k))
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
                    vectorCleanup(partition, clusters, n);
                    return EXIT_FAILURE;
                }
                assignVectorToCluster(clusters, assignedCluster, partition, vec, i, k);
            }
            else
            {
                assignVectorToCluster(clusters, assignedCluster, partition, NULL, i, k);
            }
        }
        convergences = updateCenrtoids(centroids, clusters);
    }
    vectorCleanup(partition, clusters, n);
    return EXIT_SUCCESS;
}

int wamC(int n, int d, matrix *vectors, matrix *w)
/*
calculates w, the weighed adjacency matrix of the graph created from data points in vectors.
arguments:
(1) n - number of vectors.
(2) d - dimention of each vector.
(3) vectors - array representing a matrix of dim n*d of all data points that were observed.
(4) w - array representing a matrix of dim n*n.
return : 0 if operation is successful.
*/
{
    int i, j, l;
    double euclid_dist, w_ij, diff;
    /* init w diagonal */
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
            }
            w_ij = exp((-1.0 / 2) * euclid_dist);
            setCellValue(w, i, j, w_ij);
            setCellValue(w, j, i, w_ij);
        }
    }
    return EXIT_SUCCESS;
}

int ddgC(int n, int d, matrix *vectors, matrix *dg)
/*
calculates dg, the diagonal degree matrix of the graph created from data points in vectors.
arguments:
(1) n - number of vectors.
(2) d - dimention of each vector.
(3) vectors - array representing a matrix of dim n*d of all data points that were observed.
(4) dg - array representing a matrix of dim n*n.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    double degree;
    matrix *w;
    w = mallocMatrix(n, n); /* weighed adjacency matrix */
    if (w == NULL)
    {
        return EXIT_FAILURE;
    }
    res = wamC(n, d, vectors, w);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(w);
        return EXIT_FAILURE;
    }
    for (i = 0; i < n; i++)
    {
        degree = 0;
        for (j = 0; j < n; j++)
        {
            degree += getCellValue(w, i, j);
        }
        setCellValue(dg, i, i, degree);
    }
    matrixCleanup(w);
    return EXIT_SUCCESS;
}

int glC(int n, int d, matrix *vectors, matrix *l)
/*
calculates l, the graph laplacian matrix of the graph created from data points in vectors.
arguments:
(1) n - number of vectors.
(2) d - dimention of each vector.
(3) vectors - array representing a matrix of dim n*d of all data points that were observed.
(4) l - array representing a matrix of dim n*n.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    double diff;
    matrix *w, *dg;
    w = mallocMatrix(n, n);  /* weighed adjacency matrix */
    dg = mallocMatrix(n, n); /* diagonal degree matrix */
    if (w == NULL || dg == NULL)
    {
        matrixCleanup(w);
        matrixCleanup(dg);
        return EXIT_FAILURE;
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
    res = wamC(n, d, vectors, w);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(w);
        matrixCleanup(dg);
        return EXIT_FAILURE;
    }
    res = ddgC(n, d, vectors, dg);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(w);
        matrixCleanup(dg);
        return EXIT_FAILURE;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            diff = getCellValue(dg, i, j) - getCellValue(w, i, j);
            setCellValue(l, i, j, diff);
        }
    }
    matrixCleanup(w);
    matrixCleanup(dg);
    return EXIT_SUCCESS;
}

int jacobiC(int n, matrix *a, matrix *eigenvals, matrix *eigenvecs)
/*
the Jacobi eigenvalue algorithm is an iterative method
for the calculation of eigenvalues and eigenvectors of a real symmetric matrix.
arguments:
(1) n - size of matrix a.
(2) a - array representing a symmetric matrix of dim n*n,
for which to calculate eigenvalues and eigenvectors.
(3) eigenvals - array representing a matrix of dim n*n,
in which the diagonal values are a's eigenvalues.
(4) eigenvecs - array representing a matrix of dim n*n,
in which the columns are a's eigenvectors.
return : 0 if run was successful, 1 otherwise.
*/
{
    int isDiagonal, res, count, pivotIndxes[2];
    double prevOff, currOff, t, c, s;
    matrix *p;
    pivotIndxes[0] = 0;
    pivotIndxes[1] = 0;
    isDiagonal = 0;
    prevOff = 0;
    prevOff = getSumSquares(a);
    count = 0;
    p = mallocMatrix(n, n);
    if (p == NULL)
    {
        return EXIT_FAILURE;
    }
    while (!isBreakConditionJacobi(count, isDiagonal))
    {
        count += 1;
        res = getPivot(eigenvals, pivotIndxes);
        t = getT(eigenvals, pivotIndxes[0], pivotIndxes[1]);
        c = getC(t);
        s = getS(c, t);
        buildRotationMatrix(p, s, c, pivotIndxes[0], pivotIndxes[1]);
        res = calculateEigenvectors(eigenvecs, p);
        if (res == EXIT_FAILURE)
        {
            matrixCleanup(p);
            return EXIT_FAILURE;
        }
        res = calculateEigenvalues(eigenvals, s, c, pivotIndxes[0], pivotIndxes[1]);
        if (res == EXIT_FAILURE)
        {
            matrixCleanup(p);
            return EXIT_FAILURE;
        }
        currOff = getSumSquares(eigenvals);
        if ((prevOff - currOff) <= EPS) /* the difference should always be a positive double */
        {
            isDiagonal = 1;
        }
        prevOff = currOff;
    }
    correctMinusZero(eigenvecs, eigenvals);
    matrixCleanup(p);
    return EXIT_SUCCESS;
}

static double getSumSquares(matrix *a)
/*
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

static int isBreakConditionJacobi(int count, int isDiagonal)
{
    if (count == MAX_ROTATIONS || isDiagonal)
    {
        return 1;
    }
    return 0;
}

static int getPivot(matrix *a, int *idx)
/*
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
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(getCellValue(a, i, j)) > pivot)
            {
                pivot = fabs(getCellValue(a, i, j));
                idx[0] = i;
                idx[1] = j;
            }
        }
    }
    return EXIT_SUCCESS;
}

static void buildRotationMatrix(matrix *p, double s, double c, int idx1, int idx2)
/*
calculates rotation matrix p.
*/
{
    int i, j;
    /* initiate p to be a unit matrix: */
    for (i = 0; i < p->rowsNum; i++)
    {
        for (j = 0; j < p->colsNum; j++)
        {
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
    /* factoring in the values of s, c: */
    setCellValue(p, idx1, idx1, c);
    setCellValue(p, idx2, idx2, c);
    setCellValue(p, idx1, idx2, s);
    setCellValue(p, idx2, idx1, (-1) * s);
}

static double getT(matrix *a, int i, int j)
/*
calculates value t out of matrix a.
*/
{
    double theta, denominator, t;
    theta = 0;
    theta = (getCellValue(a, j, j) - getCellValue(a, i, i)) / (2 * getCellValue(a, i, j));
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

static int calculateEigenvectors(matrix *eigenvecs, matrix *p)
/*
updates matrix eigenvecs by performing matrix multiplication with p.
*/
{
    int i, j, n;
    matrix *prod;
    n = eigenvecs->rowsNum;
    prod = mallocMatrix(n, n);
    if (prod == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            setCellValue(prod, i, j, 0);
        }
    }
    matrixMultiplication(eigenvecs, p, prod);
    free(eigenvecs->array);
    eigenvecs->array = prod->array;
    free(prod);
    return EXIT_SUCCESS;
}

static int matrixMultiplication(matrix *a, matrix *b, matrix *prod)
/*
calculates matrix multiplication between a and b.
pre : a and b are square matrices of size n.
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
    return EXIT_SUCCESS;
}

static int calculateEigenvalues(matrix *eigenvals, double s, double c, int idx1, int idx2)
/*
updates matrix eigenvals.
*/
{
    int i, j, n;
    matrix *temp;
    n = eigenvals->rowsNum;
    /* init a copy of eigenvals */
    temp = mallocMatrix(n, n);
    if (temp == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            setCellValue(temp, i, j, getCellValue(eigenvals, i, j));
        }
    }
    /* update eigenvals */
    rotate(eigenvals, temp, s, c, idx1, idx2);
    free(eigenvals->array);
    eigenvals->array = temp->array;
    free(temp);
    return EXIT_SUCCESS;
}

static int rotate(matrix *a, matrix *b, double s, double c, int i, int j)
/*
rotates matrix b to make b_ij = 0.
*/
{
    int r, n;
    double val;
    n = a->rowsNum;
    val = (c * c * getCellValue(a, i, i)) + (s * s * getCellValue(a, j, j)) - (2 * s * c * getCellValue(a, i, j));
    setCellValue(b, i, i, val);
    val = (s * s * getCellValue(a, i, i)) + (c * c * getCellValue(a, j, j)) + (2 * s * c * getCellValue(a, i, j));
    setCellValue(b, j, j, val);
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
    return EXIT_SUCCESS;
}

static int correctMinusZero(matrix *eigenvecs, matrix *eigenvals)
/*
this function loops over the diagonal of matrix eigenvals.
in the edge case that a minus zero (double) eigenvalue was calculated due to floating point precision,
this function will correct its corresponding eigenvector by multiplying the corresponding column of matrix eigenvecs by -1.
*/
{
    int i, j;
    for (i = 0; i < eigenvals->rowsNum; i++)
    {
        if (getCellValue(eigenvals, i, i) == 0 && (1 / (getCellValue(eigenvals, i, i))) != (1 / 0.0))
        /* || getCellValue(eigenvals, i, i) == -DBL_EPSILON */
        {
            for (j = 0; j < eigenvals->rowsNum; j++)
            {
                setCellValue(eigenvecs, j, i, (-1) * getCellValue(eigenvecs, j, i));
            }
        }
    }
    return EXIT_SUCCESS;
}