#include "spkmeans.h"

/*
compile using:
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans_main.c -o spkmeans_main -lm
run using:
make && ./spkmeans wam tests/test_batch/test1.txt
make && ./spkmeans jacobi tests/test_batch/test1_j.txt

test with asan:
make spkmeans_asan && ./spkmeans_asan jacobi tests/test_batch/test1_j.txt

test for case -0:
make && ./spkmeans jacobi tests/test_batch/test0_j.txt
*/

static matrix *parseData(char *fileName);
static int operation(char *goal, matrix *data);
static int wam(matrix *data);
static int ddg(matrix *data);
static int gl(matrix *data);
static int jacobi(matrix *data);

static matrix *parseData(char *fileName)
/*
reads data from file with directory fileName.
pre : input is valid - valid range of double variables, same amount of columns in each row,
all given data points are different.
return : matrix of given data if run was successful, NULL otherwise.
*/
{
    int i, j, c, rowsNum, colsNum, error;
    double coord;
    matrix *data;
    FILE *ifp;
    c = 0;
    rowsNum = 0;
    colsNum = 1;
    coord = 0;
    ifp = fopen(fileName, "r");
    if (ifp == NULL)
    {
        return NULL;
    }
    /* count length of input rows and columns: */
    while ((c = fgetc(ifp)) != EOF)
    {
        if (c == '\n')
        {
            rowsNum++;
        }
        if (c == ',' && rowsNum == 0)
        {
            colsNum++;
        }
    }
    rewind(ifp);
    data = mallocMatrix(rowsNum, colsNum);
    if (data == NULL)
    {
        fclose(ifp);
        return NULL;
    }
    for (i = 0; i < rowsNum; i++)
    {
        for (j = 0; j < colsNum; j++)
        {
            error = fscanf(ifp, "%lf", &coord);
            if (error == EOF)
            {
                matrixCleanup(data);
                fclose(ifp);
                return NULL;
            }
            setCellValue(data, i, j, coord);
            if (fgetc(ifp) == EOF)
            {
                break;
            }
        }
    }
    fclose(ifp);
    return data;
}

static int operation(char *goal, matrix *data)
{
    if (strcmp(goal, "wam") == 0)
    {
        return wam(data);
    }
    if (strcmp(goal, "ddg") == 0)
    {
        return ddg(data);
    }
    if (strcmp(goal, "gl") == 0)
    {
        return gl(data);
    }
    if (strcmp(goal, "jacobi") == 0)
    {
        return jacobi(data);
    }
    return EXIT_SUCCESS;
}

static int wam(matrix *data)
/*
calculates weighed adjacency matrix.
the weighed adjacency matrix represents an undirected graph,
which represents n datapoints given in matrix data,
and each datapoint is viewed as a vertex.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    matrix *w;
    w = mallocMatrix(data->rowsNum, data->rowsNum); /* weighed adjacency matrix */
    if (w == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < data->rowsNum; i++)
    {
        for (j = 0; j < data->rowsNum; j++)
        {
            setCellValue(w, i, j, 0);
        }
    }
    res = wamC(data->rowsNum, data->colsNum, data, w);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(w);
        return EXIT_FAILURE;
    }
    printMatrix(w);
    matrixCleanup(w);
    return EXIT_SUCCESS;
}

static int ddg(matrix *data)
/*
calculates diagonal degree matrix of an undirected graph,
which represents n datapoints given in matrix data,
and each datapoint is viewed as a vertex.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    matrix *dg;
    dg = mallocMatrix(data->rowsNum, data->rowsNum); /* diagonal degree matrix */
    if (dg == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < data->rowsNum; i++)
    {
        for (j = 0; j < data->rowsNum; j++)
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
    res = ddgC(data->rowsNum, data->colsNum, data, dg);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(dg);
        return EXIT_FAILURE;
    }
    printMatrix(dg);
    matrixCleanup(dg);
    return EXIT_SUCCESS;
}

static int gl(matrix *data)
/*
calculates the graph laplacian matrix of an undirected graph,
which represents n datapoints given in matrix data,
and each datapoint is viewed as a vertex.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    matrix *l;
    l = mallocMatrix(data->rowsNum, data->rowsNum); /* graph laplacian matrix */
    if (l == NULL)
    {
        return EXIT_FAILURE;
    }
    for (i = 0; i < data->rowsNum; i++)
    {
        for (j = 0; j < data->rowsNum; j++)
        {
            setCellValue(l, i, j, 0);
        }
    }
    res = glC(data->rowsNum, data->colsNum, data, l);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(l);
        return EXIT_FAILURE;
    }
    printMatrix(l);
    matrixCleanup(l);
    return EXIT_SUCCESS;
}

static int jacobi(matrix *data)
/*
calculates eigenvalues and eigenvectors of symmetric matrix data.
return : 0 if run was successful, 1 otherwise.
*/
{
    int i, j, res;
    matrix *eigenvals, *eigenvecs;
    eigenvals = mallocMatrix(data->rowsNum, data->rowsNum);
    eigenvecs = mallocMatrix(data->rowsNum, data->rowsNum);
    if (eigenvals == NULL || eigenvecs == NULL)
    {
        matrixCleanup(eigenvals);
        matrixCleanup(eigenvecs);
        return EXIT_FAILURE;
    }
    for (i = 0; i < data->rowsNum; i++)
    {
        for (j = 0; j < data->rowsNum; j++)
        {
            setCellValue(eigenvals, i, j, getCellValue(data, i, j));
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
    res = jacobiC(data->rowsNum, data, eigenvals, eigenvecs);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(eigenvals);
        matrixCleanup(eigenvecs);
        return EXIT_FAILURE;
    }
    printJacobi(eigenvals, eigenvecs);
    matrixCleanup(eigenvals);
    matrixCleanup(eigenvecs);
    return EXIT_SUCCESS;
}

void printMatrix(matrix *a)
/*
prints given matrix a.
*/
{
    int i, j;
    for (i = 0; i < a->rowsNum; i++)
    {
        for (j = 0; j < a->colsNum - 1; j++)
        {
            printf("%.4f,", getCellValue(a, i, j));
        }
        printf("%.4f\n", getCellValue(a, i, j));
    }
}

void printJacobi(matrix *eigenvals, matrix *eigenvecs)
/*
eigenvalues are on the diagonal of matrix eigenvalues.
eigenvectors are the columns of matrix eigenvectors.
*/
{
    int i, j, n;
    n = eigenvals->rowsNum;
    for (i = 0; i < n - 1; i++)
    {
        printf("%.4f,", getCellValue(eigenvals, i, i));
    }
    printf("%.4f\n", getCellValue(eigenvals, n - 1, n - 1));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
        {
            printf("%.4f,", getCellValue(eigenvecs, i, j));
        }
        printf("%.4f\n", getCellValue(eigenvecs, i, j));
    }
}

int main(int argc, char *argv[])
{
    int res;
    matrix *data;
    if (argc != ARGS_NUM)
    {
        fprintf(stderr, GENERAL_ERROR_MESSAGE);
        return EXIT_FAILURE;
    }
    data = parseData(argv[FILENAME_IDX]);
    if (data == NULL)
    {
        fprintf(stderr, GENERAL_ERROR_MESSAGE);
        return EXIT_FAILURE;
    }
    res = operation(argv[GOAL_IDX], data);
    if (res == EXIT_FAILURE)
    {
        matrixCleanup(data);
        fprintf(stderr, GENERAL_ERROR_MESSAGE);
        return EXIT_FAILURE;
    }
    matrixCleanup(data);
    return res;
}
