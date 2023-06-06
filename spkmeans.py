import sys
import os
import numpy as np
# import pandas as pd
import mykmeanssp as km
# import mykmeanssp_python as km


"""
build using:
python3 setup.py build_ext --inplace

run using:
python3 spkmeans.py 3 spk tests/input_1.txt
python3 spkmeans.py 3 wam tests/input_smallest.txt
python3 spkmeans.py 3 jacobi tests/jacobi_0.txt
"""


# seed
MIN_ARGS = 2
MAX_ARGS = 3
KEY = 0
ERROR_OUT_OF_MEMORY = 1
ERROR = {
    'GENERAL_ERROR_MESSAGE': 'An Error Has Occurred'
}


def parse_data(file_name):
    """
    returns matrix of given data as Numpy array.
    """
    data = np.loadtxt(file_name, delimiter=',')
    return data


def default_k():
    pass


def spk(data):
    return km.spk(data)


def wam(data):
    """
    return: weighed adjacency matrix, or None if operation failed.
    """
    w = km.wam(data)
    # if w == ERROR_OUT_OF_MEMORY:
    #     print(ERROR['GENERAL_ERROR_MESSAGE'])
    #     return
    print_matrix(w)
    return w


def ddg(data):
    """
    return: diagonal degree matrix, or None if operation failed.
    """
    d = km.ddg(data)
    # if d == ERROR_OUT_OF_MEMORY:
    #     print(ERROR['GENERAL_ERROR_MESSAGE'])
    #     return
    print_matrix(d)
    return d


def gl(data):
    """
    return: the graph laplacian matrix, or None if operation failed.
    """
    l = km.gl(data)
    # if l == ERROR_OUT_OF_MEMORY:
    #     print(ERROR['GENERAL_ERROR_MESSAGE'])
    #     return
    print_matrix(l)
    return l


def jacobi(data):
    jac = km.jacobi(data)
    # if eigenvalues == ERROR_OUT_OF_MEMORY or eigenvectors == ERROR_OUT_OF_MEMORY:
    #     print(ERROR['GENERAL_ERROR_MESSAGE'])
    #     return
    print_jacobi(jac, jac)
    return jac


def print_centroids(lst, idx_lst, k):
    d = len(lst[0])
    for i in range(k - 1):
        print(idx_lst[i], end=',')
    print(idx_lst[k - 1])
    for i in range(k):
        for j in range(d - 1):
            print('{:.4f}'.format(lst[i][j]), end=',')
        print('{:.4f}'.format(lst[i][d - 1]))


def print_matrix(x):
    print(x)


def print_jacobi(eigenvalues, eigenvectors):
    # n = len(eigenvalues)
    # for i in range(n - 1):
    #     print('{:.4f}'.format(eigenvalues[i, i]), end=',')
    # print('{:.4f}'.format(eigenvalues[n - 1, n - 1]))

    print(eigenvalues)
    print(eigenvectors)


def operation(goal, data):
    if goal == "spk":
        return spk(data)
    if goal == "wam":
        return wam(data)
    if goal == "ddg":
        return ddg(data)
    if goal == "gl":
        return gl(data)
    if goal == "jacobi":
        return jacobi(data)


def main():
    args = sys.argv[1:]
    if len(args) == MAX_ARGS:
        k, goal, file_name = args
        k = int(k)
    else:
        goal, file_name = args
        k = default_k()
    data = parse_data(file_name)
    data = [data[i].tolist()
                     for i in range(len(data))]
    res = operation(goal, data)


if __name__ == "__main__":
    main()
