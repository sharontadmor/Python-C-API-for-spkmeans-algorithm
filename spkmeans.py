import sys
import os
import numpy as np
import math
import mykmeanssp as km


np.random.seed(0)
MIN_ARGS = 2
MAX_ARGS = 3
KEY = 0
DEFAULT_K = -1
OPERATION_SUCCESSFL = 0
GENERAL_ERROR = 1
ERROR = {
    'GENERAL_ERROR_MESSAGE': 'An Error Has Occurred'
}


def parse_data(file_name):
    """
    return : matrix of given data as Numpy array.
    pre : input is valid - valid range of float variables, same amount of columns in each row,
    all given data points are different.  
    """
    data = np.loadtxt(file_name, delimiter=',')
    return data


def eigengap_heuristic(eigenvalues):
    """
    eigenvalues : sorted array of doubles in ascending order.
    return : k, the number of clusters.
    """
    n = len(eigenvalues)
    vals1 = np.array(eigenvalues[0 : n - 1])
    vals2 = np.array(eigenvalues[1 : n])
    eigengap = np.abs(vals1 - vals2)
    k = np.argmax(eigengap[:math.floor(n / 2)])
    return k + 1


def first_k_eigenvectors(data, k):
    """
    computes the first k eigenvectors of the laplacian matrix of an undirected graph,
    which represents n datapoints given in matrix data,
    the first k eigenvectors are eigenvectors corresponding to the k smallest eigenvalues.
    pre : num of eigenvectors > k
    return : u - matrix whose columns are the first k eigenvectors; k - number of required clusters.
    """
    # get eigenvalues and eigenvectors of the graph laplacian matrix.
    # the eigenvectors are the columns of matrix eigenvectors,
    # the eigenvalues are on the diagonal of matrix eigenvalues.
    data = [data[i].tolist() for i in range(len(data))]    
    eigenvalues, eigenvectors = np.array(km.jacobi(km.gl(data)))
    eigenvalues = np.diag(eigenvalues)
    eigenvectors = np.array(eigenvectors).T
    # sort eigenvalues in ascending order, while rearranging corresponding eigenvectors,
    idx = np.argsort(eigenvalues)  # get the sorted indices
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[idx]
    # get k:
    if k == DEFAULT_K:
        k = eigengap_heuristic(eigenvalues)
    # create a matrix which columns are the first k eigenvectors:
    u = np.array(eigenvectors[:k]).T
    return u, k


def init_centroids(vectors, k):
    """
    initializes k centroids out of given array of vectors.
    choice of centroids is based on weighted probability distribution.
    """
    n = len(vectors)
    centroids = np.array([None] * k)
    initial_idx = np.array([None] * k)
    distances = [None] * n
    # choose first centroid uniformly at random:
    random_vec(vectors, centroids, initial_idx, 0, None)
    # choose the rest of the centroids:
    for i in range(1, k):
        for j in range(n):
            distances[j] = min_dist(centroids, vectors[j])
        # choose next centroid using weighted probability distribution:
        pr = distances / np.sum(distances)
        random_vec(vectors, centroids, initial_idx, i, pr)
    return initial_idx, centroids


def random_vec(vectors, centroids, initial_idx, i, pr):
    """
    chooses a vector from vectors at random, with specified probability.
    """
    if not (pr is None):
        rnd_idx = np.random.choice(vectors.shape[0], p=pr)
    else:
        rnd_idx = np.random.choice(vectors.shape[0])
    centroids[i] = vectors[rnd_idx]
    initial_idx[i] = rnd_idx


def min_dist(lst, vec):
    """
    calculates the distance between given vector and the nearest centroid from given list.
    """
    distances = list(filter(lambda x: x is not None, [
                     eclidean_dist(lst[i], vec) for i in range(len(lst))]))
    return np.amin(distances)


def eclidean_dist(vec1, vec2):
    """
    calculates eclidean distance between vec1 and vec2.
    eclidean distance is exactly the norm of (vec1 - vec2).
    """
    if not (vec1 is None or vec2 is None):
        return np.linalg.norm(vec1 - vec2)


def spk(data, k):
    """
    The Unnormalized Spectral K-means clustering algorithm.
    """
    u, k = first_k_eigenvectors(data, k)
    initial_idx, initial_centroids = init_centroids(u, k)
    u = [u[i].tolist() for i in range(len(u))]
    initial_centroids = [initial_centroids[i].tolist()
                         for i in range(len(initial_centroids))]
    
    final_centroids = km.spk(u, initial_centroids)
    if final_centroids == GENERAL_ERROR:
        return GENERAL_ERROR
    print_centroids(final_centroids, initial_idx, k)
    return OPERATION_SUCCESSFL


def wam(data):
    """
    calculates weighed adjacency matrix.
    the weighed adjacency matrix represents an undirected graph,
    which represents n datapoints given in matrix data,
    and each datapoint is viewed as a vertex.
    return : 0 if run was successful, 1 otherwise.
    """
    data = [data[i].tolist()
                     for i in range(len(data))]
    w = km.wam(data)
    if w == GENERAL_ERROR:
        return GENERAL_ERROR
    print_matrix(w)
    return OPERATION_SUCCESSFL


def ddg(data):
    """
    calculates diagonal degree matrix of an undirected graph,
    which represents n datapoints given in matrix data,
    and each datapoint is viewed as a vertex.
    return : 0 if run was successful, 1 otherwise.
    """
    data = [data[i].tolist()
                     for i in range(len(data))]
    d = km.ddg(data)
    if d == GENERAL_ERROR:
        return GENERAL_ERROR
    print_matrix(d)
    return OPERATION_SUCCESSFL


def gl(data):
    """
    calculates the graph laplacian matrix of an undirected graph,
    which represents n datapoints given in matrix data,
    and each datapoint is viewed as a vertex.
    return : 0 if run was successful, 1 otherwise.
    """
    data = [data[i].tolist()
                     for i in range(len(data))]
    l = km.gl(data)
    if l == GENERAL_ERROR:
        return GENERAL_ERROR
    print_matrix(l)
    return OPERATION_SUCCESSFL


def jacobi(data):
    """
    calculates eigenvalues and eigenvectors of symmetric matrix data.
    return : 0 if run was successful, 1 otherwise.
    """
    data = [data[i].tolist()
                     for i in range(len(data))]
    eigenvalues, eigenvectors = km.jacobi(data)
    if eigenvalues == GENERAL_ERROR or eigenvectors == GENERAL_ERROR:
        return GENERAL_ERROR
    print_jacobi(eigenvalues, eigenvectors)
    return OPERATION_SUCCESSFL


def print_centroids(lst, idx_lst, k):
    """
    print centroids such that the first line has the indices
    of the observations chosen by the K-means++ algorithm as the initial centroids.
    the second line onward has the calculated final centroids.
    """
    d = len(lst[0])
    for i in range(k - 1):
        print(idx_lst[i], end=',')
    print(idx_lst[k - 1])
    for i in range(k):
        for j in range(d - 1):
            print('{:.4f}'.format(lst[i][j]), end=',')
        print('{:.4f}'.format(lst[i][d - 1]))


def print_matrix(x):
    """
    print matrix x of dim n*d.
    """
    n = len(x)
    d = len(x[0])
    for i in range(n):
        for j in range(d - 1):
            print('{:.4f}'.format(x[i][j]), end=',')
        print('{:.4f}'.format(x[i][d - 1]))


def print_jacobi(eigenvalues, eigenvectors):
    """
    eigenvalues are on the diagonal of matrix eigenvalues.
    eigenvectors are the columns of matrix eigenvectors.
    """
    n = len(eigenvalues)
    for i in range(n - 1):
        print('{:.4f}'.format(eigenvalues[i][i]), end=',')
    print('{:.4f}'.format(eigenvalues[n - 1][n - 1]))
    for i in range(n):
        for j in range(n - 1):
            print('{:.4f}'.format(eigenvectors[i][j]), end=',')
        print('{:.4f}'.format(eigenvectors[i][n - 1]))


def operation(goal, data, k):
    if goal == "spk":
        return spk(data, k)
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
        k = DEFAULT_K
    if not os.path.isfile(file_name):
        print(ERROR['GENERAL_ERROR_MESSAGE'])
        return GENERAL_ERROR
    data = parse_data(file_name)        
    res = operation(goal, data, k)
    if res == GENERAL_ERROR:
        print(ERROR['GENERAL_ERROR_MESSAGE'])
        return GENERAL_ERROR
    return OPERATION_SUCCESSFL

if __name__ == "__main__":
    main()
