import sys
import os
import numpy as np
import pandas as pd
# import mykmeanssp as km
import mykmeanssp_python as km


"""
run using:
python3 spkmeans.py 3 spk tests/input_1.txt
"""


# seed
MIN_ARGS = 2
MAX_ARGS = 3
KEY = 0


def parse_data(file_name):
    """
    returns matrix of given data as Numpy array.
    """
    data = np.loadtxt(file_name, delimiter=',')
    return data


def default_k():
    pass


def spk():
    return km.spk()


def wam():
    return km.wam()


def ddg():
    return km.ddg()


def gl():
    return km.gl()


def jacobi():
    return km.jacobi()


def operation(goal, data):
    if goal == "spk":
        return spk()
    if goal == "wam":
        return wam()
    if goal == "ddg":
        return ddg()
    if goal == "gl":
        return gl()
    if goal == "jacobi":
        return jacobi()


def main():
    args = sys.argv[1:]
    if len(args) == MAX_ARGS:
        k, goal, file_name = args
        k = int(k)
    else:
        goal, file_name = args
        k = default_k()
    data = parse_data(file_name)
    res = operation(goal, data)
    print(data)


if __name__ == "__main__":
    main()
