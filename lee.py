import random

import usefull as use
import numpy as np
import echelon as ech


def isExtract(H, s, n, r):
    """
    Initial part from Prange
    """
    limitWhile = min(r*n, 200)
    while True:
        if limitWhile == 0:
            return -1, 0, 0, 0
        Q = use.RandomPermutation(n)
        new_H = np.matmul(H, Q)
        P, V = ech.ToReducedRowEchelonForm(new_H.tolist())
        if use.isIr(np.array(V), r, n):
            break
        limitWhile -= 1
    return 0, Q, V, use.array_base2(np.matmul(P, s))


def sum_columns_from_V(I, VIr):
    """
    Add columns from V with indexes j in I
    :param I: column indexes to sum
    :param VIr: matrix from where we add columns
    :return: sum of j in I columns in base 2
    """
    columns_sum = use.extract_column(VIr, I[0])
    for i in range(1, len(I)):
        columns_sum = use.add_column_arrays(columns_sum, use.extract_column(VIr, I[i]))
    return use.array_base2(columns_sum)


def LeeAlgorithm(H, s, t, r, n):
    """
    Lee and Brickell algorithm
    :param H: Band matrix
    :param s: syndrome a.i. s = Hx
    :param t: weight
    :param r: number of rows
    :param n: number of columns
    :return: reconstructed initial vector x, or [] if fails
    """
    k = n - r
    limit = n
    while True:
        if limit == 0:
            return []
        ok, Q, VIr, s_p = isExtract(H, s, n, r)
        p = 2
        new_e = []
        if ok == 0:
            for i in range(0, k - 1):
                for j in range(i+1, k):
                    I = [i, j]
                    sums = sum_columns_from_V(I, VIr)
                    s_b = use.add_column_arrays(s_p, sums)
                    if np.count_nonzero(s_b) == t - p:
                        new_e = np.concatenate((np.zeros(shape=(1, k)), s_b), axis=None)
                        for i in I:
                            new_e[i] = 1
                        if (np.matmul(VIr, np.array(new_e).T) == np.array(s_p)).all():
                            break
        if np.count_nonzero(new_e) == t:
            break
        limit = limit - 1
    return use.from_float_to_int_array(np.matmul(new_e, Q.T))
