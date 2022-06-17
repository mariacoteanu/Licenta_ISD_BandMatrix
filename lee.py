import random

import usefull as use
import numpy as np
import echelon as ech


def isExtract(H, s, n, r):
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
    columns_sum = use.extract_column(VIr, I[0])
    for i in range(1, len(I)):
        columns_sum = use.add_column_arrays(columns_sum, use.extract_column(VIr, I[i]))
    return use.array_base2(columns_sum)


def integerToCombinations(j, k):
    values = []
    while j:
        value = random.choice(range(0, k))
        while value in values:
            value = random.choice(range(0, k))
        else:
            values.append(value)
            j = j - 1
    return values


def LeeAlgorithm(H, s, t, r, n):
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
