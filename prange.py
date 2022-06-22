import numpy as np
import echelon as ech
import usefull as use


def Prange(s, H, t, r, n):
    """
    Prange's algorithm to find error-vector e
    :param s: syndrome a.i. s = Hx
    :param H: Band matrix
    :param t: weight
    :param r: number of rows
    :param n: number of columns
    :return: reconstructed initial vector x, or [] if fails
    """

    limit = n
    while True:
        if limit == 0:
            return []
        limitWhile = min(r*n, 200)
        ok = 1
        while True:
            if limitWhile == 0:
                ok = 0
                break
            Q = use.RandomPermutation(n)
            new_H = np.matmul(H, Q)

            P, V = ech.ToReducedRowEchelonForm(new_H.tolist())
            V = np.array(V)

            if use.isIr(V, r, n):
                break
            limitWhile = limitWhile - 1

        if ok == 1:
            new_s = np.array(use.array_base2(np.matmul(P, s)))
            k = n - r
            zeros = np.zeros(shape=(1, k))
            new_s_T = new_s.T
            new_e = np.concatenate((zeros, new_s_T), axis=None)
            weight = np.count_nonzero(new_e)

            if weight == t:
                break
        limit = limit - 1
    return use.from_float_to_int_array(np.matmul(new_e, Q.T))
