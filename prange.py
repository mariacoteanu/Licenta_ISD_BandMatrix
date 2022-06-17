import numpy as np
import random
import echelon as ech
import usefull as use


def matriciSpeciale(r, n):
    rnd = random.randint(1, 2)
    mijloc = int(r / 2)
    banda = np.zeros((r, n))
    ilen_ones = random.randint(1, r - 1)  # sub
    jlen_ones = random.randint(1, r - 1)  # deasupra
    print(f'----------tip = {rnd} && ilen = {ilen_ones} && jlen = {jlen_ones} --------')
    if rnd == 1:
        for lin in range(r):
            for col in range(mijloc):
                if lin == col:
                    banda[lin][col] = 1
                elif col < lin <= col + ilen_ones:
                    banda[lin][col] = 1
                elif lin < col <= lin + jlen_ones:
                    banda[lin][col] = 1
    else:
        for lin in range(r):
            for pas, col in zip(reversed(range(1, mijloc + 1)), range(mijloc)):
                if mijloc + pas - jlen_ones <= lin <= mijloc + pas:
                    banda[lin][col] = 1
                elif mijloc + pas <= lin <= mijloc + pas + ilen_ones:
                    banda[lin][col] = 1
    for lin in range(max(0, mijloc - jlen_ones), min(r, mijloc + ilen_ones+1)):
        for col in range(mijloc, n - mijloc):
            banda[lin][col] = 1
    if rnd == 1:
        for lin in range(max(0, mijloc - jlen_ones + 1), r):
            for pas, col in zip(range(1, mijloc+1), range(n - mijloc, n)):
                if mijloc >= lin >= max(0, mijloc - jlen_ones + pas):
                    banda[lin][col] = 1
                elif mijloc < lin <= min(r, mijloc + ilen_ones+pas):
                    banda[lin][col] = 1
    else:
        for lin in range(0, min(r, mijloc + ilen_ones + 1)):
            for pas, col in zip(range(1, mijloc+1), range(n - mijloc, n)):
                if mijloc - pas - jlen_ones <= lin <= mijloc - pas:
                    banda[lin][col] = 1
                elif mijloc - pas <= lin <= mijloc - pas + ilen_ones:
                    banda[lin][col] = 1

    return banda


def generareMatriceBanda(r, n):
    rnd = random.randint(1, 4)
    banda = np.zeros((r, n))
    print(f'---------- tip = {rnd} --------')
    if rnd == 1:
        # generez matrice banda pe orizontala
        istart = random.randint(0, r - 1)
        print(f'---------- istart = {istart} --------')
        len_ones = random.randint(1, r - istart)
        print(f'---------- len = {len_ones} --------')
        for i in range(istart, istart + len_ones):
            for j in range(n):
                banda[i][j] = 1
    elif rnd == 2:
        # generez matrice banda pe verticala
        jstart = random.randint(0, n - 1)
        print(f'---------- jstart = {jstart} --------')
        len_ones = random.randint(1, n - jstart)
        print(f'---------- len = {len_ones} --------')
        for i in range(r):
            for j in range(jstart, jstart + len_ones):
                banda[i][j] = 1
    elif rnd == 3:
        # generez matrice banda diagonala principala
        ilen_ones = random.randint(1, r - 1)
        jlen_ones = random.randint(1, n - 1)

        print(f'---------- ilen = {ilen_ones} && jlen = {jlen_ones} --------')
        for i in range(r):
            for j in range(n):
                if i == j:
                    banda[i][j] = 1
                elif j < i <= j + ilen_ones:
                    banda[i][j] = 1
                elif i < j <= i + jlen_ones:
                    banda[i][j] = 1
    else:
        # generez matrice banda diagonala principala
        ilen_ones = random.randint(1, r - 1)
        jlen_ones = random.randint(1, n - 1)

        print(f'---------- ilen = {ilen_ones} && jlen = {jlen_ones} --------')
        for i in range(r):
            for j in range(n):
                if i + j == n - 1:
                    banda[i][j] = 1
                elif n - 1 < i + j <= n - 1 + ilen_ones:
                    banda[i][j] = 1
                elif n - 1 > i + j >= n - jlen_ones - 1:
                    banda[i][j] = 1
    return banda


def Prange(s, H, t, r, n):

    limit = n
    while True:
        if limit == 0:
            return []
        limitWhile = min(r*n, 500)
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
