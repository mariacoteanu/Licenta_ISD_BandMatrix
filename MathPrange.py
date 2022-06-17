import itertools
import numpy as np
import math
from sympy import Matrix, zeros
from random import randint


def RandomPermutation(n):
    identity = np.eye(n, dtype=int)
    factorial = math.factorial(n)
    randomValue = randint(0, factorial)
    i = 0
    for m in itertools.permutations(identity):
        if i == randomValue:
            return Matrix(m)
        else:
            i = i + 1


def Prange(s, H, t, n):

    while True:
        while True:
            Q = RandomPermutation(n)
            print("\nQ : {} ".format(Q))
            new_H = H*Q
            print("\nH^ : {} ".format(new_H))
            U, W = new_H.rref() # ???
            print("-----U------\n", U)
            print("-----W------\n", W)

            break
            # if W != np.eye(r, dtype=int):
            #     break

        new_s = U*s
        zero = zeros(1, 3)
        transpose = new_s.T
        new_e = transpose.row_insert(0, zero)
        weight = len(new_e.values()) # cate valori dif de 0 sunt
        if weight == t:
            break

    return new_e * Q.T


if __name__ == '__main__':
    t = 3
    print(' t = ', t)

    s = Matrix([1, 0, 0, 1, 1])
    print("\ns : {} ".format(s))

    r = len(s.col(0))
    print('\n r = ', r)

    n = 5

    H = Matrix(np.random.randint(2, size=(r, n)))
    print("\nH : {} ".format(H))

    print(Prange(s, H, t, n))
