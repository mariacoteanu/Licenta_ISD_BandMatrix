import random
import numpy as np
from matplotlib import pyplot as plt


def display_plot(alg, r, n, xs, medie):
    types = np.array(["True", "False", "Undefined"])
    values = np.array(xs)
    means = np.mean(np.array(medie))

    plt.barh(types, values, label=f"Mean time = {means} ")
    plt.legend()
    if alg == 1:
        plt.savefig(f'images/prange/{r}+{n}.png', bbox_inches='tight')
    else:
        plt.savefig(f'images/lee/{r}+{n}.png', bbox_inches='tight')
    plt.clf()


def append_to_file(file, r, n, t, tip, rez):
    with open(file, "a") as file_object:
        file_object.write("\n")
        file_object.write(f"{r} {n} {t} {tip} {rez}")


def isOkInput(r, n):
    if r > n:
        return False
    if r <= 0 or n <= 0:
        return False
    return True


def generateCodeword(t, n):
    v = [0] * n
    i = 0
    while i < t:
        elem = random.choice(np.arange(0, n-1))
        while v[elem] == 1:
            elem = random.choice(np.arange(0, n-1))
        else:
            v[elem] = 1
            i += 1

    return v


def RandomPermutation(n):
    identity = np.zeros((n, n))
    pozitions = np.arange(0, n)
    random.shuffle(pozitions)
    for i in range(n):
        identity[i][pozitions[i]] = 1

    return identity


def isIr(U, r, n):
    sub = U[0:r, n - r:n]
    m2 = np.eye(r)

    for i in range(0, r):
        for j in range(0, r):
            if sub[i][j] != m2[i][j]:
                return False

    return True


def altfelMatriciBanda(r, n):
    rnd = random.randint(1, 2)
    M = np.zeros((r, n))

    if rnd == 1:
        for i in range(r):
            for j in range(n-r+1):
                M[i][i+j] = 1
    else:
        for i in range(r):
            for j in range(n-r+1):
                M[i][n-i-j-1] = 1

    return M, rnd


def isOk(P, H, V):
    intermediar = np.matmul(P, H)
    for i in range(len(intermediar)):
        for j in range(len(intermediar[0])):
            intermediar[i][j] = intermediar[i][j] % 2
    if (V == intermediar).all():
        return True

    return False


def matrix_base2(P):
    for i in range(len(P)):
        for j in range(len(P[0])):
            P[i][j] = P[i][j] % 2

    return P


def array_base2(A):
    return [elem % 2 for elem in A]


def write_to_file(fd, M, m):
    fd.write(f"\n----- {m} ----\n")
    for line in M.tolist():
        fd.write(" ".join(str(line)) + "\n")


def from_float_to_int_array(P):
    return [int(elem) for elem in P]


def extract_column(M, i):
    return [row[i] for row in M]


def add_column_line(v_i, v_j):
    return array_base2([v_i[k][0] + v_j[k] for k in range(0, len(v_i))])


def add_column_arrays(v_i, v_j):
    return array_base2([v_i[k] + v_j[k] for k in range(0, len(v_i))])