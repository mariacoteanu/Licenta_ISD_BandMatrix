import numpy as np
import usefull as use


def ToReducedRowEchelonForm(M):
    if not M:
        return

    rowCount = len(M)
    columnCount = len(M[0])
    lead = columnCount - rowCount

    p = np.eye(rowCount)

    for r in range(rowCount):
        new_P = np.eye(rowCount)
        i = r
        if lead >= columnCount:
            return use.matrix_base2(p), M
        while M[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return use.matrix_base2(p), M
        M[i], M[r] = M[r], M[i]

        if i != r:
            new_P[i][r] = 1
            new_P[r][i] = 1
            new_P[i][i] = 0
            new_P[r][r] = 0

        p = np.matmul(new_P, p)
        lv = M[r][lead]
        M[r] = [mrx / float(lv) for mrx in M[r]]

        for i in range(rowCount):
            if i != r:
                lv = M[i][lead]

                M[i] = [(2 - (iv - lv * rv)) % 2 for rv, iv in zip(M[r], M[i])]
                if lv != 0:
                    new_P = np.eye(rowCount)
                    new_P[i][r] = 1

                    p = np.matmul(new_P, p)
        lead += 1
    return use.matrix_base2(p), M
