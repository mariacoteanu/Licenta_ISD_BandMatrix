import sys
import time
import numpy as np
import random
import usefull as use
import prange
import lee

PRANGE_RESULT_FILE = "rezultatePrange.txt"
LEE_RESULT_FILE = "rezultateLee.txt"

if __name__ == '__main__':
    # r = int(input("number of rows = "))
    # n = int(input("number of columns = "))
    # while not use.isOkInput(r, n):
    #     print("Give positive numbers and number of rows must be lower that numer of columns!\n")
    #     r = int(input("number of rows = "))
    #     n = int(input("number of columns = "))
    n = 50
    for r in range(7, 8):
        print(f"PAS == {r}")
        results = {1: [0, 0, 0], 2: [0, 0, 0]}
        times = {1: [], 2: []}

        for i in range(0, 10):
            print("\rTest {}/{}".format(i, 10), end="")
            t = random.randint(1, r - 1)

            while True:
                x = use.generateCodeword(t, n)
                H, tip = use.altfelMatriciBanda(r, n)
                s = use.array_base2(np.matmul(H, x))

                if np.count_nonzero(s) != 0:
                    break
            w_s = np.count_nonzero(s)

            start = time.time()
            print('------- PRANGE ---------\n')
            prangeResult = prange.Prange(s, H, t, r, n)
            final = time.time()
            times[1].append(final - start)

            if not prangeResult:
                correct = "Undefined"
                results[1][2] += 1
            elif x == prangeResult:
                correct = "True"
                results[1][0] += 1
            else:
                correct = "False"
                results[1][1] += 1
            use.append_to_file(PRANGE_RESULT_FILE, r, n, t, tip, correct)

            start = time.time()
            print('------- LEE-BRICKELL ---------\n')
            leeResult = lee.LeeAlgorithm(H, s, t, r, n)
            final = time.time()
            times[2].append(final - start)

            if not leeResult:
                correct = "Undefined"
                results[2][2] += 1
            elif x == leeResult:
                correct = "True"
                results[2][0] += 1
            else:
                correct = "False"
                results[2][1] += 1
            use.append_to_file(LEE_RESULT_FILE, r, n, t, tip, correct)
            sys.stdout.flush()

        use.display_plot(1, r, n, results[1], times[1])
        use.display_plot(2, r, n, results[2], times[2])
        time.sleep(5)
