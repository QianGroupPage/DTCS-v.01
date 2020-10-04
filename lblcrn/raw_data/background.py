import numpy as np


def shirley_bg(df):
    max_iter = 50
    max_rdiff = 1e-6
    n = df.size()
    ya = df.index[0]
    yb = df.index[-1]
    dy = yb - ya
    B = np.array([ya for _ in range(n)])
    PA = np.array([0 for _ in range(n)])
    old_A = 0
    for iter in range(max_iter):
        Y = df["data"] - B
        for i in range(1, n):
            PA[i] = PA[i-1] + (Y[i] + Y[i-1]) / 2 * (df.index[i] - df.index[i-1])
        rel_diff = abs(PA[n-1] - old_A) / old_A if old_A != 0. else 1.
        if rel_diff < max_rdiff:
            break
        old_A = PA[n-1]
        for i in range(n):
            B[i] = ya + dy / PA[n-1] * PA[i]

    df["shirley"] = B