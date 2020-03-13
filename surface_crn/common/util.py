import numpy as np


def weave(*threads):
    n = len(threads)
    woven = np.zeros(np.size(threads))
    for i, t in enumerate(threads):
        woven[i::n] = t
    return woven


def unweave(woven, n):
    return [woven[i::n] for i in range(n)]