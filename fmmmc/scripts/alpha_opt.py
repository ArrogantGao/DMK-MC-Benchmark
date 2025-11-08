import time
import numpy as np

from coulomb_mc import mc_fmm_mm, mc_fmm_lm


from ppmd import *

import math
import ctypes

INT64 = ctypes.c_int64
REAL = ctypes.c_double


import os
import gc

import sys

from timing_func import timing_1


def step_time(l):
    return 2.718281828459045 * l[0] + l[1]


if __name__ == "__main__":

    stdout = open("./stdout", "a")

    method = sys.argv[1]
    bc = sys.argv[2]
    npoint = 10
    nset = np.logspace(4, 5.5, npoint)

    R_set = (
        4,
        5,
    )

    if method == "lm":
        mx = mc_fmm_lm.MCFMM_LM
    elif method == "mm":
        mx = mc_fmm_mm.MCFMM_MM
    else:
        raise RuntimeError()

    h = ""
    for N in nset:
        N = int(N)

        R_times = []
        for rx in R_set:
            gc.disable()
            time.sleep(2)
            l = timing_1(bc, mx, N, R=rx)
            R_times.append(step_time(l))
            gc.enable()
            gc.collect()

        print_str = "{:8.2e} || ".format(N) + " | ".join(
            ["{:d} {:8.2e}".format(ax, bx) for ax, bx in zip(R_set, R_times)]
        )
        print(print_str)
        stdout.write(print_str + '\n')
        stdout.flush()
