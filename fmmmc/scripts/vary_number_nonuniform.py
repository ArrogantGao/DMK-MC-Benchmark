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

from timing_func import timing_nonuniform



if __name__ == '__main__':
    

    name = sys.argv[1]
    bc = sys.argv[2]
    npoint = 30
    nset = np.logspace(4, 7, npoint)
    
    tset_combined = []

    for N in nset:
        N = int(N)
        
        gc.disable()
        time.sleep(2)
        l = timing_nonuniform(bc, mc_fmm_lm.MCFMM_LM, N)
        gc.enable()
        gc.collect()


        tset_local_results = (l[0], l[1], l[2], l[3])

        gc.disable()
        time.sleep(2)
        m = timing_nonuniform(bc, mc_fmm_mm.MCFMM_MM, N)
        gc.enable()
        gc.collect()        

        tset_multipole_results = (m[0], m[1], m[2], m[3])
        
        tset_combined.append((N,) + (l[4],) + tset_local_results + tset_multipole_results)
        print(N, l[4], "|", l[0], l[1], l[2], l[3], "|", m[0], m[1], m[2], m[3])
    


    combined_name = name + '_nonuniform.csv'

    tset_combined = np.array(tset_combined)

    np.savetxt(combined_name, tset_combined, delimiter=',', header='N,e,local_propose,local_accept,local_R,local_time,multipole_propose,multipole_accept,multipole_R,multipole_time')

