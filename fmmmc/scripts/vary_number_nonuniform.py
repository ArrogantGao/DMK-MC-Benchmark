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
    
    combined_name = name + '_nonuniform.csv'
    with open(combined_name, 'w') as f:
        f.write('N,e,local_propose,local_accept,local_R,local_time\n')
        for N in nset:
            N = int(N)
            
            gc.disable()
            time.sleep(2)
            l = timing_nonuniform(bc, mc_fmm_lm.MCFMM_LM, N)
            gc.enable()
            gc.collect()


            tset_local_results = (l[0], l[1], l[2], l[3])

            # gc.disable()
            # time.sleep(2)
            # m = timing_nonuniform(bc, mc_fmm_mm.MCFMM_MM, N)
            # gc.enable()
            # gc.collect()        

            # tset_multipole_results = (m[0], m_results[1], m[2], m[3])
            
            result_line_data = (N,) + (l[4],) + tset_local_results
            f.write(','.join(map(str, result_line_data)) + '\n')
            f.flush()
            print(N, l[4], "|", l[0], l[1], l[2], l[3])

