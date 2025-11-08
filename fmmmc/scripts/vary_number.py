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



if __name__ == '__main__':
    

    name = sys.argv[1]
    bc = sys.argv[2]
    npoint = 30
    nset = np.logspace(3, 6, npoint)
    
    tset_local = []
    tset_multipole = []

    for N in nset:
        N = int(N)
        
        gc.disable()
        time.sleep(2)
        l = timing_1(bc, mc_fmm_lm.MCFMM_LM, N)
        gc.enable()
        gc.collect()


        tset_local.append((N, l[0], l[1], l[2], l[3]))

        gc.disable()
        time.sleep(2)
        m = timing_1(bc, mc_fmm_mm.MCFMM_MM, N)
        gc.enable()
        gc.collect()        

        tset_multipole.append((N, m[0], m[1], m[2], m[3]))
        print(N, "|", l[0], l[1], l[2], l[3], "|", m[0], m[1], m[2], m[3])
    


    lname = name + '_tset_local.npy'
    mname = name + '_tset_multipole.npy'

    tset_local = np.array(tset_local)
    if os.path.exists(lname):
        to = np.load(lname)
        tset_local = np.minimum(tset_local, to)

    np.save(lname, tset_local)

    tset_multipole = np.array(tset_multipole)
    if os.path.exists(mname):
        to = np.load(mname)
        tset_multipole = np.minimum(tset_multipole, to)

    np.save(mname, tset_multipole)

