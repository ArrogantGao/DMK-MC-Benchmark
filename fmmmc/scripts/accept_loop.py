import time
import numpy as np

from coulomb_mc import mc_fmm_lm, mc_fmm_mm


from ppmd import *

import math
import ctypes

INT64 = ctypes.c_int64
REAL = ctypes.c_double

import sys
import os
import gc

import cProfile


def accept_loop(bc, method, N, N_accept=1000, R=None):

    e = 3. * (N**(1./3.))

    if R is None:
        R = max(3, int(math.log(2*N, 8)))

    L = 12


    rng = np.random.RandomState()

    A = state.State()
    A.domain = domain.BaseDomainHalo(extent=(e, e, e))
    A.domain.boundary_condition = domain.BoundaryTypePeriodic()


    A.P = data.PositionDat()
    A.Q = data.ParticleDat(ncomp=1)
    A.G = data.ParticleDat(ncomp=1, dtype=INT64)

    
    pi = np.array(rng.uniform(low=-0.5*e, high=0.5*e, size=(N, 3)), REAL)
    qi = np.array(rng.uniform(low=-1, high=1, size=(N, 1)), REAL)
    gi = np.arange(N).reshape((N, 1))

    with A.modify() as m:
        m.add({
            A.P: pi,
            A.Q: qi,
            A.G: gi
        })


    MC = method(A.P, A.Q, A.domain, bc, R, L)
    
    ti0 = time.time()
    MC.initialise()
    ti1 = time.time()
    

    
    t0 = time.time()
    
    if N_accept is not None:
        prof = cProfile.Profile()
        prof.enable()
        for testx in range(N_accept):
            lid = rng.randint(0, N)
            pos = rng.uniform(-0.5*e, 0.5*e, (3,))
            MC.accept((lid, pos), 0.0)
        prof.disable()
        prof.dump_stats('accept_loop.prof')
    else:
        while True:
            lid = rng.randint(0, N)
            pos = rng.uniform(-0.5*e, 0.5*e, (3,))
            MC.accept((lid, pos))

    t1 = time.time()
    
    t_accept = t1 - t0
    t_accept /= N_accept



    MC.free()

    return

if __name__ == '__main__':
    
    bc = sys.argv[1]
        
    m = sys.argv[2]
    if m == 'lm':
        method = mc_fmm_lm.MCFMM_LM
    elif m == 'mm':
        method = mc_fmm_mm.MCFMM_MM

    N = int(sys.argv[3])

    Np = None
    if len(sys.argv) > 4:
        Np = int(sys.argv[4])
    
    R = None
    if len(sys.argv) > 5:
        R = int(sys.argv[5])

    p = accept_loop(bc, method, N, Np, R)



