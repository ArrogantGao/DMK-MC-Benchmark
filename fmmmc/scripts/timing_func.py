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


def timing_1(bc, method, N, N_propose=1000, N_accept=1000, R=None, local=False):

    e = 3. * (N**(1./3.))
    
    if method == mc_fmm_lm.MCFMM_LM:
        alpha = 0.327
    elif method == mc_fmm_mm.MCFMM_MM:
        alpha = 0.138
    else:
        raise RuntimeError()


    if R is None:
        R = max(3, int(math.log(alpha*N, 8)))

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
    
    if not np.isfinite(MC.energy):
        raise RuntimeError('Energy was not a normal float')

    current_hop_distance = 0.25
    gc.collect()
    gc.disable()

    t0 = time.time()
    for testx in range(N_propose):

        lid = rng.randint(0, N)
        if local:
            pos = A.P[lid, :].copy().reshape((3,)) + rng.uniform(-current_hop_distance, current_hop_distance, (3,))
            for dx in (0, 1, 2):
                if pos[dx] < e * -0.5: pos[dx] += e
                elif pos[dx] > e * 0.5: pos[dx] -= e
                assert pos[dx] >= e * -0.5
                assert pos[dx] <= e *  0.5
        else:
            pos = rng.uniform(-0.5*e, 0.5*e, (3,))


        MC.propose((lid, pos))
    t1 = time.time()
    
    t_propose = t1 - t0
    t_propose /= N_propose

    t0 = time.time()
    for testx in range(N_accept):

        lid = rng.randint(0, N)
        if local:
            pos = A.P[lid, :].copy().reshape((3,)) + rng.uniform(-current_hop_distance, current_hop_distance, (3,))
            for dx in (0, 1, 2):
                if pos[dx] < e * -0.5: pos[dx] += e
                elif pos[dx] > e * 0.5: pos[dx] -= e
                assert pos[dx] >= e * -0.5
                assert pos[dx] <= e *  0.5
        else:
            pos = rng.uniform(-0.5*e, 0.5*e, (3,))

        MC.accept((lid, pos), 0.0)

    t1 = time.time()
    
    gc.enable()
    gc.collect()

    t_accept = t1 - t0
    t_accept /= N_accept

    MC.free()

    return t_propose, t_accept, R, ti1 - ti0, e


