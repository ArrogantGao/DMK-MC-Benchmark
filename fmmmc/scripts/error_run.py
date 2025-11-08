import os
import json
import numpy as np
import sys
import ctypes
from itertools import product
import tempfile

from ppmd.utility.lattice import cubic_lattice
from make_dl_monte_inputs import get_dlmonte_energy


from ppmd import *

import ppmd.coulomb.ewald_half
import ppmd.coulomb.fmm

import ppmd.lib as lib


from coulomb_mc.mc_fmm_lm import MCFMM_LM
from coulomb_mc.mc_fmm_mm import MCFMM_MM
from coulomb_mc.mc_short_range import NonBondedDiff

from scipy import constants as spc

PairLoop = pairloop.CellByCellOMP
ParticleLoop = loop.ParticleLoopOMP
State = state.State
PositionDat = data.PositionDat
ParticleDat = data.ParticleDat
ScalarArray = data.ScalarArray
Kernel = kernel.Kernel
GlobalArray = data.GlobalArray
Constant = kernel.Constant




if __name__ == "__main__":

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "input.json"

    with open(filename) as fh:
        data = json.loads(fh.read())
    print(data)
    L_set = data["L_set"]
    rng = np.random.RandomState(data.get('seed', None))

    perturb = float(data.get("perturb", 0.01))

    Ns = int(data["side_lens"][0])
    assert Ns % 2 == 0
    N = Ns * Ns * Ns
    Nt = N // 2
    a = float(data["lattice_spacing"])
    E = a * Ns
    l = cubic_lattice((Ns, Ns, Ns), (E, E, E))

    i0 = []
    i1 = []

    o = np.zeros(3, "float64") + 0.5 * (E - a)
    oc = np.zeros(3, "float64") + 0.5 * E
    
    q = np.zeros((N, 1), "float64")
    gid = np.zeros((N, 1), ctypes.c_int)

    for li, lx in enumerate(l):
        s = lx + o

        xk = round(s[0] / a)
        yk = round(s[1] / a)
        zk = round(s[2] / a)

        tb = (xk == yk) ^ (zk % 2 == 0)
        gid[li, 0] = li

        if tb:
            i0.append(li)
            q[li, 0] = -1.0
        else:
            i1.append(li)
            q[li, 0] = 1.0


    r = perturb
    l += rng.uniform((-r) * a, r * a, l.shape)

    # l = rng.uniform(low=-0.5*E, high=0.5*E, size=l.shape)

    
    DLM_CURRENT_ENERGY = get_dlmonte_energy(N, data, E, i0, i1, l)
    
    A = [(lx, State()) for lx in product(L_set, (MCFMM_LM, MCFMM_MM))]
    for si, sx in enumerate(A):
        sx[1].domain = domain.BaseDomainHalo(extent=E)
        sx[1].domain.boundary_condition = domain.BoundaryTypePeriodic()
        sx[1].pos = PositionDat(ncomp=3)
        sx[1].charge = ParticleDat(ncomp=1)
        sx[1].gid = ParticleDat(ncomp=1, dtype=ctypes.c_int)

        with sx[1].modify() as mv:
            mv.add({
                sx[1].pos: l,
                sx[1].charge: q,
                sx[1].gid: gid
            })

    LM_MM = [sx[0][1](sx[1].pos, sx[1].charge, sx[1].domain, 'pbc', r=3, l=sx[0][0]) for sx in A]

    for lx in LM_MM:
        lx.initialise()


    B = State()
    B.domain = domain.BaseDomainHalo(extent=E)
    B.domain.boundary_condition = domain.BoundaryTypePeriodic()
    B.pos = PositionDat(ncomp=3)
    B.charge = ParticleDat(ncomp=1)
    B.gid = ParticleDat(ncomp=1, dtype=ctypes.c_int)

    with B.modify() as mv:
        mv.add({
            B.pos: l,
            B.charge: q,
            B.gid: gid
        })

    LM_B = MCFMM_LM(B.pos, B.charge, B.domain, 'pbc', r=3, l=data["L_correct"])
    LM_B.initialise()




    internal_to_ev = ppmd.coulomb.fmm.internal_to_ev()
    current_hop_distance = 0.25

    

    err_data = {
            'method' : [],
            'p': [],
            'diff': [],
            'correct_diff': [],
            'rel_err_raw': [],
    }
    
    energy_data = {
        'step': [],
        'correct_diff': [],
        'total_energy': [],
    }


    def get_method(lm_mm):
        if lm_mm == MCFMM_LM:
            return 'lm'
        elif lm_mm == MCFMM_MM:
            return 'mm'
        else:
            raise RuntimeError()


    for stepx in range(int(data['steps'])):
        print("STEP", stepx)

        direction = rng.uniform(low=-1.0*current_hop_distance, high=current_hop_distance, size=3)
        particle_id = rng.randint(N)
        
        original_pos = l[particle_id, :].copy()
        prop_pos = original_pos + direction
    
        for dx in (0, 1, 2):
            if prop_pos[dx] < E * -0.5: prop_pos[dx] += E
            elif prop_pos[dx] > E * 0.5: prop_pos[dx] -= E

            assert prop_pos[dx] >= E * -0.5
            assert prop_pos[dx] <= E *  0.5       
        
        print(stepx, particle_id, prop_pos)


        l[particle_id, :] = prop_pos
        DLM_PROP_ENERGY = get_dlmonte_energy(N, data, E, i0, i1, l)
        l[particle_id, :] = original_pos
        
        DLM_DIFF = DLM_PROP_ENERGY - DLM_CURRENT_ENERGY

        local_id = np.where(B.gid.view == particle_id)[0]
        move = (particle_id, prop_pos)
        correct_pe = LM_B.propose(move) * internal_to_ev

        print("CORRECT", correct_pe)

        energy_data['step'].append(stepx)
        energy_data['correct_diff'].append(correct_pe)
        energy_data['total_energy'].append(LM_B.energy)

        DLM_ERR = abs(DLM_DIFF - correct_pe) / abs(correct_pe)
        print("DLM_ERR", DLM_ERR)

        err_data['method'].append('dlm')
        err_data['p'].append(-1)
        err_data['diff'].append(DLM_DIFF)
        err_data['correct_diff'].append(correct_pe)
        err_data['rel_err_raw'].append(DLM_ERR)

        for si, sx in enumerate(A):
            ss = sx[1]
            lm_mm = LM_MM[si]

            local_id = np.where(ss.gid.view == particle_id)[0]
            move = (particle_id, prop_pos)
            pe = lm_mm.propose(move) * internal_to_ev
            err = abs(pe - correct_pe) / abs(correct_pe)
            print(sx[0], err)

            err_data['method'].append(get_method(sx[0][1]))
            err_data['p'].append(sx[0][0])
            err_data['diff'].append(pe)
            err_data['correct_diff'].append(correct_pe)
            err_data['rel_err_raw'].append(err)
    
    open('err_data.json', 'w').write(json.dumps(err_data))
    open('energy_data.json', 'w').write(json.dumps(energy_data))




        


















