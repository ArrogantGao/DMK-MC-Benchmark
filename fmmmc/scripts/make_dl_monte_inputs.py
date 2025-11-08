import shutil
import os
import json
import numpy as np
import sys
from itertools import product
import tempfile
from ppmd.utility.lattice import cubic_lattice
import subprocess
import yaml

# /build-2.02.bsh SRL dir gfortran

CONTROL = """GENERATED CONTROL FILE
use ortho
finish
temperature       273.0
maxatmdist 1.0
acceptatmmoveupdate 100
acceptmolmoveupdate  100        # Period (in moves) at which the maximum move size is recalculated
steps             {STEPS}        # Number of moves to perform in simulation
equilibration     0        # Equilibration time before statistics are gathered (in moves)
print              {STEPS}        # Information is output every 'print' moves     
stack              {PRINT_STEPS}        # Size of blocks (in moves) for block averaging
yamldata {NEVER}
move atom 2 100  
Na core
Cl core
check {NEVER}
start
"""

FIELD = """GENERATED FIELD FILE
CUTOFF {CUTOFF}
UNITS EV
NCONFIGS 1
ATOMS 2
Na core 23.0 1.0
Cl core 35.0 -1.0
MOLTYPES 1
salt
MAXATOM {N}
FINISH
VDW       3
Na core       Na core       slj    0.005637 2.35
Cl core       Cl core       slj    0.004336 4.45
Cl core       Na core       slj    0.004943 3.40
CLOSE
"""

# Na Na LJ 0.005637 2.35
# Cl Cl LJ 0.004336 4.45
# Na Cl LJ 0.004943 3.40

OLD_LJ = """
Na core       Na core       slj    1.00000000 3
Cl core       Na core       slj    1.00000000 3
Cl core       Cl core       slj    1.00000000 3
"""

CONFIG = """NCONFIG 1 
         0         1
        {EXTENT_X}         0.0000000000        0.0000000000
        0.0000000000       {EXTENT_Y}          0.0000000000
        0.0000000000       0.0000000000        {EXTENT_Z}
NUMMOL 1 1 
molecule salt     {N}     {N}"""


def write_dlm_config(dir_name, N, data, E, i0, i1, l):

    with open(os.path.join(dir_name, "FIELD"), "w") as fh:
        fh.write(FIELD.format(CUTOFF=data["cutoff"], N=N))

    with open(os.path.join(dir_name, "CONTROL"), "w") as fh:
        fh.write(CONTROL.format(STEPS=data["steps"], NEVER=data["steps"] + 2, PRINT_STEPS=data["print_steps"]))

    with open(os.path.join(dir_name, "CONFIG"), "w") as fh:
        fh.write(CONFIG.format(N=N, EXTENT_X=E, EXTENT_Y=E, EXTENT_Z=E,))
        for ax in i0:
            fh.write(
                """
Cl core
{EX}    {EY}    {EZ}""".format(
                    EX=l[ax, 0], EY=l[ax, 1], EZ=l[ax, 2],
                )
            )
        for ax in i1:
            fh.write(
                """
Na core
{EX}    {EY}    {EZ}""".format(
                    EX=l[ax, 0], EY=l[ax, 1], EZ=l[ax, 2],
                )
            )


def get_dlm_yaml_energy(filename):

    with open(filename) as fh:
        data = yaml.safe_load(fh)
    
    e_lr = data[0]['energyrcp']
    e_sr = data[0]['energyreal']
    
    return e_lr + e_sr


def get_dlmonte_energy(N, data, E, i0, i1, l, DLM=None):
        
    if DLM is None:
        DLM = shutil.which('DLMONTE-SRL.X')
    if DLM is None:
        raise RuntimeError('could not find DLMONTE-SRL.X')


    tdir = tempfile.TemporaryDirectory()
    dir_name = tdir.name
    write_dlm_config(dir_name, N, data, E, i0, i1, l)
    subprocess.check_call(DLM, cwd=dir_name)
    yaml_file = os.path.join(dir_name, 'YAMLDATA.000')
    return get_dlm_yaml_energy(yaml_file)
    


if __name__ == "__main__":

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "input.json"

    with open(filename) as fh:
        data = json.loads(fh.read())
    print(data)

    rng = np.random.RandomState(seed=data["seed"])

    n_sample = int(data["n_sample"])
    perturb = float(data.get("perturb", 0.01))

    for Ns, samplex in product(data["side_lens"], range(n_sample)):
        Ns = int(Ns)
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

        for li, lx in enumerate(l):
            s = lx + o

            xk = round(s[0] / a)
            yk = round(s[1] / a)
            zk = round(s[2] / a)

            tb = (xk == yk) ^ (zk % 2 == 0)

            if tb:
                i0.append(li)
            else:
                i1.append(li)

        r = perturb
        l += rng.uniform((-r) * a, r * a, l.shape)


        dir_name = "DLM_" + str(Ns) + "_" + str(samplex)
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)

        L_set_dict = {"L_set": data["L_set"], "sample_index": samplex}
        with open(os.path.join(dir_name, "L_set.json"), "w") as fh:
            fh.write(json.dumps(L_set_dict))

        write_dlm_config(dir_name, N, data, E, i0, i1, l)
