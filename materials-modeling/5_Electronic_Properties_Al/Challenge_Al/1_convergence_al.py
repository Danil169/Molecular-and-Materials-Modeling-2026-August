#!/usr/bin/env python
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import numpy as np
import os, sys

# 1. Experimental Al structure (a=4.049 A, 1 atom primitive cell)
atoms = bulk('Al', 'fcc', a=4.049)

# 2. Setup
pseudopotentials = {'Al': 'Al.upf'}
base_input_data = {
    'control': {'calculation': 'scf', 'prefix': 'al_conv', 'outdir': './tmp', 'pseudo_dir': '../../', 'verbosity': 'low'},
    'system': {'occupations': 'smearing', 'smearing': 'mv', 'degauss': 0.02, 'ecutwfc': 30},
    'electrons': {'conv_thr': 1.0e-8}
}

os.environ['OMP_NUM_THREADS'] = '1'
profile = EspressoProfile(command='mpirun -np 4 pw.x', pseudo_dir='../../')

# 3. Test values
ecut_values = np.arange(20, 80, 10)
kpoints_values = [(k,k,k) for k in range(4, 16, 2)]

print("--- 1. Ecutwfc Convergence ---")
k_mid = (8,8,8)
for ecut in ecut_values:
    inp = base_input_data.copy()
    inp['system']['ecutwfc'] = int(ecut)
    atoms.calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=inp, kpts=k_mid)
    print(f"ecutwfc = {ecut:2d} Ry | Energy = {atoms.get_potential_energy():.5f} eV")

print("\n--- 2. K-points Convergence (at ecut=50) ---")
for k in kpoints_values:
    inp = base_input_data.copy()
    inp['system']['ecutwfc'] = 50
    atoms.calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=inp, kpts=k)
    print(f"kpts = {k} | Energy = {atoms.get_potential_energy():.5f} eV")
