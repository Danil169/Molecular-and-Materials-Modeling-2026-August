#!/usr/bin/env python3
"""Challenge II.8: relax C8 and recompute H adsorption at a higher cutoff."""
import os

import numpy as np
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import read, write
from ase.optimize import BFGS

os.environ['OMP_NUM_THREADS'] = '1'
profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir='../..')
ECUT = 100

def calculator(prefix, species, nat):
    return Espresso(
        profile=profile, pseudopotentials=species,
        input_data={
            'control': {'calculation': 'scf', 'prefix': prefix, 'outdir': './tmp', 'verbosity': 'low',
                        'tstress': True, 'tprnfor': True},
            'system': {'ecutwfc': ECUT, 'ibrav': 0, 'occupations': 'smearing', 'smearing': 'gauss',
                       'degauss': 0.01, 'nspin': 2, 'vdw_corr': 'DFT-D3', 'dftd3_version': 4,
                       'nat': nat, 'ntyp': len(species)},
            'electrons': {'conv_thr': 1.0e-10},
        }, kpts=(1, 1, 1),
    )

def relax(atoms, label, species):
    atoms.calc = calculator(label, species, len(atoms))
    opt = BFGS(atoms, trajectory=f'{label}.traj', logfile=f'{label}.log')
    opt.run(fmax=0.01, steps=100)
    energy = atoms.get_potential_energy()
    fmax = np.abs(atoms.get_forces()).max()
    write(f'{label}_final.vasp', atoms, direct=False)
    return energy, fmax

c8 = read('../c.vasp', format='vasp')
e_c8, f_c8 = relax(c8, 'c8_100Ry', {'C': 'C.upf'})

hc8 = read('../ch.vasp', format='vasp')
e_hc8, f_hc8 = relax(hc8, 'hc8_100Ry', {'C': 'C.upf', 'H': 'H.upf'})

h = read('../h.vasp', format='vasp')
h.calc = calculator('h_100Ry', {'H': 'H.upf'}, 1)
e_h = h.get_potential_energy()
write('h_100Ry_final.vasp', h, direct=False)

binding = e_c8 + e_h - e_hc8
print('Challenge II.8 — C8 relaxation and higher-cutoff adsorption energy')
print(f'ecutwfc: {ECUT} Ry (original workflow: 80 Ry)')
print(f'E(C8 relaxed): {e_c8:.6f} eV; max force: {f_c8:.6f} eV/Å')
print(f'E(H@C8 relaxed): {e_hc8:.6f} eV; max force: {f_hc8:.6f} eV/Å')
print(f'E(H): {e_h:.6f} eV')
print(f'E_bind = E(C8)+E(H)-E(H@C8): {binding:.6f} eV')
