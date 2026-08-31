#!/usr/bin/env python
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.filters import UnitCellFilter
from ase.optimize import BFGS
from ase.io import write
import os

# Start with experimental guess
atoms = bulk('Al', 'fcc', a=4.049)

pseudopotentials = {'Al': 'Al.upf'}
input_data = {
    'control': {'calculation': 'scf', 'prefix': 'al_relax', 'outdir': './tmp', 'pseudo_dir': '../../', 'tstress': True, 'tprnfor': True},
    'system': {'ecutwfc': 50, 'occupations': 'smearing', 'smearing': 'mv', 'degauss': 0.02},
    'electrons': {'conv_thr': 1.0e-8}
}

os.environ['OMP_NUM_THREADS'] = '1'
profile = EspressoProfile(command='mpirun -np 4 pw.x', pseudo_dir='../../')
atoms.calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=input_data, kpts=(12,12,12))

print("Starting cell relaxation...")
ucf = UnitCellFilter(atoms)
opt = BFGS(ucf, logfile='relax.log')
opt.run(fmax=0.005)

write('al_relaxed.vasp', atoms, format='vasp')
print("Relaxed structure saved to al_relaxed.vasp")
print(f"Relaxed Cell Vectors:\n{atoms.get_cell()}")
