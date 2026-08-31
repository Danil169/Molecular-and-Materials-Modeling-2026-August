#!/usr/bin/env python
from ase.io import read
from ase.calculators.espresso import Espresso, EspressoProfile
import subprocess, os, shutil, sys

# 1. Read relaxed structure
if not os.path.exists('al_relaxed.vasp'):
    print("Run 2_relax_al.py first!")
    sys.exit(1)
atoms = read('al_relaxed.vasp')

# 2. Setup
pseudopotentials = {'Al': 'Al.upf'}
base_input_data = {
    'control': {'prefix': 'Al_Elec', 'outdir': './tmp', 'pseudo_dir': '../../', 'wf_collect': True},
    'system': {'ecutwfc': 50, 'occupations': 'smearing', 'smearing': 'mv', 'degauss': 0.02},
    'electrons': {'conv_thr': 1.0e-8}
}

os.environ['OMP_NUM_THREADS'] = '1'
profile = EspressoProfile(command='mpirun -np 4 pw.x', pseudo_dir='../../')

# SCF
print("Running SCF...")
scf_inp = base_input_data.copy()
scf_inp['control']['calculation'] = 'scf'
atoms.calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=scf_inp, kpts=(12,12,12))
print(f"Energy: {atoms.get_potential_energy():.4f} eV")

# Move outdir so we don't clobber it (ASE trick for NSCF)
# Actually, ASE NSCF directly uses the same profile, we just change to NSCF
print("Running NSCF for DOS...")
nscf_inp = base_input_data.copy()
nscf_inp['control']['calculation'] = 'nscf'
nscf_calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=nscf_inp, kpts=(24,24,24))
atoms.calc = nscf_calc
nscf_calc.calculate(atoms=atoms, properties=[], system_changes=['positions', 'cell', 'numbers', 'pbc'])

# DOS tool
print("Calculating DOS...")
with open('dos.in', 'w') as f:
    f.write(f"&DOS\n  prefix='Al_Elec',\n  outdir='./tmp',\n  Emin=-10.0, Emax=20.0, DeltaE=0.05,\n  fildos='al_total_dos.dat'\n/\n")
subprocess.run("dos.x < dos.in > dos.out 2>&1", shell=True)
print("Done! DOS saved to al_total_dos.dat")
