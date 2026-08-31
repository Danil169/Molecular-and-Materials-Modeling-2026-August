#!/usr/bin/env python
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import read
import os

# ==============================================
# 1. Quantum Espresso input parameters
# ==============================================
pseudopotentials = {
    'Hg': 'Hg.upf',
}

input_data = {
    'control': {
        'prefix': 'energy_hg',
        'outdir': './tmp',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 80,
        'occupations': 'smearing',
        'smearing': 'gauss',
        'degauss': 0.01,
        'ibrav': 0,
        'nat': 1,
        'ntyp': 1,
        'vdw_corr': 'DFT-D3',
        'dftd3_version': 4
    },
    'electrons': {
        'conv_thr': 1.0e-10
    }
}

# ==============================================
# 2. Atomic Structure 
# ==============================================
vasp_file = 'hg.vasp'  

try:
    atoms = read(vasp_file, format='vasp')
    print(f"Successfully read {len(atoms)} atoms from {vasp_file}")
except FileNotFoundError:
    print(f"Error: {vasp_file} not found!")
    sys.exit(1)
except Exception as e:
    print(f"Error reading {vasp_file}: {e}")
    sys.exit(1)

# ==============================================
# 3. Calculator configuration
# ==============================================
# Set QE bin directory 
qe_bin = ""

# Parallel calculation 
os.environ['OMP_NUM_THREADS'] = '1'
pw_command = f'mpirun -np 4 pw.x'

pw_profile = EspressoProfile(
    command=pw_command,
    pseudo_dir='./'
)

# Set k-grids
scf_kpts = (1, 1, 1)

# ==============================================
# 4. SCF calculation
# ==============================================
scf_calc = Espresso(
    profile=pw_profile,
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=scf_kpts
)

print("\nRunning SCF calculation...", flush=True)
os.makedirs('tmp', exist_ok=True)

atoms.calc = scf_calc
total_energy = atoms.get_potential_energy()
print(f"  Total energy: {total_energy:.6f} eV", flush=True)
fermi_ev = atoms.calc.get_fermi_level()
print(f"  Fermi level: {fermi_ev:.6f} eV", flush=True)
