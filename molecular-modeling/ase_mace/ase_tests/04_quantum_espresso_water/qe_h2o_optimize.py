#
# ASE driven Quantum Espresso geometry optimization of water
#
import os
from ase import Atoms
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import BFGS
from ase.units import Ry, Bohr

os.environ['OMP_NUM_THREADS'] = '1'

# Set the path to your pseudopotential directory
PSEUDO_DIR = '/usr/share/espresso/pseudo/'

# 1. Create the initial water molecule structure
h2o = molecule('H2O')
h2o.set_cell([12, 12, 12])  # Create a vacuum box
h2o.center()               # Center the molecule in the box

# 2. Configure the Quantum ESPRESSO Calculator
# Update pseudopotential names to match what's in /usr/share/espresso/pseudo/
pseudopotentials = {
    'O': 'O.pbe-kjpaw.UPF',
    'H': 'H.pbe-kjpaw.UPF'
}

input_data = {
    'control': {
        'calculation': 'scf', 
        'prefix': 'h2o',
        'outdir': './outdir',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 46.0,    # Plane wave cutoff (Ry)
   #     'ecutrho': 350.0,   # Charge density cutoff (Ry)
        'ibrav': 0,         # Use 0 for ASE to handle cell parameters
        'nosym': True,      # No symmetry for molecules
        'noinv': True,
        'occupations': 'smearing',  # Better convergence for molecules
        'smearing': 'gaussian',
        'degauss': 0.02,
    },
    'electrons': {
        'conv_thr': 1e-8,
        'mixing_beta': 0.7,
        'diagonalization': 'david'
    }
}

command = 'mpirun -np 4 pw.x'
profile = EspressoProfile(command, pseudo_dir=PSEUDO_DIR)

calc = Espresso(profile=profile, 
                pseudopotentials=pseudopotentials,
                input_data=input_data,
                kpts=(1, 1, 1))  # Gamma-point only for molecules

h2o.calc = calc

# 3. Run the Geometry Optimization
dyn = BFGS(h2o, trajectory='h2o_opt.traj', logfile='h2o_opt.log')
dyn.run(fmax=0.01)

# 4. Print optimized geometry
print("Optimized Bond Lengths (A):")
print(f"O-H1: {h2o.get_distance(0, 1):.4f}")
print(f"O-H2: {h2o.get_distance(0, 2):.4f}")
print(f"H-O-H Angle (deg): {h2o.get_angle(1, 0, 2):.2f}")
print(f"Final Energy (eV): {h2o.get_total_energy():.4f}")
