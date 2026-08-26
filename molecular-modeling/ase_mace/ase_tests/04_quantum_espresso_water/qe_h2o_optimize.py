#
# ASE driven Quantum Espresso geometry optimization of water
#
import os
from ase import Atoms
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import BFGS
from ase.units import Ry, Bohr


# Path to the directory containing your .UPF files
#os.environ['ESPRESSO_PSEUDO'] = '/usr/share/espresso/pseudo/'
# Path to the directory containing the pw.x executable
#os.environ['ASE_ESPRESSO_COMMAND'] = 'mpirun -np 4 /opt/espresso/7.5/pw.x -in PREFIX.pwi > PREFIX.pwo'

# 1. Create the initial water molecule structure
# H-O-H structure, approx bond length 0.96 A, angle 104.5 deg
h2o = molecule('H2O')
h2o.set_cell([10, 10, 10])  # Create a 10x10x10 Angstrom vacuum box
h2o.center()               # Center the molecule in the box

# 2. Configure the Quantum ESPRESSO Calculator
# Update pseudo_dir to your actual pseudopotential folder
pseudopotentials = {
    'O': 'O.pbe-kjpaw.UPF',
    'H': 'H.pbe-kjpaw.UPF'
}

input_data = {
    'control': {
        'calculation': 'scf', 
        'prefix': 'h2o',
        'pseudo_dir': '/usr/share/espresso/pseudo/',  # CHANGE THIS
        'outdir': './outdir',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 40.0,    # Plane wave cutoff (Ry)
        'ecutrho': 320.0,   # Charge density cutoff (Ry)
        'ibrav': 0,         # 
        'nosym': True,      # No symmetry for molecules
        'noinv': True,
    },
    'electrons': {
        'conv_thr': 1e-8,
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    }
}

#command='mpirun -np 4 /opt/espresso/7.5/pw.x'
command='mpirun -np 2 /opt/espresso/7.5/pw.x'
profile = EspressoProfile(command,pseudo_dir='/usr/share/espresso/pseudo/')

calc = Espresso(profile=profile, pseudopotentials=pseudopotentials,
                input_data=input_data,
                kpts=(1, 1, 1) ) # Gamma-point only for molecules )

h2o.calc = calc

# 3. Run the Geometry Optimization
# fmax: Max force threshold for convergence (eV/Angstrom)
dyn = BFGS(h2o, trajectory='h2o_opt.traj', logfile='h2o_opt.log')
dyn.run(fmax=0.01)

# 4. Print optimized geometry
print("Optimized Bond Lengths (A):")
print(f"O-H1: {h2o.get_distance(0, 1):.4f}")
print(f"O-H2: {h2o.get_distance(0, 2):.4f}")
print(f"H-O-H Angle (deg): {h2o.get_angle(1, 0, 2):.2f}")
print(f"Final Energy (eV): {h2o.get_total_energy():.4f}")

