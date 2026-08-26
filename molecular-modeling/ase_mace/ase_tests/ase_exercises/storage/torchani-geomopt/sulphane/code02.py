import torch
import torchani
from ase import Atoms
from ase.optimize import BFGS

# 1. Define H2S Molecule (Sulfur at index 0, Hydrogens at 1 and 2)
pos = [[0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [0.0, -1.0, 1.0]]
atoms = Atoms('SH2', positions=pos)

# 2. Setup ANI-2x (Supports Sulfur)
calculator = torchani.models.ANI2x().ase()
atoms.set_calculator(calculator)

# 3. Optimize
opt = BFGS(atoms, logfile=None) # Set logfile=None for cleaner output
opt.run(fmax=0.01)

# 4. Extract Geometric Parameters
# get_distance(atom1, atom2)
# get_angle(atom1, vertex, atom2)
d_sh1 = atoms.get_distance(0, 1)
d_sh2 = atoms.get_distance(0, 2)
angle_hsh = atoms.get_angle(1, 0, 2)

# 5. Experimental Data (Reference)
# Experimental H2S: Bond ~1.336 A, Angle ~92.1 deg
exp_bond = 1.3356
exp_angle = 92.11

print("-" * 30)
print(f"{'Parameter':<15} | {'ANI-2x':<10} | {'Exp.':<10} | {'Diff.'}")
print("-" * 30)
print(f"{'S-H1 (A)':<15} | {d_sh1:<10.4f} | {exp_bond:<10.4f} | {abs(d_sh1-exp_bond):.4f}")
print(f"{'S-H2 (A)':<15} | {d_sh2:<10.4f} | {exp_bond:<10.4f} | {abs(d_sh2-exp_bond):.4f}")
print(f"{'H-S-H (deg)':<15} | {angle_hsh:<10.2f} | {exp_angle:<10.2f} | {abs(angle_hsh-exp_angle):.2f}")
print("-" * 30)

