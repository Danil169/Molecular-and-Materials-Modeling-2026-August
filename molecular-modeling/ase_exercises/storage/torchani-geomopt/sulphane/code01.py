import torch
import torchani
from ase import Atoms
from ase.optimize import BFGS
from ase.visualize import view

# 1. Define the sulphane molecule (Hydrogen Sulfide)
# Approximate initial positions (Sulfur at center, Hydrogens nearby)
pos = [[0.0, 0.0, 0.0],          # Sulfur
       [0.0, 0.9, 1.2],          # Hydrogen 1
       [0.0, -0.9, 1.2]]         # Hydrogen 2
atoms = Atoms('SH2', positions=pos)

# 2. Set up the TorchANI calculator
# Note: ANI-2x supports Sulfur (S), while ANI-1x/1ccx does not.
calculator = torchani.models.ANI2x().ase()
atoms.set_calculator(calculator)

# 3. Perform Geometry Optimization
# Using BFGS optimizer to find the energy minimum
opt = BFGS(atoms, trajectory='sulphane_opt.traj')
opt.run(fmax=0.01)  # Stop when forces are below 0.01 eV/Angstrom

print("Optimized Positions:\n", atoms.get_positions())
print("Final Energy:", atoms.get_potential_energy())

