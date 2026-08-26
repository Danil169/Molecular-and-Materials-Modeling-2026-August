import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping  # Updated import path
from ase.units import kB

# 1. Setup: Create 7 Argon atoms in a random starting configuration
n_atoms = 7
np.random.seed(42)
# Starting atoms in a small 3x3x3 A box to ensure interaction
positions = np.random.rand(n_atoms, 3) * 3.0 
atoms = Atoms('Ar' * n_atoms, positions=positions)

# 2. Assign Calculator: Use CamelCase LennardJones
atoms.calc = LennardJones()

# 3. Initialize Global Optimizer:
# - temperature: Controls structure acceptance. Using units like kB is common.
# - dr: Maximum displacement per step
bh = BasinHopping(atoms, 
                  temperature=0.05, 
                  dr=0.5, 
                  optimizer_logfile=None)

# 4. Execute the search
print("Starting search for Ar7 global minimum...")
bh.run(steps=100)

# 5. Output Results
final_energy = atoms.get_potential_energy()
print("-" * 30)
print(f"Lowest Energy Found: {final_energy:.6f} eV")
print(f"Target Global Minimum: ~ -16.505 eV")

# To visualize the result if you have a display:
# from ase.visualize import view
# view(atoms)

