import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import read

# 1. Setup: 7 Argon atoms
n_atoms = 7
np.random.seed(42)
positions = np.random.rand(n_atoms, 3) * 1.5 
atoms = Atoms('Ar' * n_atoms, positions=positions)
atoms.calc = LennardJones()

# 2. Optimizer with Trajectory
# Adding 'trajectory="ar7_search.traj"' saves the progress
bh = BasinHopping(atoms, 
                  temperature=0.05, 
                  dr=0.6, 
                  optimizer_logfile=None,
                  trajectory="ar7_search.traj") # <--- This saves the file

# 3. Run search
print("Searching and saving to ar7_search.traj...")
bh.run(steps=50)

# 4. Final Energy
final_energy = atoms.get_potential_energy()
print("-" * 40)
print(f"Lowest Energy Found: {final_energy:.6f} eV")

# 5. How to check the trajectory
# You can read the last (best) structure from the file like this:
# best_structure = read("ar7_search.traj", index="-1")

