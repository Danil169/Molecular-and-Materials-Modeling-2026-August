import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import read

# 1. Setup
n_atoms = 7
np.random.seed(42)
positions = np.random.rand(n_atoms, 3) * 2.0 
atoms = Atoms('Ar' * n_atoms, positions=positions)
atoms.calc = LennardJones()

# 2. Optimizer - Increased dr to 1.5 to cross energy barriers
bh = BasinHopping(atoms, 
                  temperature=0.1, # Slightly higher temp to escape traps
                  dr=1.5,          # Much larger jumps
                  optimizer_logfile=None,
                  trajectory="ar7_search.traj")

# 3. Run
print("Searching for the pentagonal bipyramid...")
bh.run(steps=200)

# 4. IMPORTANT: Get the absolute best structure found, not the current one
# The trajectory file stores every 'accepted' step.
traj = read("ar7_search.traj", index=":")
energies = [s.get_potential_energy() for s in traj]
best_idx = np.argmin(energies)
best_atoms = traj[best_idx]

print("-" * 40)
print(f"Absolute Lowest Energy Found: {best_atoms.get_potential_energy():.6f} eV")
print(f"Target: -16.505385 eV")

# Save only the winner
best_atoms.write("ar7_global_min.xyz")

