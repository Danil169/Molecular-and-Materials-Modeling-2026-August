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

# 2. Optimizer
# dr=1.5 is the "sweet spot" for Ar7 to jump between different shapes
bh = BasinHopping(atoms, 
                  temperature=0.1, 
                  dr=1.5, 
                  optimizer_logfile=None,
                  trajectory="ar7_search.traj")

# 3. Run search
print("Searching for the global minimum...")
bh.run(steps=300)

# 4. Analyze Trajectory to find the "Winner" and the "When"
traj = read("ar7_search.traj", index=":")
energies = [s.get_potential_energy() for s in traj]

# Find the best energy and the step (index) it occurred at
best_idx = np.argmin(energies)
best_energy = energies[best_idx]
best_atoms = traj[best_idx]

print("-" * 45)
print(f"Absolute Lowest Energy Found: {best_energy:.6f} eV")
print(f"Found at Step:                {best_idx}") # This matches the trajectory frame
print(f"Target Global Minimum:        -16.505385 eV")
print("-" * 45)

# Save the best structure
best_atoms.write("ar7_best.xyz")
print("Best structure saved to 'ar7_best.xyz'")

