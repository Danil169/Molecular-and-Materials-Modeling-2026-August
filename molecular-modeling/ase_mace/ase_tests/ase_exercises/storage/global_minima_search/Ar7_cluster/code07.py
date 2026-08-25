import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import read

# 1. Physics Parameters
sigma_ar = 3.405
epsilon_ar = 0.010323
r_min_ar = sigma_ar * (2**(1/6)) # ~3.822 Å

# 2. Setup - Start them closer together (3x3x3 box)
n_atoms = 7
np.random.seed(42)
positions = np.random.rand(n_atoms, 3) * 3.0 
atoms = Atoms('Ar' * n_atoms, positions=positions)
atoms.calc = LennardJones(sigma=sigma_ar, epsilon=epsilon_ar)

# 3. Search - 500 steps ensures we find the global minimum
bh = BasinHopping(atoms, 
                  temperature=0.01, # Keep it cool enough to settle
                  dr=2.0,           # Big jumps to escape local traps
                  optimizer_logfile=None,
                  trajectory="ar7_final.traj")

print("Searching for the 15-bond global minimum...")
bh.run(steps=500)

# 4. Analysis of the BEST structure found
traj = read("ar7_final.traj", index=":")
energies = [s.get_potential_energy() for s in traj]
best_atoms = traj[np.argmin(energies)]

dist_matrix = best_atoms.get_all_distances()
# Use a slightly more generous threshold for physical Ar (~4.3 Å)
threshold = 4.3 
bond_lengths = [dist_matrix[i, j] for i in range(n_atoms) for j in range(i+1, n_atoms) if dist_matrix[i,j] < threshold]

print("\n--- Final Physical Ar7 Analysis ---")
print(f"Number of bonds:      {len(bond_lengths)} (Target: 15)")
print(f"Mean Bond Length:     {np.mean(bond_lengths):.4f} Å")
print(f"Reference r_min:      {r_min_ar:.4f} Å")
print(f"Total Potential Energy: {best_atoms.get_potential_energy():.6f} eV")

