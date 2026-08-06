import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import read

# 1. Physical Parameters for Argon
sigma_ar = 3.405   # Angstroms
epsilon_ar = 0.010323 # eV (approx 120 K)
r_min_ar = sigma_ar * (2**(1/6)) # Theoretical min (~3.822 Å)

# 2. Setup (Note: Starting positions scaled to Argon size)
n_atoms = 7
np.random.seed(42)
atoms = Atoms('Ar' * n_atoms, positions=np.random.rand(n_atoms, 3) * 5.0)

# 3. Assign Calculator with Argon Units
atoms.calc = LennardJones(sigma=sigma_ar, epsilon=epsilon_ar)

# 4. Search for Global Minimum
# We use a larger dr (0.5 * sigma) for physical units
bh = BasinHopping(atoms, temperature=0.005, dr=1.5, trajectory="ar7_physical.traj")
print("Searching for physical Ar7 minimum...")
bh.run(steps=200)

# 5. Analysis
best_atoms = read("ar7_physical.traj", index="-1")
dist_matrix = best_atoms.get_all_distances()

# We set a threshold relative to sigma (e.g., 1.2 * sigma)
threshold = 1.2 * sigma_ar 
bond_lengths = [dist_matrix[i, j] for i in range(n_atoms) for j in range(i+1, n_atoms) if dist_matrix[i,j] < threshold]

print("\n--- Physical Ar7 Analysis ---")
print(f"Number of bonds:      {len(bond_lengths)} (Should be 15)")
print(f"Mean Bond Length:     {np.mean(bond_lengths):.4f} Å")
print(f"Reference r_min:      {r_min_ar:.4f} Å")
print(f"Deviation:            {((np.mean(bond_lengths)-r_min_ar)/r_min_ar)*100:.2f}%")

