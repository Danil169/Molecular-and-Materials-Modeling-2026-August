import numpy as np
from ase.io import read

# 1. Load the best structure
atoms = read("ar7_best.xyz")
n_atoms = len(atoms)

# 2. Calculate the distance matrix
# This gives the distance between every pair of atoms
dist_matrix = atoms.get_all_distances()

# 3. Identify "Bonds"
# For Argon (LJ), the minimum of the potential is at ~3.82 Å.
# We consider any pair within 4.2 Å to be bonded.
bond_lengths = []
for i in range(n_atoms):
    for j in range(i + 1, n_atoms):
        d = dist_matrix[i, j]
        if d < 4.2:  # Threshold to capture the first shell only
            bond_lengths.append(d)

mean_bond_length = np.mean(bond_lengths)
num_bonds = len(bond_lengths)

# 4. Print Results
print("--- Ar7 Geometric Analysis ---")
print(f"Number of bonds found:      {num_bonds} (Expected for PBP: 15)")
print(f"Mean Bond Length:           {mean_bond_length:.4f} Å")

# 5. Experimental/Reference Comparison
# The Lennard-Jones sigma for Argon is ~3.405 Å. 
# The equilibrium distance (r_m) is sigma * 2^(1/6) ≈ 3.822 Å.
ref_rm = 3.822
deviation = ((mean_bond_length - ref_rm) / ref_rm) * 100

print(f"Reference LJ r_min:         {ref_rm:.4f} Å")
print(f"Deviation from Reference:    {deviation:.2f}%")
print("-" * 30)

if num_bonds == 15:
    print("Structure confirmed: Pentagonal Bipyramid geometry.")
else:
    print(f"Structure check: Found {num_bonds} bonds. Expected 15 for PBP.")

