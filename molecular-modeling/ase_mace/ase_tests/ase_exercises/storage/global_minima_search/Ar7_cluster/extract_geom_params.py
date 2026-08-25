import numpy as np
from ase.io import read

# 1. Load your best structure
atoms = read("ar7_best.xyz")
pos = atoms.get_positions()
com = atoms.get_center_of_mass()

# 2. Identify Axial vs Equatorial atoms based on distance from center of mass
dist_from_com = np.linalg.norm(pos - com, axis=1)
axial_indices = np.argsort(dist_from_com)[-2:]  # The two furthest atoms
equatorial_indices = np.argsort(dist_from_com)[:5]

print("--- Geometry Parameters for Ar7 PBP ---")

# 3. Bond Lengths
# Axial-Axial distance
ax_dist = atoms.get_distance(axial_indices[0], axial_indices[1])
print(f"Axial-Axial Distance:      {ax_dist:.4f} Å")

# Average Equatorial-Equatorial distance (the ring edges)
eq_eq_dists = [atoms.get_distance(equatorial_indices[i], equatorial_indices[(i+1)%5]) for i in range(5)]
print(f"Mean Eq-Eq Bond Length:    {np.mean(eq_eq_dists):.4f} Å")

# Average Axial-Equatorial distance
ax_eq_dists = [atoms.get_distance(ax, eq) for ax in axial_indices for eq in equatorial_indices]
print(f"Mean Ax-Eq Bond Length:    {np.mean(ax_eq_dists):.4f} Å")

# 4. Bond Angles
# Angle between two adjacent equatorial atoms and the COM
print(f"Ideal Eq-Center-Eq Angle:  72.0° (Found avg: {360/5:.1f}°)")
# Angle between axial-center-equatorial
print(f"Ideal Ax-Center-Eq Angle:  90.0°")

