from ase.build import molecule
from ase.optimize import BFGS
from mace.calculators import mace_mp

# 1. Take water from the internal ASE database
water = molecule("H2O")

# 2. Attach the MACE calculator
calc = mace_mp(model="medium", device="cpu")
water.calc = calc

# 3. Set up and run the geometry optimization
opt = BFGS(water, trajectory="water_opt.traj")
opt.run(fmax=0.01)

# 4. Print optimized results
print("\nOptimization Complete!")
print(f"Final Potential Energy: {water.get_potential_energy():.4f} eV")
print("Optimized Positions (Å):")
print(water.get_positions())

# Calculate and print the final O-H bond length
bond_length = water.get_distance(0, 1)
print(f"Optimized O-H Bond Length: {bond_length:.4f} Å")

