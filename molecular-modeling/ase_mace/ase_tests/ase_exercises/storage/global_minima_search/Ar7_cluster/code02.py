import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping

# 1. Setup: 7 Argon atoms
n_atoms = 7
np.random.seed(42)
# Start them in a very tight cluster (1.5A box) to ensure they feel the LJ potential immediately
positions = np.random.rand(n_atoms, 3) * 1.5 
atoms = Atoms('Ar' * n_atoms, positions=positions)

# 2. Calculator
atoms.calc = LennardJones()

# 3. Optimizer Settings
# Increased dr to 0.6 to help "push" atoms out of local traps
# temperature 0.05 is usually good for LJ clusters
bh = BasinHopping(atoms, 
                  temperature=0.05, 
                  dr=0.6, 
                  optimizer_logfile=None)

# 4. Run more steps
# 500 steps is overkill for Ar7 but guarantees the global minimum
print("Searching for the global minimum (Target: -16.50538)...")
bh.run(steps=50)

# 5. Results
final_energy = atoms.get_potential_energy()
print("-" * 40)
print(f"Lowest Energy Found: {final_energy:.6f} eV")

if final_energy < -16.50:
    print("Success! You found the global minimum.")
else:
    print("Still in a local minimum. Try increasing 'dr' or 'steps'.")

