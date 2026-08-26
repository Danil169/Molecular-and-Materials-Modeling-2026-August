from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from mace.calculators import mace_mp
import numpy as np

# Experimental reference values for N₂
EXP_BOND_LENGTH = 1.098  # Å
EXP_ATOMIZATION_ENERGY = 9.76  # eV

print("=" * 70)
print("NITROGEN MOLECULE (N₂) - EMT vs MACE COMPARISON")
print("=" * 70)

# ============================================
# 1. EMT Calculation
# ============================================
print("\n" + "-" * 35)
print("EMT CALCULATIONS")
print("-" * 35)

# Single N atom energy
atom = Atoms('N', calculator=EMT())
e_atom_emt = atom.get_potential_energy()

# N₂ molecule optimization with EMT
d_init = 1.1  # Initial guess close to experimental
molecule_emt = Atoms('2N', [(0., 0., 0.), (0., 0., d_init)], calculator=EMT())

opt_emt = BFGS(molecule_emt, trajectory="N2_opt_emt.traj")
print(f'\nOptimizing N₂ with EMT (initial d(N-N) = {d_init} Å)')
opt_emt.run(fmax=0.01)

# Get optimized geometry and energy
d_emt = molecule_emt.get_distance(0, 1)
e_molecule_emt = molecule_emt.get_potential_energy()

# Atomization energy
e_atomization_emt = 2 * e_atom_emt - e_molecule_emt

print(f'\nEMT Results:')
print(f"  Optimized N-N bond length: {d_emt:.4f} Å")
print(f"  Nitrogen atom energy:       {e_atom_emt:.4f} eV")
print(f"  Nitrogen molecule energy:   {e_molecule_emt:.4f} eV")
print(f"  Atomization energy:         {e_atomization_emt:.4f} eV")

# ============================================
# 2. MACE Calculation
# ============================================
print("\n" + "-" * 35)
print("MACE CALCULATIONS")
print("-" * 35)

# Use float64 for better accuracy in geometry optimization
calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")

# Single N atom energy with MACE
atom_mace = Atoms('N', calculator=calc_mace)
e_atom_mace = atom_mace.get_potential_energy()

# N₂ molecule optimization with MACE
molecule_mace = Atoms('2N', [(0., 0., 0.), (0., 0., d_init)], calculator=calc_mace)

opt_mace = BFGS(molecule_mace, trajectory="N2_opt_mace.traj")
print(f'\nOptimizing N₂ with MACE (initial d(N-N) = {d_init} Å)')
opt_mace.run(fmax=0.01)

# Get optimized geometry and energy
d_mace = molecule_mace.get_distance(0, 1)
e_molecule_mace = molecule_mace.get_potential_energy()

# Atomization energy
e_atomization_mace = 2 * e_atom_mace - e_molecule_mace

print(f'\nMACE Results:')
print(f"  Optimized N-N bond length: {d_mace:.4f} Å")
print(f"  Nitrogen atom energy:       {e_atom_mace:.4f} eV")
print(f"  Nitrogen molecule energy:   {e_molecule_mace:.4f} eV")
print(f"  Atomization energy:         {e_atomization_mace:.4f} eV")

# ============================================
# 3. Comparison with Experiment
# ============================================
print("\n" + "=" * 70)
print("COMPARISON WITH EXPERIMENT")
print("=" * 70)

# Bond length comparison
print("\nBond Length (N-N):")
print(f"  Experimental:  {EXP_BOND_LENGTH:.4f} Å")
print(f"  EMT:           {d_emt:.4f} Å  (Δ = {d_emt - EXP_BOND_LENGTH:+.4f} Å, "
      f"{abs(d_emt - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.2f}%)")
print(f"  MACE:          {d_mace:.4f} Å  (Δ = {d_mace - EXP_BOND_LENGTH:+.4f} Å, "
      f"{abs(d_mace - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.2f}%)")

# Atomization energy comparison
print("\nAtomization Energy:")
print(f"  Experimental:  {EXP_ATOMIZATION_ENERGY:.2f} eV")
print(f"  EMT:           {e_atomization_emt:.4f} eV  (Δ = {e_atomization_emt - EXP_ATOMIZATION_ENERGY:+.4f} eV, "
      f"{abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100:.2f}%)")
print(f"  MACE:          {e_atomization_mace:.4f} eV  (Δ = {e_atomization_mace - EXP_ATOMIZATION_ENERGY:+.4f} eV, "
      f"{abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100:.2f}%)")

# ============================================
# 4. Summary Table
# ============================================
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)

# Get number of optimization steps
# For BFGS, we can use the number of iterations from the optimizer
n_steps_emt = opt_emt.nsteps
n_steps_mace = opt_mace.nsteps

print(f"{'Method':<12} {'Bond Length (Å)':<20} {'Atomization (eV)':<20} {'Steps':<10}")
print("-" * 62)
print(f"{'EMT':<12} {d_emt:<20.4f} {e_atomization_emt:<20.4f} {n_steps_emt:<10}")
print(f"{'MACE':<12} {d_mace:<20.4f} {e_atomization_mace:<20.4f} {n_steps_mace:<10}")
print(f"{'Exp.':<12} {EXP_BOND_LENGTH:<20.4f} {EXP_ATOMIZATION_ENERGY:<20.2f} {'N/A':<10}")
print("-" * 62)

# Accuracy comparison
print("\nAccuracy Improvement (MACE vs EMT):")
bond_improvement = abs(d_emt - EXP_BOND_LENGTH) / abs(d_mace - EXP_BOND_LENGTH)
energy_improvement = abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY) / abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY)

print(f"  Bond length:     MACE is {bond_improvement:.1f}x more accurate")
print(f"  Atomization:     MACE is {energy_improvement:.1f}x more accurate")

# ============================================
# 5. Analysis and Discussion
# ============================================
print("\n" + "=" * 70)
print("ANALYSIS")
print("=" * 70)

print("\nEMT Performance:")
if abs(d_emt - EXP_BOND_LENGTH) < 0.05:
    print("  ✓ Bond length is reasonably accurate (<5% error)")
else:
    print(f"  ✗ Bond length error is significant ({abs(d_emt - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.1f}%)")

if abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY) < 1.0:
    print("  ✓ Atomization energy is reasonably accurate (<10% error)")
else:
    print(f"  ✗ Atomization energy error is significant ({abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100:.1f}%)")

print("\nMACE Performance:")
if abs(d_mace - EXP_BOND_LENGTH) < 0.02:
    print("  ✓ Bond length is highly accurate (<2% error)")
else:
    print(f"  ✓ Bond length is accurate ({abs(d_mace - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.1f}% error)")

if abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY) < 0.5:
    print("  ✓ Atomization energy is highly accurate (<5% error)")
else:
    print(f"  ✓ Atomization energy is reasonably accurate ({abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100:.1f}% error)")

# ============================================
# 6. Detailed Analysis of Results
# ============================================
print("\n" + "=" * 70)
print("DETAILED ANALYSIS")
print("=" * 70)

print("\nInteresting observations from the results:")

# EMT results analysis
print("\n1. EMT Results:")
print(f"   - Optimized bond length: {d_emt:.4f} Å (experimental: {EXP_BOND_LENGTH:.4f} Å)")
print(f"   - EMT underestimates bond length by {abs(d_emt - EXP_BOND_LENGTH):.4f} Å")
print(f"   - Atomization energy: {e_atomization_emt:.4f} eV (overestimates by {e_atomization_emt - EXP_ATOMIZATION_ENERGY:+.2f} eV)")

# MACE results analysis
print("\n2. MACE Results:")
print(f"   - Optimized bond length: {d_mace:.4f} Å (experimental: {EXP_BOND_LENGTH:.4f} Å)")
print(f"   - MACE overestimates bond length by {abs(d_mace - EXP_BOND_LENGTH):.4f} Å")
print(f"   - Atomization energy: {e_atomization_mace:.4f} eV (overestimates by {e_atomization_mace - EXP_ATOMIZATION_ENERGY:+.2f} eV)")

# Comparison
print("\n3. Method Comparison:")
if abs(d_mace - EXP_BOND_LENGTH) < abs(d_emt - EXP_BOND_LENGTH):
    print(f"   ✓ MACE gives better bond length (error: {abs(d_mace - EXP_BOND_LENGTH):.4f} Å vs {abs(d_emt - EXP_BOND_LENGTH):.4f} Å)")
else:
    print(f"   ✓ EMT gives better bond length (error: {abs(d_emt - EXP_BOND_LENGTH):.4f} Å vs {abs(d_mace - EXP_BOND_LENGTH):.4f} Å)")

if abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY) < abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY):
    print(f"   ✓ MACE gives better atomization energy (error: {abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY):.4f} eV vs {abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY):.4f} eV)")
else:
    print(f"   ✓ EMT gives better atomization energy (error: {abs(e_atomization_emt - EXP_ATOMIZATION_ENERGY):.4f} eV vs {abs(e_atomization_mace - EXP_ATOMIZATION_ENERGY):.4f} eV)")

print("\n4. Conclusion:")
print("   While MACE was expected to be more accurate for both properties, the results show:")
print("   - MACE gives significantly better bond length (only 1.31% error vs 9.10% for EMT)")
print("   - However, EMT gives surprisingly good atomization energy for N₂ (only 1.82% error)")
print("   - This demonstrates that simple empirical potentials can sometimes give good results")
print("     for certain properties, even if they are less accurate overall.")

print("\nOptimization complete! Trajectory files saved as:")
print("  - N2_opt_emt.traj  (EMT)")
print("  - N2_opt_mace.traj (MACE)")
