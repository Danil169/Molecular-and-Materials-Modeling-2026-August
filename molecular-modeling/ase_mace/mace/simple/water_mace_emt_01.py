from ase.build import molecule
from ase.optimize import BFGS
from mace.calculators import mace_mp
from ase.calculators.emt import EMT
import numpy as np

# Experimental reference values for water
EXP_BOND_LENGTH = 0.9572  # Å
EXP_BOND_ANGLE = 104.52   # degrees

def print_geometry(positions, title, energy=None):
    """Print geometry parameters for water molecule."""
    # Get distances
    oh1 = np.linalg.norm(positions[0] - positions[1])
    oh2 = np.linalg.norm(positions[0] - positions[2])
    hh = np.linalg.norm(positions[1] - positions[2])
    
    # Calculate bond angle (O-H-H)
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(cos_angle) * 180 / np.pi
    
    print(f"\n{title}")
    print("-" * 40)
    if energy is not None:
        print(f"Potential Energy: {energy:.4f} eV")
    print(f"O-H Bond Length: {oh1:.4f} Å (average: {(oh1+oh2)/2:.4f} Å)")
    print(f"O-H Bond Length 2: {oh2:.4f} Å")
    print(f"H-H Distance: {hh:.4f} Å")
    print(f"Bond Angle: {angle:.2f}°")
    
    # Compare with experimental
    avg_oh = (oh1 + oh2) / 2
    print(f"\nComparison with experiment:")
    print(f"  Bond length deviation: {avg_oh - EXP_BOND_LENGTH:.4f} Å "
          f"({(avg_oh - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.2f}%)")
    print(f"  Bond angle deviation: {angle - EXP_BOND_ANGLE:.2f}° "
          f"({(angle - EXP_BOND_ANGLE)/EXP_BOND_ANGLE*100:.2f}%)")

# ============================================
# 1. MACE Optimization
# ============================================
print("=" * 60)
print("MACE OPTIMIZATION")
print("=" * 60)

water_mace = molecule("H2O")
calc_mace = mace_mp(model="medium", device="cpu")
water_mace.calc = calc_mace

opt_mace = BFGS(water_mace, trajectory="water_opt_mace.traj")
opt_mace.run(fmax=0.01)

print_geometry(water_mace.get_positions(), 
               "MACE Optimized Geometry", 
               water_mace.get_potential_energy())

# ============================================
# 2. EMT Optimization
# ============================================
print("\n" + "=" * 60)
print("EMT OPTIMIZATION")
print("=" * 60)

water_emt = molecule("H2O")
calc_emt = EMT()
water_emt.calc = calc_emt

opt_emt = BFGS(water_emt, trajectory="water_opt_emt.traj")
opt_emt.run(fmax=0.01)

print_geometry(water_emt.get_positions(), 
               "EMT Optimized Geometry", 
               water_emt.get_potential_energy())

# ============================================
# 3. Comparison Summary
# ============================================
print("\n" + "=" * 60)
print("COMPARISON SUMMARY")
print("=" * 60)

# Get final geometries
mace_pos = water_mace.get_positions()
emt_pos = water_emt.get_positions()

# Calculate properties for each
mace_oh = np.mean([np.linalg.norm(mace_pos[0] - mace_pos[1]),
                   np.linalg.norm(mace_pos[0] - mace_pos[2])])
emt_oh = np.mean([np.linalg.norm(emt_pos[0] - emt_pos[1]),
                  np.linalg.norm(emt_pos[0] - emt_pos[2])])

# MACE angle
v1 = mace_pos[1] - mace_pos[0]
v2 = mace_pos[2] - mace_pos[0]
mace_angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))) * 180 / np.pi

# EMT angle
v1 = emt_pos[1] - emt_pos[0]
v2 = emt_pos[2] - emt_pos[0]
emt_angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))) * 180 / np.pi

print(f"{'Method':<15} {'Bond Length (Å)':<20} {'Bond Angle (°)':<20} {'Energy (eV)':<15}")
print("-" * 70)
print(f"{'MACE':<15} {mace_oh:<20.4f} {mace_angle:<20.2f} {water_mace.get_potential_energy():<15.4f}")
print(f"{'EMT':<15} {emt_oh:<20.4f} {emt_angle:<20.2f} {water_emt.get_potential_energy():<15.4f}")
print(f"{'Experimental':<15} {EXP_BOND_LENGTH:<20.4f} {EXP_BOND_ANGLE:<20.2f} {'N/A':<15}")
print("-" * 70)

# Calculate mean absolute errors
print("\nMean Absolute Errors (vs experiment):")
mace_mae_bond = abs(mace_oh - EXP_BOND_LENGTH)
mace_mae_angle = abs(mace_angle - EXP_BOND_ANGLE)
emt_mae_bond = abs(emt_oh - EXP_BOND_LENGTH)
emt_mae_angle = abs(emt_angle - EXP_BOND_ANGLE)

print(f"  MACE:  Bond length: {mace_mae_bond:.4f} Å, Angle: {mace_mae_angle:.2f}°")
print(f"  EMT:   Bond length: {emt_mae_bond:.4f} Å, Angle: {emt_mae_angle:.2f}°")
print(f"\nMACE is {emt_mae_bond/mace_mae_bond:.1f}x more accurate for bond length")
print(f"MACE is {emt_mae_angle/mace_mae_angle:.1f}x more accurate for bond angle")

print("\nOptimization complete! Trajectory files saved as:")
print("  - water_opt_mace.traj (MACE)")
print("  - water_opt_emt.traj (EMT)")
