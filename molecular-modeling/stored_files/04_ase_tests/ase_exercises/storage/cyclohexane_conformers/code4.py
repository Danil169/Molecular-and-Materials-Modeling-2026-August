import numpy as np
import warnings
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from chgnet.model.dynamics import CHGNetCalculator

# Suppress PyTorch warnings for a cleaner terminal output
warnings.filterwarnings("ignore", message="Converting a tensor with requires_grad=True")

def get_cyclohexane_chair():
    """
    Creates cyclohexane chair coordinates with a 15A vacuum box.
    CHGNet requires a non-singular unit cell for its internal calculations.
    """
    coords = [
        # Carbons (Indices 0-5)
        [0.000, 1.454, 0.228],   # C0
        [1.259, 0.727, -0.228],  # C1
        [1.259, -0.727, 0.228],  # C2
        [0.000, -1.454, -0.228], # C3 (The 'footrest' to flip)
        [-1.259, -0.727, 0.228], # C4
        [-1.259, 0.727, -0.228], # C5
        # Hydrogens (Indices 6-17)
        [0.000, 2.499, -0.108], [0.000, 1.454, 1.328],   # H on C0
        [2.164, 1.250, 0.155],  [1.259, 0.727, -1.328],  # H on C1
        [2.164, -1.250, -0.155], [1.259, -0.727, 1.328], # H on C2
        [0.000, -2.499, 0.108], [0.000, -1.454, -1.328], # H on C3
        [-2.164, -1.250, 0.155], [-1.259, -0.727, 1.328],# H on C4
        [-2.164, 1.250, -0.155], [-1.259, 0.727, -1.328] # H on C5
    ]
    mol = Atoms('C6H12', positions=coords)
    mol.set_cell([15, 15, 15]) 
    mol.center()
    return mol

def analyze_geometry(atoms, title):
    """Prints energy and ring dihedrals for structural verification."""
    print(f"\n--- {title} Analysis ---")
    # Ring indices: 0-1-2-3-4-5 and wrapping back for dihedrals
    idx = [0, 1, 2, 3, 4, 5, 0, 1, 2]
    print("Dihedral Angles (C-C-C-C):")
    for i in range(6):
        # The '*' unpacks the list into the 4 required arguments for get_dihedral
        d = atoms.get_dihedral(*idx[i:i+4])
        if d > 180: d -= 360
        print(f"  C{idx[i]}-C{idx[i+1]}-C{idx[i+2]}-C{idx[i+3]}: {d:6.1f}°")

def run_study():
    # 1. Initialize Calculator
    calc = CHGNetCalculator(use_device='cpu')

    # 2. Process Chair
    chair = get_cyclohexane_chair()
    chair.calc = calc
    print("Optimizing Chair Conformation...")
    opt_chair = BFGS(chair, logfile=None)
    opt_chair.run(fmax=0.05)
    e_chair = chair.get_potential_energy()
    write('chair_optimized.xyz', chair)
    analyze_geometry(chair, "CHAIR")

    # 3. Process Boat
    # Start from optimized chair and flip C3 + its Hydrogens (12, 13)
    boat = chair.copy()
    pos = boat.get_positions()
    flip_indices = [3, 12, 13]
    pos[flip_indices, 2] += 2.2 
    boat.set_positions(pos)
    boat.calc = calc
    
    print("\nOptimizing Boat/Twist-Boat Conformation...")
    opt_boat = BFGS(boat, logfile=None)
    opt_boat.run(fmax=0.05)
    e_boat = boat.get_potential_energy()
    write('boat_optimized.xyz', boat)
    analyze_geometry(boat, "BOAT / TWIST-BOAT")

    # 4. Results
    diff_kcal = (e_boat - e_chair) * 23.061
    print("\n" + "="*45)
    print(f"Chair Energy: {e_chair:.4f} eV")
    print(f"Boat Energy:  {e_boat:.4f} eV")
    print(f"Energy Gap:   {diff_kcal:.2f} kcal/mol")
    print("-" * 45)
    print("EXPERIMENTAL VALUES:")
    print("Chair -> Twist-Boat: ~5.5 kcal/mol")
    print("Chair -> Pure Boat:  ~7.0 kcal/mol")
    print("="*45)

if __name__ == "__main__":
    run_study()

