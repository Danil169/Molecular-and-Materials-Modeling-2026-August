import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from chgnet.model.dynamics import CHGNetCalculator

def get_cyclohexane_chair():
    """
    Manual coordinates for a cyclohexane chair.
    Includes a 15A vacuum box to satisfy CHGNet's periodic requirements.
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
        [0.000, -2.499, 0.108], [0.000, -1.454, -1.328], # H on C3 (Move these)
        [-2.164, -1.250, 0.155], [-1.259, -0.727, 1.328],# H on C4
        [-2.164, 1.250, -0.155], [-1.259, 0.727, -1.328] # H on C5
    ]
    mol = Atoms('C6H12', positions=coords)
    mol.set_cell([15, 15, 15]) # Define box to avoid 'Singular Matrix'
    mol.center()
    return mol

def run_cyclohexane_study():
    # 1. Setup Calculator (CHGNet)
    # This captures organic chemistry much better than EMT.
    calc = CHGNetCalculator(use_device='cpu')

    # 2. Optimize Chair
    chair = get_cyclohexane_chair()
    chair.calc = calc
    
    print("--- Optimizing Chair Conformation ---")
    opt_chair = BFGS(chair, logfile=None)
    opt_chair.run(fmax=0.05)
    e_chair = chair.get_potential_energy()
    write('cyclohexane_chair_optimized.xyz', chair)

    # 3. Build a Better Boat
    # We take the optimized chair and flip Carbon 3 and its Hydrogens up.
    boat = chair.copy()
    pos = boat.get_positions()
    
    # Indices for C3 and its two hydrogens (12, 13)
    # Flipping these ~2.2A converts the chair into a boat.
    target_indices = [3, 12, 13]
    pos[target_indices, 2] += 2.2 
    
    boat.set_positions(pos)
    boat.calc = calc
    
    print("--- Optimizing Boat (Twist-Boat) Conformation ---")
    opt_boat = BFGS(boat, logfile=None)
    opt_boat.run(fmax=0.05)
    e_boat = boat.get_potential_energy()
    write('cyclohexane_boat_optimized.xyz', boat)

    # 4. Results and Comparison
    # Conversion: 1 eV = 23.061 kcal/mol
    diff_ev = e_boat - e_chair
    diff_kcal = diff_ev * 23.061
    
    print("\n" + "="*45)
    print(f"Chair Energy: {e_chair:.4f} eV")
    print(f"Boat Energy:  {e_boat:.4f} eV")
    print(f"Computed Gap: {diff_kcal:.2f} kcal/mol")
    print("-" * 45)
    print("EXPERIMENTAL VALUES:")
    print("Chair -> Twist-Boat (Min): ~5.5 kcal/mol")
    print("Chair -> Pure Boat (TS):   ~7.0 kcal/mol")
    print("="*45)
    print("Geometries saved to .xyz files.")

if __name__ == "__main__":
    run_cyclohexane_study()

