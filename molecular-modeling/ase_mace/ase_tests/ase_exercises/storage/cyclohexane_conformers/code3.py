import numpy as np
import warnings
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from chgnet.model.dynamics import CHGNetCalculator

# Suppress the PyTorch UserWarning for a cleaner output
warnings.filterwarnings("ignore", message="Converting a tensor with requires_grad=True")

def get_cyclohexane_chair():
    coords = [
        [0.000, 1.454, 0.228],   [1.259, 0.727, -0.228],
        [1.259, -0.727, 0.228],  [0.000, -1.454, -0.228],
        [-1.259, -0.727, 0.228], [-1.259, 0.727, -0.228],
        [0.000, 2.499, -0.108],  [0.000, 1.454, 1.328],
        [2.164, 1.250, 0.155],   [1.259, 0.727, -1.328],
        [2.164, -1.250, -0.155], [1.259, -0.727, 1.328],
        [0.000, -2.499, 0.108],  [0.000, -1.454, -1.328],
        [-2.164, -1.250, 0.155], [-1.259, -0.727, 1.328],
        [-2.164, 1.250, -0.155], [-1.259, 0.727, -1.328]
    ]
    mol = Atoms('C6H12', positions=coords)
    mol.set_cell([15, 15, 15])
    mol.center()
    return mol

def print_dihedrals(atoms, title):
    """Prints the 6 endocyclic C-C-C-C dihedral angles."""
    print(f"\n--- {title} Dihedral Angles (Degrees) ---")
    # Ring indices: 0-1-2-3-4-5 and wrap back
    idx = [0, 1, 2, 3, 4, 5, 0, 1, 2] 
    for i in range(6):
        # Use * to unpack the four indices into separate arguments
        d = atoms.get_dihedral(*idx[i:i+4])
        if d > 180: d -= 360
        print(f"C{idx[i]}-C{idx[i+1]}-C{idx[i+2]}-C{idx[i+3]}: {d:6.1f}°")

def run_cyclohexane_study():
    calc = CHGNetCalculator(use_device='cpu')

    # 1. Chair
    chair = get_cyclohexane_chair()
    chair.calc = calc
    print("Optimizing Chair...")
    opt_chair = BFGS(chair, logfile=None)
    opt_chair.run(fmax=0.05)
    e_chair = chair.get_potential_energy()
    write('chair_final.xyz', chair)
    print_dihedrals(chair, "Chair")

    # 2. Boat
    boat = chair.copy()
    pos = boat.get_positions()
    # Target C3 (index 3) and its hydrogens (indices 12, 13)
    target = [3, 12, 13]
    pos[target, 2] += 2.2 
    boat.set_positions(pos)
    boat.calc = calc
    
    print("\nOptimizing Boat...")
    opt_boat = BFGS(boat, logfile=None)
    opt_boat.run(fmax=0.05)
    e_boat = boat.get_potential_energy()
    write('boat_final.xyz', boat)
    print_dihedrals(boat, "Twist-Boat")

    # 3. Final Summary
    diff_kcal = (e_boat - e_chair) * 23.061
    print("\n" + "="*45)
    print(f"Energy Gap: {diff_kcal:.2f} kcal/mol")
    print(f"Target (Exp): ~5.5 kcal/mol")
    print("="*45)

if __name__ == "__main__":
    run_cyclohexane_study()

