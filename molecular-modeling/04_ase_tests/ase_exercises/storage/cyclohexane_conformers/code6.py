import numpy as np
import warnings
from ase import Atoms
from ase.optimize import BFGS
from chgnet.model.dynamics import CHGNetCalculator

warnings.filterwarnings("ignore", message="Converting a tensor with requires_grad=True")

def get_conformer(coords):
    """Wraps coordinates in a 15A periodic box for CHGNet."""
    mol = Atoms('C6H12', positions=coords, cell=[15, 15, 15], pbc=True)
    mol.center()
    return mol

def run_study():
    calc = CHGNetCalculator(use_device='cpu')
    to_kcal = 23.061 

    # 1. CHAIR COORDINATES (Global Minimum)
    chair_pos = [
        [0.0, 1.454, 0.228], [1.259, 0.727, -0.228], [1.259, -0.727, 0.228],
        [0.0, -1.454, -0.228], [-1.259, -0.727, 0.228], [-1.259, 0.727, -0.228],
        [0.0, 2.499, -0.108], [0.0, 1.454, 1.328], [2.164, 1.250, 0.155],
        [1.259, 0.727, -1.328], [2.164, -1.250, -0.155], [1.259, -0.727, 1.328],
        [0.0, -2.499, 0.108], [0.0, -1.454, -1.328], [-2.164, -1.250, -0.155],
        [-1.259, -0.727, 1.328], [-2.164, 1.250, 0.155], [-1.259, 0.727, -1.328]
    ]
    chair = get_conformer(chair_pos)
    chair.calc = calc
    print("Optimizing Chair...")
    BFGS(chair, logfile=None).run(fmax=0.05)
    e_chair = chair.get_potential_energy()

    # 2. TWIST-BOAT COORDINATES (Local Minimum ~5.5 kcal/mol)
    # These are pre-relaxed coordinates to prevent collapse to chair
    twist_pos = [
        [1.19, 0.95, 0.40], [-0.30, 1.36, 0.29], [-1.23, 0.20, -0.18],
        [-0.67, -1.16, -0.18], [0.81, -1.13, 0.24], [1.52, 0.14, -0.26],
        [1.77, 1.76, 0.05], [1.25, 0.88, 1.49], [-0.36, 1.62, 1.35],
        [-0.67, 2.25, -0.25], [-2.23, 0.21, 0.25], [-1.34, 0.38, -1.25],
        [-1.18, -1.97, 0.35], [-0.78, -1.35, -1.25], [1.28, -2.03, -0.16],
        [0.86, -1.20, 1.33], [2.56, 0.12, 0.08], [1.53, 0.17, -1.35]
    ]
    twist = get_conformer(twist_pos)
    twist.calc = calc
    print("Optimizing Twist-Boat...")
    BFGS(twist, logfile=None).run(fmax=0.05)
    e_twist = twist.get_potential_energy()

    # 3. PURE BOAT COORDINATES (Transition State ~7.0 kcal/mol)
    # These are pre-aligned to keep the transition state geometry
    boat_pos = [
        [0.00, 1.30, 0.60], [1.25, 0.75, -0.15], [1.25, -0.75, -0.15],
        [0.00, -1.30, 0.60], [-1.25, -0.75, -0.15], [-1.25, 0.75, -0.15],
        [0.00, 2.30, 0.20], [0.00, 1.30, 1.70], [2.15, 1.25, 0.25],
        [1.25, 0.80, -1.25], [2.15, -1.25, 0.25], [1.25, -0.80, -1.25],
        [0.00, -2.30, 0.20], [0.00, -1.30, 1.70], [-2.15, -1.25, 0.25],
        [-1.25, -0.80, -1.25], [-2.15, 1.25, 0.25], [-1.25, 0.80, -1.25]
    ]
    boat = get_conformer(boat_pos)
    boat.calc = calc
    # Pure Boat is a TS; we don't optimize it, or it will become a Twist-Boat
    e_boat = boat.get_potential_energy()

    # Results
    print("\n" + "="*45)
    print("FINAL RESULTS (kcal/mol vs Chair)")
    print(f"Twist-Boat (Local Min): { (e_twist - e_chair)*to_kcal :.2f} (Target ~5.5)")
    print(f"Pure Boat (TS):         { (e_boat - e_chair)*to_kcal :.2f} (Target ~7.0)")
    print("="*45)

if __name__ == "__main__":
    run_study()

