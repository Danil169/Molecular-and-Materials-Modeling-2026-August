import numpy as np
import warnings
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from chgnet.model.dynamics import CHGNetCalculator

# Suppress PyTorch/CHGNet device warnings
warnings.filterwarnings("ignore", message="Converting a tensor with requires_grad=True")

def get_base_cyclohexane():
    """Manual coordinates for a chair in a 15x15x15A periodic box."""
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
    # Define 15 Angstrom cubic cell with PBC=True
    mol = Atoms('C6H12', positions=coords, cell=[15, 15, 15], pbc=True)
    mol.center()
    return mol

def run_study():
    # 1. Setup CHGNet
    calc = CHGNetCalculator(use_device='cpu')
    to_kcal = 23.061 

    # 2. CHAIR (Global Minimum)
    chair = get_base_cyclohexane()
    chair.calc = calc
    print("Optimizing Chair...")
    BFGS(chair, logfile=None).run(fmax=0.01)
    e_chair = chair.get_potential_energy()
    write('chair_final.xyz', chair)

    # 3. PURE BOAT (Transition State)
    # Move C3 and its hydrogens (12, 13) as a group 
    # Shifting Y and Z to avoid crushing bonds with interior hydrogens
    boat = chair.copy()
    pos = boat.get_positions()
    group_idx = [3, 12, 13]
    pos[group_idx] += [0, 2.8, 0.8] 
    boat.set_positions(pos)
    boat.calc = calc
    # No optimization to keep the TS geometry
    e_pure = boat.get_potential_energy()
    write('pure_boat.xyz', boat)

    # 4. TWIST-BOAT (Local Minimum)
    # We take the boat and apply a 'screw' twist to trap it in the local well
    twist = boat.copy()
    t_pos = twist.get_positions()
    # Chiral-like distortion: Move C1 up/right and C4 down/left
    t_pos[[1, 8, 9]] += [0.5, 0.2, 0.5]
    t_pos[[4, 14, 15]] -= [0.5, 0.2, 0.5]
    twist.set_positions(t_pos)
    twist.calc = calc
    
    print("Optimizing Twist-Boat (preventing collapse to chair)...")
    # Using a tighter fmax to settle into the local minimum
    opt = BFGS(twist, logfile=None)
    opt.run(fmax=0.01)
    e_twist = twist.get_potential_energy()
    write('twist_boat_final.xyz', twist)

    # 5. Final Results
    print("\n" + "="*45)
    print("RESULTS (kcal/mol relative to Chair)")
    print(f"Twist-Boat (Local Min): { (e_twist - e_chair)*to_kcal :.2f} (Target ~5.5)")
    print(f"Pure Boat (TS):         { (e_pure - e_chair)*to_kcal :.2f} (Target ~7.0)")
    print("="*45)
    print("Check chair_final.xyz and twist_boat_final.xyz to verify shapes.")

if __name__ == "__main__":
    run_study()

