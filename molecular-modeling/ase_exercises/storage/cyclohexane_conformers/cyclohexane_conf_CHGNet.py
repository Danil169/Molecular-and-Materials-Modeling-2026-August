import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
# This is the missing import causing your NameError:
from chgnet.model.dynamics import CHGNetCalculator

def get_cyclohexane_chair():
    """Manual coordinates with a defined unit cell to avoid Singular Matrix error."""
    coords = [
        [0.000, 1.454, 0.228],  [1.259, 0.727, -0.228],
        [1.259, -0.727, 0.228], [0.000, -1.454, -0.228],
        [-1.259, -0.727, 0.228], [-1.259, 0.727, -0.228],
        [0.000, 2.499, -0.108], [0.000, 1.454, 1.328],
        [2.164, 1.250, 0.155],  [1.259, 0.727, -1.328],
        [2.164, -1.250, -0.155], [1.259, -0.727, 1.328],
        [0.000, -2.499, 0.108], [0.000, -1.454, -1.328],
        [-2.164, -1.250, 0.155], [-1.259, -0.727, 1.328],
        [-2.164, 1.250, -0.155], [-1.259, 0.727, -1.328]
    ]
    # Use a large vacuum box (15A) to satisfy CHGNet's periodic requirements
    mol = Atoms('C6H12', positions=coords)
    mol.set_cell([15, 15, 15])
    mol.center()
    return mol

def run_chgnet_calculation():
    # 1. Initialize Calculator
    calc = CHGNetCalculator(use_device='cpu') 

    # 2. Optimize Chair
    chair = get_cyclohexane_chair()
    chair.calc = calc
    
    print("--- Optimizing Chair with CHGNet ---")
    opt_chair = BFGS(chair, logfile=None)
    opt_chair.run(fmax=0.05)
    e_chair = chair.get_potential_energy()
    
    # Save Chair Geometry
    write('cyclohexane_chair_optimized.xyz', chair)
    print("Saved: cyclohexane_chair_optimized.xyz")

    # 3. Create and Optimize Boat
    boat = chair.copy()
    pos = boat.get_positions()
    
    # Manually 'flip' the footrest carbon to reach the boat basin
    # We move the carbon at index 3 and its hydrogens up
    pos[3, 2] += 2.0  # C4
    pos[12, 2] += 2.0 # H on C4
    pos[13, 2] += 2.0 # H on C4
    
    boat.set_positions(pos)
    boat.calc = calc
    
    print("--- Optimizing Boat with CHGNet ---")
    opt_boat = BFGS(boat, logfile=None)
    opt_boat.run(fmax=0.05)
    e_boat = boat.get_potential_energy()
    
    # Save Boat Geometry
    write('cyclohexane_boat_optimized.xyz', boat)
    print("Saved: cyclohexane_boat_optimized.xyz")

    # 4. Results & Experimental Comparison
    diff_ev = e_boat - e_chair
    diff_kcal = diff_ev * 23.061
    
    print("\n" + "="*45)
    print(f"CHGNet Chair Energy: {e_chair:.4f} eV")
    print(f"CHGNet Boat Energy:  {e_boat:.4f} eV")
    print(f"Computed Difference: {diff_kcal:.2f} kcal/mol")
    print("-" * 45)
    print("EXPERIMENTAL VALUES:")
    print("Chair -> Twist-Boat: ~5.5 kcal/mol")
    print("Chair -> Boat (TS):  ~7.0 kcal/mol")
    print("="*45)

if __name__ == "__main__":
    run_chgnet_calculation()

