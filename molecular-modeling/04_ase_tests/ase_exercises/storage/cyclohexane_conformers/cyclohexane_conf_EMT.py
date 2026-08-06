import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT

def get_cyclohexane_chair():
    """Manual coordinates for cyclohexane chair (C6H12) 
    to avoid KeyError in ASE databases."""
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
    return Atoms('C6H12', positions=coords)

def run_calculation():
    # Experimental values (approximate)
    exp_twist_boat = 5.5  # kcal/mol (Local Minimum)
    exp_pure_boat = 7.0   # kcal/mol (Transition State)

    # 1. Initialize Structure
    chair = get_cyclohexane_chair()
    chair.calc = EMT()

    print("--- Attempting Chair Optimization ---")
    try:
        # EMT calculation will trigger NotImplementedError here because of Hydrogen
        opt = BFGS(chair, logfile=None)
        opt.run(fmax=0.05)
        e_chair = chair.get_potential_energy()

        # 2. Create Boat (flip one carbon up)
        boat = chair.copy()
        pos = boat.get_positions()
        pos[3, 2] += 2.0  # Flip the 'foot' up
        boat.set_positions(pos)
        boat.calc = EMT()

        print("--- Attempting Boat Optimization ---")
        opt_b = BFGS(boat, logfile=None)
        opt_b.run(fmax=0.05)
        e_boat = boat.get_potential_energy()

        # 3. Results
        diff_kcal = (e_boat - e_chair) * 23.061
        print(f"\nComputed Difference: {diff_kcal:.2f} kcal/mol")

    except NotImplementedError:
        print("\n[STOP] EMT calculator cannot handle Hydrogen (H).")
        print("To get a real result, install xtb: 'pip install xtb-python'")
        print("And change 'EMT()' to 'XTB()' in the script.")
    
    # 4. Experimental Printout
    print("\n" + "="*35)
    print("EXPERIMENTAL ENERGY DIFFERENCES")
    print(f"Chair -> Twist-Boat: ~{exp_twist_boat} kcal/mol")
    print(f"Chair -> Boat (TS):  ~{exp_pure_boat} kcal/mol")
    print("="*35)

if __name__ == "__main__":
    run_calculation()

