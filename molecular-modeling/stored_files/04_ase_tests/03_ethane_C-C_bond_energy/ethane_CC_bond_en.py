#
#
#
from ase.build import molecule
from ase.optimize import BFGS
# Replace EMT with a higher-level calculator (e.g., GPAW, ORCA) for accuracy
from ase.calculators.emt import EMT
import numpy as np
#from ase.units import Å

def calculate_ethane_cc_bond():
    # 1. Setup Ethane (C2H6)
    c2h6 = molecule('C2H6')
    c2h6.calc = EMT()
    
    # Relax ethane geometry
    dyn_e = BFGS(c2h6, trajectory='c2h6_geom_opt.traj', logfile=None)
    dyn_e.run(fmax=0.01)
    e_ethane = c2h6.get_potential_energy()
    print(f"Ethane Energy: {e_ethane:.4f} eV")

    #  Compute C-C bond length
    # In ase.build.molecule('C2H6'), atoms 0 and 1 are the carbons
    c1 = c2h6[0]
    c2 = c2h6[1]
    bond_length = c2h6.get_distance(0, 1)

    print(f"Optimized C-C bond length in ethane: {bond_length:.4f} Å")


    # 2. Setup Methyl Radical (CH3)
    # The methyl radical is generally planar-like (D3h symmetry)
    ch3 = molecule('CH3')
    ch3.calc = EMT()
    
    # Relax methyl radical (crucial for accurate BDE)
    # Using spin-polarized calculator is necessary for radicals
    # ch3.set_initial_magnetic_moments([1]) # Uncomment if using DFT
    
    dyn_m = BFGS(ch3, trajectory='ch3_geom_opt.traj', logfile=None)
    dyn_m.run(fmax=0.01)
    e_methyl = ch3.get_potential_energy()
    print(f"Methyl Radical Energy: {e_methyl:.4f} eV")

    # 3. Calculate Bond Dissociation Energy (BDE)
    # BDE = 2 * E(CH3) - E(C2H6)
    # Note: If C2H6 splits into 2 radicals, 
    # BDE = E(fragment1) + E(fragment2) - E(parent)
    
    bond_energy = (2 * e_methyl) - e_ethane
    
    print("-" * 30)
    print(f"Calculated C-C Bond Energy: {bond_energy:.4f} eV")
    print(f"Calculated C-C Bond Energy: {bond_energy * 96.485:.2f} kJ/mol")
    print("-" * 30)

if __name__ == "__main__":
    calculate_ethane_cc_bond()

