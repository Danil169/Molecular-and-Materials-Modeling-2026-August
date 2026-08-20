import numpy as np
from ase import Atoms
from ase.calculators.mopac import MOPAC
from ase.calculators.xtb import XTB

def calculate_ie(calculator_name, calculator_obj):
    """Calculates the vertical ionization energy for a given ASE calculator."""
    # 1. Define the neutral Hg atom
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = calculator_obj
    
    # 2. Define the cationic Hg+ atom
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = calculator_obj
    
    # Apply charging parameters based on the calculator type
    if calculator_name == 'MOPAC':
        # MOPAC keywords for a radical cation (Charge = 1, Doublet spin multiplicity)
        calculator_obj.set(keywords='1SCF CHARGE=1 DOUBLET C.I.=2')
    elif calculator_name == 'xTB':
        # xTB accepts explicit charge and unpaired electron counts
        calculator_obj.set(charge=1, uhf=1)
        
    try:
        # 3. Compute potential energies (eV)
        e_neutral = neutral_atom.get_potential_energy()
        e_cation = cation_atom.get_potential_energy()
        
        # IE = E(cation) - E(neutral)
        ie = e_cation - e_neutral
        
        print(f"\n--- {calculator_name} Results ---")
        print(f"Neutral Energy: {e_neutral:.4f} eV")
        print(f"Cation Energy:  {e_cation:.4f} eV")
        print(f"Vertical IE:    {ie:.4f} eV")
        return ie
    except Exception as e:
        print(f"\n❌ {calculator_name} calculation failed: {e}")
        return None

if __name__ == "__main__":
    print("Starting Ionization Energy calculations for Hg atom...")
    
    # Initialize separate calculator instances for the neutral run
    # MOPAC default keywords: 1SCF (single point), PM7 parameterization
    mopac_calc = MOPAC(keywords='1SCF PM7')
    xtb_calc = XTB(method='GFN2-xTB')
    
    # Run benchmarks
    ie_mopac = calculate_ie('MOPAC', mopac_calc)
    ie_xtb = calculate_ie('xTB', xtb_calc)

