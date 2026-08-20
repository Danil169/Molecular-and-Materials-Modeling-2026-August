from ase import Atoms
from ase.calculators.mopac import MOPAC
from xtb.ase.calculator import XTB

def run_mopac_ie():
    print("\n--- Running MOPAC (PM7) ---")
    
    # 1. Neutral Calculation 
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    # Charge=0 (Singlet is implicit for even electrons, but can be specified via keywords)
    neutral_atom.calc = MOPAC(label='mopac_neutral/mop', task='1SCF GRADIENTS', keywords='PM7')
    e_neutral = neutral_atom.get_potential_energy()
    
    # 2. Cation Calculation
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    # Use explicit charge parameter + specify doublet spin state via keywords 
    cation_atom.calc = MOPAC(label='mopac_cation/mop', task='1SCF GRADIENTS', charge=1, keywords='PM7 DOUBLET')
    e_cation = cation_atom.get_potential_energy()
    
    ie = e_cation - e_neutral
    print(f"Neutral Energy: {e_neutral:.4f} eV")
    print(f"Cation Energy:  {e_cation:.4f} eV")
    print(f"Vertical IE:    {ie:.4f} eV")

def run_xtb_ie():
    print("\n--- Running xTB (GFN2-xTB) ---")
    
    # 1. Neutral Calculation
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = XTB(method='GFN2-xTB', directory='xtb_neutral')
    e_neutral = neutral_atom.get_potential_energy()
    
    # 2. Cation Calculation
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    # Pass charge and unpaired electrons (uhf=1 for doublet)
    cation_atom.calc = XTB(method='GFN2-xTB', charge=1, uhf=1, directory='xtb_cation')
    e_cation = cation_atom.get_potential_energy()
    
    ie = e_cation - e_neutral
    print(f"Neutral Energy: {e_neutral:.4f} eV")
    print(f"Cation Energy:  {e_cation:.4f} eV")
    print(f"Vertical IE:    {ie:.4f} eV")

if __name__ == "__main__":
    print("Starting Ionization Energy calculations for Hg atom...")
    run_mopac_ie()
    run_xtb_ie()

