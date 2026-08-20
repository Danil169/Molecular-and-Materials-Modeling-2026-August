#!/usr/bin/env python
"""
Pure HF-SCF calculation for Hg ionization energy using ASE-NWChem
- Neutral Hg: RHF (closed-shell, singlet)
- Cation Hg+: UHF (open-shell, doublet)
"""

from ase import Atoms
from ase.calculators.nwchem import NWChem

def run_hf_calculation():
    print("="*70)
    print("NWChem HF-SCF Calculation for Hg")
    print("="*70)
    
    # === Neutral Hg (RHF, singlet) ===
    print("\n--- Neutral Hg (RHF, Singlet) ---")
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        theory='scf',
        task='energy',
        basis='def2-tzvp',      # Triple-zeta valence with polarization
        # singlet is default for closed-shell
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Neutral energy: {e_neutral:.6f} eV")
    
    # === Cation Hg+ (UHF, doublet) ===
    print("\n--- Cation Hg+ (UHF, Doublet) ---")
    hg_cation = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation.calc = NWChem(
        theory='scf',
        task='energy',
        basis='def2-tzvp',
        charge=1,               # Net charge +1
        uhf=True,               # Unrestricted HF for open-shell
        mult=2,                 # Doublet state
    )
    
    e_cation = hg_cation.get_potential_energy()
    print(f"Cation energy:   {e_cation:.6f} eV")
    
    # === Ionization Energy ===
    ie = e_cation - e_neutral
    exp_ie = 10.44  # Experimental value
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"Neutral energy (RHF):  {e_neutral:>12.6f} eV")
    print(f"Cation energy (UHF):   {e_cation:>12.6f} eV")
    print(f"Vertical IE:           {ie:>12.6f} eV")
    print(f"Experimental IE:       {exp_ie:>12.4f} eV")
    print(f"Difference:            {ie-exp_ie:>+12.4f} eV")
    print("="*70)
    
    return {
        'neutral': e_neutral,
        'cation': e_cation,
        'ie': ie,
    }

if __name__ == "__main__":
    results = run_hf_calculation()
