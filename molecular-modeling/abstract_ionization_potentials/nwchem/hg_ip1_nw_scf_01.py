#!/usr/bin/env python
"""
NWChem HF-SCF Calculation for Hg
- Neutral Hg: RHF (closed-shell, singlet)
- Cation Hg+: UHF (open-shell, doublet)
"""

from ase import Atoms
from ase.calculators.nwchem import NWChem
import os

def run_hf_calculation():
    print("="*70)
    print("NWChem HF-SCF Calculation for Hg")
    print("="*70)
    
    # Create directory for output
    os.makedirs('nwchem_scf', exist_ok=True)
    
    # === Neutral Hg (RHF, singlet) ===
    print("\n--- Neutral Hg (RHF, Singlet) ---")
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label='nwchem_scf/hg_neutral',
        theory='scf',
        task='energy',
        basis={'Hg': 'def2-svp'},
        ecp={'Hg': 'def2-svp'},  # Use ECP for Hg
        scf__thresh=1.0e-7,
        scf__maxiter=100,
        title="Hg Neutral RHF Energy Calculation"
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Neutral energy: {e_neutral:.6f} eV")
    
    # === Cation Hg+ (UHF, doublet) ===
    print("\n--- Cation Hg+ (UHF, Doublet) ---")
    hg_cation = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation.calc = NWChem(
        label='nwchem_scf/hg_cation',
        theory='scf',
        task='energy',
        basis={'Hg': 'def2-svp'},
        ecp={'Hg': 'def2-svp'},  # Use ECP for Hg
        charge=1,                  # Net charge +1
        scf__uhf=True,             # UHF flag (becomes "uhf" in input)
        scf__doublet=True,         # Doublet state (becomes "doublet" in input)
        scf__thresh=1.0e-7,
        scf__maxiter=100,
        title="Hg+ Cation UHF Energy Calculation"
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
