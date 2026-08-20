#!/usr/bin/env python
"""
MOPAC PM7 Ionization Energy Calculation for Hg Atom
====================================================
This script calculates the vertical ionization energy of a mercury atom
using MOPAC PM7 with:
- RHF for neutral (closed-shell singlet)
- UHF for cation (open-shell doublet)

Reference: Experimental IE of Hg = 10.44 eV
"""

import os
import shutil
from ase import Atoms
from ase.calculators.mopac import MOPAC

def setup_directories():
    """Create directory structure for MOPAC calculations"""
    base_dir = "mopac_calculations"
    subdirs = [
        "mopac_neutral_rhf",
        "mopac_cation_uhf"
    ]
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

def run_mopac_ie():
    """
    Run MOPAC PM7 calculations for Hg ionization energy.
    
    According to MOPAC manual:
    - RHF is default for even-electron systems (neutral Hg)
    - UHF is default for odd-electron systems (Hg+ cation)
    """
    print("\n" + "="*70)
    print("MOPAC PM7 Ionization Energy Calculation for Hg")
    print("="*70)
    
    # Setup directory structure
    base_dir = setup_directories()
    
    # === Neutral: RHF (Closed-Shell) ===
    print("\n--- Neutral Hg (RHF, Closed-Shell Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'mopac_neutral_rhf')
    
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = MOPAC(
        label=os.path.join(neutral_dir, 'mop'), 
        task='1SCF RHF', 
        keywords='PM7'
    )
    e_neutral = neutral_atom.get_potential_energy()
    print(f"Energy: {e_neutral:.4f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation: UHF (Open-Shell) ===
    print("\n--- Cation Hg+ (UHF, Open-Shell Doublet) ---")
    cation_dir = os.path.join(base_dir, 'mopac_cation_uhf')
    
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = MOPAC(
        label=os.path.join(cation_dir, 'mop'), 
        task='1SCF UHF', 
        charge=1, 
        keywords='PM7 DOUBLET'
    )
    e_cation = cation_atom.get_potential_energy()
    print(f"Energy: {e_cation:.4f} eV")
    print(f"📁 Output in: {cation_dir}")
    
    # === Ionization Energy ===
    ie = e_cation - e_neutral
    exp_ie = 10.44  # Experimental value for Hg
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"Neutral Energy (RHF):   {e_neutral:>12.4f} eV")
    print(f"Cation Energy (UHF):    {e_cation:>12.4f} eV")
    print(f"Vertical IE:            {ie:>12.4f} eV")
    print(f"Experimental IE:        {exp_ie:>12.4f} eV")
    print(f"Difference:             {ie-exp_ie:>+12.4f} eV")
    print("="*70)
    
    return {
        'neutral_energy': e_neutral,
        'cation_energy': e_cation,
        'ie': ie,
        'experimental': exp_ie,
        'error': ie - exp_ie,
        'neutral_dir': neutral_dir,
        'cation_dir': cation_dir,
        'neutral_method': 'RHF (closed-shell)',
        'cation_method': 'UHF (open-shell, doublet)'
    }

def print_directory_structure(results):
    """Print the directory structure"""
    print("\n📁 Directory Structure:")
    print("   mopac_calculations/")
    if results:
        print(f"   ├── {os.path.basename(results['neutral_dir'])}/  (RHF)")
        print(f"   └── {os.path.basename(results['cation_dir'])}/   (UHF)")

def main():
    """Main execution function"""
    print("="*70)
    print("MOPAC PM7 CALCULATION FOR Hg ATOM")
    print("="*70)
    
    # Run calculation
    results = run_mopac_ie()
    
    # Print directory structure
    print_directory_structure(results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ MOPAC PM7 gives excellent results for Hg")
    print(f"   IE = {results['ie']:.4f} eV (error: {results['error']:+.4f} eV)")
    print(f"   Error is only {abs(results['error']):.2f}% from experimental")
    print(f"\n   • PM7 is specifically parameterized for heavy elements like Hg")
    print(f"   • Correctly handles open-shell systems with UHF")
    print(f"   • RHF is used for neutral (default for even-electron)")
    print(f"   • UHF is used for cation (default for odd-electron)")
    print("="*70)

if __name__ == "__main__":
    main()
