#!/usr/bin/env python
"""
NWChem HF-SCF Calculation for Hg using raw input
- Neutral Hg: RHF (closed-shell, singlet)
- Cation Hg+: UHF (open-shell, doublet)
"""

from ase import Atoms
from ase.calculators.nwchem import NWChem
import os
import shutil

def setup_directories():
    """Create directory structure for NWChem calculations"""
    base_dir = "nwchem_calculations"
    subdirs = [
        "nwchem_neutral_rhf",
        "nwchem_cation_uhf"
    ]
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

def run_hf_calculation():
    """Run NWChem HF-SCF calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem HF-SCF Calculation for Hg (Raw Input)")
    print("="*70)
    
    # Setup directory structure
    base_dir = setup_directories()
    
    # === Neutral Hg (RHF, singlet) ===
    print("\n--- Neutral Hg (RHF, Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'nwchem_neutral_rhf')
    
    neutral_input = '''title "Hg Neutral RHF Energy Calculation"
geometry units angstroms nocenter noautoz
  Hg 0.000000 0.000000 0.000000
end
basis
  * library def2-svp
end
ecp
  * library def2-svp
end
scf
  thresh 1.0e-7
  maxiter 100
end
task scf energy'''
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg_neutral'),
        input=neutral_input,
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (UHF, doublet) ===
    print("\n--- Cation Hg+ (UHF, Doublet) ---")
    cation_dir = os.path.join(base_dir, 'nwchem_cation_uhf')
    
    cation_input = '''title "Hg+ Cation UHF Energy Calculation"
charge 1
geometry units angstroms nocenter noautoz
  Hg 0.000000 0.000000 0.000000
end
basis
  * library def2-svp
end
ecp
  * library def2-svp
end
scf
  uhf
  doublet
  thresh 1.0e-7
  maxiter 100
end
task scf energy'''
    
    hg_cation = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation.calc = NWChem(
        label=os.path.join(cation_dir, 'hg_cation'),
        input=cation_input,
    )
    
    e_cation = hg_cation.get_potential_energy()
    print(f"Energy: {e_cation:.6f} eV")
    print(f"📁 Output in: {cation_dir}")
    
    # === Ionization Energy ===
    ie = e_cation - e_neutral
    exp_ie = 10.44  # Experimental value for Hg
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"Neutral Energy (RHF):  {e_neutral:>12.6f} eV")
    print(f"Cation Energy (UHF):   {e_cation:>12.6f} eV")
    print(f"Vertical IE:           {ie:>12.6f} eV")
    print(f"Experimental IE:       {exp_ie:>12.4f} eV")
    print(f"Difference:            {ie-exp_ie:>+12.4f} eV")
    print("="*70)
    
    # Print directory structure
    print("\n📁 Directory Structure:")
    print("   nwchem_calculations/")
    print(f"   ├── {os.path.basename(neutral_dir)}/  (RHF)")
    print(f"   └── {os.path.basename(cation_dir)}/   (UHF)")
    
    return {
        'neutral': e_neutral,
        'cation': e_cation,
        'ie': ie,
        'neutral_dir': neutral_dir,
        'cation_dir': cation_dir,
        'neutral_method': 'RHF (closed-shell)',
        'cation_method': 'UHF (open-shell, doublet)'
    }

def print_file_summary(results):
    """Print summary of generated files"""
    print("\n📄 Generated Files:")
    
    if results:
        # Neutral files
        neutral_dir = results['neutral_dir']
        if os.path.exists(neutral_dir):
            files = [f for f in os.listdir(neutral_dir) if os.path.isfile(os.path.join(neutral_dir, f))]
            print(f"\n   {os.path.basename(neutral_dir)}/ (RHF):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(neutral_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation files
        cation_dir = results['cation_dir']
        if os.path.exists(cation_dir):
            files = [f for f in os.listdir(cation_dir) if os.path.isfile(os.path.join(cation_dir, f))]
            print(f"\n   {os.path.basename(cation_dir)}/ (UHF):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_dir, f))
                print(f"      - {f} ({size} bytes)")

def main():
    """Main execution function"""
    print("="*70)
    print("NWChem SCF CALCULATION FOR Hg ATOM")
    print("Using raw input strings for reliability")
    print("="*70)
    
    # Run calculation
    results = run_hf_calculation()
    
    # Print file summary
    print_file_summary(results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem HF-SCF calculation completed")
    print(f"   IE = {results['ie']:.6f} eV (error: {results['ie']-10.44:+.6f} eV)")
    print(f"\n   • Neutral uses RHF (closed-shell, singlet)")
    print(f"   • Cation uses UHF (open-shell, doublet)")
    print(f"   • ECP: def2-svp for both states")
    print(f"   • Basis set: def2-svp with ECP")
    print("="*70)

if __name__ == "__main__":
    main()
