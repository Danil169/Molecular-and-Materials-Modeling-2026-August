#!/usr/bin/env python
"""
NWChem SCF-ECP Calculation for Hg using ASE parameters only
- Neutral Hg: RHF (closed-shell, singlet)
- Cation Hg+: UHF (open-shell, doublet)
- No raw input strings - all parameters through ASE interface
"""

from ase import Atoms
from ase.calculators.nwchem import NWChem
import os
import shutil

# ============================================================
# USER-DEFINED PARAMETERS - Modify these as needed
# ============================================================

# Basis set and ECP
BASIS_SET = "def2-svp"
ECP_SET = "def2-svp"

# SCF settings
SCF_THRESH = 1.0e-7
SCF_MAXITER = 100

# Memory (must include units: mb, gb, etc.)
MEMORY = "2000 mb"  # <- Changed from 2000 to "2000 mb"

# Experimental reference
EXP_IE = 10.437  # eV

# ============================================================
# END OF USER PARAMETERS
# ============================================================

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

def run_scf_calculation():
    """Run NWChem SCF-ECP calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem SCF-ECP Calculation for Hg")
    print("="*70)
    print(f"Basis set:     {BASIS_SET}")
    print(f"ECP:           {ECP_SET}")
    print(f"SCF threshold: {SCF_THRESH}")
    print(f"SCF maxiter:   {SCF_MAXITER}")
    print(f"Memory:        {MEMORY}")
    print("="*70)
    
    # Setup directories
    base_dir = setup_directories()
    
    # === Neutral Hg (RHF, singlet) ===
    print("\n--- Neutral Hg (RHF, Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'nwchem_neutral_rhf')
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg_neutral'),
        theory='scf',
        task='energy',
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'thresh':SCF_THRESH, 'maxiter': SCF_MAXITER},
        memory=MEMORY  # Now "2000 mb" with units
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (UHF, doublet) ===
    print("\n--- Cation Hg+ (UHF, Doublet) ---")
    cation_dir = os.path.join(base_dir, 'nwchem_cation_uhf')
    
    hg_cation = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation.calc = NWChem(
        label=os.path.join(cation_dir, 'hg_cation'),
        theory='scf',
        task='energy',
        charge=1,
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'thresh':SCF_THRESH, 'maxiter':SCF_MAXITER,'uhf':True,'doublet':True},
        memory=MEMORY  # Now "2000 mb" with units
    )
    
    e_cation = hg_cation.get_potential_energy()
    print(f"Energy: {e_cation:.6f} eV")
    print(f"📁 Output in: {cation_dir}")
    
    # === Ionization Energy ===
    ie = e_cation - e_neutral
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"Neutral Energy (RHF):  {e_neutral:>12.6f} eV")
    print(f"Cation Energy (UHF):   {e_cation:>12.6f} eV")
    print(f"Vertical IE:           {ie:>12.6f} eV")
    print(f"Experimental IE:       {EXP_IE:>12.4f} eV")
    print(f"Difference:            {ie-EXP_IE:>+12.4f} eV")
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
        'basis': BASIS_SET,
        'ecp': ECP_SET
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
    print("NWChem SCF-ECP CALCULATION FOR Hg ATOM")
    print("Using ASE parameters only (no raw input)")
    print("="*70)
    
    # Run calculation
    results = run_scf_calculation()
    
    # Print file summary
    print_file_summary(results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem SCF-ECP calculation completed")
    print(f"   IE = {results['ie']:.6f} eV (error: {results['ie']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: SCF (Hartree-Fock)")
    print(f"   • Neutral uses RHF (closed-shell, singlet)")
    print(f"   • Cation uses UHF (open-shell, doublet)")
    print(f"   • Basis set: {results['basis']}")
    print(f"   • ECP: {results['ecp']}")
    print("="*70)

if __name__ == "__main__":
    main()
