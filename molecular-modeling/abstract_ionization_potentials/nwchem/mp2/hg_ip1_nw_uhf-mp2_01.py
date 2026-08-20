#!/usr/bin/env python
"""
NWChem MP2-ECP Calculation for Hg using ASE parameters with MPI
- Neutral Hg: RHF + MP2 (closed-shell, singlet)
- Cation Hg+: UHF + MP2 (open-shell, doublet)
- Uses mpirun with number of processors from an input file
- Only ASE NWChem flags, no raw input strings
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

# Memory (MP2 needs more memory)
MEMORY = "4000 mb"

# Experimental reference
EXP_IE = 10.437  # eV

# File containing number of processors
NPROC_FILE = "nproc.txt"

# ============================================================
# END OF USER PARAMETERS
# ============================================================

def read_nproc_from_file(filename="nproc.txt"):
    """
    Read number of processors from a text file.
    The file should contain a single integer (e.g., "4").
    """
    nproc = 1  # Default to 1 if file doesn't exist
    
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as f:
                content = f.read().strip()
                nproc = int(content)
                print(f"✅ Read nproc = {nproc} from {filename}")
                return nproc
        except (ValueError, IOError) as e:
            print(f"⚠️  Error reading {filename}: {e}")
            print(f"   Using default nproc = 1")
            return 1
    else:
        print(f"⚠️  {filename} not found. Creating default file with nproc=1")
        with open(filename, 'w') as f:
            f.write("1\n")
        return 1

def setup_directories():
    """Create directory structure for NWChem calculations"""
    base_dir = "nwchem_mp2_calculations"
    subdirs = [
        "nwchem_neutral_rhf_mp2",
        "nwchem_cation_uhf_mp2"
    ]
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

def run_mp2_calculation(nproc):
    """Run NWChem MP2-ECP calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem MP2-ECP Calculation for Hg (with MPI)")
    print("="*70)
    print(f"Number of processors: {nproc}")
    print(f"Basis set:            {BASIS_SET}")
    print(f"ECP:                  {ECP_SET}")
    print(f"SCF threshold:        {SCF_THRESH}")
    print(f"SCF maxiter:          {SCF_MAXITER}")
    print(f"Memory:               {MEMORY}")
    print("="*70)
    
    # Setup directories
    base_dir = setup_directories()
    
    # === Neutral Hg (RHF + MP2) ===
    print("\n--- Neutral Hg (RHF + MP2, Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'nwchem_neutral_rhf_mp2')
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg_neutral'),
        theory='mp2',
        task='energy',
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (UHF + MP2) ===
    print("\n--- Cation Hg+ (UHF + MP2, Open-Shell Doublet) ---")
    cation_uhf_dir = os.path.join(base_dir, 'nwchem_cation_uhf_mp2')
    
    hg_cation_uhf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_uhf.calc = NWChem(
        label=os.path.join(cation_uhf_dir, 'hg_cation_uhf'),
        theory='mp2',
        task='energy',
        charge=1,
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'uhf': True, 'doublet': True, 'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_cation_uhf = hg_cation_uhf.get_potential_energy()
    print(f"Energy: {e_cation_uhf:.6f} eV")
    print(f"📁 Output in: {cation_uhf_dir}")
    
    # === Ionization Energies ===
    ie_uhf = e_cation_uhf - e_neutral
    
    print("\n" + "="*70)
    print("RESULTS (MP2 Level)")
    print("="*70)
    print(f"Neutral Energy (RHF+MP2):      {e_neutral:>12.6f} eV")
    print(f"Cation Energy (UHF+MP2):       {e_cation_uhf:>12.6f} eV")
    print("-"*70)
    print(f"Vertical IE (UHF+MP2 cation):  {ie_uhf:>12.6f} eV")
    print(f"Experimental IE:               {EXP_IE:>12.4f} eV")
    print("-"*70)
    print(f"Error (UHF+MP2 cation):        {ie_uhf-EXP_IE:>+12.4f} eV")
    print("="*70)
    
    # Print directory structure
    print("\n📁 Directory Structure:")
    print("   nwchem_mp2_calculations/")
    print(f"   ├── {os.path.basename(neutral_dir)}/      (RHF+MP2)")
    print(f"   └── {os.path.basename(cation_uhf_dir)}/   (UHF+MP2)")
    
    return {
        'neutral': e_neutral,
        'cation_uhf': e_cation_uhf,
        'ie_uhf': ie_uhf,
        'neutral_dir': neutral_dir,
        'cation_uhf_dir': cation_uhf_dir,
        'basis': BASIS_SET,
        'ecp': ECP_SET,
        'nproc': nproc,
        'method': 'MP2'
    }

def print_file_summary(results):
    """Print summary of generated files"""
    print("\n📄 Generated Files:")
    
    if results:
        # Neutral files
        neutral_dir = results['neutral_dir']
        if os.path.exists(neutral_dir):
            files = [f for f in os.listdir(neutral_dir) if os.path.isfile(os.path.join(neutral_dir, f))]
            print(f"\n   {os.path.basename(neutral_dir)}/ (RHF+MP2):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(neutral_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation UHF files
        cation_uhf_dir = results['cation_uhf_dir']
        if os.path.exists(cation_uhf_dir):
            files = [f for f in os.listdir(cation_uhf_dir) if os.path.isfile(os.path.join(cation_uhf_dir, f))]
            print(f"\n   {os.path.basename(cation_uhf_dir)}/ (UHF+MP2):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_uhf_dir, f))
                print(f"      - {f} ({size} bytes)")

def main():
    """Main execution function"""
    print("="*70)
    print("NWChem MP2-ECP CALCULATION FOR Hg ATOM")
    print("UHF+MP2 for cation (ROHF not supported for MP2)")
    print("="*70)
    
    # Read number of processors from file
    nproc = read_nproc_from_file(NPROC_FILE)
    
    # Run MP2 calculation
    mp2_results = run_mp2_calculation(nproc)
    
    # Print file summary
    print_file_summary(mp2_results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem MP2-ECP calculations completed")
    print(f"\n   UHF+MP2 IE = {mp2_results['ie_uhf']:.6f} eV (error: {mp2_results['ie_uhf']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: MP2 (second-order Møller-Plesset perturbation theory)")
    print(f"   • Neutral: RHF+MP2 (closed-shell, singlet)")
    print(f"   • Cation:  UHF+MP2 (unrestricted open-shell, doublet)")
    print(f"   • Note: ROHF is NOT supported for MP2 in NWChem")
    print(f"   • Basis set: {mp2_results['basis']}")
    print(f"   • ECP: {mp2_results['ecp']}")
    print(f"   • Processors: {mp2_results['nproc']}")
    print("="*70)

if __name__ == "__main__":
    main()
