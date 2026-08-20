#!/usr/bin/env python
"""
NWChem MP2-ECP Calculation for Hg using ASE parameters with MPI
- Neutral Hg: RHF + MP2 (closed-shell, singlet)
- Cation Hg+: ROHF + MP2 and UHF + MP2 (open-shell, doublet)
- Uses mpirun with number of processors from an input file
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

# MP2 settings
MP2_THRESH = 1.0e-8

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
        "nwchem_cation_rohf_mp2",
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
    print(f"MP2 threshold:        {MP2_THRESH}")
    print(f"Memory:               {MEMORY}")
    print("="*70)
    
    # Setup directories
    base_dir = setup_directories()
    
    # === Neutral Hg (RHF + MP2) ===
    print("\n--- Neutral Hg (RHF + MP2, Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'nwchem_neutral_rhf_mp2')
    
    neutral_input = f'''title "Hg Neutral RHF+MP2 Calculation"
memory {MEMORY}
geometry units angstroms nocenter noautoz
  Hg 0.000000 0.000000 0.000000
end
basis noprint
  Hg library {BASIS_SET}
end
ecp
  * library {ECP_SET}
end
scf
  thresh {SCF_THRESH}
  maxiter {SCF_MAXITER}
end
mp2
  thresh {MP2_THRESH}
end
task scf energy
task mp2 energy'''
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg_neutral'),
        input=neutral_input,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo",
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (ROHF + MP2) ===
    print("\n--- Cation Hg+ (ROHF + MP2, Open-Shell Doublet) ---")
    cation_rohf_dir = os.path.join(base_dir, 'nwchem_cation_rohf_mp2')
    
    cation_rohf_input = f'''title "Hg+ Cation ROHF+MP2 Calculation"
memory {MEMORY}
charge 1
geometry units angstroms nocenter noautoz
  Hg 0.000000 0.000000 0.000000
end
basis noprint
  Hg library {BASIS_SET}
end
ecp
  * library {ECP_SET}
end
scf
  rohf
  doublet
  thresh {SCF_THRESH}
  maxiter {SCF_MAXITER}
end
mp2
  thresh {MP2_THRESH}
end
task scf energy
task mp2 energy'''
    
    hg_cation_rohf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_rohf.calc = NWChem(
        label=os.path.join(cation_rohf_dir, 'hg_cation_rohf'),
        input=cation_rohf_input,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo",
    )
    
    e_cation_rohf = hg_cation_rohf.get_potential_energy()
    print(f"Energy: {e_cation_rohf:.6f} eV")
    print(f"📁 Output in: {cation_rohf_dir}")
    
    # === Cation Hg+ (UHF + MP2) ===
    print("\n--- Cation Hg+ (UHF + MP2, Open-Shell Doublet) ---")
    cation_uhf_dir = os.path.join(base_dir, 'nwchem_cation_uhf_mp2')
    
    cation_uhf_input = f'''title "Hg+ Cation UHF+MP2 Calculation"
memory {MEMORY}
charge 1
geometry units angstroms nocenter noautoz
  Hg 0.000000 0.000000 0.000000
end
basis noprint
  Hg library {BASIS_SET}
end
ecp
  * library {ECP_SET}
end
scf
  uhf
  doublet
  thresh {SCF_THRESH}
  maxiter {SCF_MAXITER}
end
mp2
  thresh {MP2_THRESH}
end
task scf energy
task mp2 energy'''
    
    hg_cation_uhf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_uhf.calc = NWChem(
        label=os.path.join(cation_uhf_dir, 'hg_cation_uhf'),
        input=cation_uhf_input,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo",
    )
    
    e_cation_uhf = hg_cation_uhf.get_potential_energy()
    print(f"Energy: {e_cation_uhf:.6f} eV")
    print(f"📁 Output in: {cation_uhf_dir}")
    
    # === Ionization Energies ===
    ie_rohf = e_cation_rohf - e_neutral
    ie_uhf = e_cation_uhf - e_neutral
    
    print("\n" + "="*70)
    print("RESULTS (MP2 Level)")
    print("="*70)
    print(f"Neutral Energy (RHF+MP2):      {e_neutral:>12.6f} eV")
    print(f"Cation Energy (ROHF+MP2):      {e_cation_rohf:>12.6f} eV")
    print(f"Cation Energy (UHF+MP2):       {e_cation_uhf:>12.6f} eV")
    print("-"*70)
    print(f"Vertical IE (ROHF+MP2 cation): {ie_rohf:>12.6f} eV")
    print(f"Vertical IE (UHF+MP2 cation):  {ie_uhf:>12.6f} eV")
    print(f"Experimental IE:               {EXP_IE:>12.4f} eV")
    print("-"*70)
    print(f"Error (ROHF+MP2 cation):       {ie_rohf-EXP_IE:>+12.4f} eV")
    print(f"Error (UHF+MP2 cation):        {ie_uhf-EXP_IE:>+12.4f} eV")
    print("="*70)
    
    # Print directory structure
    print("\n📁 Directory Structure:")
    print("   nwchem_mp2_calculations/")
    print(f"   ├── {os.path.basename(neutral_dir)}/      (RHF+MP2)")
    print(f"   ├── {os.path.basename(cation_rohf_dir)}/  (ROHF+MP2)")
    print(f"   └── {os.path.basename(cation_uhf_dir)}/   (UHF+MP2)")
    
    return {
        'neutral': e_neutral,
        'cation_rohf': e_cation_rohf,
        'cation_uhf': e_cation_uhf,
        'ie_rohf': ie_rohf,
        'ie_uhf': ie_uhf,
        'neutral_dir': neutral_dir,
        'cation_rohf_dir': cation_rohf_dir,
        'cation_uhf_dir': cation_uhf_dir,
        'basis': BASIS_SET,
        'ecp': ECP_SET,
        'nproc': nproc,
        'method': 'MP2'
    }

def parse_mp2_energy(output_file):
    """Parse MP2 total energy from NWChem output file"""
    if not os.path.exists(output_file):
        return None
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            # Look for MP2 total energy
            for line in content.split('\n'):
                if 'Total MP2 energy' in line or 'MP2 energy' in line:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    for num in numbers:
                        try:
                            val = float(num)
                            if abs(val) > 1.0:  # Energy should be large
                                return val * 27.211386245988  # Hartree -> eV
                        except ValueError:
                            continue
                # Look for SCF energy if MP2 not found
                if 'Total SCF energy' in line:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    for num in numbers:
                        try:
                            val = float(num)
                            if abs(val) > 1.0:
                                return val * 27.211386245988
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Error parsing {output_file}: {e}")
    return None

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
        
        # Cation ROHF files
        cation_rohf_dir = results['cation_rohf_dir']
        if os.path.exists(cation_rohf_dir):
            files = [f for f in os.listdir(cation_rohf_dir) if os.path.isfile(os.path.join(cation_rohf_dir, f))]
            print(f"\n   {os.path.basename(cation_rohf_dir)}/ (ROHF+MP2):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_rohf_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation UHF files
        cation_uhf_dir = results['cation_uhf_dir']
        if os.path.exists(cation_uhf_dir):
            files = [f for f in os.listdir(cation_uhf_dir) if os.path.isfile(os.path.join(cation_uhf_dir, f))]
            print(f"\n   {os.path.basename(cation_uhf_dir)}/ (UHF+MP2):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_uhf_dir, f))
                print(f"      - {f} ({size} bytes)")

def compare_with_scf(scf_results, mp2_results):
    """Compare SCF and MP2 results"""
    print("\n" + "="*70)
    print("COMPARISON: SCF vs MP2")
    print("="*70)
    print(f"{'Method':<20} {'IE (eV)':<15} {'Error (eV)':<15} {'Improvement':<15}")
    print("-"*70)
    
    if scf_results and mp2_results:
        # SCF results
        scf_ie = scf_results['ie_uhf']  # Using UHF from SCF
        scf_error = scf_ie - EXP_IE
        print(f"{'SCF (UHF)':<20} {scf_ie:<15.6f} {scf_error:<+15.4f} {'-':<15}")
        
        # MP2 results
        mp2_ie = mp2_results['ie_uhf']
        mp2_error = mp2_ie - EXP_IE
        improvement = scf_error - mp2_error
        print(f"{'MP2 (UHF)':<20} {mp2_ie:<15.6f} {mp2_error:<+15.4f} {improvement:<+15.4f}")
        
        # ROHF vs UHF at MP2 level
        print("-"*70)
        print(f"{'MP2 (ROHF)':<20} {mp2_results['ie_rohf']:<15.6f} {mp2_results['ie_rohf']-EXP_IE:<+15.4f} {'-':<15}")
        print(f"{'MP2 (UHF)':<20} {mp2_results['ie_uhf']:<15.6f} {mp2_results['ie_uhf']-EXP_IE:<+15.4f} {'-':<15}")
        print("="*70)

def main():
    """Main execution function"""
    print("="*70)
    print("NWChem MP2-ECP CALCULATION FOR Hg ATOM")
    print("With ROHF+MP2 and UHF+MP2 for cation")
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
    print(f"\n   ROHF+MP2 IE = {mp2_results['ie_rohf']:.6f} eV (error: {mp2_results['ie_rohf']-EXP_IE:+.6f} eV)")
    print(f"   UHF+MP2 IE  = {mp2_results['ie_uhf']:.6f} eV (error: {mp2_results['ie_uhf']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: SCF + MP2")
    print(f"   • Neutral: RHF+MP2 (closed-shell, singlet)")
    print(f"   • Cation ROHF: ROHF+MP2 (restricted open-shell, doublet)")
    print(f"   • Cation UHF:  UHF+MP2 (unrestricted open-shell, doublet)")
    print(f"   • Basis set: {mp2_results['basis']}")
    print(f"   • ECP: {mp2_results['ecp']}")
    print(f"   • Processors: {mp2_results['nproc']}")
    print("="*70)

if __name__ == "__main__":
    main()
