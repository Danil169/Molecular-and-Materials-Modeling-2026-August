#!/usr/bin/env python
"""
NWChem TCE-CCSD-ECP Calculation for Hg using ASE parameters with MPI
- Neutral Hg: RHF + TCE-CCSD (closed-shell, singlet)
- Cation Hg+: UHF + TCE-CCSD only (ROHF has file handling issues)
- Uses mpirun with number of processors from an input file
- Only ASE NWChem flags, no raw input strings
- TCE with io ga for optimized parallel memory handling
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

# TCE-CCSD settings
TCE_IO_GA = True  # Optimizes parallel memory handling

# Memory (CCSD via TCE needs more memory)
MEMORY = "8000 mb"

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
    base_dir = "nwchem_tce_ccsd_calculations"
    subdirs = [
        "nwchem_neutral_rhf_tce_ccsd",
        "nwchem_cation_uhf_tce_ccsd"
    ]
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

def get_tce_options():
    """Return TCE options dictionary with io ga if enabled"""
    tce_opts = {'ccsd': True}
    if TCE_IO_GA:
        tce_opts['io'] = 'ga'  # Optimizes parallel memory handling
    return tce_opts

def parse_tce_energy(output_file):
    """Parse TCE-CCSD total energy from NWChem output file"""
    if not os.path.exists(output_file):
        return None
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            # Look for TCE CCSD total energy
            for line in content.split('\n'):
                if 'TCE energy' in line or 'CCSD energy' in line or 'Total CCSD energy' in line:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    for num in numbers:
                        try:
                            val = float(num)
                            if abs(val) > 1.0:  # Energy should be large
                                return val * 27.211386245988  # Hartree -> eV
                        except ValueError:
                            continue
                # Look for SCF energy if TCE energy not found
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

def run_tce_ccsd_calculation(nproc):
    """Run NWChem TCE-CCSD-ECP calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem TCE-CCSD-ECP Calculation for Hg (with MPI)")
    print("="*70)
    print(f"Number of processors: {nproc}")
    print(f"Basis set:            {BASIS_SET}")
    print(f"ECP:                  {ECP_SET}")
    print(f"SCF threshold:        {SCF_THRESH}")
    print(f"SCF maxiter:          {SCF_MAXITER}")
    print(f"TCE io ga:            {TCE_IO_GA}")
    print(f"Memory:               {MEMORY}")
    print("="*70)
    
    # Setup directories
    base_dir = setup_directories()
    
    # Get TCE options
    tce_opts = get_tce_options()
    
    # === Neutral Hg (RHF + TCE-CCSD) ===
    print("\n--- Neutral Hg (RHF + TCE-CCSD, Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'nwchem_neutral_rhf_tce_ccsd')
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg_neutral'),
        theory='tce',           # Use TCE module
        task='energy',          # Single point energy
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        tce=tce_opts,           # TCE options with ccsd and io ga
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (UHF + TCE-CCSD only) ===
    # ROHF+CCSD has file handling issues in NWChem
    print("\n--- Cation Hg+ (UHF + TCE-CCSD, Open-Shell Doublet) ---")
    cation_uhf_dir = os.path.join(base_dir, 'nwchem_cation_uhf_tce_ccsd')
    
    hg_cation_uhf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_uhf.calc = NWChem(
        label=os.path.join(cation_uhf_dir, 'hg_cation_uhf'),
        theory='tce',
        task='energy',
        charge=1,
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'uhf': True, 'doublet': True, 'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        tce=tce_opts,           # TCE options with ccsd and io ga
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_cation_uhf = hg_cation_uhf.get_potential_energy()
    print(f"Energy: {e_cation_uhf:.6f} eV")
    print(f"📁 Output in: {cation_uhf_dir}")
    
    # === Ionization Energies ===
    ie_uhf = e_cation_uhf - e_neutral
    
    print("\n" + "="*70)
    print("RESULTS (TCE-CCSD Level)")
    print("="*70)
    print(f"Neutral Energy (RHF+TCE-CCSD):        {e_neutral:>12.6f} eV")
    print(f"Cation Energy (UHF+TCE-CCSD):         {e_cation_uhf:>12.6f} eV")
    print("-"*70)
    print(f"Vertical IE (UHF+TCE-CCSD cation):    {ie_uhf:>12.6f} eV")
    print(f"Experimental IE:                      {EXP_IE:>12.4f} eV")
    print("-"*70)
    print(f"Error (UHF+TCE-CCSD cation):          {ie_uhf-EXP_IE:>+12.4f} eV")
    print("="*70)
    
    # Print directory structure
    print("\n📁 Directory Structure:")
    print("   nwchem_tce_ccsd_calculations/")
    print(f"   ├── {os.path.basename(neutral_dir)}/        (RHF+TCE-CCSD)")
    print(f"   └── {os.path.basename(cation_uhf_dir)}/     (UHF+TCE-CCSD)")
    
    return {
        'neutral': e_neutral,
        'cation_uhf': e_cation_uhf,
        'ie_uhf': ie_uhf,
        'neutral_dir': neutral_dir,
        'cation_uhf_dir': cation_uhf_dir,
        'basis': BASIS_SET,
        'ecp': ECP_SET,
        'nproc': nproc,
        'method': 'TCE-CCSD',
        'tce_io_ga': TCE_IO_GA
    }

def print_file_summary(results):
    """Print summary of generated files"""
    print("\n📄 Generated Files:")
    
    if results:
        # Neutral files
        neutral_dir = results['neutral_dir']
        if os.path.exists(neutral_dir):
            files = [f for f in os.listdir(neutral_dir) if os.path.isfile(os.path.join(neutral_dir, f))]
            print(f"\n   {os.path.basename(neutral_dir)}/ (RHF+TCE-CCSD):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(neutral_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation UHF files
        cation_uhf_dir = results['cation_uhf_dir']
        if os.path.exists(cation_uhf_dir):
            files = [f for f in os.listdir(cation_uhf_dir) if os.path.isfile(os.path.join(cation_uhf_dir, f))]
            print(f"\n   {os.path.basename(cation_uhf_dir)}/ (UHF+TCE-CCSD):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_uhf_dir, f))
                print(f"      - {f} ({size} bytes)")

def main():
    """Main execution function"""
    print("="*70)
    print("NWChem TCE-CCSD-ECP CALCULATION FOR Hg ATOM")
    print("UHF+TCE-CCSD for cation (ROHF has file handling issues)")
    print("="*70)
    
    # Read number of processors from file
    nproc = read_nproc_from_file(NPROC_FILE)
    
    # Run TCE-CCSD calculation
    tce_ccsd_results = run_tce_ccsd_calculation(nproc)
    
    # Print file summary
    print_file_summary(tce_ccsd_results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem TCE-CCSD-ECP calculations completed")
    print(f"\n   UHF+TCE-CCSD IE = {tce_ccsd_results['ie_uhf']:.6f} eV (error: {tce_ccsd_results['ie_uhf']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: TCE-CCSD (Tensor Contraction Engine - Coupled Cluster Singles and Doubles)")
    print(f"   • TCE io ga: {tce_ccsd_results['tce_io_ga']} (optimized parallel memory handling)")
    print(f"   • Neutral: RHF+TCE-CCSD (closed-shell, singlet)")
    print(f"   • Cation:  UHF+TCE-CCSD (unrestricted open-shell, doublet)")
    print(f"   • Note: ROHF+TCE-CCSD has file handling issues in NWChem")
    print(f"   • Basis set: {tce_ccsd_results['basis']}")
    print(f"   • ECP: {tce_ccsd_results['ecp']}")
    print(f"   • Processors: {tce_ccsd_results['nproc']}")
    print(f"   • Memory: {MEMORY}")
    print("="*70)

if __name__ == "__main__":
    main()
