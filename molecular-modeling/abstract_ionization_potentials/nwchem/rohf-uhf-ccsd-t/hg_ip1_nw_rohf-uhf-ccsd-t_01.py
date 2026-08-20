#!/usr/bin/env python
"""
NWChem TCE-CCSD(T)-ECP Calculation for Hg using ASE parameters with MPI
- Neutral Hg: RHF + TCE-CCSD(T) (closed-shell, singlet)
- Cation Hg+: ROHF + TCE-CCSD(T) and UHF + TCE-CCSD(T) (open-shell, doublet)
- Uses mpirun with number of processors from an input file
- Only ASE directives, no raw input strings
- No permanent_dir or scratch_dir specified (uses NWChem defaults)
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

# TCE-CCSD(T) settings
TCE_IO_GA = True  # Optimizes parallel memory handling
TCE_CCSDT = True  # Enable CCSD(T) - includes perturbative triples

# Memory (CCSD(T) needs more memory)
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
    base_dir = "nwchem_tce_ccsdt"
    subdirs = [
        "neutral",
        "cation_rohf",
        "cation_uhf"
    ]
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

def get_tce_options():
    """Return TCE options dictionary with CCSD(T) and io ga if enabled"""
    tce_opts = {'ccsd': True}
    if TCE_CCSDT:
        tce_opts['t'] = True  # Enable perturbative triples (T)
    if TCE_IO_GA:
        tce_opts['io'] = 'ga'  # Optimizes parallel memory handling
    return tce_opts

def run_tce_ccsdt_calculation(nproc):
    """Run NWChem TCE-CCSD(T)-ECP calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem TCE-CCSD(T)-ECP Calculation for Hg (with MPI)")
    print("="*70)
    print(f"Number of processors: {nproc}")
    print(f"Basis set:            {BASIS_SET}")
    print(f"ECP:                  {ECP_SET}")
    print(f"SCF threshold:        {SCF_THRESH}")
    print(f"SCF maxiter:          {SCF_MAXITER}")
    print(f"TCE CCSD(T):          {TCE_CCSDT}")
    print(f"TCE io ga:            {TCE_IO_GA}")
    print(f"Memory:               {MEMORY}")
    print("="*70)
    
    # Setup directories with short names
    base_dir = setup_directories()
    
    # Get TCE options
    tce_opts = get_tce_options()
    
    # === Neutral Hg (RHF + TCE-CCSD(T)) ===
    print("\n--- Neutral Hg (RHF + TCE-CCSD(T), Singlet) ---")
    neutral_dir = os.path.join(base_dir, 'neutral')
    
    hg_neutral = Atoms('Hg', positions=[(0, 0, 0)])
    hg_neutral.calc = NWChem(
        label=os.path.join(neutral_dir, 'hg'),
        theory='tce',
        task='energy',
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        tce=tce_opts,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_neutral = hg_neutral.get_potential_energy()
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # === Cation Hg+ (ROHF + TCE-CCSD(T)) ===
    print("\n--- Cation Hg+ (ROHF + TCE-CCSD(T), Open-Shell Doublet) ---")
    cation_rohf_dir = os.path.join(base_dir, 'cation_rohf')
    
    hg_cation_rohf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_rohf.calc = NWChem(
        label=os.path.join(cation_rohf_dir, 'hg'),
        theory='tce',
        task='energy',
        charge=1,
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'rohf': True, 'doublet': True, 'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        tce=tce_opts,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_cation_rohf = hg_cation_rohf.get_potential_energy()
    print(f"Energy: {e_cation_rohf:.6f} eV")
    print(f"📁 Output in: {cation_rohf_dir}")
    
    # === Cation Hg+ (UHF + TCE-CCSD(T)) ===
    print("\n--- Cation Hg+ (UHF + TCE-CCSD(T), Open-Shell Doublet) ---")
    cation_uhf_dir = os.path.join(base_dir, 'cation_uhf')
    
    hg_cation_uhf = Atoms('Hg', positions=[(0, 0, 0)])
    hg_cation_uhf.calc = NWChem(
        label=os.path.join(cation_uhf_dir, 'hg'),
        theory='tce',
        task='energy',
        charge=1,
        basis={'Hg': BASIS_SET},
        ecp={'Hg library': ECP_SET},
        scf={'uhf': True, 'doublet': True, 'thresh': SCF_THRESH, 'maxiter': SCF_MAXITER},
        tce=tce_opts,
        memory=MEMORY,
        command=f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
    )
    
    e_cation_uhf = hg_cation_uhf.get_potential_energy()
    print(f"Energy: {e_cation_uhf:.6f} eV")
    print(f"📁 Output in: {cation_uhf_dir}")
    
    # === Ionization Energies ===
    ie_rohf = e_cation_rohf - e_neutral
    ie_uhf = e_cation_uhf - e_neutral
    
    print("\n" + "="*70)
    print("RESULTS (TCE-CCSD(T) Level)")
    print("="*70)
    print(f"Neutral Energy (RHF+TCE-CCSD(T)):     {e_neutral:>12.6f} eV")
    print(f"Cation Energy (ROHF+TCE-CCSD(T)):     {e_cation_rohf:>12.6f} eV")
    print(f"Cation Energy (UHF+TCE-CCSD(T)):      {e_cation_uhf:>12.6f} eV")
    print("-"*70)
    print(f"Vertical IE (ROHF+TCE-CCSD(T)):       {ie_rohf:>12.6f} eV")
    print(f"Vertical IE (UHF+TCE-CCSD(T)):        {ie_uhf:>12.6f} eV")
    print(f"Experimental IE:                      {EXP_IE:>12.4f} eV")
    print("-"*70)
    print(f"Error (ROHF+TCE-CCSD(T)):             {ie_rohf-EXP_IE:>+12.4f} eV")
    print(f"Error (UHF+TCE-CCSD(T)):              {ie_uhf-EXP_IE:>+12.4f} eV")
    print("="*70)
    
    print("\n📁 Directory Structure:")
    print("   nwchem_tce_ccsdt/")
    print(f"   ├── {os.path.basename(neutral_dir)}/        (RHF+TCE-CCSD(T))")
    print(f"   ├── {os.path.basename(cation_rohf_dir)}/    (ROHF+TCE-CCSD(T))")
    print(f"   └── {os.path.basename(cation_uhf_dir)}/     (UHF+TCE-CCSD(T))")
    
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
        'method': 'TCE-CCSD(T)',
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
            print(f"\n   {os.path.basename(neutral_dir)}/ (RHF+TCE-CCSD(T)):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(neutral_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation ROHF files
        cation_rohf_dir = results['cation_rohf_dir']
        if os.path.exists(cation_rohf_dir):
            files = [f for f in os.listdir(cation_rohf_dir) if os.path.isfile(os.path.join(cation_rohf_dir, f))]
            print(f"\n   {os.path.basename(cation_rohf_dir)}/ (ROHF+TCE-CCSD(T)):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_rohf_dir, f))
                print(f"      - {f} ({size} bytes)")
        
        # Cation UHF files
        cation_uhf_dir = results['cation_uhf_dir']
        if os.path.exists(cation_uhf_dir):
            files = [f for f in os.listdir(cation_uhf_dir) if os.path.isfile(os.path.join(cation_uhf_dir, f))]
            print(f"\n   {os.path.basename(cation_uhf_dir)}/ (UHF+TCE-CCSD(T)):")
            for f in sorted(files):
                size = os.path.getsize(os.path.join(cation_uhf_dir, f))
                print(f"      - {f} ({size} bytes)")

def main():
    """Main execution function"""
    print("="*70)
    print("NWChem TCE-CCSD(T)-ECP CALCULATION FOR Hg ATOM")
    print("With ROHF+TCE-CCSD(T) and UHF+TCE-CCSD(T) for cation")
    print("Using ASE directives with short paths")
    print("="*70)
    
    # Read number of processors from file
    nproc = read_nproc_from_file(NPROC_FILE)
    
    # Run TCE-CCSD(T) calculation
    tce_ccsdt_results = run_tce_ccsdt_calculation(nproc)
    
    # Print file summary
    print_file_summary(tce_ccsdt_results)
    
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem TCE-CCSD(T)-ECP calculations completed")
    print(f"\n   ROHF+TCE-CCSD(T) IE = {tce_ccsdt_results['ie_rohf']:.6f} eV (error: {tce_ccsdt_results['ie_rohf']-EXP_IE:+.6f} eV)")
    print(f"   UHF+TCE-CCSD(T) IE  = {tce_ccsdt_results['ie_uhf']:.6f} eV (error: {tce_ccsdt_results['ie_uhf']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: TCE-CCSD(T) (Tensor Contraction Engine - Coupled Cluster with Singles, Doubles, and perturbative Triples)")
    print(f"   • TCE io ga: {tce_ccsdt_results['tce_io_ga']} (optimized parallel memory handling)")
    print(f"   • Neutral: RHF+TCE-CCSD(T) (closed-shell, singlet)")
    print(f"   • Cation ROHF: ROHF+TCE-CCSD(T) (restricted open-shell, doublet)")
    print(f"   • Cation UHF:  UHF+TCE-CCSD(T) (unrestricted open-shell, doublet)")
    print(f"   • Note: No permanent_dir or scratch_dir specified (uses NWChem defaults)")
    print(f"   • Basis set: {tce_ccsdt_results['basis']}")
    print(f"   • ECP: {tce_ccsdt_results['ecp']}")
    print(f"   • Processors: {tce_ccsdt_results['nproc']}")
    print(f"   • Memory: {MEMORY}")
    print("="*70)

if __name__ == "__main__":
    main()
