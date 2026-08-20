#!/usr/bin/env python
"""
NWChem TCE-CCSD(T)-ECP Calculation for Hg using ASE parameters with MPI
- Neutral Hg: RHF + TCE-CCSD(T) (closed-shell, singlet)
- Cation Hg+: ROHF + TCE-CCSD(T) and UHF + TCE-CCSD(T) (open-shell, doublet)
- Uses larger basis set def2-qzvp for improved accuracy
- Uses mpirun with number of processors from an input file
- Includes comprehensive timing information
"""

from ase import Atoms
from ase.calculators.nwchem import NWChem
import os
import shutil
import time
import datetime
from collections import defaultdict

# ============================================================
# USER-DEFINED PARAMETERS - Modify these as needed
# ============================================================

# Basis set and ECP - UPGRADED to Quadruple-Zeta
BASIS_SET = "def2-qzvp"      # Quadruple-zeta valence with polarization
ECP_SET = "def2-qzvp"        # Quadruple-zeta ECP for Hg

# SCF settings (tighter for larger basis)
SCF_THRESH = 1.0e-9          # Tighter threshold for quadruple-zeta
SCF_MAXITER = 300            # More iterations for larger basis

# TCE-CCSD(T) settings
TCE_IO_GA = True             # Optimizes parallel memory handling
TCE_THRESH = 1.0e-9          # Tighter TCE threshold for accuracy

# Memory (quadruple-zeta needs much more memory)
MEMORY = "32000 mb"          # Increased memory for def2-qzvp

# Experimental reference
EXP_IE = 10.437  # eV

# File containing number of processors
NPROC_FILE = "nproc.txt"

# ============================================================
# END OF USER PARAMETERS
# ============================================================

class Timer:
    """Simple timer class for tracking calculation times"""
    def __init__(self):
        self.start_times = {}
        self.end_times = {}
        self.elapsed = defaultdict(float)
        self._running = {}
    
    def start(self, name):
        """Start timing a section"""
        self.start_times[name] = time.time()
        self._running[name] = True
    
    def stop(self, name):
        """Stop timing a section"""
        if name in self._running and self._running[name]:
            self.end_times[name] = time.time()
            elapsed = self.end_times[name] - self.start_times[name]
            self.elapsed[name] += elapsed
            self._running[name] = False
        return self.elapsed[name]
    
    def get(self, name):
        """Get elapsed time for a section"""
        return self.elapsed.get(name, 0.0)
    
    def print_summary(self):
        """Print timing summary"""
        print("\n" + "="*70)
        print("TIMING SUMMARY")
        print("="*70)
        
        total_time = self.get('total')
        
        # Calculate percentages for each section
        sections = [
            ('neutral_scf', 'Neutral SCF'),
            ('neutral_ccsd', 'Neutral CCSD'),
            ('neutral_total', 'Neutral Total'),
            ('cation_rohf_scf', 'Cation ROHF SCF'),
            ('cation_rohf_ccsd', 'Cation ROHF CCSD'),
            ('cation_rohf_total', 'Cation ROHF Total'),
            ('cation_uhf_scf', 'Cation UHF SCF'),
            ('cation_uhf_ccsd', 'Cation UHF CCSD'),
            ('cation_uhf_total', 'Cation UHF Total'),
        ]
        
        print(f"{'Section':<25} {'Time (s)':<15} {'Time (min)':<15} {'Percentage':<12}")
        print("-"*70)
        
        for key, label in sections:
            t = self.get(key)
            if t > 0:
                pct = (t / total_time * 100) if total_time > 0 else 0
                print(f"{label:<25} {t:<15.2f} {t/60:<15.2f} {pct:<11.1f}%")
        
        print("-"*70)
        print(f"{'TOTAL':<25} {total_time:<15.2f} {total_time/60:<15.2f} {'100.0':<11}%")
        print("="*70)
        
        # Print wall-clock time
        print(f"\n⏱️  Wall-clock time: {total_time/60:.2f} minutes ({total_time:.2f} seconds)")
        print(f"📅 Completed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

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
    base_dir = "nwchem_tce_ccsdt_qzvp"  # New directory for def2-qzvp
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
    """
    Return TCE options dictionary with CCSD(T) syntax.
    The correct syntax is 'ccsd(t)' for CCSD with perturbative triples.
    """
    tce_opts = {'ccsd(t)': True}  # Correct syntax for CCSD(T)
    if TCE_IO_GA:
        tce_opts['io'] = 'ga'     # Global Arrays for distributed memory
    if TCE_THRESH:
        tce_opts['thresh'] = TCE_THRESH
    return tce_opts

def parse_nwchem_timing(output_file):
    """
    Parse timing information from NWChem output file
    """
    times = {}
    if not os.path.exists(output_file):
        return times
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            
            # Look for timing summary
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if 'wall-time' in line.lower():
                    # Parse wall time
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        times['wall_time'] = float(numbers[0])
                if 'cpu-time' in line.lower():
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        times['cpu_time'] = float(numbers[0])
    except Exception as e:
        print(f"Warning: Could not parse timing from {output_file}: {e}")
    
    return times

def run_tce_ccsdt_calculation(nproc, timer):
    """Run NWChem TCE-CCSD(T)-ECP calculations for Hg ionization energy"""
    print("="*70)
    print("NWChem TCE-CCSD(T)-ECP Calculation for Hg (with MPI)")
    print("="*70)
    print(f"Number of processors: {nproc}")
    print(f"Basis set:            {BASIS_SET} (Quadruple-Zeta)")
    print(f"ECP:                  {ECP_SET}")
    print(f"SCF threshold:        {SCF_THRESH}")
    print(f"SCF maxiter:          {SCF_MAXITER}")
    print(f"TCE method:           CCSD(T)")
    print(f"TCE io ga:            {TCE_IO_GA}")
    print(f"TCE thresh:           {TCE_THRESH}")
    print(f"Memory:               {MEMORY}")
    print(f"Start time:           {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)
    
    # Start total timer
    timer.start('total')
    
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
    
    # Start timing for neutral
    timer.start('neutral_scf')
    timer.start('neutral_total')
    
    e_neutral = hg_neutral.get_potential_energy()
    
    timer.stop('neutral_scf')
    timer.stop('neutral_total')
    
    print(f"Energy: {e_neutral:.6f} eV")
    print(f"📁 Output in: {neutral_dir}")
    
    # Parse neutral timing
    neutral_timing = parse_nwchem_timing(os.path.join(neutral_dir, 'hg.nwo'))
    if neutral_timing:
        print(f"   Wall time: {neutral_timing.get('wall_time', 0):.2f} s")
    
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
    
    # Start timing for cation ROHF
    timer.start('cation_rohf_scf')
    timer.start('cation_rohf_total')
    
    e_cation_rohf = hg_cation_rohf.get_potential_energy()
    
    timer.stop('cation_rohf_scf')
    timer.stop('cation_rohf_total')
    
    print(f"Energy: {e_cation_rohf:.6f} eV")
    print(f"📁 Output in: {cation_rohf_dir}")
    
    # Parse cation ROHF timing
    cation_rohf_timing = parse_nwchem_timing(os.path.join(cation_rohf_dir, 'hg.nwo'))
    if cation_rohf_timing:
        print(f"   Wall time: {cation_rohf_timing.get('wall_time', 0):.2f} s")
    
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
    
    # Start timing for cation UHF
    timer.start('cation_uhf_scf')
    timer.start('cation_uhf_total')
    
    e_cation_uhf = hg_cation_uhf.get_potential_energy()
    
    timer.stop('cation_uhf_scf')
    timer.stop('cation_uhf_total')
    
    print(f"Energy: {e_cation_uhf:.6f} eV")
    print(f"📁 Output in: {cation_uhf_dir}")
    
    # Parse cation UHF timing
    cation_uhf_timing = parse_nwchem_timing(os.path.join(cation_uhf_dir, 'hg.nwo'))
    if cation_uhf_timing:
        print(f"   Wall time: {cation_uhf_timing.get('wall_time', 0):.2f} s")
    
    # === Ionization Energies ===
    ie_rohf = e_cation_rohf - e_neutral
    ie_uhf = e_cation_uhf - e_neutral
    
    # Stop total timer
    timer.stop('total')
    
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
    print("   nwchem_tce_ccsdt_qzvp/")
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
        'tce_io_ga': TCE_IO_GA,
        'tce_thresh': TCE_THRESH
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

def print_basis_set_comparison(prev_results, current_results):
    """Compare results from different basis sets"""
    if prev_results and current_results:
        print("\n" + "="*70)
        print("BASIS SET COMPARISON")
        print("="*70)
        print(f"{'Basis Set':<15} {'IE (eV)':<15} {'Error (eV)':<15} {'Improvement':<15}")
        print("-"*70)
        
        # Previous basis
        prev_ie = prev_results['ie_uhf']
        prev_error = prev_ie - EXP_IE
        print(f"{prev_results['basis']:<15} {prev_ie:<15.6f} {prev_error:<+15.4f} {'-':<15}")
        
        # Current basis
        curr_ie = current_results['ie_uhf']
        curr_error = curr_ie - EXP_IE
        improvement = prev_error - curr_error
        print(f"{current_results['basis']:<15} {curr_ie:<15.6f} {curr_error:<+15.4f} {improvement:<+15.4f}")
        
        print("-"*70)
        print(f"{'ROHF':<15} {current_results['ie_rohf']:<15.6f} {current_results['ie_rohf']-EXP_IE:<+15.4f} {'-':<15}")
        print(f"{'Experimental':<15} {'10.4370':<15} {'0.0000':<15} {'-':<15}")
        print("="*70)

def main():
    """Main execution function"""
    # Initialize timer
    timer = Timer()
    timer.start('total_script')
    
    print("="*70)
    print("NWChem TCE-CCSD(T)-ECP CALCULATION FOR Hg ATOM")
    print(f"Using {BASIS_SET} basis set (Quadruple-Zeta)")
    print("With ROHF+TCE-CCSD(T) and UHF+TCE-CCSD(T) for cation")
    print(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)
    
    # Read number of processors from file
    nproc = read_nproc_from_file(NPROC_FILE)
    
    # Run TCE-CCSD(T) calculation
    tce_ccsdt_results = run_tce_ccsdt_calculation(nproc, timer)
    
    # Print file summary
    print_file_summary(tce_ccsdt_results)
    
    # Print timing summary
    timer.stop('total_script')
    timer.print_summary()
    
    # Print recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print(f"✅ NWChem TCE-CCSD(T)-ECP calculations completed with {BASIS_SET}")
    print(f"\n   ROHF+TCE-CCSD(T) IE = {tce_ccsdt_results['ie_rohf']:.6f} eV (error: {tce_ccsdt_results['ie_rohf']-EXP_IE:+.6f} eV)")
    print(f"   UHF+TCE-CCSD(T) IE  = {tce_ccsdt_results['ie_uhf']:.6f} eV (error: {tce_ccsdt_results['ie_uhf']-EXP_IE:+.6f} eV)")
    print(f"\n   • Method: TCE-CCSD(T) (Tensor Contraction Engine - CCSD with perturbative triples)")
    print(f"   • Basis set: {BASIS_SET} (Quadruple-Zeta Valence with Polarization)")
    print(f"   • ECP: {ECP_SET}")
    print(f"   • TCE syntax: ccsd(t) with io ga and thresh {tce_ccsdt_results['tce_thresh']}")
    print(f"   • TCE io ga: {tce_ccsdt_results['tce_io_ga']} (optimized parallel memory handling)")
    print(f"   • Processors: {tce_ccsdt_results['nproc']}")
    print(f"   • Memory: {MEMORY}")
    print(f"   • Total time: {timer.get('total_script')/60:.2f} minutes")
    
    print("\n📈 Expected Convergence:")
    print(f"   def2-svp:  10.018 eV (error: -0.419 eV)")
    print(f"   def2-tzvp: 10.139 eV (error: -0.298 eV)")
    print(f"   def2-qzvp: ~10.22-10.26 eV (error: ~-0.18 to -0.22 eV)")
    print(f"   CBS limit: ~10.34-10.38 eV (error: ~-0.06 to -0.10 eV)")
    print(f"   Experimental: 10.437 eV")
    print("="*70)

if __name__ == "__main__":
    main()
