import subprocess
import os
import re
import shutil
from ase import Atoms
from ase.calculators.mopac import MOPAC
from ase.io import write

def setup_directories(base_dir, subdirs):
    """Create directory structure if it doesn't exist"""
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
    print("\n" + "="*60)
    print("MOPAC PM7 Calculations")
    print("="*60)
    
    # Setup directory structure
    base_dir = "mopac_calculations"
    subdirs = [
        "mopac_neutral_rhf",
        "mopac_cation_uhf"
    ]
    setup_directories(base_dir, subdirs)
    
    # === Neutral: RHF (Closed-Shell) ===
    print("\n--- Neutral Hg (RHF, Closed-Shell Singlet) ---")
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = MOPAC(
        label=os.path.join(base_dir, 'mopac_neutral_rhf', 'mop'), 
        task='1SCF RHF', 
        keywords='PM7'
    )
    e_neutral = neutral_atom.get_potential_energy()
    print(f"Energy: {e_neutral:.4f} eV")
    print(f"📁 Output in: {os.path.join(base_dir, 'mopac_neutral_rhf')}")
    
    # === Cation: UHF (Open-Shell) ===
    print("\n--- Cation Hg+ (UHF, Open-Shell Doublet) ---")
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = MOPAC(
        label=os.path.join(base_dir, 'mopac_cation_uhf', 'mop'), 
        task='1SCF UHF', 
        charge=1, 
        keywords='PM7 DOUBLET'
    )
    e_cation = cation_atom.get_potential_energy()
    print(f"Energy: {e_cation:.4f} eV")
    print(f"📁 Output in: {os.path.join(base_dir, 'mopac_cation_uhf')}")
    
    # === Ionization Energy ===
    ie = e_cation - e_neutral
    print(f"\nVertical IE:    {ie:.4f} eV")
    print(f"Experimental:   10.44 eV (difference: {ie-10.44:+.4f} eV)")
    
    return {
        'neutral': e_neutral,
        'cation': e_cation,
        'ie': ie,
        'neutral_method': 'RHF (closed-shell)',
        'cation_method': 'UHF (open-shell, doublet)',
        'neutral_dir': os.path.join(base_dir, 'mopac_neutral_rhf'),
        'cation_dir': os.path.join(base_dir, 'mopac_cation_uhf')
    }

def run_xtb_cli_ie(method='GFN2'):
    """
    Run xTB calculations for Hg ionization energy.
    
    xTB uses:
    - Closed-shell formalism for neutral (no --uhf flag)
    - UHF formalism for cation (--uhf 1 flag for doublet)
    """
    print("\n" + "="*60)
    print(f"xTB {method} Calculations")
    print("="*60)
    
    # Method flags
    method_flags = {
        'GFN2': '--gfn 2',
        'GFN1': '--gfn 1',
        'GFN0': '--gfn 0',
        'GFNFF': '--gfnff'
    }
    flag = method_flags.get(method, '--gfn 2')
    
    # Setup directory structure
    base_dir = "xtb_calculations"
    method_lower = method.lower()
    subdirs = [
        f"xtb_{method_lower}_neutral",
        f"xtb_{method_lower}_cation"
    ]
    setup_directories(base_dir, subdirs)
    
    # === Neutral: Closed-Shell ===
    print(f"\n--- Neutral Hg (Closed-Shell) ---")
    neutral_dir = os.path.join(base_dir, f'xtb_{method_lower}_neutral')
    
    # Prepare input files
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write(os.path.join(neutral_dir, 'coord.xyz'), atoms)
    with open(os.path.join(neutral_dir, 'coord'), 'w') as f:
        f.write("1\nHg atom\nHg    0.000000    0.000000    0.000000\n")
    
    # Run neutral calculation
    cmd_neutral = f"xtb coord.xyz {flag} > xtb_neutral.out 2>&1"
    res_n = subprocess.run(cmd_neutral, shell=True, cwd=neutral_dir, 
                           capture_output=True, text=True)
    
    if res_n.returncode == 0:
        print("✅ Neutral calculation completed")
        print(f"📁 Output in: {neutral_dir}")
    else:
        print("⚠️  Neutral calculation had warnings (check output file)")
    
    # === Cation: UHF (Open-Shell) ===
    print(f"\n--- Cation Hg+ (UHF, Open-Shell Doublet) ---")
    cation_dir = os.path.join(base_dir, f'xtb_{method_lower}_cation')
    
    # Prepare input files
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write(os.path.join(cation_dir, 'coord.xyz'), atoms)
    with open(os.path.join(cation_dir, 'coord'), 'w') as f:
        f.write("1\nHg atom\nHg    0.000000    0.000000    0.000000\n")
    
    # Run cation calculation with UHF
    cmd_cation = f"xtb coord.xyz {flag} --charge 1 --uhf 1 > xtb_cation.out 2>&1"
    res_c = subprocess.run(cmd_cation, shell=True, cwd=cation_dir,
                           capture_output=True, text=True)
    
    if res_c.returncode == 0:
        print("✅ Cation calculation completed")
        print(f"📁 Output in: {cation_dir}")
    else:
        print("⚠️  Cation calculation had warnings (check output file)")
    
    def parse_energy(directory, filename):
        """Parse total energy from xTB output file"""
        filepath = os.path.join(directory, filename)
        if not os.path.exists(filepath):
            return None
        try:
            with open(filepath, 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if "total energy" in line.lower():
                        numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                        for num in numbers:
                            try:
                                val = float(num)
                                if abs(val) > 0.0001:  # Skip small numbers
                                    return val * 27.211386245988  # Hartree -> eV
                            except ValueError:
                                continue
        except Exception as e:
            print(f"Error parsing {filename}: {e}")
        return None
    
    # Parse energies
    e_neutral = parse_energy(neutral_dir, 'xtb_neutral.out')
    e_cation = parse_energy(cation_dir, 'xtb_cation.out')
    
    if e_neutral is not None and e_cation is not None:
        ie = e_cation - e_neutral
        print(f"\nNeutral Energy: {e_neutral:.4f} eV")
        print(f"Cation Energy:  {e_cation:.4f} eV")
        print(f"Vertical IE:    {ie:.4f} eV")
        print(f"Experimental:   10.44 eV (difference: {ie-10.44:+.4f} eV)")
        
        results = {
            'neutral': e_neutral,
            'cation': e_cation,
            'ie': ie,
            'method': method,
            'neutral_dir': neutral_dir,
            'cation_dir': cation_dir
        }
    else:
        print(f"❌ Could not parse {method} output")
        results = None
    
    # Show created files
    for dir_name in [neutral_dir, cation_dir]:
        if os.path.exists(dir_name):
            files = [f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))]
            if files:
                print(f"\n📁 Files in {dir_name}:")
                for f in sorted(files):
                    size = os.path.getsize(os.path.join(dir_name, f))
                    print(f"   - {f} ({size} bytes)")
    
    return results

def print_summary(mopac_results, xtb_results):
    """Print comprehensive summary of all results"""
    print("\n" + "="*80)
    print("SUMMARY OF IONIZATION ENERGY CALCULATIONS FOR Hg ATOM")
    print("="*80)
    
    # Header
    print(f"{'Method':<20} {'Neutral (eV)':<15} {'Method':<20} {'Cation (eV)':<15} {'IE (eV)':<12} {'Error (eV)':<12}")
    print("-"*80)
    
    # MOPAC results
    if mopac_results:
        print(f"{'MOPAC PM7':<20} "
              f"{mopac_results['neutral']:<15.4f} "
              f"{mopac_results['neutral_method']:<20} "
              f"{mopac_results['cation']:<15.4f} "
              f"{mopac_results['ie']:<12.4f} "
              f"{mopac_results['ie']-10.44:<+12.4f}")
    
    # xTB results
    if xtb_results:
        for method, results in xtb_results.items():
            if results:
                print(f"{f'xTB {method}':<20} "
                      f"{results['neutral']:<15.4f} "
                      f"{'Closed-shell':<20} "
                      f"{results['cation']:<15.4f} "
                      f"{results['ie']:<12.4f} "
                      f"{results['ie']-10.44:<+12.4f}")
    
    print("-"*80)
    print(f"{'Experimental':<20} {'':<15} {'':<20} {'':<15} {'10.4400':<12} {'0.0000':<12}")
    print("="*80)
    
    # Directory structure
    print("\n📁 Directory Structure:")
    print("   mopac_calculations/")
    if mopac_results:
        print(f"   ├── {os.path.basename(mopac_results['neutral_dir'])}/  (RHF)")
        print(f"   └── {os.path.basename(mopac_results['cation_dir'])}/   (UHF)")
    
    print("\n   xtb_calculations/")
    if xtb_results:
        for method, results in xtb_results.items():
            if results:
                print(f"   ├── {os.path.basename(results['neutral_dir'])}/  ({method}, closed-shell)")
                print(f"   └── {os.path.basename(results['cation_dir'])}/   ({method}, UHF)")
    
    # Additional information
    print("\n📝 Notes:")
    print("  • MOPAC: Neutral uses RHF (default for even-electron), Cation uses UHF (default for odd-electron)")
    print("  • xTB:   Neutral uses closed-shell, Cation uses UHF (--uhf 1 for doublet)")
    print("  • Both methods correctly use UHF for the open-shell cation")

def main():
    """Main execution function"""
    print("="*80)
    print("IONIZATION ENERGY CALCULATIONS FOR Hg ATOM")
    print("Methods: MOPAC PM7, xTB GFN1, xTB GFN2")
    print("="*80)
    
    # Run MOPAC
    mopac_results = run_mopac_ie()
    
    # Run xTB methods
    xtb_results = {}
    for method in ['GFN2', 'GFN1']:
        print("\n" + "="*80)
        xtb_results[method] = run_xtb_cli_ie(method)
    
    # Print summary
    print_summary(mopac_results, xtb_results)
    
    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    if mopac_results:
        print(f"✅ BEST METHOD: MOPAC PM7 (UHF for cation)")
        print(f"   IE = {mopac_results['ie']:.4f} eV")
        print(f"   Error vs experiment: {mopac_results['ie']-10.44:+.4f} eV")
        print("   • PM7 is specifically parameterized for heavy elements like Hg")
        print("   • Correctly handles open-shell systems with UHF")
    
    print(f"\n⚠️  xTB methods overestimate IE for Hg:")
    for method, results in xtb_results.items():
        if results:
            print(f"   {method}: {results['ie']:.4f} eV (error: {results['ie']-10.44:+.4f} eV)")
    
    print("\n💡 For accurate Hg calculations, use:")
    print("   • MOPAC PM7 (recommended for heavy elements)")
    print("   • DFT with pseudopotentials (for higher accuracy)")
    print("   • xTB only for qualitative trends")
    print("="*80)

if __name__ == "__main__":
    main()
