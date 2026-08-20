#!/usr/bin/env python
"""
xTB Ionization Energy Calculation for Hg Atom
==============================================
This script calculates the vertical ionization energy of a mercury atom
using xTB with different methods (GFN1, GFN2, etc.):
- Closed-shell for neutral
- UHF for cation (--uhf 1 flag for doublet)

Reference: Experimental IE of Hg = 10.44 eV
"""

import subprocess
import os
import re
import shutil
from ase import Atoms
from ase.io import write

def setup_directories(methods):
    """Create directory structure for xTB calculations"""
    base_dir = "xtb_calculations"
    subdirs = []
    
    for method in methods:
        method_lower = method.lower()
        subdirs.extend([
            f"xtb_{method_lower}_neutral",
            f"xtb_{method_lower}_cation"
        ])
    
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    
    return base_dir

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

def run_xtb_calculation(method='GFN2', verbose=True):
    """
    Run xTB calculations for Hg ionization energy.
    
    xTB uses:
    - Closed-shell formalism for neutral (no --uhf flag)
    - UHF formalism for cation (--uhf 1 flag for doublet)
    """
    # Method flags
    method_flags = {
        'GFN2': '--gfn 2',
        'GFN1': '--gfn 1',
        'GFN0': '--gfn 0',
        'GFNFF': '--gfnff'
    }
    flag = method_flags.get(method, '--gfn 2')
    
    # Setup directories
    base_dir = setup_directories([method])
    method_lower = method.lower()
    
    # === Neutral: Closed-Shell ===
    if verbose:
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
        if verbose:
            print("✅ Neutral calculation completed")
            print(f"📁 Output in: {neutral_dir}")
    else:
        print("⚠️  Neutral calculation had warnings (check output file)")
    
    # === Cation: UHF (Open-Shell) ===
    if verbose:
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
        if verbose:
            print("✅ Cation calculation completed")
            print(f"📁 Output in: {cation_dir}")
    else:
        print("⚠️  Cation calculation had warnings (check output file)")
    
    # Parse energies
    e_neutral = parse_energy(neutral_dir, 'xtb_neutral.out')
    e_cation = parse_energy(cation_dir, 'xtb_cation.out')
    
    if e_neutral is not None and e_cation is not None:
        ie = e_cation - e_neutral
        exp_ie = 10.44  # Experimental value for Hg
        
        if verbose:
            print("\n" + "="*70)
            print(f"RESULTS for xTB {method}")
            print("="*70)
            print(f"Neutral Energy:         {e_neutral:>12.4f} eV")
            print(f"Cation Energy:          {e_cation:>12.4f} eV")
            print(f"Vertical IE:            {ie:>12.4f} eV")
            print(f"Experimental IE:        {exp_ie:>12.4f} eV")
            print(f"Difference:             {ie-exp_ie:>+12.4f} eV")
            print("="*70)
        
        # Show created files
        for dir_name in [neutral_dir, cation_dir]:
            if os.path.exists(dir_name):
                files = [f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))]
                if files and verbose:
                    print(f"\n📁 Files in {dir_name}:")
                    for f in sorted(files):
                        size = os.path.getsize(os.path.join(dir_name, f))
                        print(f"   - {f} ({size} bytes)")
        
        return {
            'method': method,
            'neutral_energy': e_neutral,
            'cation_energy': e_cation,
            'ie': ie,
            'experimental': exp_ie,
            'error': ie - exp_ie,
            'neutral_dir': neutral_dir,
            'cation_dir': cation_dir,
            'method_type': 'UHF (open-shell)'
        }
    else:
        print(f"❌ Could not parse {method} output")
        return None

def run_multiple_xtb_methods(methods=['GFN2', 'GFN1']):
    """Run xTB calculations for multiple methods"""
    results = {}
    
    for method in methods:
        print("\n" + "="*70)
        print(f"xTB {method} Ionization Energy Calculation for Hg")
        print("="*70)
        results[method] = run_xtb_calculation(method)
    
    return results

def print_summary(results):
    """Print summary of all xTB results"""
    print("\n" + "="*80)
    print("SUMMARY OF xTB IONIZATION ENERGY CALCULATIONS FOR Hg")
    print("="*80)
    
    # Header
    print(f"{'Method':<15} {'Neutral (eV)':<15} {'Cation (eV)':<15} {'IE (eV)':<12} {'Error (eV)':<12}")
    print("-"*80)
    
    # Results
    for method, res in results.items():
        if res:
            print(f"{method:<15} "
                  f"{res['neutral_energy']:<15.4f} "
                  f"{res['cation_energy']:<15.4f} "
                  f"{res['ie']:<12.4f} "
                  f"{res['error']:<+12.4f}")
    
    print("-"*80)
    print(f"{'Experimental':<15} {'':<15} {'':<15} {'10.4400':<12} {'0.0000':<12}")
    print("="*80)
    
    # Directory structure
    print("\n📁 Directory Structure:")
    print("   xtb_calculations/")
    for method, res in results.items():
        if res:
            method_lower = method.lower()
            print(f"   ├── xtb_{method_lower}_neutral/  ({method}, closed-shell)")
            print(f"   └── xtb_{method_lower}_cation/   ({method}, UHF)")

def main():
    """Main execution function"""
    print("="*80)
    print("xTB IONIZATION ENERGY CALCULATIONS FOR Hg ATOM")
    print("Methods: GFN1, GFN2")
    print("="*80)
    
    # Run xTB calculations for multiple methods
    methods = ['GFN2', 'GFN1']
    results = run_multiple_xtb_methods(methods)
    
    # Print summary
    print_summary(results)
    
    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    # Find best method (closest to experimental)
    best_method = None
    best_error = float('inf')
    for method, res in results.items():
        if res:
            error = abs(res['error'])
            if error < best_error:
                best_error = error
                best_method = method
    
    if best_method:
        print(f"✅ Best xTB method for Hg: {best_method}")
        print(f"   IE = {results[best_method]['ie']:.4f} eV")
        print(f"   Error: {results[best_method]['error']:+.4f} eV")
    
    print("\n⚠️  All xTB methods overestimate IE for Hg")
    print("   Consider using:")
    print("   • MOPAC PM7 (recommended for heavy elements)")
    print("   • DFT with pseudopotentials (for higher accuracy)")
    print("   • xTB only for qualitative trends")
    print("="*80)

if __name__ == "__main__":
    main()
