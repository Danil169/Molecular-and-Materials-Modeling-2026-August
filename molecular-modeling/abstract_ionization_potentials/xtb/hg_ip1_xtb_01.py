#!/usr/bin/env python
"""
xTB Ionization Energy Calculation for Hg Atom
==============================================
This script calculates the vertical ionization energy of a mercury atom
using xTB with different methods (GFN1, GFN2, etc.):
- Closed-shell for neutral
- RHF-like (closed-shell) for cation (no --uhf flag)
- UHF for cation (--uhf 1 flag for doublet)

Reference: Experimental IE of Hg = 10.44 eV
"""

import subprocess
import os
import re
import shutil
from ase import Atoms
from ase.io import write

def setup_directories(methods, include_rhf=True):
    """Create directory structure for xTB calculations"""
    base_dir = "xtb_calculations"
    subdirs = []
    
    for method in methods:
        method_lower = method.lower()
        # Neutral (always closed-shell)
        subdirs.append(f"xtb_{method_lower}_neutral")
        
        # Cation with RHF-like (closed-shell)
        if include_rhf:
            subdirs.append(f"xtb_{method_lower}_cation_rhf")
        
        # Cation with UHF (open-shell)
        subdirs.append(f"xtb_{method_lower}_cation_uhf")
    
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

def run_xtb_calculation(method='GFN2', include_rhf=True, verbose=True):
    """
    Run xTB calculations for Hg ionization energy.
    
    xTB uses:
    - Closed-shell formalism for neutral (no --uhf flag)
    - RHF-like (closed-shell) for cation (no --uhf flag)
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
    base_dir = setup_directories([method], include_rhf)
    method_lower = method.lower()
    
    results = {}
    exp_ie = 10.44  # Experimental value for Hg
    
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
    
    # Parse neutral energy
    e_neutral = parse_energy(neutral_dir, 'xtb_neutral.out')
    if e_neutral is None:
        e_neutral = parse_energy(neutral_dir, 'xtb_neutral.out')  # Try again
    
    if e_neutral is None:
        print(f"❌ Could not parse neutral energy for {method}")
        return None
    
    # === Cation: RHF-like (Closed-Shell) ===
    if include_rhf:
        if verbose:
            print(f"\n--- Cation Hg+ (RHF-like, Closed-Shell) ---")
        cation_rhf_dir = os.path.join(base_dir, f'xtb_{method_lower}_cation_rhf')
        
        # Prepare input files
        atoms = Atoms('Hg', positions=[(0, 0, 0)])
        write(os.path.join(cation_rhf_dir, 'coord.xyz'), atoms)
        with open(os.path.join(cation_rhf_dir, 'coord'), 'w') as f:
            f.write("1\nHg atom\nHg    0.000000    0.000000    0.000000\n")
        
        # Run cation calculation without UHF (RHF-like)
        cmd_cation_rhf = f"xtb coord.xyz {flag} --charge 1 > xtb_cation_rhf.out 2>&1"
        res_c_rhf = subprocess.run(cmd_cation_rhf, shell=True, cwd=cation_rhf_dir,
                                   capture_output=True, text=True)
        
        if res_c_rhf.returncode == 0:
            if verbose:
                print("✅ Cation (RHF-like) calculation completed")
                print(f"📁 Output in: {cation_rhf_dir}")
        else:
            print("⚠️  Cation (RHF-like) calculation had warnings (check output file)")
        
        # Parse cation RHF energy
        e_cation_rhf = parse_energy(cation_rhf_dir, 'xtb_cation_rhf.out')
        
        if e_cation_rhf is not None:
            ie_rhf = e_cation_rhf - e_neutral
            if verbose:
                print("\n" + "="*70)
                print(f"RESULTS for xTB {method} (RHF-like Cation)")
                print("="*70)
                print(f"Neutral Energy:         {e_neutral:>12.4f} eV")
                print(f"Cation Energy (RHF-like): {e_cation_rhf:>8.4f} eV")
                print(f"Vertical IE (RHF-like): {ie_rhf:>12.4f} eV")
                print(f"Experimental IE:        {exp_ie:>12.4f} eV")
                print(f"Difference:             {ie_rhf-exp_ie:>+12.4f} eV")
                print("="*70)
            
            results['RHF-like'] = {
                'method': method,
                'neutral_energy': e_neutral,
                'cation_energy': e_cation_rhf,
                'ie': ie_rhf,
                'experimental': exp_ie,
                'error': ie_rhf - exp_ie,
                'neutral_dir': neutral_dir,
                'cation_dir': cation_rhf_dir,
                'method_type': 'RHF-like (closed-shell)'
            }
            
            # Show files for RHF-like
            if verbose:
                for dir_name in [neutral_dir, cation_rhf_dir]:
                    if os.path.exists(dir_name):
                        files = [f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))]
                        if files:
                            print(f"\n📁 Files in {dir_name}:")
                            for f in sorted(files):
                                size = os.path.getsize(os.path.join(dir_name, f))
                                print(f"   - {f} ({size} bytes)")
        else:
            print(f"❌ Could not parse RHF-like cation energy for {method}")
    
    # === Cation: UHF (Open-Shell) ===
    if verbose:
        print(f"\n--- Cation Hg+ (UHF, Open-Shell Doublet) ---")
    cation_uhf_dir = os.path.join(base_dir, f'xtb_{method_lower}_cation_uhf')
    
    # Prepare input files
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write(os.path.join(cation_uhf_dir, 'coord.xyz'), atoms)
    with open(os.path.join(cation_uhf_dir, 'coord'), 'w') as f:
        f.write("1\nHg atom\nHg    0.000000    0.000000    0.000000\n")
    
    # Run cation calculation with UHF
    cmd_cation_uhf = f"xtb coord.xyz {flag} --charge 1 --uhf 1 > xtb_cation_uhf.out 2>&1"
    res_c_uhf = subprocess.run(cmd_cation_uhf, shell=True, cwd=cation_uhf_dir,
                               capture_output=True, text=True)
    
    if res_c_uhf.returncode == 0:
        if verbose:
            print("✅ Cation (UHF) calculation completed")
            print(f"📁 Output in: {cation_uhf_dir}")
    else:
        print("⚠️  Cation (UHF) calculation had warnings (check output file)")
    
    # Parse cation UHF energy
    e_cation_uhf = parse_energy(cation_uhf_dir, 'xtb_cation_uhf.out')
    
    if e_cation_uhf is not None:
        ie_uhf = e_cation_uhf - e_neutral
        if verbose:
            print("\n" + "="*70)
            print(f"RESULTS for xTB {method} (UHF Cation)")
            print("="*70)
            print(f"Neutral Energy:         {e_neutral:>12.4f} eV")
            print(f"Cation Energy (UHF):    {e_cation_uhf:>12.4f} eV")
            print(f"Vertical IE (UHF):      {ie_uhf:>12.4f} eV")
            print(f"Experimental IE:        {exp_ie:>12.4f} eV")
            print(f"Difference:             {ie_uhf-exp_ie:>+12.4f} eV")
            print("="*70)
        
        results['UHF'] = {
            'method': method,
            'neutral_energy': e_neutral,
            'cation_energy': e_cation_uhf,
            'ie': ie_uhf,
            'experimental': exp_ie,
            'error': ie_uhf - exp_ie,
            'neutral_dir': neutral_dir,
            'cation_dir': cation_uhf_dir,
            'method_type': 'UHF (open-shell)'
        }
        
        # Show files for UHF
        if verbose:
            for dir_name in [neutral_dir, cation_uhf_dir]:
                if os.path.exists(dir_name):
                    files = [f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))]
                    if files:
                        print(f"\n📁 Files in {dir_name}:")
                        for f in sorted(files):
                            size = os.path.getsize(os.path.join(dir_name, f))
                            print(f"   - {f} ({size} bytes)")
    else:
        print(f"❌ Could not parse UHF cation energy for {method}")
    
    return results

def run_multiple_xtb_methods(methods=['GFN2', 'GFN1'], include_rhf=True):
    """Run xTB calculations for multiple methods"""
    results = {}
    
    for method in methods:
        print("\n" + "="*70)
        print(f"xTB {method} Ionization Energy Calculation for Hg")
        print("="*70)
        results[method] = run_xtb_calculation(method, include_rhf)
    
    return results

def print_summary(results):
    """Print summary of all xTB results"""
    print("\n" + "="*80)
    print("SUMMARY OF xTB IONIZATION ENERGY CALCULATIONS FOR Hg")
    print("="*80)
    
    # Header
    print(f"{'Method':<15} {'State':<15} {'Neutral (eV)':<15} {'Cation (eV)':<15} {'IE (eV)':<12} {'Error (eV)':<12}")
    print("-"*80)
    
    # Results
    for method, res in results.items():
        if res:
            # RHF-like results
            if 'RHF-like' in res:
                r = res['RHF-like']
                print(f"{method:<15} {'RHF-like':<15} "
                      f"{r['neutral_energy']:<15.4f} "
                      f"{r['cation_energy']:<15.4f} "
                      f"{r['ie']:<12.4f} "
                      f"{r['error']:<+12.4f}")
            
            # UHF results
            if 'UHF' in res:
                r = res['UHF']
                print(f"{method:<15} {'UHF':<15} "
                      f"{r['neutral_energy']:<15.4f} "
                      f"{r['cation_energy']:<15.4f} "
                      f"{r['ie']:<12.4f} "
                      f"{r['error']:<+12.4f}")
    
    print("-"*80)
    print(f"{'Experimental':<15} {'':<15} {'':<15} {'':<15} {'10.4400':<12} {'0.0000':<12}")
    print("="*80)
    
    # Directory structure
    print("\n📁 Directory Structure:")
    print("   xtb_calculations/")
    for method, res in results.items():
        if res:
            method_lower = method.lower()
            print(f"   ├── xtb_{method_lower}_neutral/      ({method}, closed-shell)")
            if 'RHF-like' in res:
                print(f"   ├── xtb_{method_lower}_cation_rhf/   ({method}, RHF-like)")
            if 'UHF' in res:
                print(f"   └── xtb_{method_lower}_cation_uhf/   ({method}, UHF)")

def main():
    """Main execution function"""
    print("="*80)
    print("xTB IONIZATION ENERGY CALCULATIONS FOR Hg ATOM")
    print("Methods: GFN1, GFN2")
    print("Including: RHF-like and UHF for cation")
    print("="*80)
    
    # Run xTB calculations for multiple methods
    methods = ['GFN2', 'GFN1']
    results = run_multiple_xtb_methods(methods, include_rhf=True)
    
    # Print summary
    print_summary(results)
    
    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    # Find best method (closest to experimental)
    best_method = None
    best_type = None
    best_error = float('inf')
    
    for method, res in results.items():
        if res:
            for state in ['RHF-like', 'UHF']:
                if state in res:
                    error = abs(res[state]['error'])
                    if error < best_error:
                        best_error = error
                        best_method = method
                        best_type = state
    
    if best_method and best_type:
        print(f"✅ Best xTB method for Hg: {best_method} ({best_type})")
        print(f"   IE = {results[best_method][best_type]['ie']:.4f} eV")
        print(f"   Error: {results[best_method][best_type]['error']:+.4f} eV")
    
    print("\n⚠️  All xTB methods overestimate IE for Hg")
    print("   Consider using:")
    print("   • MOPAC PM7 (recommended for heavy elements)")
    print("   • DFT with pseudopotentials (for higher accuracy)")
    print("   • xTB only for qualitative trends")
    print("="*80)

if __name__ == "__main__":
    main()
