import subprocess
import os
import re
import shutil
from ase import Atoms
from ase.calculators.mopac import MOPAC
from ase.io import write

def run_mopac_ie():
    """Run MOPAC calculations for both RHF (closed-shell) and UHF (open-shell)"""
    print("\n" + "="*60)
    print("MOPAC PM7 Calculations")
    print("="*60)
    
    results = {}
    
    # === RHF Calculation (Closed-Shell) ===
    print("\n--- MOPAC PM7 (RHF - Closed Shell) ---")
    
    # Neutral (RHF, singlet)
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = MOPAC(label='mopac_neutral_rhf/mop', 
                              task='1SCF RHF', 
                              keywords='PM7')
    e_neutral_rhf = neutral_atom.get_potential_energy()
    
    # Cation (RHF forced)
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = MOPAC(label='mopac_cation_rhf/mop', 
                             task='1SCF RHF', 
                             charge=1, 
                             keywords='PM7 CLOSED SHELL')
    e_cation_rhf = cation_atom.get_potential_energy()
    
    ie_rhf = e_cation_rhf - e_neutral_rhf
    print(f"Neutral Energy: {e_neutral_rhf:.4f} eV")
    print(f"Cation Energy:  {e_cation_rhf:.4f} eV")
    print(f"Vertical IE:    {ie_rhf:.4f} eV")
    print(f"Experimental:   10.44 eV (difference: {ie_rhf-10.44:+.4f} eV)")
    
    results['RHF'] = {
        'neutral': e_neutral_rhf,
        'cation': e_cation_rhf,
        'ie': ie_rhf
    }
    
    # === UHF Calculation (Open-Shell) ===
    print("\n--- MOPAC PM7 (UHF - Open Shell) ---")
    
    # Neutral (UHF, singlet - can also do RHF for singlet)
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = MOPAC(label='mopac_neutral_uhf/mop', 
                              task='1SCF UHF', 
                              keywords='PM7')
    e_neutral_uhf = neutral_atom.get_potential_energy()
    
    # Cation (UHF, doublet)
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = MOPAC(label='mopac_cation_uhf/mop', 
                             task='1SCF UHF', 
                             charge=1, 
                             keywords='PM7 DOUBLET')
    e_cation_uhf = cation_atom.get_potential_energy()
    
    ie_uhf = e_cation_uhf - e_neutral_uhf
    print(f"Neutral Energy: {e_neutral_uhf:.4f} eV")
    print(f"Cation Energy:  {e_cation_uhf:.4f} eV")
    print(f"Vertical IE:    {ie_uhf:.4f} eV")
    print(f"Experimental:   10.44 eV (difference: {ie_uhf-10.44:+.4f} eV)")
    
    results['UHF'] = {
        'neutral': e_neutral_uhf,
        'cation': e_cation_uhf,
        'ie': ie_uhf
    }
    
    return results

def run_xtb_cli_ie(method='GFN2'):
    """Run xTB calculations for both RHF (closed-shell) and UHF (open-shell)
       with full output files saved"""
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
    
    results = {}
    
    # === RHF-like Calculation (Closed-Shell) ===
    print(f"\n--- xTB {method} (RHF-like - Closed Shell) ---")
    
    # Clean directory for RHF
    xtb_dir_rhf = f'xtb_run_{method.lower()}_rhf'
    if os.path.exists(xtb_dir_rhf):
        shutil.rmtree(xtb_dir_rhf)
    os.makedirs(xtb_dir_rhf, exist_ok=True)
    
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write(f'{xtb_dir_rhf}/coord.xyz', atoms)
    
    # Also create a .xyz file with proper format
    with open(f'{xtb_dir_rhf}/coord', 'w') as f:
        f.write("1\n")
        f.write("Hg atom\n")
        f.write("Hg    0.000000    0.000000    0.000000\n")
    
    # Neutral (closed-shell, singlet) - Save full output
    cmd_neutral = f"xtb coord.xyz {flag} --opt  > xtb_neutral.out 2>&1"
    res_n_rhf = subprocess.run(cmd_neutral, shell=True, cwd=xtb_dir_rhf, 
                               capture_output=True, text=True)
    
    if res_n_rhf.returncode != 0:
        print(f"❌ xTB neutral ({method}) failed!")
        print(res_n_rhf.stderr)
        # Still try to read output files
    else:
        print("✅ Neutral calculation completed")
    
    # Cation (forced closed-shell)
    cmd_cation = f"xtb coord.xyz {flag} --charge 1 --opt > xtb_cation.out 2>&1"
    res_c_rhf = subprocess.run(cmd_cation, shell=True, cwd=xtb_dir_rhf,
                               capture_output=True, text=True)
    
    if res_c_rhf.returncode != 0:
        print(f"❌ xTB cation ({method}) failed!")
        print(res_c_rhf.stderr)
    else:
        print("✅ Cation calculation completed")
    
    def get_energy_from_file(directory, filename):
        """Parse energy from xTB output file"""
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
                                if abs(val) > 0.0001:
                                    return val * 27.211386245988
                            except ValueError:
                                continue
        except:
            return None
        return None
    
    def get_energy_from_stdout(output_text):
        """Parse energy from stdout"""
        for line in output_text.split('\n'):
            if "total energy" in line.lower():
                numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                for num in numbers:
                    try:
                        val = float(num)
                        if abs(val) > 0.0001:
                            return val * 27.211386245988
                    except ValueError:
                        continue
        return None
    
    # Try to get energies from output files first
    e_neutral_rhf = get_energy_from_file(xtb_dir_rhf, 'xtb_neutral.out')
    e_cation_rhf = get_energy_from_file(xtb_dir_rhf, 'xtb_cation.out')
    
    # If not found in files, try stdout
    if e_neutral_rhf is None:
        e_neutral_rhf = get_energy_from_stdout(res_n_rhf.stdout)
    if e_cation_rhf is None:
        e_cation_rhf = get_energy_from_stdout(res_c_rhf.stdout)
    
    if e_neutral_rhf is not None and e_cation_rhf is not None:
        ie_rhf = e_cation_rhf - e_neutral_rhf
        print(f"Neutral Energy: {e_neutral_rhf:.4f} eV")
        print(f"Cation Energy:  {e_cation_rhf:.4f} eV")
        print(f"Vertical IE:    {ie_rhf:.4f} eV")
        print(f"Experimental:   10.44 eV (difference: {ie_rhf-10.44:+.4f} eV)")
        
        results['RHF-like'] = {
            'neutral': e_neutral_rhf,
            'cation': e_cation_rhf,
            'ie': ie_rhf
        }
    else:
        print(f"❌ Could not parse {method} RHF output")
        print(f"Neutral found: {e_neutral_rhf is not None}, Cation found: {e_cation_rhf is not None}")
        results['RHF-like'] = None
    
    # === UHF Calculation (Open-Shell) ===
    print(f"\n--- xTB {method} (UHF - Open Shell) ---")
    
    # Clean directory for UHF
    xtb_dir_uhf = f'xtb_run_{method.lower()}_uhf'
    if os.path.exists(xtb_dir_uhf):
        shutil.rmtree(xtb_dir_uhf)
    os.makedirs(xtb_dir_uhf, exist_ok=True)
    
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write(f'{xtb_dir_uhf}/coord.xyz', atoms)
    
    # Also create a .xyz file with proper format
    with open(f'{xtb_dir_uhf}/coord', 'w') as f:
        f.write("1\n")
        f.write("Hg atom\n")
        f.write("Hg    0.000000    0.000000    0.000000\n")
    
    # Neutral (closed-shell, singlet)
    cmd_neutral = f"xtb coord.xyz {flag} --opt > xtb_neutral.out 2>&1"
    res_n_uhf = subprocess.run(cmd_neutral, shell=True, cwd=xtb_dir_uhf, 
                               capture_output=True, text=True)
    
    if res_n_uhf.returncode != 0:
        print(f"❌ xTB neutral ({method}) failed!")
        print(res_n_uhf.stderr)
    else:
        print("✅ Neutral calculation completed")
    
    # Cation (open-shell, doublet with UHF)
    cmd_cation = f"xtb coord.xyz {flag} --charge 1 --uhf 1 --opt > xtb_cation.out 2>&1"
    res_c_uhf = subprocess.run(cmd_cation, shell=True, cwd=xtb_dir_uhf,
                               capture_output=True, text=True)
    
    if res_c_uhf.returncode != 0:
        print(f"❌ xTB cation ({method}) failed!")
        print(res_c_uhf.stderr)
    else:
        print("✅ Cation calculation completed")
    
    # Try to get energies from output files first
    e_neutral_uhf = get_energy_from_file(xtb_dir_uhf, 'xtb_neutral.out')
    e_cation_uhf = get_energy_from_file(xtb_dir_uhf, 'xtb_cation.out')
    
    # If not found in files, try stdout
    if e_neutral_uhf is None:
        e_neutral_uhf = get_energy_from_stdout(res_n_uhf.stdout)
    if e_cation_uhf is None:
        e_cation_uhf = get_energy_from_stdout(res_c_uhf.stdout)
    
    if e_neutral_uhf is not None and e_cation_uhf is not None:
        ie_uhf = e_cation_uhf - e_neutral_uhf
        print(f"Neutral Energy: {e_neutral_uhf:.4f} eV")
        print(f"Cation Energy:  {e_cation_uhf:.4f} eV")
        print(f"Vertical IE:    {ie_uhf:.4f} eV")
        print(f"Experimental:   10.44 eV (difference: {ie_uhf-10.44:+.4f} eV)")
        
        results['UHF'] = {
            'neutral': e_neutral_uhf,
            'cation': e_cation_uhf,
            'ie': ie_uhf
        }
    else:
        print(f"❌ Could not parse {method} UHF output")
        print(f"Neutral found: {e_neutral_uhf is not None}, Cation found: {e_cation_uhf is not None}")
        results['UHF'] = None
    
    # List files created
    print(f"\n📁 Files created in {xtb_dir_rhf}:")
    if os.path.exists(xtb_dir_rhf):
        files = os.listdir(xtb_dir_rhf)
        for f in sorted(files):
            size = os.path.getsize(os.path.join(xtb_dir_rhf, f))
            print(f"   - {f} ({size} bytes)")
    
    print(f"\n📁 Files created in {xtb_dir_uhf}:")
    if os.path.exists(xtb_dir_uhf):
        files = os.listdir(xtb_dir_uhf)
        for f in sorted(files):
            size = os.path.getsize(os.path.join(xtb_dir_uhf, f))
            print(f"   - {f} ({size} bytes)")
    
    return results

def print_summary(mopac_results, xtb_results):
    """Print a summary of all results"""
    print("\n" + "="*70)
    print("SUMMARY OF IONIZATION ENERGY CALCULATIONS FOR Hg")
    print("="*70)
    print(f"{'Method':<20} {'Neutral (eV)':<15} {'Cation (eV)':<15} {'IE (eV)':<12} {'Error (eV)':<12}")
    print("-"*70)
    
    # MOPAC results
    if mopac_results:
        for method_type in ['RHF', 'UHF']:
            if method_type in mopac_results:
                print(f"{f'MOPAC PM7 {method_type}':<20} "
                      f"{mopac_results[method_type]['neutral']:<15.4f} "
                      f"{mopac_results[method_type]['cation']:<15.4f} "
                      f"{mopac_results[method_type]['ie']:<12.4f} "
                      f"{mopac_results[method_type]['ie']-10.44:<+12.4f}")
    
    # xTB results
    if xtb_results:
        for method, results in xtb_results.items():
            if results and results.get('RHF-like'):
                print(f"{f'xTB {method} RHF-like':<20} "
                      f"{results['RHF-like']['neutral']:<15.4f} "
                      f"{results['RHF-like']['cation']:<15.4f} "
                      f"{results['RHF-like']['ie']:<12.4f} "
                      f"{results['RHF-like']['ie']-10.44:<+12.4f}")
            
            if results and results.get('UHF'):
                print(f"{f'xTB {method} UHF':<20} "
                      f"{results['UHF']['neutral']:<15.4f} "
                      f"{results['UHF']['cation']:<15.4f} "
                      f"{results['UHF']['ie']:<12.4f} "
                      f"{results['UHF']['ie']-10.44:<+12.4f}")
    
    print("-"*70)
    print(f"{'Experimental':<20} {'':<15} {'':<15} {'10.4400':<12} {'0.0000':<12}")
    print("="*70)

if __name__ == "__main__":
    print("Starting Ionization Energy calculations for Hg atom...")
    print("="*70)
    
    # Run MOPAC with both RHF and UHF
    mopac_results = run_mopac_ie()
    
    # Run xTB with both RHF-like and UHF for different methods
    xtb_results = {}
    
    # Try different xTB methods
    for method in ['GFN2', 'GFN1']:  # Add 'GFN0', 'GFNFF' if needed
        print("\n" + "="*70)
        xtb_results[method] = run_xtb_cli_ie(method)
    
    # Print summary
    print_summary(mopac_results, xtb_results)
    
    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS:")
    print("="*70)
    print("✅ For Hg atom, MOPAC PM7 with UHF gives the best results")
    print(f"   IE = {mopac_results['UHF']['ie']:.4f} eV (error: {mopac_results['UHF']['ie']-10.44:+.4f} eV)")
    print("\n⚠️  xTB GFN2 significantly overestimates IE for Hg")
    print("   Consider using:")
    print("   - MOPAC PM7 for heavy elements (recommended)")
    print("   - DFT with pseudopotentials for more accurate results")
    print("   - xTB GFN1/GFN-FF for qualitative trends only")
    print("\n📁 All xTB output files are saved in their respective directories:")
    print("   - xtb_run_gfn1_rhf/, xtb_run_gfn1_uhf/")
    print("   - xtb_run_gfn2_rhf/, xtb_run_gfn2_uhf/")
    print("="*70)
