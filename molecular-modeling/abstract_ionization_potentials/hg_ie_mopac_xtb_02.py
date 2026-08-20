import subprocess
import os
import re
from ase import Atoms
from ase.calculators.mopac import MOPAC
from ase.io import write

def run_mopac_ie():
    print("\n--- Running MOPAC (PM7) ---")
    neutral_atom = Atoms('Hg', positions=[(0, 0, 0)])
    neutral_atom.calc = MOPAC(label='mopac_neutral/mop', task='1SCF GRADIENTS', keywords='PM7')
    e_neutral = neutral_atom.get_potential_energy()
    
    cation_atom = Atoms('Hg', positions=[(0, 0, 0)])
    cation_atom.calc = MOPAC(label='mopac_cation/mop', task='1SCF GRADIENTS', charge=1, keywords='PM7 DOUBLET')
    e_cation = cation_atom.get_potential_energy()
    
    ie = e_cation - e_neutral
    print(f"Neutral Energy: {e_neutral:.4f} eV")
    print(f"Cation Energy:  {e_cation:.4f} eV")
    print(f"Vertical IE:    {ie:.4f} eV")
    print(f"Experimental:   10.44 eV (difference: {ie-10.44:.4f} eV)")

def run_xtb_cli_ie():
    print("\n--- Running xTB (GFN2) via CLI ---")
    
    # Create fresh directory
    import shutil
    if os.path.exists('xtb_run'):
        shutil.rmtree('xtb_run')
    os.makedirs('xtb_run', exist_ok=True)
    
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write('xtb_run/coord.xyz', atoms)
    
    # Run neutral (singlet)
    cmd_neutral = "xtb coord.xyz --gfn 2"
    res_n = subprocess.run(cmd_neutral, shell=True, cwd='xtb_run', 
                          capture_output=True, text=True)
    
    if res_n.returncode != 0:
        print("❌ xTB neutral calculation failed!")
        print(res_n.stderr)
        return
    
    # Run cation (doublet)
    cmd_cation = "xtb coord.xyz --gfn 2 --charge 1 --uhf 1"
    res_c = subprocess.run(cmd_cation, shell=True, cwd='xtb_run',
                          capture_output=True, text=True)
    
    if res_c.returncode != 0:
        print("❌ xTB cation calculation failed!")
        print(res_c.stderr)
        return
    
    def get_energy(output_text):
        for line in output_text.split('\n'):
            if "total energy" in line.lower():
                # Find numbers in the line
                numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                for num in numbers:
                    try:
                        val = float(num)
                        # Skip small numbers (like 0.000)
                        if abs(val) > 0.0001:
                            # Return in eV
                            return val * 27.211386245988
                    except ValueError:
                        continue
        return None
    
    e_neutral = get_energy(res_n.stdout)
    e_cation = get_energy(res_c.stdout)
    
    if e_neutral is not None and e_cation is not None:
        ie = e_cation - e_neutral
        print(f"Neutral Energy: {e_neutral:.4f} eV")
        print(f"Cation Energy:  {e_cation:.4f} eV")
        print(f"Vertical IE:    {ie:.4f} eV")
        print(f"Experimental:   10.44 eV (difference: {ie-10.44:.4f} eV)")
    else:
        print("❌ Could not parse xTB output files.")
        print("Neutral output (first 500 chars):")
        print(res_n.stdout[:500])
        print("Cation output (first 500 chars):")
        print(res_c.stdout[:500])

if __name__ == "__main__":
    print("Starting Ionization Energy calculations for Hg atom...")
    run_mopac_ie()
    run_xtb_cli_ie()
