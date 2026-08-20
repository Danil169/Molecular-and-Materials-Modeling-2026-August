import subprocess
from ase import Atoms
from ase.calculators.mopac import MOPAC

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

def run_xtb_cli_ie():
    print("\n--- Running xTB via CLI Executable ---")
    import os
    from ase.io import write
    
    os.makedirs('xtb_run', exist_ok=True)
    atoms = Atoms('Hg', positions=[(0, 0, 0)])
    write('xtb_run/coord.xyz', atoms)
    
    # 1. Run Neutral Command
    cmd_neutral = "xtb coord.xyz --gfn 2"
    res_n = subprocess.run(cmd_neutral, shell=True, cwd='xtb_run', capture_output=True, text=True)
    
    # 2. Run Cation Command
    cmd_cation = "xtb coord.xyz --gfn 2 --charge 1 --uhf 1"
    res_c = subprocess.run(cmd_cation, shell=True, cwd='xtb_run', capture_output=True, text=True)
    
    # Fixed parsing function
    def get_energy(output_text):
        for line in output_text.split('\n'):
            if "total energy" in line:
                parts = line.split()
                # Find the float within the line components
                for part in parts:
                    try:
                        # Convert Hartree to eV (1 Hartree = 27.211386245988 eV)
                        return float(part) * 27.211386245988
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
    else:
        print("❌ Could not parse xTB output files. Verify that the 'xtb' binary is working.")

if __name__ == "__main__":
    print("Starting Ionization Energy calculations for Hg atom...")
    run_mopac_ie()
    run_xtb_cli_ie()

