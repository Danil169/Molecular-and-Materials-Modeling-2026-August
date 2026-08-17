import os
import subprocess
from ase.build import molecule
from ase.calculators.mopac import MOPAC
from ase.calculators.xtb import XTB
from ase.calculators.nwchem import NWChem
from ase.calculators.espresso import Espresso
# ... (Assume PySCFCalculator class defined here for brevity)

# --- Environment Setup ---
conda_prefix = os.environ.get('CONDA_PREFIX')
if conda_prefix:
    os.environ['ASE_MOPAC_COMMAND'] = f"{os.path.join(conda_prefix, 'bin', 'mopac')} PREFIX.mop"
    os.environ['ASE_NWCHEM_COMMAND'] = os.path.join(conda_prefix, 'bin', 'nwchem')
    os.environ['ASE_ESPRESSO_COMMAND'] = f"{os.path.join(conda_prefix, 'bin', 'pw.x')} -in PREFIX.pwi > PREFIX.pwo"
    pseudo_dir = os.path.join(os.getcwd(), "pseudo_mock")
    os.makedirs(pseudo_dir, exist_ok=True)
    with open(os.path.join(pseudo_dir, "H.upf"), "w") as f: f.write("<!-- Mock -->")
else:
    pseudo_dir = os.getcwd()

# --- Execution Pipeline ---
if __name__ == "__main__":
    h2o = molecule('H2O')
    h2 = molecule('H2')
    h2.center(vacuum=5.0)

    # 1. MOPAC
    try:
        h2o.calc = MOPAC(method='PM7')
        print(f"MOPAC: {h2o.get_potential_energy():.4f} eV")
    except Exception as e: print(f"MOPAC FAILED: {e}")

    # 2. xTB
    try:
        h2o.calc = XTB(method='GFN2-xTB')
        print(f"xTB: {h2o.get_potential_energy():.4f} eV")
    except Exception as e: print(f"xTB FAILED: {e}")

    # 3. CREST (CLI Check)
    try:
        crest_path = os.path.join(conda_prefix, 'bin', 'crest') if conda_prefix else 'crest'
        result = subprocess.run([crest_path, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"CREST Active: {result.returncode == 0}")
    except Exception as e: print(f"CREST FAILED: {e}")

    # 4. NWChem
    try:
        h2o.calc = NWChem(xc='PBE', basis='3-21G')
        print(f"NWChem: {h2o.get_potential_energy():.4f} eV")
    except Exception as e: print(f"NWChem FAILED: {e}")

    # 5. Quantum Espresso
    try:
        input_data = {'control': {'calculation': 'scf', 'pseudo_dir': pseudo_dir}, 'system': {'ecutwfc': 20.0}, 'electrons': {'conv_thr': 1e-4}}
        h2.calc = Espresso(input_data=input_data, pseudopotentials={'H': 'H.upf'})
        print(f"QE: {h2.get_potential_energy():.4f} eV")
    except Exception as e: print(f"QE FAILED/Mocked: {e}")

