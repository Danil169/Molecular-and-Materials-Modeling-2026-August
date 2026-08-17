import os
import sys
import time
import subprocess
import shutil
from ase import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator
from ase.calculators.mopac import MOPAC

# =====================================================================
# 1. Environment & Path Resolution
# =====================================================================
conda_prefix = os.environ.get('CONDA_PREFIX')
if conda_prefix:
    os.environ['ASE_MOPAC_COMMAND'] = f"{os.path.join(conda_prefix, 'bin', 'mopac')} PREFIX.mop"
    xtb_bin = os.path.join(conda_prefix, 'bin', 'xtb')
    nwchem_bin = os.path.join(conda_prefix, 'bin', 'nwchem')
else:
    xtb_bin = "xtb"
    nwchem_bin = "nwchem"

# =====================================================================
# 2. Custom Monolithic Calculator Wrappers
# =====================================================================
class CustomPySCFCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='RHF', basis='3-21g', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        import pyscf
        xyz_coords = [f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}" 
                      for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions())]
        atom_str = "; ".join(xyz_coords)
        mol = pyscf.gto.Mole()
        mol.atom = atom_str
        mol.basis = self.basis
        mol.verbose = 0
        mol.build()
        mf = pyscf.scf.RHF(mol) if self.method.upper() == 'RHF' else pyscf.scf.KS(mol)
        if self.method.upper() != 'RHF': mf.xc = self.method
        self.results['energy'] = mf.kernel() * 27.211386245988

class CustomXTBCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='2', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        xyz_file = "tmp_bench_xtb.xyz"
        with open(xyz_file, "w") as f:
            f.write(f"{len(self.atoms)}\n\n")
            for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
                f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        try:
            cmd = [xtb_bin, xyz_file, "--gfn", self.method]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            for line in res.stdout.splitlines():
                if "total energy" in line.lower():
                    for token in line.split():
                        try:
                            self.results['energy'] = float(token.strip('|').strip()) * 27.211386245988
                            break
                        except ValueError: continue
                    break
        finally:
            for f in [xyz_file, "xtbopt.xyz", "xtbopt.log", "wbo", "charges", "xtbrestart", "gfnff_topo"]:
                if os.path.exists(f): os.remove(f)

class CustomNWChemCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='dft', basis='3-21g', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        nw_file = "tmp_bench_nwchem.nw"
        
        # ✅ FIX: Simplified geometry signature block to avoid syntax aborts
        with open(nw_file, "w") as f:
            f.write("geometry\n")
            for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
                f.write(f"  {sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
            f.write(f"end\nbasis\n  * library {self.basis}\nend\ntask {self.method} energy\n")
            
        try:
            res = subprocess.run([nwchem_bin, nw_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            energy = None
            for line in res.stdout.splitlines():
                if "total dft energy" in line.lower() or "total energy" in line.lower():
                    for token in line.split():
                        try:
                            energy = float(token) * 27.211386245988
                            break
                        except ValueError: continue
                    if energy is not None: break
            
            if energy is not None:
                self.results['energy'] = energy
            else:
                raise RuntimeError(f"Could not parse energy header. NWChem exit status: {res.returncode}")
        finally:
            if os.path.exists(nw_file): os.remove(nw_file)

# =====================================================================
# 3. Cluster Generation Routine
# =====================================================================
def generate_water_cluster(n_molecules):
    cluster = Atoms()
    for i in range(n_molecules):
        h2o = molecule('H2O')
        h2o.translate([i * 3.5, 0.0, 0.0])
        cluster += h2o
    return cluster

# =====================================================================
# 4. Main Benchmark Setup
# =====================================================================
if __name__ == "__main__":
    cluster_sizes = [1, 2, 4, 8]
    methods = {
        'MOPAC (PM7)': lambda: MOPAC(method='PM7'),
        'xTB (GFN2)': lambda: CustomXTBCalculator(method='2'),
        'PySCF (RHF/3-21G)': lambda: CustomPySCFCalculator(method='RHF', basis='3-21g'),
        'NWChem (DFT/3-21G)': lambda: CustomNWChemCalculator(method='dft', basis='3-21g')
    }

    results = {m: [] for m in methods}

    print("==========================================================")
    print("      COMPUTATIONAL CHEMISTRY ENGINE CPU WALL BENCHMARK   ")
    print("==========================================================")

    for n in cluster_sizes:
        atoms = generate_water_cluster(n)
        n_atoms = len(atoms)
        print(f"\n🚀 Cluster Size: {n} H2O molecules ({n_atoms} atoms)")
        print("-" * 58)

        for name, calc_init in methods.items():
            try:
                atoms.calc = calc_init()
                
                start_time = time.perf_counter()
                energy = atoms.get_potential_energy()
                end_time = time.perf_counter()
                
                wall_time = end_time - start_time
                results[name].append(f"{wall_time:.4f}s")
                print(f"  |-- {name:<20} : {wall_time:>8.4f} seconds (E = {energy:.2f} eV)")
            except Exception as e:
                results[name].append("FAILED")
                # ✅ diagnostic error tracking visible to avoid ghost bugs
                print(f"  |-- {name:<20} : ❌ FAILED (Reason: {e})")

            for junk in ['mopac.out', 'mopac.arc', 'mopac.mop']:
                if os.path.exists(junk): os.remove(junk)

    # =====================================================================
    # 5. Summary Generation
    # =====================================================================
    print("\n\n==========================================================")
    print("                    BENCHMARK SUMMARY                     ")
    print("==========================================================")
    
    header_str = f"| {'Method':<20} "
    for n in cluster_sizes:
        header_str += f"| {n} H2O ({n*3} atoms) "
    print(header_str + "|")
    print("|" + "----------------------|" * (len(cluster_sizes) + 1))
    
    for name in methods:
        row_str = f"| {name:<20} "
        for idx, n in enumerate(cluster_sizes):
            row_str += f"| {results[name][idx]:^15} "
        print(row_str + "|")
    print("==========================================================")

