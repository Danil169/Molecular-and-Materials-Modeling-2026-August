#!/usr/bin/env python
"""
PySCF SCF + CCSD(T) Calculation for Hg Atom
- Neutral Hg: RHF + CCSD(T) (closed-shell, singlet)
- Cation Hg+: UHF + CCSD(T) (open-shell, doublet)
- Uses ECP and appropriate basis sets for heavy elements
- Supports MPI parallelization with nproc from file
- Minimal output for cleaner results
"""

import time
import datetime
import numpy as np
import os
from pyscf import gto, scf, cc

# ============================================================
# USER-DEFINED PARAMETERS - Modify these as needed
# ============================================================

# Basis set and ECP
BASIS_SET = "def2-tzvp"      # Can also use "def2-qzvp" for higher accuracy
ECP_SET = "def2-tzvp"        # ECP for Hg

# SCF settings
SCF_CONV_TOL = 1e-9
SCF_MAX_CYCLE = 200

# CCSD settings
CCSD_CONV_TOL = 1e-8
CCSD_MAX_CYCLE = 100
CCSD_DIIS_SPACE = 10

# Memory (per processor)
MEMORY_PER_PROC = 16000  # MB per processor

# Verbosity (0 = silent, 1 = minimal, 2 = normal, 3 = verbose)
VERBOSE = 1  # Set to 0 for minimal output

# Experimental reference
EXP_IE = 10.437  # eV

# File containing number of processors
NPROC_FILE = "nproc.txt"

# ============================================================
# END OF USER PARAMETERS
# ============================================================

def read_nproc_from_file(filename="nproc.txt"):
    """Read number of processors from a text file."""
    nproc = 1
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as f:
                nproc = int(f.read().strip())
                if VERBOSE >= 1:
                    print(f"✅ Read nproc = {nproc} from {filename}")
                return nproc
        except (ValueError, IOError):
            return 1
    else:
        with open(filename, 'w') as f:
            f.write("1\n")
        return 1

def setup_molecule(charge=0, spin=0, nproc=1):
    """Set up the Hg atom/molecule with appropriate basis and ECP."""
    mol = gto.Mole()
    mol.atom = [['Hg', (0, 0, 0)]]
    mol.basis = {'Hg': BASIS_SET}
    mol.ecp = {'Hg': ECP_SET}
    mol.charge = charge
    mol.spin = spin
    mol.verbose = VERBOSE  # Control output verbosity
    mol.max_memory = MEMORY_PER_PROC * nproc
    mol.build()
    
    if VERBOSE >= 1:
        print(f"\n📋 {['Neutral', 'Cation'][charge]} Hg ({'RHF' if spin==0 else 'UHF'})")
        print(f"   Electrons: {mol.nelectron}, Basis: {mol.nao}, Memory: {mol.max_memory} MB")
    
    return mol

def run_scf(mol):
    """Run Hartree-Fock calculation with minimal output."""
    if VERBOSE >= 1:
        print("   SCF...", end=" ", flush=True)
    
    mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = SCF_MAX_CYCLE
    mf.diis_space = 10
    mf.verbose = 0  # Suppress SCF output
    
    start_time = time.time()
    mf.kernel()
    scf_time = time.time() - start_time
    
    if VERBOSE >= 1:
        print(f"done ({scf_time:.2f}s)")
    
    return mf, scf_time

def get_t1_diagnostic(mycc):
    """Get T1 diagnostic for both RHF and UHF."""
    try:
        return mycc.get_t1_diagnostic()
    except AttributeError:
        try:
            t1 = mycc.t1
            nocc = mycc.nocc
            if isinstance(nocc, tuple):
                nocc_alpha, nocc_beta = nocc
                t1_alpha_norm = np.linalg.norm(t1[0].ravel())
                t1_beta_norm = np.linalg.norm(t1[1].ravel())
                t1_norm = (t1_alpha_norm + t1_beta_norm) / 2
                return t1_norm / np.sqrt(nocc_alpha + nocc_beta)
            else:
                return np.linalg.norm(t1.ravel()) / np.sqrt(nocc)
        except (AttributeError, TypeError):
            return 0.0

def run_ccsd_t(mf):
    """Run CCSD(T) calculation with minimal output."""
    if VERBOSE >= 1:
        print("   CCSD(T)...", end=" ", flush=True)
    
    mycc = cc.CCSD(mf)
    mycc.conv_tol = CCSD_CONV_TOL
    mycc.max_cycle = CCSD_MAX_CYCLE
    mycc.diis_space = CCSD_DIIS_SPACE
    mycc.verbose = 0  # Suppress CCSD output
    
    start_time = time.time()
    e_corr, t1, t2 = mycc.kernel()
    ccsd_time = time.time() - start_time
    
    start_time_t = time.time()
    e_t = mycc.ccsd_t()
    t_time = time.time() - start_time_t
    
    e_ccsd = mf.e_tot + e_corr
    e_ccsdt = e_ccsd + e_t
    
    if VERBOSE >= 1:
        print(f"done ({ccsd_time + t_time:.2f}s)")
    
    return {
        'e_ccsdt': e_ccsdt,
        'ccsd_time': ccsd_time,
        't_time': t_time,
        't1_diagnostic': get_t1_diagnostic(mycc)
    }

def calculate_ionization_potential(nproc):
    """Calculate the vertical ionization potential of Hg."""
    print("="*60)
    print(f"Hg IP: {BASIS_SET} basis, {nproc} procs")
    print(f"Started: {datetime.datetime.now().strftime('%H:%M:%S')}")
    print("="*60)
    
    total_start = time.time()
    
    # Neutral Hg (RHF + CCSD(T))
    print("\n🔹 Neutral Hg (RHF)")
    mol_neutral = setup_molecule(charge=0, spin=0, nproc=nproc)
    mf_neutral, scf_time_neutral = run_scf(mol_neutral)
    ccsdt_neutral = run_ccsd_t(mf_neutral)
    
    # Cation Hg+ (UHF + CCSD(T))
    print("\n🔹 Cation Hg+ (UHF)")
    mol_cation = setup_molecule(charge=1, spin=1, nproc=nproc)
    mf_cation, scf_time_cation = run_scf(mol_cation)
    ccsdt_cation = run_ccsd_t(mf_cation)
    
    # Ionization Potential
    ie_hartree = ccsdt_cation['e_ccsdt'] - ccsdt_neutral['e_ccsdt']
    ie_ev = ie_hartree * 27.211386245988
    total_time = time.time() - total_start
    
    # Results
    print("\n" + "="*60)
    print("RESULTS")
    print("="*60)
    print(f"Neutral CCSD(T):  {ccsdt_neutral['e_ccsdt']:.8f} Eh")
    print(f"Cation CCSD(T):   {ccsdt_cation['e_ccsdt']:.8f} Eh")
    print("-"*60)
    print(f"IE:               {ie_ev:.4f} eV")
    print(f"Experiment:       {EXP_IE:.4f} eV")
    print(f"Error:            {ie_ev-EXP_IE:+.4f} eV")
    print("="*60)
    
    # Timing summary
    print(f"\n⏱️  SCF:     {scf_time_neutral + scf_time_cation:.2f}s")
    print(f"   CCSD(T): {ccsdt_neutral['ccsd_time'] + ccsdt_neutral['t_time'] + ccsdt_cation['ccsd_time'] + ccsdt_cation['t_time']:.2f}s")
    print(f"   Total:   {total_time:.2f}s ({total_time/60:.2f} min)")
    
    # T1 diagnostics
    print(f"\n📊 T1 diag: Neutral={ccsdt_neutral['t1_diagnostic']:.4f}, Cation={ccsdt_cation['t1_diagnostic']:.4f}")
    
    return {
        'ie_ev': ie_ev,
        'total_time': total_time,
        'nproc': nproc
    }

def main():
    """Main execution function."""
    nproc = read_nproc_from_file(NPROC_FILE)
    results = calculate_ionization_potential(nproc)
    
    print("\n" + "="*60)
    print(f"✅ IE = {results['ie_ev']:.4f} eV  (error: {results['ie_ev']-EXP_IE:+.4f} eV)")
    print(f"   Processors: {results['nproc']}, Time: {results['total_time']/60:.2f} min")
    print("="*60)

if __name__ == "__main__":
    main()
