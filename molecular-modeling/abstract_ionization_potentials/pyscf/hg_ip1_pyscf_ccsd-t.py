#!/usr/bin/env python
"""
PySCF SCF + CCSD(T) Calculation for Hg Atom
- Neutral Hg: RHF + CCSD(T) (closed-shell, singlet)
- Cation Hg+: UHF + CCSD(T) (open-shell, doublet)
- Uses ECP and appropriate basis sets for heavy elements
"""

import time
import datetime
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

# Experimental reference
EXP_IE = 10.437  # eV

# ============================================================
# END OF USER PARAMETERS
# ============================================================

def setup_molecule(charge=0, spin=0, basis=BASIS_SET, ecp=ECP_SET):
    """
    Set up the Hg atom/molecule with appropriate basis and ECP.
    
    Args:
        charge: Net charge (0 for neutral, 1 for cation)
        spin: 2S (0 for singlet, 1 for doublet)
        basis: Basis set name
        ecp: ECP name
    
    Returns:
        pyscf.gto.Mole object
    """
    mol = gto.Mole()
    mol.atom = [['Hg', (0, 0, 0)]]
    mol.basis = {'Hg': basis}
    mol.ecp = {'Hg': ecp}
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 4
    mol.max_memory = 32000  # MB
    mol.build()
    
    print(f"\n📋 Molecule Setup:")
    print(f"   Atom: Hg (Z={mol.atom_charge(0):.0f})")
    print(f"   Charge: {charge}")
    print(f"   Spin: {spin} (2S={spin})")
    print(f"   Basis: {basis}")
    print(f"   ECP: {ecp}")
    print(f"   Number of electrons: {mol.nelectron}")
    print(f"   Number of basis functions: {mol.nao}")
    
    return mol

def run_scf(mol, conv_tol=SCF_CONV_TOL, max_cycle=SCF_MAX_CYCLE):
    """
    Run Hartree-Fock calculation.
    
    Args:
        mol: pyscf.gto.Mole object
        conv_tol: SCF convergence tolerance
        max_cycle: Maximum SCF iterations
    
    Returns:
        SCF object
    """
    print("\n" + "="*70)
    print("SCF CALCULATION")
    print("="*70)
    
    # Determine if RHF or UHF based on spin
    if mol.spin == 0:
        mf = scf.RHF(mol)
        method = "RHF"
    else:
        mf = scf.UHF(mol)
        method = "UHF"
    
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    mf.diis_space = 10
    mf.chkfile = f'hg_{method.lower()}.chk'
    
    print(f"Method: {method}")
    print(f"Conv. tol: {conv_tol}")
    print(f"Max cycles: {max_cycle}")
    print("-"*70)
    
    start_time = time.time()
    mf.kernel()
    scf_time = time.time() - start_time
    
    print(f"\n✅ SCF converged in {mf.converged} cycles")
    print(f"   Energy: {mf.e_tot:.10f} Hartree")
    print(f"   Energy: {mf.e_tot * 27.211386245988:.6f} eV")
    print(f"   Time: {scf_time:.2f} seconds")
    
    return mf, scf_time

def run_ccsd_t(mf, method_name="CCSD(T)"):
    """
    Run CCSD(T) calculation from SCF object.
    
    Args:
        mf: SCF object (RHF or UHF)
        method_name: Name for printing
    
    Returns:
        Dictionary with CCSD(T) results
    """
    print("\n" + "="*70)
    print(f"{method_name} CALCULATION")
    print("="*70)
    
    # Initialize CCSD
    mycc = cc.CCSD(mf)
    mycc.conv_tol = CCSD_CONV_TOL
    mycc.max_cycle = CCSD_MAX_CYCLE
    mycc.diis_space = CCSD_DIIS_SPACE
    
    print(f"Conv. tol: {CCSD_CONV_TOL}")
    print(f"Max cycles: {CCSD_MAX_CYCLE}")
    print(f"DIIS space: {CCSD_DIIS_SPACE}")
    print("-"*70)
    
    # Run CCSD
    start_time = time.time()
    e_corr, t1, t2 = mycc.kernel()
    ccsd_time = time.time() - start_time
    
    # Calculate (T) correction
    start_time_t = time.time()
    e_t = mycc.ccsd_t()
    t_time = time.time() - start_time_t
    
    # Total energies
    e_ccsd = mf.e_tot + e_corr
    e_ccsdt = e_ccsd + e_t
    
    print(f"\n✅ CCSD converged in {mycc.converged} cycles")
    print(f"   CCSD correlation energy: {e_corr:.10f} Hartree")
    print(f"   CCSD total energy: {e_ccsd:.10f} Hartree")
    print(f"   (T) correction: {e_t:.10f} Hartree")
    print(f"   CCSD(T) total energy: {e_ccsdt:.10f} Hartree")
    print(f"   CCSD(T) energy: {e_ccsdt * 27.211386245988:.6f} eV")
    print(f"   CCSD time: {ccsd_time:.2f} seconds")
    print(f"   (T) time: {t_time:.2f} seconds")
    
    # T1 diagnostic (for assessing multi-reference character)
    t1_diagnostic = mycc.get_t1_diagnostic()
    print(f"\n📊 T1 diagnostic: {t1_diagnostic:.4f}")
    if t1_diagnostic < 0.02:
        print("   ✅ Single-reference method is appropriate (T1 < 0.02)")
    elif t1_diagnostic < 0.04:
        print("   ⚠️  Moderate multi-reference character (T1 = 0.02-0.04)")
    else:
        print("   ❌ Strong multi-reference character (T1 > 0.04)")
    
    return {
        'e_corr': e_corr,
        'e_ccsd': e_ccsd,
        'e_t': e_t,
        'e_ccsdt': e_ccsdt,
        't1': t1,
        't2': t2,
        'ccsd_time': ccsd_time,
        't_time': t_time,
        'converged': mycc.converged,
        't1_diagnostic': t1_diagnostic
    }

def calculate_ionization_potential():
    """
    Calculate the vertical ionization potential of Hg.
    """
    print("="*70)
    print("Hg IONIZATION POTENTIAL CALCULATION")
    print(f"Using {BASIS_SET} basis set with ECP")
    print(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)
    
    total_start = time.time()
    
    # === Neutral Hg (RHF + CCSD(T)) ===
    print("\n" + "="*70)
    print("NEUTRAL Hg (RHF + CCSD(T))")
    print("="*70)
    
    mol_neutral = setup_molecule(charge=0, spin=0)
    mf_neutral, scf_time_neutral = run_scf(mol_neutral)
    ccsdt_neutral = run_ccsd_t(mf_neutral, "CCSD(T) (Neutral)")
    
    # === Cation Hg+ (UHF + CCSD(T)) ===
    print("\n" + "="*70)
    print("CATION Hg+ (UHF + CCSD(T))")
    print("="*70)
    
    mol_cation = setup_molecule(charge=1, spin=1)
    mf_cation, scf_time_cation = run_scf(mol_cation)
    ccsdt_cation = run_ccsd_t(mf_cation, "CCSD(T) (Cation)")
    
    # === Ionization Potential ===
    ie_hartree = ccsdt_cation['e_ccsdt'] - ccsdt_neutral['e_ccsdt']
    ie_ev = ie_hartree * 27.211386245988
    
    total_time = time.time() - total_start
    
    # Print results
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"{'Neutral Energy (CCSD(T)):':<35} {ccsdt_neutral['e_ccsdt']:>15.10f} Hartree")
    print(f"{'Cation Energy (CCSD(T)):':<35} {ccsdt_cation['e_ccsdt']:>15.10f} Hartree")
    print("-"*70)
    print(f"{'Vertical IE:':<35} {ie_hartree:>15.10f} Hartree")
    print(f"{'Vertical IE:':<35} {ie_ev:>15.6f} eV")
    print(f"{'Experimental IE:':<35} {EXP_IE:>15.4f} eV")
    print(f"{'Error:':<35} {ie_ev-EXP_IE:>+15.4f} eV")
    print("="*70)
    
    # Timing summary
    print("\n" + "="*70)
    print("TIMING SUMMARY")
    print("="*70)
    print(f"{'Neutral SCF:':<35} {scf_time_neutral:>10.2f} s")
    print(f"{'Neutral CCSD(T):':<35} {ccsdt_neutral['ccsd_time'] + ccsdt_neutral['t_time']:>10.2f} s")
    print(f"{'Cation SCF:':<35} {scf_time_cation:>10.2f} s")
    print(f"{'Cation CCSD(T):':<35} {ccsdt_cation['ccsd_time'] + ccsdt_cation['t_time']:>10.2f} s")
    print("-"*70)
    print(f"{'Total time:':<35} {total_time:>10.2f} s")
    print(f"{'Total time:':<35} {total_time/60:>10.2f} min")
    print("="*70)
    
    return {
        'neutral': ccsdt_neutral,
        'cation': ccsdt_cation,
        'ie_hartree': ie_hartree,
        'ie_ev': ie_ev,
        'total_time': total_time
    }

def main():
    """Main execution function"""
    results = calculate_ionization_potential()
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"✅ PySCF CCSD(T) calculations completed")
    print(f"\n   IE = {results['ie_ev']:.6f} eV")
    print(f"   Error = {results['ie_ev']-EXP_IE:+.4f} eV")
    print(f"\n   • Method: CCSD(T) with {BASIS_SET} basis")
    print(f"   • Neutral: RHF + CCSD(T)")
    print(f"   • Cation:  UHF + CCSD(T)")
    print(f"   • T1 diagnostic (neutral): {results['neutral']['t1_diagnostic']:.4f}")
    print(f"   • T1 diagnostic (cation):  {results['cation']['t1_diagnostic']:.4f}")
    print(f"   • Total time: {results['total_time']/60:.2f} minutes")
    print("="*70)

if __name__ == "__main__":
    main()
