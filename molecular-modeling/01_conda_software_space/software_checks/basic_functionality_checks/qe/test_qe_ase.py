"""
ASE Quantum ESPRESSO Testing Script
Uses local pseudopotential files (H.upf and Si.upf) in the same directory
"""

import os
import sys
import numpy as np
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
from ase.io import write
import matplotlib.pyplot as plt
import subprocess

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Set the QE command explicitly for your environment
os.environ['ASE_ESPRESSO_COMMAND'] = '/home/miroi/miniconda3/envs/mace_env/bin/pw.x'


def check_pseudopotentials():
    """
    Check if pseudopotential files exist in the script directory.
    """
    print("\nChecking pseudopotential files...")
    print("-" * 40)
    
    pseudo_files = ['Si.upf', 'H.upf']
    all_exist = True
    
    for filename in pseudo_files:
        filepath = os.path.join(SCRIPT_DIR, filename)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            print(f"✓ Found: {filename} ({size:,} bytes)")
        else:
            print(f"✗ Missing: {filename}")
            all_exist = False
    
    if not all_exist:
        print("\n" + "="*60)
        print("ERROR: Missing pseudopotential files")
        print("="*60)
        print(f"Please place the following files in: {SCRIPT_DIR}")
        print("  - Si.upf")
        print("  - H.upf")
        print("\nYou can download them from:")
        print("  https://pseudopotentials.quantum-espresso.org/")
        print("  Or use the SSSP library from Materials Cloud")
        return False
    
    return True


def setup_calculator(calculation='scf', kpts=(4, 4, 4)):
    """
    Set up the Espresso calculator using local pseudopotentials.
    
    Args:
        calculation: Type of calculation ('scf', 'relax', 'bands', 'nscf')
        kpts: K-point grid as tuple (nx, ny, nz)
    """
    if not check_pseudopotentials():
        return None
    
    # Create profile
    try:
        profile = EspressoProfile(
            command=os.environ['ASE_ESPRESSO_COMMAND'],
            pseudo_dir=SCRIPT_DIR
        )
    except Exception as e:
        print(f"Error creating EspressoProfile: {e}")
        return None
    
    # Define pseudopotentials with local filenames
    pseudopotentials = {
        'Si': 'Si.upf',
        'H': 'H.upf',
    }
    
    # Input parameters for pw.x
    input_data = {
        'control': {
            'calculation': calculation,
            'restart_mode': 'from_scratch',
            'prefix': 'qe_test',
            'tprnfor': True,
            'tstress': True,
            'outdir': os.path.join(SCRIPT_DIR, 'tmp/'),
            'verbosity': 'high',
        },
        'system': {
            'ecutwfc': 30.0,
            'ecutrho': 120.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.01,
            'nbnd': 8,
        },
        'electrons': {
            'diagonalization': 'david',
            'conv_thr': 1e-8,
            'mixing_beta': 0.7,
            'maxstep': 100,
        },
    }
    
    return Espresso(
        profile=profile,
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=kpts,
        koffset=(0, 0, 0),
    )


def test_single_point():
    """
    Test 1: Single point SCF calculation for Silicon.
    """
    print("\n" + "="*60)
    print("TEST 1: Single Point SCF Calculation")
    print("="*60)
    
    # Create silicon crystal
    si = bulk('Si', 'diamond', a=5.43)
    print(f"Structure: {len(si)} atoms, volume={si.get_volume():.2f} A^3")
    
    # Set up calculator with local pseudos
    calc = setup_calculator(calculation='scf')
    if calc is None:
        return False
    
    si.calc = calc
    
    try:
        # Run calculation
        print("\nRunning SCF calculation...")
        energy = si.get_potential_energy()
        forces = si.get_forces()
        stress = si.get_stress()
        
        print(f"\n✓ Calculation completed successfully!")
        print(f"Total energy: {energy:.6f} eV")
        print(f"Force norms: {np.linalg.norm(forces, axis=1)}")
        print(f"Stress tensor (eV/A^3): {stress}")
        return True
        
    except Exception as e:
        print(f"✗ Error in single point calculation: {e}")
        print("\nDebugging tips:")
        print("1. Check that Si.upf is in:", SCRIPT_DIR)
        print("2. Check QE output file 'qe_test.pwo' for details")
        
        # Try to read the output file for debugging
        if os.path.exists('qe_test.pwo'):
            print("\nLast 20 lines of QE output:")
            with open('qe_test.pwo', 'r') as f:
                lines = f.readlines()
                for line in lines[-20:]:
                    print(f"  {line.strip()}")
        return False


def test_geometry_optimization():
    """
    Test 2: Geometry optimization using LBFGS.
    """
    print("\n" + "="*60)
    print("TEST 2: Geometry Optimization")
    print("="*60)
    
    # Create silicon with slightly strained lattice
    si = bulk('Si', 'diamond', a=5.50)
    print(f"Initial: a=5.500 A, volume={si.get_volume():.2f} A^3")
    
    # Setup with relaxation parameters
    calc = setup_calculator(calculation='relax')
    if calc is None:
        return False
    
    si.calc = calc
    
    try:
        print("\nRunning geometry optimization...")
        opt = LBFGS(si)
        opt.run(fmax=0.01)
        
        final_energy = si.get_potential_energy()
        final_volume = si.get_volume()
        final_a = (4 * final_volume / len(si)) ** (1/3)
        
        print(f"\n✓ Optimization completed successfully!")
        print(f"Final: a={final_a:.3f} A, volume={final_volume:.2f} A^3")
        print(f"Final energy: {final_energy:.6f} eV")
        print(f"Number of optimization steps: {opt.nsteps}")
        return True
        
    except Exception as e:
        print(f"✗ Error in geometry optimization: {e}")
        return False


def test_kpoint_convergence():
    """
    Test 3: K-point convergence test.
    """
    print("\n" + "="*60)
    print("TEST 3: K-point Convergence Test")
    print("="*60)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Test different k-point grids
    k_grids = [(2,2,2), (4,4,4), (6,6,6), (8,8,8)]
    energies = []
    grid_labels = []
    
    print("\nTesting k-point grids...")
    print("-" * 40)
    
    for kpts in k_grids:
        print(f"Testing grid: {kpts[0]}x{kpts[1]}x{kpts[2]}")
        
        try:
            calc = setup_calculator(calculation='scf', kpts=kpts)
            if calc is None:
                return False
            
            si.calc = calc
            energy = si.get_potential_energy()
            energies.append(energy)
            grid_labels.append(f"{kpts[0]}x{kpts[1]}x{kpts[2]}")
            print(f"  Energy: {energy:.6f} eV")
            
        except Exception as e:
            print(f"  Error: {e}")
            energies.append(np.nan)
            grid_labels.append(f"{kpts[0]}x{kpts[1]}x{kpts[2]}")
    
    # Plot convergence if we have valid results
    valid = ~np.isnan(energies)
    if np.any(valid) and len(valid) > 1:
        plt.figure(figsize=(8, 5))
        valid_energies = [e for e in energies if not np.isnan(e)]
        valid_labels = [g for g, v in zip(grid_labels, valid) if v]
        
        plt.plot(range(len(valid_energies)), valid_energies, 'bo-', linewidth=2, markersize=8)
        plt.xticks(range(len(valid_labels)), valid_labels, rotation=45)
        plt.xlabel('K-point grid', fontsize=12)
        plt.ylabel('Total Energy (eV)', fontsize=12)
        plt.title('K-point Convergence Test for Si', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(SCRIPT_DIR, 'kpoint_convergence.png'), dpi=150)
        print(f"\n✓ K-point convergence plot saved to: {os.path.join(SCRIPT_DIR, 'kpoint_convergence.png')}")
        return True
    elif np.any(valid) and len(valid) == 1:
        print(f"\n✓ Single k-point test completed: {energies[0]:.6f} eV")
        return True
    else:
        print("\n✗ No valid k-point convergence data")
        return False


def test_cutoff_convergence():
    """
    Test 4: Energy cutoff convergence test.
    """
    print("\n" + "="*60)
    print("TEST 4: Energy Cutoff Convergence Test")
    print("="*60)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Test different energy cutoffs
    cutoffs = [15, 20, 25, 30, 35, 40]
    energies = []
    
    print("\nTesting energy cutoffs (Ry)...")
    print("-" * 40)
    
    for ecut in cutoffs:
        print(f"Testing ecutwfc: {ecut} Ry")
        
        try:
            calc = setup_calculator(calculation='scf')
            if calc is None:
                return False
            
            calc.input_data['system']['ecutwfc'] = ecut
            calc.input_data['system']['ecutrho'] = 4 * ecut
            
            si.calc = calc
            energy = si.get_potential_energy()
            energies.append(energy)
            print(f"  Energy: {energy:.6f} eV")
            
        except Exception as e:
            print(f"  Error: {e}")
            energies.append(np.nan)
    
    # Plot convergence if we have valid results
    valid = ~np.isnan(energies)
    if np.any(valid) and len(valid) > 1:
        plt.figure(figsize=(8, 5))
        valid_cutoffs = [c for c, v in zip(cutoffs, valid) if v]
        valid_energies = [e for e in valid if not np.isnan(e)]
        
        plt.plot(valid_cutoffs, valid_energies, 'ro-', linewidth=2, markersize=8)
        plt.xlabel('ecutwfc (Ry)', fontsize=12)
        plt.ylabel('Total Energy (eV)', fontsize=12)
        plt.title('Energy Cutoff Convergence Test for Si', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(SCRIPT_DIR, 'cutoff_convergence.png'), dpi=150)
        print(f"\n✓ Energy cutoff convergence plot saved to: {os.path.join(SCRIPT_DIR, 'cutoff_convergence.png')}")
        return True
    elif np.any(valid) and len(valid) == 1:
        print(f"\n✓ Single cutoff test completed: {energies[0]:.6f} eV")
        return True
    else:
        print("\n✗ No valid cutoff convergence data")
        return False


def main():
    """
    Run all tests.
    """
    print("\n" + "="*60)
    print("ASE Quantum ESPRESSO Testing Suite")
    print("="*60)
    print(f"\nScript directory: {SCRIPT_DIR}")
    print(f"QE executable: {os.environ['ASE_ESPRESSO_COMMAND']}")
    print("Using pseudopotentials: Si.upf and H.upf (in script directory)")
    
    # Verify QE executable exists
    if not os.path.exists(os.environ['ASE_ESPRESSO_COMMAND']):
        print(f"\nERROR: QE executable not found at: {os.environ['ASE_ESPRESSO_COMMAND']}")
        print("Please check your conda environment activation.")
        return
    
    # Create necessary directories
    os.makedirs(os.path.join(SCRIPT_DIR, 'tmp/'), exist_ok=True)
    
    # Check pseudopotentials
    if not check_pseudopotentials():
        print("\nPlease place Si.upf and H.upf in:", SCRIPT_DIR)
        print("You can download them from:")
        print("  https://pseudopotentials.quantum-espresso.org/")
        return
    
    # Run selected tests (start with basic ones first)
    tests = [
        ("Single Point SCF", test_single_point),
        ("Geometry Optimization", test_geometry_optimization),
        ("K-point Convergence", test_kpoint_convergence),
        ("Cutoff Convergence", test_cutoff_convergence),
    ]
    
    results = {}
    for name, test_func in tests:
        results[name] = test_func()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, status in results.items():
        print(f"{'✓' if status else '✗'} {name}")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All tests passed! Quantum ESPRESSO is working correctly with ASE.")
    else:
        print(f"\n⚠️ {total - passed} test(s) failed.")
        print("\nTroubleshooting tips:")
        print("1. Check that Si.upf and H.upf are valid pseudopotential files")
        print("2. Check the QE output files (*.pwo) for error messages")
        print("3. Verify you have sufficient disk space in ./tmp/")


if __name__ == "__main__":
    main()
