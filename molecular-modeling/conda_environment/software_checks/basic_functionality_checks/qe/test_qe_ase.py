"""
ASE Quantum ESPRESSO Testing Script - Optimized for Speed
Uses local pseudopotential files (H.upf and Si.upf) 
Optimized cutoff parameters for fast testing
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

# Read QE command from a text file
def read_qe_command():
    """
    Read the Quantum ESPRESSO command from a text file.
    The file should be named 'qe_command.txt' in the same directory.
    """
    command_file = os.path.join(SCRIPT_DIR, 'qe_command.txt')
    
    # Default command if file doesn't exist
    default_command = 'pw.x'
    
    try:
        if os.path.exists(command_file):
            with open(command_file, 'r') as f:
                command = f.read().strip()
                if command:
                    print(f"✓ Read QE command from: {command_file}")
                    print(f"  Command: {command}")
                    return command
                else:
                    print(f"⚠️ {command_file} is empty. Using default: {default_command}")
                    return default_command
        else:
            print(f"⚠️ {command_file} not found. Using default: {default_command}")
            print(f"  To specify a custom command, create {command_file}")
            print(f"  with the command path (e.g., /path/to/pw.x)")
            return default_command
    except Exception as e:
        print(f"⚠️ Error reading {command_file}: {e}")
        print(f"  Using default: {default_command}")
        return default_command

# Set the QE command from the text file
qe_command = read_qe_command()
os.environ['ASE_ESPRESSO_COMMAND'] = qe_command


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
        return False
    
    return True


def setup_calculator(calculation='scf', kpts=(2, 2, 2), ecutwfc=20.0):
    """
    Set up the Espresso calculator using local pseudopotentials.
    OPTIMIZED FOR SPEED: Lower cutoffs and fewer k-points.
    
    Args:
        calculation: Type of calculation ('scf', 'relax', 'bands', 'nscf')
        kpts: K-point grid as tuple (nx, ny, nx)
        ecutwfc: Wavefunction cutoff in Ry (optimized: 20 Ry for testing)
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
    
    # Input parameters for pw.x - OPTIMIZED FOR SPEED
    input_data = {
        'control': {
            'calculation': calculation,
            'restart_mode': 'from_scratch',
            'prefix': 'qe_test',
            'tprnfor': True,
            'tstress': True,
            'outdir': os.path.join(SCRIPT_DIR, 'tmp/'),
            'verbosity': 'low',  # Reduce output verbosity for speed
        },
        'system': {
            'ecutwfc': ecutwfc,  # Lowered to 20 Ry for testing
            'ecutrho': 4.0 * ecutwfc,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,  # Slightly larger smearing for faster convergence
            'nbnd': 8,
        },
        'electrons': {
            'diagonalization': 'david',
            'conv_thr': 1e-6,  # Less strict convergence for testing
            'mixing_beta': 0.7,
            'electron_maxstep': 50,  # Reduced from 100 for speed
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
    OPTIMIZED: Low cutoffs and minimal k-points.
    """
    print("\n" + "="*60)
    print("TEST 1: Single Point SCF Calculation (Quick Test)")
    print("="*60)
    print("Using ecutwfc=20 Ry, kpts=(2,2,2) for speed")
    
    # Create silicon crystal
    si = bulk('Si', 'diamond', a=5.43)
    print(f"Structure: {len(si)} atoms, volume={si.get_volume():.2f} A^3")
    
    # Set up calculator with optimized parameters
    calc = setup_calculator(calculation='scf', kpts=(2, 2, 2), ecutwfc=20.0)
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
            print("\nLast 30 lines of QE output:")
            with open('qe_test.pwo', 'r') as f:
                lines = f.readlines()
                for line in lines[-30:]:
                    print(f"  {line.strip()}")
        return False


def test_geometry_optimization():
    """
    Test 2: Geometry optimization using LBFGS.
    OPTIMIZED: Fast convergence with relaxed parameters.
    """
    print("\n" + "="*60)
    print("TEST 2: Geometry Optimization (Quick Test)")
    print("="*60)
    print("Using ecutwfc=25 Ry, kpts=(2,2,2) for speed")
    
    # Create silicon with slightly strained lattice
    si = bulk('Si', 'diamond', a=5.50)
    print(f"Initial: a=5.500 A, volume={si.get_volume():.2f} A^3")
    
    # Setup with relaxation parameters
    calc = setup_calculator(calculation='relax', kpts=(2, 2, 2), ecutwfc=25.0)
    if calc is None:
        return False
    
    si.calc = calc
    
    try:
        print("\nRunning geometry optimization...")
        opt = LBFGS(si)
        opt.run(fmax=0.05)  # Looser convergence for speed
        
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


def test_quick_convergence():
    """
    Test 3: Very quick convergence test with minimal parameters.
    """
    print("\n" + "="*60)
    print("TEST 3: Quick Convergence Test")
    print("="*60)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Test just two k-point grids for speed
    k_grids = [(2,2,2), (4,4,4)]
    energies = []
    grid_labels = []
    
    print("\nTesting k-point grids (minimal)...")
    print("-" * 40)
    
    for kpts in k_grids:
        print(f"Testing grid: {kpts[0]}x{kpts[1]}x{kpts[2]}")
        
        try:
            calc = setup_calculator(calculation='scf', kpts=kpts, ecutwfc=20.0)
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
    
    # Simple summary
    valid = ~np.isnan(energies)
    if np.any(valid):
        print(f"\n✓ K-point test completed")
        for i, (label, energy) in enumerate(zip(grid_labels, energies)):
            if not np.isnan(energy):
                print(f"  {label}: {energy:.6f} eV")
        return True
    else:
        print("\n✗ No valid k-point data")
        return False


def main():
    """
    Run all tests.
    """
    print("\n" + "="*60)
    print("ASE Quantum ESPRESSO Testing Suite (Optimized for Speed)")
    print("="*60)
    print(f"\nScript directory: {SCRIPT_DIR}")
    
    # Show QE command being used
    qe_cmd = os.environ.get('ASE_ESPRESSO_COMMAND', 'Not set')
    print(f"QE executable: {qe_cmd}")
    
    # Show where the command came from
    command_file = os.path.join(SCRIPT_DIR, 'qe_command.txt')
    if os.path.exists(command_file):
        print(f"QE command loaded from: {command_file}")
    else:
        print(f"QE command using default: pw.x (no qe_command.txt found)")
    
    print("\n⚠️  OPTIMIZED FOR SPEED: Using low cutoffs (20-25 Ry) and coarse k-points (2x2x2)")
    print("   These parameters are NOT suitable for production calculations!")
    print("   For accurate results, increase ecutwfc to 30-40 Ry and kpts to 4x4x4 or higher")
    
    # Verify QE executable exists
    if not os.path.exists(qe_cmd) and qe_cmd != 'pw.x':
        print(f"\nERROR: QE executable not found at: {qe_cmd}")
        print("Please check your qe_command.txt file.")
        return
    
    # Create necessary directories
    os.makedirs(os.path.join(SCRIPT_DIR, 'tmp/'), exist_ok=True)
    
    # Check pseudopotentials
    if not check_pseudopotentials():
        print("\nPlease place Si.upf and H.upf in:", SCRIPT_DIR)
        return
    
    # Run tests (only quick ones)
    tests = [
        ("Single Point SCF (Quick)", test_single_point),
        ("Geometry Optimization (Quick)", test_geometry_optimization),
        ("Quick Convergence Test", test_quick_convergence),
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
        print("\n🎉 All tests passed!")
        print("\n📝 Next steps for production calculations:")
        print("   - Increase ecutwfc to 30-40 Ry")
        print("   - Increase kpts to (4,4,4) or higher")
        print("   - Use conv_thr = 1e-8")
        print("   - Use smaller smearing (degauss = 0.01)")
    else:
        print(f"\n⚠️ {total - passed} test(s) failed.")


if __name__ == "__main__":
    main()
