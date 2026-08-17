"""
ASE Quantum ESPRESSO Testing Script - Fully Automated with MPI
Tests Silicon (diamond cubic structure) using Si.upf pseudopotential
Automatically detects available CPUs using /proc/cpuinfo and uses mpirun
"""

import os
import sys
import numpy as np
import multiprocessing
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
from ase.io import write
import matplotlib.pyplot as plt
import subprocess
import time
import re

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def get_optimal_cpu_count():
    """
    Detect the optimal number of CPU cores/threads to use.
    Uses /proc/cpuinfo for accurate detection on Linux systems.
    Returns the number of cores for MPI.
    """
    try:
        # Method 1: Use /proc/cpuinfo (most reliable on Linux)
        total_cpus = 0
        try:
            with open('/proc/cpuinfo', 'r') as f:
                content = f.read()
                # Count processor entries
                total_cpus = len(re.findall(r'^processor\s+:', content, re.MULTILINE))
        except (FileNotFoundError, IOError):
            # Fallback to multiprocessing if /proc/cpuinfo not available
            total_cpus = multiprocessing.cpu_count()
        
        # Check environment variables for HPC/scheduler settings
        if 'SLURM_CPUS_PER_TASK' in os.environ:
            total_cpus = min(total_cpus, int(os.environ['SLURM_CPUS_PER_TASK']))
        elif 'OMP_NUM_THREADS' in os.environ:
            total_cpus = min(total_cpus, int(os.environ['OMP_NUM_THREADS']))
        
        # Determine optimal number of cores to use
        # Leave some cores for system and other processes
        if total_cpus <= 2:
            n_cpus = 1
        elif total_cpus <= 4:
            n_cpus = 2
        elif total_cpus <= 8:
            n_cpus = 4
        elif total_cpus <= 16:
            n_cpus = 6
        else:
            n_cpus = 8
        
        return n_cpus, total_cpus
        
    except Exception as e:
        print(f"Warning: Could not detect CPU count: {e}")
        return 2, 2

def test_mpi_availability():
    """Check if MPI is available and working."""
    try:
        result = subprocess.run(['which', 'mpirun'], capture_output=True, text=True)
        if result.returncode != 0:
            return False
        
        test_result = subprocess.run(['mpirun', '--version'], 
                                   capture_output=True, text=True, timeout=5)
        return test_result.returncode == 0
        
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return False

def find_pw_x():
    """Find the pw.x executable in common locations."""
    try:
        result = subprocess.run(['which', 'pw.x'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    common_paths = [
        '/usr/local/bin/pw.x',
        '/opt/quantum_espresso/bin/pw.x',
        os.path.expanduser('~/anaconda3/bin/pw.x'),
        os.path.expanduser('~/miniconda3/bin/pw.x'),
        os.path.expanduser('~/anaconda3/envs/mace_env/bin/pw.x'),
        os.path.expanduser('~/miniconda3/envs/mace_env/bin/pw.x'),
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            return path
    
    return None

def read_qe_command():
    """Read QE command from file or auto-detect."""
    command_file = os.path.join(SCRIPT_DIR, 'qe_command.txt')
    
    if os.path.exists(command_file):
        with open(command_file, 'r') as f:
            command = f.read().strip()
            if command:
                print(f"✓ Using custom command from: {command_file}")
                print(f"  Command: {command}")
                return command
    
    print("🔍 Auto-detecting system configuration...")
    
    pw_path = find_pw_x()
    if not pw_path:
        print("  ⚠️ pw.x not found")
        return 'pw.x'
    
    print(f"  ✓ Found pw.x at: {pw_path}")
    
    mpi_available = test_mpi_availability()
    
    if mpi_available:
        n_cpus, total_cpus = get_optimal_cpu_count()
        print(f"  ✓ MPI available")
        print(f"  ✓ Total available CPUs: {total_cpus}")
        print(f"  ✓ Using {n_cpus} CPU cores for MPI")
        command = f'mpirun -np {n_cpus} {pw_path}'
        print(f"  ✓ Command: {command}")
        return command
    else:
        print("  ⚠️ MPI not available, using serial mode")
        return pw_path

# Set the QE command
qe_command = read_qe_command()
os.environ['ASE_ESPRESSO_COMMAND'] = qe_command


def check_pseudopotentials():
    """Check if Si.upf exists in the script directory."""
    print("\nChecking pseudopotential files...")
    print("-" * 40)
    
    filename = 'Si.upf'
    filepath = os.path.join(SCRIPT_DIR, filename)
    
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"✓ Found: {filename} ({size:,} bytes)")
        return True
    else:
        print(f"✗ Missing: {filename}")
        print("\n" + "="*60)
        print("ERROR: Missing pseudopotential file")
        print("="*60)
        print(f"Please place Si.upf in: {SCRIPT_DIR}")
        print("\nYou can download it from the Pseudo-Dojo website:")
        print("  https://www.pseudo-dojo.org/")
        print("  (Select Silicon and download the 'UPF' file)")
        return False


def setup_calculator(calculation='scf', kpts=(2, 2, 2), ecutwfc=20.0):
    """
    Set up the Espresso calculator for Silicon using Si.upf.
    """
    if not check_pseudopotentials():
        return None
    
    try:
        profile = EspressoProfile(
            command=os.environ['ASE_ESPRESSO_COMMAND'],
            pseudo_dir=SCRIPT_DIR
        )
    except Exception as e:
        print(f"Error creating EspressoProfile: {e}")
        return None
    
    # Only Si pseudopotential needed
    pseudopotentials = {
        'Si': 'Si.upf',
    }
    
    input_data = {
        'control': {
            'calculation': calculation,
            'restart_mode': 'from_scratch',
            'prefix': 'si_test',
            'tprnfor': True,
            'tstress': True,
            'outdir': os.path.join(SCRIPT_DIR, 'tmp/'),
            'verbosity': 'low',
        },
        'system': {
            'ecutwfc': ecutwfc,
            'ecutrho': 4.0 * ecutwfc,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'nbnd': 8,
        },
        'electrons': {
            'diagonalization': 'david',
            'conv_thr': 1e-6,
            'mixing_beta': 0.7,
            'electron_maxstep': 50,
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
    """Test 1: Single point SCF calculation for Silicon (diamond cubic)."""
    print("\n" + "="*60)
    print("TEST 1: Single Point SCF Calculation")
    print("="*60)
    print("Silicon crystal: Diamond cubic (Fd3̄m)")
    print("Lattice constant: a = 5.43 Å")
    print("Using ecutwfc=20 Ry, kpts=(2,2,2) for speed")
    
    # Create silicon in diamond cubic structure
    # ASE's bulk('Si', 'diamond', a=5.43) creates the correct structure
    si = bulk('Si', 'diamond', a=5.43)
    
    print(f"\nStructure information:")
    print(f"  Number of atoms: {len(si)} (primitive cell)")
    print(f"  Volume: {si.get_volume():.2f} Å³")
    print(f"  Space group: Fd3̄m (diamond cubic)")
    
    # Display atomic positions
    print("\n  Atomic positions (fractional coordinates):")
    for i, atom in enumerate(si):
        pos = atom.position / si.cell.cellpar()[0]  # Normalize by lattice constant
        print(f"    Si({i+1}): ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
    
    calc = setup_calculator(calculation='scf', kpts=(2, 2, 2), ecutwfc=20.0)
    if calc is None:
        return False
    
    si.calc = calc
    
    try:
        print("\nRunning SCF calculation...")
        start_time = time.time()
        energy = si.get_potential_energy()
        forces = si.get_forces()
        stress = si.get_stress()
        elapsed = time.time() - start_time
        
        print(f"\n✓ Calculation completed successfully!")
        print(f"  Time elapsed: {elapsed:.2f} seconds")
        print(f"  Total energy: {energy:.6f} eV")
        print(f"  Force norms: {np.linalg.norm(forces, axis=1)}")
        print(f"  Stress tensor (eV/Å³): {stress}")
        return True
        
    except Exception as e:
        print(f"✗ Error in single point calculation: {e}")
        if os.path.exists('si_test.pwo'):
            print("\nLast 30 lines of QE output:")
            with open('si_test.pwo', 'r') as f:
                lines = f.readlines()
                for line in lines[-30:]:
                    print(f"  {line.strip()}")
        return False


def test_geometry_optimization():
    """Test 2: Geometry optimization for Silicon (diamond cubic)."""
    print("\n" + "="*60)
    print("TEST 2: Geometry Optimization")
    print("="*60)
    print("Silicon crystal: Diamond cubic (Fd3̄m)")
    print("Initial lattice constant: a = 5.50 Å (strained)")
    print("Using ecutwfc=25 Ry, kpts=(2,2,2) for speed")
    
    # Create silicon with slightly strained lattice
    si = bulk('Si', 'diamond', a=5.50)
    print(f"\nInitial: a=5.500 Å, volume={si.get_volume():.2f} Å³")
    
    calc = setup_calculator(calculation='relax', kpts=(2, 2, 2), ecutwfc=25.0)
    if calc is None:
        return False
    
    si.calc = calc
    
    try:
        print("\nRunning geometry optimization...")
        opt = LBFGS(si)
        opt.run(fmax=0.05)
        
        final_energy = si.get_potential_energy()
        final_volume = si.get_volume()
        final_a = (4 * final_volume / len(si)) ** (1/3)  # For primitive cell
        
        print(f"\n✓ Optimization completed successfully!")
        print(f"  Final lattice constant: a={final_a:.3f} Å")
        print(f"  Final volume: {final_volume:.2f} Å³")
        print(f"  Final energy: {final_energy:.6f} eV")
        print(f"  Number of optimization steps: {opt.nsteps}")
        
        # Compare with experimental value
        exp_a = 5.431
        diff = abs(final_a - exp_a) / exp_a * 100
        print(f"\n  Experimental a = 5.431 Å")
        print(f"  Difference: {diff:.2f}%")
        return True
        
    except Exception as e:
        print(f"✗ Error in geometry optimization: {e}")
        return False


def test_kpoint_convergence():
    """Test 3: K-point convergence test for Silicon."""
    print("\n" + "="*60)
    print("TEST 3: K-point Convergence Test")
    print("="*60)
    print("Silicon crystal: Diamond cubic (Fd3̄m)")
    
    si = bulk('Si', 'diamond', a=5.43)
    
    k_grids = [(2,2,2), (4,4,4)]
    energies = []
    grid_labels = []
    
    print("\nTesting k-point grids...")
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
    
    valid = ~np.isnan(energies)
    if np.any(valid):
        print(f"\n✓ K-point test completed")
        for i, (label, energy) in enumerate(zip(grid_labels, energies)):
            if not np.isnan(energy):
                print(f"  {label}: {energy:.6f} eV")
        
        if len([e for e in energies if not np.isnan(e)]) > 1:
            try:
                plt.figure(figsize=(8, 5))
                valid_energies = [e for e in energies if not np.isnan(e)]
                valid_labels = [g for g, v in zip(grid_labels, valid) if v]
                plt.plot(range(len(valid_energies)), valid_energies, 'bo-', linewidth=2, markersize=8)
                plt.xticks(range(len(valid_labels)), valid_labels, rotation=45)
                plt.xlabel('K-point grid', fontsize=12)
                plt.ylabel('Total Energy (eV)', fontsize=12)
                plt.title('K-point Convergence Test for Si (Diamond Cubic)', fontsize=14)
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(SCRIPT_DIR, 'kpoint_convergence.png'), dpi=150)
                print(f"  ✓ Plot saved to: kpoint_convergence.png")
            except Exception as e:
                print(f"  Could not create plot: {e}")
        
        return True
    else:
        print("\n✗ No valid k-point data")
        return False


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("ASE Quantum ESPRESSO Testing Suite")
    print("="*60)
    print("\n📐 Silicon Crystal Structure:")
    print("   - Structure: Diamond cubic")
    print("   - Space group: Fd3̄m (No. 227)")
    print("   - Lattice constant: a = 5.431 Å (experimental)")
    print("   - Atoms per primitive cell: 2")
    print("   - Atomic positions: (0,0,0) and (1/4,1/4,1/4)")
    print("   - Coordination: 4-fold tetrahedral")
    
    print(f"\nScript directory: {SCRIPT_DIR}")
    
    print("\n📊 System Configuration:")
    print("-" * 40)
    print(f"  QE command: {os.environ.get('ASE_ESPRESSO_COMMAND', 'Not set')}")
    
    # Show CPU information
    try:
        result = subprocess.run(['cat', '/proc/cpuinfo'], capture_output=True, text=True)
        if result.returncode == 0:
            processor_lines = [line for line in result.stdout.split('\n') 
                             if line.startswith('processor')]
            total_cpus = len(processor_lines)
            print(f"  Total available CPUs: {total_cpus}")
        else:
            total_cpus = multiprocessing.cpu_count()
            print(f"  Total available CPUs: {total_cpus}")
    except:
        try:
            total_cpus = multiprocessing.cpu_count()
            print(f"  Total available CPUs: {total_cpus}")
        except:
            print("  Total available CPUs: Unknown")
    
    qe_cmd = os.environ.get('ASE_ESPRESSO_COMMAND', '')
    if 'mpirun' in qe_cmd:
        match = re.search(r'-np\s+(\d+)', qe_cmd)
        if match:
            n_procs = int(match.group(1))
            print(f"  MPI processes: {n_procs}")
    
    print("  Pseudopotential: Si.upf (from Pseudo-Dojo)")
    
    command_file = os.path.join(SCRIPT_DIR, 'qe_command.txt')
    if os.path.exists(command_file):
        print(f"  Using custom command from: qe_command.txt")
    else:
        print("  Using auto-detected MPI configuration")
    
    print("\n⚠️  OPTIMIZED FOR SPEED: Using low cutoffs (20-25 Ry) and coarse k-points (2x2x2)")
    print("   These parameters are NOT suitable for production calculations!")
    print("   For accurate results, increase ecutwfc to 30-40 Ry and kpts to (4,4,4) or higher")
    
    # Verify QE executable
    qe_cmd = os.environ.get('ASE_ESPRESSO_COMMAND', 'pw.x')
    qe_base = qe_cmd.split()[-1] if 'mpirun' in qe_cmd else qe_cmd
    
    try:
        subprocess.run(['which', qe_base], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        if not os.path.exists(qe_base):
            print(f"\n❌ ERROR: '{qe_base}' not found")
            print("Please ensure Quantum ESPRESSO is installed and in your PATH")
            print("\nYou can create a 'qe_command.txt' file with the correct command:")
            print("  echo 'pw.x' > qe_command.txt")
            print("  echo 'mpirun -np 2 pw.x' > qe_command.txt")
            return
    
    os.makedirs(os.path.join(SCRIPT_DIR, 'tmp/'), exist_ok=True)
    
    if not check_pseudopotentials():
        return
    
    print("\n🚀 Running tests...\n")
    tests = [
        ("Single Point SCF", test_single_point),
        ("Geometry Optimization", test_geometry_optimization),
        ("K-point Convergence", test_kpoint_convergence),
    ]
    
    results = {}
    for name, test_func in tests:
        results[name] = test_func()
    
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, status in results.items():
        print(f"{'✅' if status else '❌'} {name}")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All tests passed!")
        print("\n📝 Next steps for production calculations:")
        print("   - Increase ecutwfc to 30-40 Ry")
        print("   - Increase kpts to (4,4,4) or higher")
        print("   - Use conv_thr = 1e-8")
        print("   - Use smaller smearing (degauss = 0.01)")
        print("\n💻 Current configuration:")
        print(f"   {os.environ.get('ASE_ESPRESSO_COMMAND', 'Not set')}")
    else:
        print(f"\n⚠️ {total - passed} test(s) failed.")
        print("\n🔧 Troubleshooting suggestions:")
        print("  1. Create a qe_command.txt file with:")
        print("     echo 'pw.x' > qe_command.txt")
        print("     echo 'mpirun -np 2 pw.x' > qe_command.txt")
        print("  2. Check if Quantum ESPRESSO is installed:")
        print("     which pw.x")
        print("  3. Check if MPI is working:")
        print("     mpirun --version")
        print("  4. Verify Si.upf is a valid pseudopotential file")


if __name__ == "__main__":
    main()
