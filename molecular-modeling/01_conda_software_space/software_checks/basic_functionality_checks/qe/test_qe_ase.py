"""
ASE Quantum ESPRESSO Testing Script
Comprehensive test suite for QE calculations via ASE
"""

import os
import numpy as np
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
from ase.io import write, read
import matplotlib.pyplot as plt


def setup_calculator(pseudo_dir="./pseudos", command="pw.x"):
    """
    Set up the Espresso calculator with basic parameters.
    
    Args:
        pseudo_dir: Directory containing pseudopotential files
        command: Command to run pw.x (can include mpirun)
    
    Returns:
        Espresso calculator instance
    """
    # Create profile for the calculator
    # If using MPI: command = "mpirun -np 4 pw.x"
    profile = EspressoProfile(
        command=command,
        pseudo_dir=pseudo_dir
    )
    
    # Define pseudopotentials (SSSP Efficiency v1.3.0 recommended)
    # Download from: https://www.materialscloud.org/discover/sssp/table/efficiency
    pseudopotentials = {
        'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
        'H': 'H.pbe-n-rrkjus_psl.1.0.0.UPF',
    }
    
    # Input parameters for pw.x
    input_data = {
        'control': {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'prefix': 'si_test',
            'tprnfor': True,      # Print forces
            'tstress': True,      # Print stress
            'outdir': './tmp/',
        },
        'system': {
            'ecutwfc': 30.0,      # Wavefunction cutoff in Ry
            'ecutrho': 120.0,     # Density cutoff in Ry (4x ecutwfc typical)
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.01,
            'nbnd': 8,
        },
        'electrons': {
            'diagonalization': 'david',
            'conv_thr': 1e-8,
            'mixing_beta': 0.7,
        },
    }
    
    # Create calculator
    calc = Espresso(
        profile=profile,
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=(4, 4, 4),          # Monkhorst-Pack grid
        koffset=(0, 0, 0),       # No offset (gamma-centered grid)
    )
    
    return calc


def test_single_point():
    """
    Test 1: Single point SCF calculation for Silicon.
    This tests the basic functionality of the QE calculator.
    """
    print("\n" + "="*50)
    print("TEST 1: Single Point SCF Calculation")
    print("="*50)
    
    # Create silicon crystal
    si = bulk('Si', 'diamond', a=5.43)
    print(f"Initial structure: {len(si)} atoms, volume={si.get_volume():.2f} A^3")
    
    # Set up calculator
    calc = setup_calculator()
    si.calc = calc
    
    try:
        # Run calculation
        energy = si.get_potential_energy()
        forces = si.get_forces()
        stress = si.get_stress()
        
        print(f"Total energy: {energy:.6f} eV")
        print(f"Force norms: {np.linalg.norm(forces, axis=1)}")
        print(f"Stress tensor (eV/A^3): {stress}")
        print("✓ Single point calculation successful")
        return True
        
    except Exception as e:
        print(f"✗ Error in single point calculation: {e}")
        return False


def test_geometry_optimization():
    """
    Test 2: Geometry optimization using LBFGS.
    Tests the relaxation capabilities of the calculator.
    """
    print("\n" + "="*50)
    print("TEST 2: Geometry Optimization")
    print("="*50)
    
    # Create silicon with slightly strained lattice
    si = bulk('Si', 'diamond', a=5.50)  # 1.3% larger than equilibrium
    print(f"Initial: a={5.50:.3f} A, volume={si.get_volume():.2f} A^3")
    
    # Setup with relaxation parameters
    calc = setup_calculator()
    # Override calculation type for relaxation
    calc.input_data['control']['calculation'] = 'relax'
    si.calc = calc
    
    try:
        # Run optimization
        opt = LBFGS(si)
        opt.run(fmax=0.01)  # Converge forces below 0.01 eV/A
        
        final_energy = si.get_potential_energy()
        final_volume = si.get_volume()
        final_a = (4 * final_volume / len(si)) ** (1/3)
        
        print(f"Final: a={final_a:.3f} A, volume={final_volume:.2f} A^3")
        print(f"Final energy: {final_energy:.6f} eV")
        print(f"Number of optimization steps: {opt.nsteps}")
        print("✓ Geometry optimization successful")
        return True
        
    except Exception as e:
        print(f"✗ Error in geometry optimization: {e}")
        return False


def test_band_structure_setup():
    """
    Test 3: Set up band structure calculation.
    Tests the ability to create input for band structure.
    """
    print("\n" + "="*50)
    print("TEST 3: Band Structure Setup")
    print("="*50)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Define k-point path for band structure
    # High symmetry points for diamond cubic Si
    from ase.dft.kpoints import bandpath
    
    # Band path in reduced coordinates
    kpts = bandpath(
        path='GXWKGLUWLK',  # Standard path for FCC
        cell=si.cell,
        npoints=100
    )
    
    # Setup calculator for band structure
    calc = setup_calculator()
    calc.input_data['control']['calculation'] = 'bands'
    calc.input_data['system']['occupations'] = 'fixed'
    # Keep kpts as path for band structure
    calc.kpts = kpts
    si.calc = calc
    
    try:
        # Write input file only (don't run)
        from ase.io import write
        write('band_structure_input.in', si, format='espresso-in')
        print(f"Band structure input written with {len(kpts)} k-points")
        print("✓ Band structure setup successful")
        return True
        
    except Exception as e:
        print(f"✗ Error in band structure setup: {e}")
        return False


def test_dos_setup():
    """
    Test 4: Set up density of states (DOS) calculation.
    Tests the ability to create input for DOS.
    """
    print("\n" + "="*50)
    print("TEST 4: Density of States Setup")
    print("="*50)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Setup for DOS with dense k-point grid
    calc = setup_calculator()
    calc.input_data['control']['calculation'] = 'nscf'
    calc.kpts = (8, 8, 8)  # Denser grid for DOS
    
    # Add DOS-specific parameters
    calc.input_data['system']['occupations'] = 'tetrahedra'
    si.calc = calc
    
    try:
        # Write input file
        from ase.io import write
        write('dos_input.in', si, format='espresso-in')
        print("DOS input written with 8x8x8 k-point grid")
        print("✓ DOS setup successful")
        return True
        
    except Exception as e:
        print(f"✗ Error in DOS setup: {e}")
        return False


def test_kpoint_convergence():
    """
    Test 5: K-point convergence test.
    Demonstrates how to systematically test convergence.
    """
    print("\n" + "="*50)
    print("TEST 5: K-point Convergence Test")
    print("="*50)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Test different k-point grids
    k_grids = [(2,2,2), (4,4,4), (6,6,6), (8,8,8)]
    energies = []
    grid_labels = []
    
    print("Testing k-point grids...")
    print("-" * 40)
    
    for kpts in k_grids:
        print(f"Testing grid: {kpts[0]}x{kpts[1]}x{kpts[2]}")
        
        try:
            # Create calculator for this grid
            calc = setup_calculator()
            calc.kpts = kpts
            si.calc = calc
            
            # Run SCF
            energy = si.get_potential_energy()
            energies.append(energy)
            grid_labels.append(f"{kpts[0]}x{kpts[1]}x{kpts[2]}")
            print(f"  Energy: {energy:.6f} eV")
            
        except Exception as e:
            print(f"  Error: {e}")
            energies.append(np.nan)
            grid_labels.append(f"{kpts[0]}x{kpts[1]}x{kpts[2]}")
    
    # Plot convergence
    valid = ~np.isnan(energies)
    if np.any(valid):
        plt.figure(figsize=(8, 5))
        plt.plot(range(sum(valid)), [e for e in energies if not np.isnan(e)], 'bo-')
        plt.xticks(range(sum(valid)), [g for g, v in zip(grid_labels, valid) if v], rotation=45)
        plt.xlabel('K-point grid')
        plt.ylabel('Total Energy (eV)')
        plt.title('K-point Convergence Test for Si')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('kpoint_convergence.png', dpi=150)
        print("✓ K-point convergence plot saved to kpoint_convergence.png")
        print("✓ K-point convergence test complete")
        return True
    
    return False


def test_cutoff_convergence():
    """
    Test 6: Energy cutoff convergence test.
    Demonstrates how to test convergence with respect to ecutwfc.
    """
    print("\n" + "="*50)
    print("TEST 6: Energy Cutoff Convergence Test")
    print("="*50)
    
    si = bulk('Si', 'diamond', a=5.43)
    
    # Test different energy cutoffs
    cutoffs = [15, 20, 25, 30, 35, 40, 50]
    energies = []
    
    print("Testing energy cutoffs (Ry)...")
    print("-" * 40)
    
    for ecut in cutoffs:
        print(f"Testing ecutwfc: {ecut} Ry")
        
        try:
            calc = setup_calculator()
            calc.input_data['system']['ecutwfc'] = ecut
            calc.input_data['system']['ecutrho'] = 4 * ecut  # 4x is typical
            
            si.calc = calc
            energy = si.get_potential_energy()
            energies.append(energy)
            print(f"  Energy: {energy:.6f} eV")
            
        except Exception as e:
            print(f"  Error: {e}")
            energies.append(np.nan)
    
    # Plot convergence
    valid = ~np.isnan(energies)
    if np.any(valid):
        plt.figure(figsize=(8, 5))
        plt.plot([c for c, v in zip(cutoffs, valid) if v], 
                 [e for e in valid if not np.isnan(e)], 'ro-')
        plt.xlabel('ecutwfc (Ry)')
        plt.ylabel('Total Energy (eV)')
        plt.title('Energy Cutoff Convergence Test for Si')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('cutoff_convergence.png', dpi=150)
        print("✓ Energy cutoff convergence plot saved to cutoff_convergence.png")
        print("✓ Cutoff convergence test complete")
        return True
    
    return False


def test_read_results():
    """
    Test 7: Read results from existing QE output.
    Tests the ability to parse QE output files.
    """
    print("\n" + "="*50)
    print("TEST 7: Read Results from Output")
    print("="*50)
    
    # First run a calculation to generate output
    si = bulk('Si', 'diamond', a=5.43)
    calc = setup_calculator()
    si.calc = calc
    
    try:
        # Run calculation
        si.get_potential_energy()
        
        # Read results from output file
        from ase.io import read
        atoms = read('espresso.pwo', format='espresso-out')
        
        print(f"Read {len(atoms)} atoms")
        print(f"Energy from read results: {atoms.get_potential_energy():.6f} eV")
        
        # Check if forces were parsed
        if hasattr(atoms.calc, 'results') and 'forces' in atoms.calc.results:
            forces = atoms.calc.results['forces']
            print(f"Force norms from read results: {np.linalg.norm(forces, axis=1)}")
        
        print("✓ Reading results successful")
        return True
        
    except Exception as e:
        print(f"✗ Error reading results: {e}")
        return False


def main():
    """
    Run all tests.
    """
    print("\n" + "="*50)
    print("ASE Quantum ESPRESSO Testing Suite")
    print("="*50)
    print("\nNote: Make sure pseudopotentials are available in ./pseudos/")
    print("Set environment variable: export ASE_ESPRESSO_COMMAND='pw.x'")
    print("or modify the 'command' parameter in setup_calculator()")
    print()
    
    # Check environment
    if 'ASE_ESPRESSO_COMMAND' not in os.environ:
        print("Warning: ASE_ESPRESSO_COMMAND not set.")
        print("You may need to export ASE_ESPRESSO_COMMAND='pw.x'")
        print("Or modify the command parameter in setup_calculator()")
    
    # Run tests
    tests = [
        ("Single Point SCF", test_single_point),
        ("Geometry Optimization", test_geometry_optimization),
        ("Band Structure Setup", test_band_structure_setup),
        ("DOS Setup", test_dos_setup),
        ("K-point Convergence", test_kpoint_convergence),
        ("Cutoff Convergence", test_cutoff_convergence),
        ("Read Results", test_read_results),
    ]
    
    results = {}
    for name, test_func in tests:
        results[name] = test_func()
    
    # Summary
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, status in results.items():
        print(f"{'✓' if status else '✗'} {name}")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All tests passed!")
    else:
        print(f"\n⚠️ {total - passed} test(s) failed. Check the output above for details.")


if __name__ == "__main__":
    main()
