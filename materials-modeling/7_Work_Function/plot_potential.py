import numpy as np
import matplotlib.pyplot as plt

try:
    # 1. Get Fermi level from out file
    fermi_ev = 0.0
    with open('work_function.out', 'r') as f:
        for line in f:
            if 'Fermi level:' in line:
                fermi_ev = float(line.split()[2])

    # 2. Read data
    # Column 0: z in Bohr
    # Column 1: planar avg in Ry
    # Column 2: macroscopic avg in Ry
    data = np.loadtxt('potential_results/avg.dat')
    
    # Constants
    bohr_to_ang = 0.529177249
    ry_to_ev = 13.6056980659
    
    z_ang = data[:, 0] * bohr_to_ang
    planar_ev = data[:, 1] * ry_to_ev
    macro_ev = data[:, 2] * ry_to_ev

    # Calculate Work Function
    vacuum_level = macro_ev[0]  # The potential at z=0 (in the vacuum region)
    work_function = vacuum_level - fermi_ev

    # 3. Plot
    plt.figure(figsize=(12, 8))
    plt.plot(z_ang, planar_ev, color='red', label='Planar Average', linewidth=1.5, alpha=0.7)
    plt.plot(z_ang, macro_ev, color='blue', label='Macroscopic Average', linewidth=3)
    
    plt.axhline(y=fermi_ev, color='green', linestyle='--', linewidth=2, label=f'Fermi Level: {fermi_ev:.3f} eV')
    plt.axhline(y=vacuum_level, color='purple', linestyle=':', linewidth=2, label=f'Vacuum Level: {vacuum_level:.3f} eV')
    
    # Formatting
    plt.xlim(0, max(z_ang))
    plt.xlabel('z Coordinate (Å)', fontsize=14)
    plt.ylabel('Electrostatic Potential (eV)', fontsize=14)
    plt.title(f'Graphene Work Function = {work_function:.3f} eV', fontsize=16)
    plt.legend(fontsize=12, loc='upper right')
    plt.grid(True, linestyle=':', alpha=0.7)
    
    plt.savefig('potential_plot_python.png', dpi=300)
    print(f"Plot saved as potential_plot_python.png. Work Function: {work_function:.3f} eV")

except Exception as e:
    print(f"Error: {e}")
