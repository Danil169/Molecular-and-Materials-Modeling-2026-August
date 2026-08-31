import numpy as np
import matplotlib.pyplot as plt

# Read total DOS
try:
    data = np.loadtxt('total_dos.dat')
    energy = data[:, 0]
    dos = data[:, 1]
    
    # Try to find Fermi level from electronic_properties_al.out
    fermi_level = 0.0
    import os
    if os.path.exists('electronic_properties_al.out'):
        with open('electronic_properties_al.out', 'r') as f:
            for line in f:
                if 'Fermi level from NSCF calculation' in line:
                    fermi_level = float(line.split()[-2])
    
    # Shift energy so Fermi level is at 0
    energy_shifted = energy - fermi_level
    
    plt.figure(figsize=(10, 6))
    plt.plot(energy_shifted, dos, label='Al Total DOS', color='darkorange')
    plt.fill_between(energy_shifted, 0, dos, where=(energy_shifted <= 0), color='bisque', alpha=0.7)
    
    plt.axvline(x=0, color='red', linestyle='--', label=f'Уровень Ферми (0 эВ)')
    
    plt.xlim(-15, 10)
    plt.ylim(0, max(dos)*1.1)
    plt.xlabel('Энергия (ЭВ)', fontsize=14)
    plt.ylabel('Плотность состояний', fontsize=14)
    plt.title('Плотность состояний Алюминия (Металл)', fontsize=16)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(fontsize=12)
    
    plt.savefig('Al_DOS.png', dpi=300)
    print("Plot saved as Al_DOS.png")
except Exception as e:
    print(f"Error: {e}")
