import numpy as np
import matplotlib.pyplot as plt

# Read total DOS
data = np.loadtxt('total_dos.dat')
energy = data[:, 0]
dos = data[:, 1]

# Fermi level from the output
fermi_level = 6.3215

# Shift energy so Fermi level is at 0
energy_shifted = energy - fermi_level

plt.figure(figsize=(10, 6))
plt.plot(energy_shifted, dos, label='Total DOS', color='blue')
plt.fill_between(energy_shifted, 0, dos, where=(energy_shifted <= 0), color='lightblue', alpha=0.5)

plt.axvline(x=0, color='red', linestyle='--', label='Fermi Level (0 eV)')

plt.xlim(-12, 10)
plt.ylim(0, max(dos)*1.1)
plt.xlabel('Energy - E_F (eV)', fontsize=14)
plt.ylabel('Density of States', fontsize=14)
plt.title('Silicon Density of States (DOS)', fontsize=16)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=12)

plt.savefig('Silicon_DOS.png', dpi=300)
print("Plot saved as Silicon_DOS.png")
