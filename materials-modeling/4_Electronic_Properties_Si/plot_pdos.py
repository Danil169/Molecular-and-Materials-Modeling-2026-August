import numpy as np
import matplotlib.pyplot as plt

# Fermi level from Task 4
fermi_level = 6.3215

# Load Total DOS for background reference
tdos_data = np.loadtxt('total_dos.dat')
e_tdos = tdos_data[:, 0] - fermi_level
tdos = tdos_data[:, 1]

# Load PDOS files
s_data = np.loadtxt('pdos_results/pdos.pdos_atm#1(Si)_wfc#1(s)')
p_data = np.loadtxt('pdos_results/pdos.pdos_atm#1(Si)_wfc#2(p)')

e_pdos = s_data[:, 0] - fermi_level
s_dos = s_data[:, 1]
p_dos = p_data[:, 1]

plt.figure(figsize=(10, 6))

# Plot Total DOS as light grey background
plt.fill_between(e_tdos, 0, tdos, color='grey', alpha=0.2, label='Total DOS')

# Plot s and p orbitals
plt.plot(e_pdos, s_dos, label='s orbital', color='red', linewidth=2)
plt.plot(e_pdos, p_dos, label='p orbital', color='blue', linewidth=2)

plt.axvline(x=0, color='black', linestyle='--', label='Fermi Level (0 eV)')

plt.xlim(-15, 10)
plt.ylim(0, max(tdos)*1.1)
plt.xlabel('Energy (eV)', fontsize=14)
plt.ylabel('Density of States (PDOS)', fontsize=14)
plt.title('Silicon Projected Density of States (PDOS)', fontsize=16)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=12)

plt.savefig('Silicon_PDOS.png', dpi=300)
print("Saved Silicon_PDOS.png")
