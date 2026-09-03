import numpy as np
import matplotlib.pyplot as plt

try:
    # Read Fermi level
    fermi_level = 0.0
    with open('electronic_properties_graphene.out', 'r') as f:
        for line in f:
            if 'Fermi level from NSCF calculation' in line:
                fermi_level = float(line.split()[-2])

    # Read Total DOS
    tdos_data = np.loadtxt('total_dos.dat')
    e_tdos = tdos_data[:, 0] - fermi_level
    tdos = tdos_data[:, 1]

    # Read p-orbital PDOS file
    p_data = np.loadtxt('pdos_results/pdos.pdos_atm#1(C)_wfc#2(p)')
    e_pdos = p_data[:, 0] - fermi_level
    
    # Column 1 (index 1) is the total LDOS for all p-orbitals (px + py + pz)
    p_total_dos = p_data[:, 1]

    plt.figure(figsize=(10, 6))
    
    # Plot Total DOS in background
    plt.fill_between(e_tdos, 0, tdos, color='grey', alpha=0.2, label='Total DOS')
    
    # Plot Total p-orbital
    plt.plot(e_pdos, p_total_dos, label='Total p orbital (px + py + pz)', color='purple', linewidth=2.5)

    plt.axvline(x=0, color='black', linestyle='--', label='Fermi Level (0 eV)')

    plt.xlim(-15, 10)
    plt.ylim(0, max(p_total_dos[ (e_pdos > -10) & (e_pdos < 10) ]) * 1.5)
    plt.xlabel('Energy (eV)', fontsize=14)
    plt.ylabel('Density of States', fontsize=14)
    plt.title('Graphene: Total p Orbital', fontsize=16)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(fontsize=12)

    plt.savefig('Graphene_Total_p_PDOS.png', dpi=300)
    print("Plot saved as Graphene_Total_p_PDOS.png")

except Exception as e:
    print(f"Error: {e}")
