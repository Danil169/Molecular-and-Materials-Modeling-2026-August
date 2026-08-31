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

    # Read PDOS files (we only need the p-orbitals)
    # The file is 'pdos.pdos_atm#1(C)_wfc#2(p)'
    p_data = np.loadtxt('pdos_results/pdos.pdos_atm#1(C)_wfc#2(p)')
    e_pdos = p_data[:, 0] - fermi_level
    
    # According to readme: columns 3-5 describe pz, px, py in that order
    # (In Python indices: col 2 -> pz, col 3 -> px, col 4 -> py)
    pz_dos = p_data[:, 2]
    px_dos = p_data[:, 3]
    py_dos = p_data[:, 4]

    plt.figure(figsize=(10, 6))
    
    # Plot Total DOS in background
    plt.fill_between(e_tdos, 0, tdos, color='grey', alpha=0.2, label='Total DOS')
    
    # Plot p-orbitals
    plt.plot(e_pdos, pz_dos, label='pz (out-of-plane)', color='blue', linewidth=2)
    plt.plot(e_pdos, px_dos, label='px (in-plane)', color='red', linewidth=1.5, linestyle='--')
    plt.plot(e_pdos, py_dos, label='py (in-plane)', color='green', linewidth=1.5, linestyle='--')

    plt.axvline(x=0, color='black', linestyle='--', label='Fermi Level (0 eV)')

    plt.xlim(-10, 10)
    plt.ylim(0, max(pz_dos[ (e_pdos > -10) & (e_pdos < 10) ]) * 1.5)
    plt.xlabel('Energy - E_F (eV)', fontsize=14)
    plt.ylabel('Density of States', fontsize=14)
    plt.title('Graphene PDOS (Dirac Cone Formation)', fontsize=16)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(fontsize=12)

    plt.savefig('Graphene_PDOS.png', dpi=300)
    print("Plot saved as Graphene_PDOS.png")

except Exception as e:
    print(f"Waiting for files to be generated... (Error: {e})")
