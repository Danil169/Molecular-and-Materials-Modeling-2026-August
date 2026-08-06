from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.calculators.emt import EMT
from ase.optimize import BFGS

# --- 1. SET UP CALCULATIONS ---
# Isolated Oxygen atom (Reference)
adsorbate = Atoms('O')
adsorbate.calc = EMT()
e_o = adsorbate.get_potential_energy()

# Clean Gold surface (Reference)
slab = fcc111('Au', size=(2, 2, 3), vacuum=10.0)
slab.calc = EMT()
e_slab = slab.get_potential_energy()

# Oxygen on Gold (Adsorption System)
# We use 'fcc' as it is the most stable site in literature
add_adsorbate(slab, 'O', height=1.5, position='fcc')

# --- 2. OPTIMIZE STRUCTURE ---
# Let the Oxygen atom find its optimal height
dyn = BFGS(slab, trajectory='adsorption.traj', logfile=None)
dyn.run(fmax=0.05)

# --- 3. RESULTS & COMPARISON ---
e_total = slab.get_potential_energy()
e_bind_calc = e_total - (e_slab + e_o)

# Published data for O/Au(111) at fcc site:
# Reference: J. Phys. Chem. C 2009, 113, 2, 635–644
lit_val_rpbe = -2.69  # RPBE functional
lit_val_pw91 = -3.08  # PW91 functional

print(f"--- ASE ADSORPTION REPORT ---")
print(f"Calculated Binding Energy (EMT): {e_bind_calc:8.3f} eV")
print(f"Literature Value (DFT-RPBE):     {lit_val_rpbe:8.3f} eV")
print(f"Literature Value (DFT-PW91):     {lit_val_pw91:8.3f} eV")
print(f"Difference (vs RPBE):           {abs(e_bind_calc - lit_val_rpbe):8.3f} eV")
print("-" * 29)

if e_bind_calc < 0:
    print("Status: Adsorption is EXOTHERMIC (stable).")
else:
    print("Status: Adsorption is ENDOTHERMIC (unstable).")

