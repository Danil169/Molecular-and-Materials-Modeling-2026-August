import numpy as np
import matplotlib.pyplot as plt
from ase.build import molecule
from ase.md.langevin import Langevin
from ase import units
from scipy.signal import find_peaks
from chgnet.model.dynamics import CHGNetCalculator

# 1. Setup System with a Unit Cell (The Fix)
atoms = molecule('H2O')
atoms.center(vacuum=5.0)  # Adds 5 Angstroms of empty space around the molecule
atoms.pbc = True          # CHGNet requires Periodic Boundary Conditions

# 2. Initialize Calculator
try:
    atoms.calc = CHGNetCalculator()
except Exception as e:
    print("Ensure chgnet is installed: pip install chgnet")
    raise e

# 3. Parameters
temp_k = 300
dt_fs = 0.5
equil_steps = 1000
prod_steps = 10000

# 4. Equilibration
print(f"--- Equilibration at {temp_k} K ---")
dyn = Langevin(atoms, dt_fs * units.fs, temperature_K=temp_k, friction=0.01)
dyn.run(equil_steps)

# 5. Production
print(f"--- Production Run ---")
traj_velocities = []
def collect_velocities():
    traj_velocities.append(atoms.get_velocities())

dyn.attach(collect_velocities, interval=1)
dyn.run(prod_steps)

# 6. Processing and Comparison
vels = np.array(traj_velocities)
n_fft = 2**15 
freq = np.fft.rfftfreq(n_fft, d=dt_fs)
fft_data = np.fft.rfft(vels, n=n_fft, axis=0)
power_spectrum = np.sum(np.abs(fft_data)**2, axis=(1, 2))

# Convert to cm^-1
freq_cm = freq * (1e15 / 29979245800.0)

# Extract and Compare
mask = freq_cm > 500
peaks, _ = find_peaks(power_spectrum[mask], prominence=np.max(power_spectrum)*0.05)
found_freqs = sorted(freq_cm[mask][peaks])

experimental = {"Bending": 1595, "Sym-Stretch": 3657, "Asym-Stretch": 3756}

print(f"\nResults at Temperature: {temp_k} K")
print(f"{'-'*50}")
for label, expt in experimental.items():
    if found_freqs:
        # Match peak closest to experimental value
        calc = min(found_freqs, key=lambda x: abs(x - expt))
        error = abs(calc - expt) / expt * 100
        print(f"{label:<15}: Calc {calc:>7.1f} | Expt {expt:>7.1f} | Error {error:>5.2f}%")

plt.plot(freq_cm, power_spectrum)
plt.xlim(500, 4500)
plt.xlabel("cm$^{-1}$")
plt.title(f"H2O Power Spectrum ({temp_k} K)")
plt.show()

