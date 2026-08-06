import numpy as np
import matplotlib.pyplot as plt
from ase.build import molecule
from ase.md.langevin import Langevin
from ase import units
from scipy.signal import find_peaks
from chgnet.model.dynamics import CHGNetCalculator

# 1. Setup the System & Calculator
atoms = molecule('H2O')
try:
    atoms.calc = CHGNetCalculator()
except Exception as e:
    print("Ensure chgnet is installed: pip install chgnet")
    raise e

# 2. Simulation Parameters
temp_k = 300
dt_fs = 0.5 
equil_steps = 1000  # Warm up
prod_steps = 10000   # Longer production for better resolution

# 3. Equilibration (Crucial for ML Force Fields)
print(f"Equilibrating at {temp_k} K...")
dyn = Langevin(atoms, dt_fs * units.fs, temperature_K=temp_k, friction=0.01)
dyn.run(equil_steps)

# 4. Production Run (Data Collection)
print(f"Production Run starting...")
traj_velocities = []

def collect_velocities():
    traj_velocities.append(atoms.get_velocities())

dyn.attach(collect_velocities, interval=1)
dyn.run(prod_steps)

# 5. Signal Processing
vels = np.array(traj_velocities)
# Zero-padding the signal can improve frequency resolution
n_fft = 2**14 
freq = np.fft.rfftfreq(n_fft, d=dt_fs)
fft_data = np.fft.rfft(vels, n=n_fft, axis=0)
power_spectrum = np.sum(np.abs(fft_data)**2, axis=(1, 2))

# Convert to cm^-1
freq_cm = freq * (1e15 / 29979245800.0)

# 6. Peak Extraction
# We search for peaks specifically in the bending (1200-1800) and stretching (3000+) regions
peaks, _ = find_peaks(power_spectrum, prominence=np.max(power_spectrum)*0.02)
extracted_freqs = sorted(freq_cm[peaks])

# Clean up extracted peaks (filtering out low-freq rotation/translation)
real_modes = [f for f in extracted_freqs if f > 500]

# 7. Comparison & Printing
experimental = {"Bending": 1595, "Sym-Stretch": 3657, "Asym-Stretch": 3756}
labels = list(experimental.keys())

print(f"\n{'='*60}")
print(f"RESULTS FOR CHGNet AT {temp_k} K")
print(f"{'='*60}")
print(f"{'Mode':<20} | {'CHGNet (cm^-1)':<15} | {'Expt (cm^-1)':<12} | {'Error (%)'}")
print(f"{'-'*60}")

# Match extracted peaks to experimental labels
for i, label in enumerate(labels):
    if i < len(real_modes):
        calc = real_modes[i]
        expt = experimental[label]
        error = abs(calc - expt) / expt * 100
        print(f"{label:<20} | {calc:15.1f} | {expt:12.1f} | {error:.2f}%")
    else:
        print(f"{label:<20} | {'Not Found':<15} | {expt:12.1f} | N/A")

# Plot
plt.plot(freq_cm, power_spectrum)
plt.xlim(500, 4500)
plt.title(f"CHGNet Anharmonic Spectrum ({temp_k} K)")
plt.xlabel("cm$^{-1}$")
plt.show()

