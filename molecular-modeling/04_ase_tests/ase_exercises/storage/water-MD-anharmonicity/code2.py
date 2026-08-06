import numpy as np
import matplotlib.pyplot as plt
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.md.langevin import Langevin
from ase import units
from scipy.signal import find_peaks

# 1. Setup and Simulation Parameters
temp_k = 300
dt_fs = 0.5
total_steps = 5000  # Increased for better resolution

atoms = molecule('H2O')
atoms.calc = EMT()

# 2. Run MD
print(f"Starting MD simulation at Temperature: {temp_k} K")
traj_velocities = []
dyn = Langevin(atoms, dt_fs * units.fs, temperature_K=temp_k, friction=0.01)

def collect_velocities():
    traj_velocities.append(atoms.get_velocities())

dyn.attach(collect_velocities, interval=1)
dyn.run(total_steps)

# 3. Process Power Spectrum
vels = np.array(traj_velocities)
n_steps = vels.shape[0]

freq = np.fft.rfftfreq(n_steps, d=dt_fs)
fft_data = np.fft.rfft(vels, axis=0)
power_spectrum = np.sum(np.abs(fft_data)**2, axis=(1, 2))

# Convert to cm^-1
freq_cm = freq * (1e15 / 29979245800.0)

# 4. Extract Peaks
# 'prominence' helps ignore small noise peaks
peaks, props = find_peaks(power_spectrum, prominence=np.max(power_spectrum)*0.05)
extracted_freqs = sorted(freq_cm[peaks])

# 5. Experimental Data (Gas Phase H2O)
# Source: NIST CCCBDB
experimental = {
    "Bending": 1595,
    "Symmetric Stretch": 3657,
    "Asymmetric Stretch": 3756
}

print(f"\n--- Results (Temp: {temp_k} K) ---")
print(f"{'Mode':<20} | {'MD Peak (cm^-1)':<15} | {'Expt (cm^-1)':<12} | {'Error (%)'}")
print("-" * 65)

# Mapping logic: simplest version assumes the 3 peaks correspond to 3 modes
labels = ["Bending", "Symmetric Stretch", "Asymmetric Stretch"]
for i, (label, calc) in enumerate(zip(labels, extracted_freqs)):
    expt = experimental[label]
    error = abs(calc - expt) / expt * 100
    print(f"{label:<20} | {calc:15.1f} | {expt:12.1f} | {error:.2f}%")

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(freq_cm, power_spectrum, label='MD Power Spectrum')
plt.plot(freq_cm[peaks], power_spectrum[peaks], "x", color='red', label='Extracted Peaks')
plt.xlabel('Frequency (cm$^{-1}$)')
plt.ylabel('Intensity')
plt.title(f'Vibrational Spectrum at {temp_k}K')
plt.xlim(0, 4500)
plt.legend()
plt.show()

