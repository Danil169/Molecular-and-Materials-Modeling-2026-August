import numpy as np
import matplotlib.pyplot as plt
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.md.langevin import Langevin
from ase import units

# 1. Setup the System
atoms = molecule('H2O')
atoms.calc = EMT() # Replace with your calculator (VASP, Gaussian, etc.)

# 2. Run MD to collect velocities
# We use a small timestep (0.5 fs) to capture high-frequency vibrations
traj_velocities = []
dyn = Langevin(atoms, 0.5 * units.fs, temperature_K=300, friction=0.01)

def collect_velocities():
    traj_velocities.append(atoms.get_velocities())

dyn.attach(collect_velocities, interval=1)
dyn.run(2000) # Run for 1000 fs (1 ps)

# 3. Calculate Power Spectrum (Fourier Transform of Velocities)
vels = np.array(traj_velocities)  # Shape: (steps, atoms, 3)
n_steps = vels.shape[0]
dt = 0.5 # fs

# Compute the FFT of the velocities
# We sum across atoms and dimensions to get the total VDOS
freq = np.fft.rfftfreq(n_steps, d=dt) 
fft_data = np.fft.rfft(vels, axis=0)
power_spectrum = np.sum(np.abs(fft_data)**2, axis=(1, 2))

# 4. Convert frequency to cm^-1
# freq is in 1/fs. 1/fs -> Hz (1e15) -> cm^-1 (divide by speed of light in cm/s)
freq_cm = freq * (1e15 / 29979245800.0)

# Plotting
plt.plot(freq_cm, power_spectrum)
plt.xlabel('Frequency (cm$^{-1}$)')
plt.ylabel('Intensity (arb. units)')
plt.title('Anharmonic Vibrational Spectrum (MD at 300K)')
plt.xlim(0, 4000) 
plt.show()

