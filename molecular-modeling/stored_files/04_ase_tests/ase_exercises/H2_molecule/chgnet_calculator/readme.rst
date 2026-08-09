==============================
H2 potential curve with CHGNet
==============================

@deepseek:  replace EMT with CHGNet calculator, do not forget to use periodic boxes

fix error...

python code01.py

Data saved to H2_potential_curve_CHGNet.txt

Summary of Results (CHGNet Machine Learning Potential)
------------------------------------------------------
Equilibrium Bond Length:              0.7590 Å
Experimental bond length:             0.7414 Å
Bond length error:                    0.0176 Å
Bond length relative error:           2.37%

Simulated Dissociation Energy (De):   4.3663 eV
Experimental De:                      4.5200 eV
Absolute Error in De:                 0.1537 eV
Relative Error in De:                 3.40%

============================================================
Periodic Boundary Conditions Verification
============================================================
✓ Periodic cell size: 10.0 Å (sufficient vacuum to isolate molecules)
✓ H2 bond vector aligned along x-axis
✓ Minimum image convention applied automatically by ASE
✓ Minimum distance between periodic images: 9.2 Å
  -> Sufficient vacuum (well-converged isolated molecule)

============================================================
About CHGNet
============================================================
CHGNet = Crystal Hamiltonian Graph Neural Network
Pre-trained on ~1.5M DFT calculations from Materials Project
Designed for universal interatomic potentials in materials
Expected accuracy for H2: De ~ 4.3-4.5 eV (vs experiment 4.52 eV)


morse, EMT does not work..
