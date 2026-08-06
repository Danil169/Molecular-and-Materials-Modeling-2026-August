H2 potential curve with MACE calculator
=======================================

@deepseek  : code04.py .. . replace CHGNet with MACE calculator, write the full code

python code01.py
cuequivariance or cuequivariance_torch is not available. Cuequivariance acceleration will be disabled.
✓ mace-models imported successfully

============================================================
H2 Potential Energy Curve with MACE Machine Learning Potential
============================================================

1. Loading MACE model...
   Trying MACE-MP-0_small...

        You're using the MACE-MP-0_small model. The model is released under the MIT license.
        Note:
        If you are using this model, please cite the relevant paper for the Materials Project,
        any paper associated with the MACE model, and also the following:
        - MACE-Universal by Yuan Chiang, 2023, Hugging Face, Revision e5ebd9b,
            DOI: 10.57967/hf/1202, URL: https://huggingface.co/cyrusyc/mace-universal
        - Matbench Discovery by Janosh Riebesell, Rhys EA Goodall, Philipp Benner, Yuan Chiang,
            Alpha A Lee, Anubhav Jain, Kristin A Persson, 2023, arXiv:2308.14920
        - https://arxiv.org/abs/2401.00096

   ✓ Loaded MACE-MP-0_small

2. Calculating isolated H atom energy...
   Energy of single H atom: -1.21355788 eV
   Energy of two isolated H atoms: -2.42711576 eV

3. Calculating potential energy curve...
--------------------------------------------------
   r = 0.500 Å | E = -4.301223 eV
   r = 1.114 Å | E = -5.351814 eV
   r = 1.795 Å | E = -2.996450 eV
   r = 2.477 Å | E = -2.370749 eV
   r = 3.159 Å | E = -2.405785 eV
   r = 3.500 Å | E = -2.414478 eV
--------------------------------------------------

4. Finding equilibrium bond length...
   Quadratic fit: r_e = 0.7540 Å, E_min = -6.554578 eV

5. Dissociation Energy Calculation:
   Simulated De = 4.1275 eV
   Experimental De = 4.5200 eV
   Error = 0.3925 eV

6. Data saved to H2_potential_curve_MACE-MP-0_small.txt

============================================================
✓ Calculation completed successfully!
============================================================

downloaded model
----------------
(venv) milias@DESKTOP-7OTLCGO:~/work/projects/Molecular-and-Materials-Modeling-2026-May/exercises/H2_potential_curve/mace_calculator/.ls data/
2023-12-10-mace-128-L0_energy_epoch-249.model
