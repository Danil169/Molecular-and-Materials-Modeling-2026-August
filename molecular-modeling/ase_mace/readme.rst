==============
ASE with MLIPs
==============

Working in the Atomic Simulation Environment (ASE) with Machine-Learned
Interatomic Potentials (MLIPs), specifically MACE.

Installation
------------

::

    pip install mace-torch ase

or via conda::

    conda install -c conda-forge ase pymace

Exercises
---------

Exercise I.3 — ASE + MACE: atomization energies of diatomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Located in ``ase_tests/01_atomization_energy/``.

Compute the H₂ potential-energy curve and atomization energies of N₂ and O₂
with the MACE-MP-0 universal force field.  Key results:

- H₂ potential-energy curve: minimum at r ≈ 0.74 Å (see ``../report_figures_part1/h2_potential_curve_mace.png``)
- N₂ atomization energy: MACE vs experiment comparison
- O₂ (triplet): r(O–O) = 1.2367 Å (exp. 1.2070 Å), ν = 1419 cm⁻¹ (exp. 1580 cm⁻¹)

Also in ``ase_tests/``:

- ``02_thermodyn_prop/`` — thermodynamic properties via ASE thermochemistry module
- ``03_ethane_C-C_bond_energy/`` — C–C and C–H bond dissociation energies (EMT vs MACE)
- ``04_quantum_espresso_water/`` — QE water geometry optimisation driven from ASE (Exercise I.7)
- ``05_openbabel_water/`` — OpenBabel structure conversion
- ``ase_exercises/`` — additional ASE exercises

Exercise I.4 — MACE installation test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Located in ``mace/``.

Run the official MACE test suite to verify installation::

    python -m mace.calculators.tests.test_mace_calculator

22 of 25 tests pass (3 failures related to optional CUDA/MPS backends not
present on the test machine).

Challenge I.3
~~~~~~~~~~~~~
Modify the atomization-energy script (using an AI tool) to compute energies
of other diatomics (e.g. NO, CO, Cl₂).

Challenge I.4
~~~~~~~~~~~~~
Run all testing scripts shown in the video by yourself.


Installing collected packages: python-hostlist, nvidia-cusparselt-cu13, mpmath, cuda-toolkit, wcwidth, triton, tqdm, sympy, smmap, setuptools, pyYAML, orjson, opt_einsum, nvidia-nvtx, nvidia-nvshmem-cu13, nvidia-nvjitlink, nvidia-nccl-cu13, nvidia-curand, nvidia-cufile, nvidia-cuda-runtime, nvidia-cuda-nvrtc, nvidia-cuda-cupti, networkx, MarkupSafe, lmdb, lightning-utilities, h5py, fsspec, filelock, cuda-pathfinder, configargparse, prettytable, pandas, nvidia-cusparse, nvidia-cufft, nvidia-cublas, jinja2, gitdb, cuda-bindings, nvidia-cusolver, nvidia-cudnn-cu13, GitPython, matscipy, torch, torchmetrics, torch-ema, opt-einsum-fx, e3nn, mace-torch
Successfully installed GitPython-3.1.58 MarkupSafe-3.0.3 configargparse-1.7.5 cuda-bindings-13.3.1 cuda-pathfinder-1.6.0 cuda-toolkit-13.0.3.0 e3nn-0.4.4 filelock-3.32.2 fsspec-2026.7.0 gitdb-4.0.12 h5py-3.16.0 jinja2-3.1.6 lightning-utilities-0.15.3 lmdb-2.3.0 mace-torch-0.3.16 matscipy-1.2.0 mpmath-1.3.0 networkx-3.6.1 nvidia-cublas-13.1.1.3 nvidia-cuda-cupti-13.0.85 nvidia-cuda-nvrtc-13.0.88 nvidia-cuda-runtime-13.0.96 nvidia-cudnn-cu13-9.20.0.48 nvidia-cufft-12.0.0.61 nvidia-cufile-1.15.1.6 nvidia-curand-10.4.0.35 nvidia-cusolver-12.0.4.66 nvidia-cusparse-12.6.3.3 nvidia-cusparselt-cu13-0.8.1 nvidia-nccl-cu13-2.29.7 nvidia-nvjitlink-13.3.33 nvidia-nvshmem-cu13-3.4.5 nvidia-nvtx-13.0.85 opt-einsum-fx-0.1.4 opt_einsum-3.4.0 orjson-3.11.9 pandas-3.0.5 prettytable-3.18.0 pyYAML-6.0.3 python-hostlist-2.3.0 setuptools-84.0.0 smmap-5.0.3 sympy-1.14.0 torch-2.13.0 torch-ema-0.3 torchmetrics-1.9.0 tqdm-4.70.0 triton-3.7.1 wcwidth-0.8.2





