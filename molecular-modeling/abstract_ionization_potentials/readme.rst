=============================================================================================================================
Ionization potentials of Zn, Cd, Hg and dipole polarizabilities of Zn+, Cd+, Hg+: correlation and relativistic effects (v0.1)
=============================================================================================================================

https://doi.org/10.1016/S0009-2614(99)00665-X

Recalculating published data using computationally cheaper methods.

In this computational exercise we focus on the heaviest atom of mercury.

The first ionization potential represents the energy for removal of the electron from the neutral atom.
It is calculated as the difference between energy of ionized atom (M+) and the neutral atom (M):

IP = E(Hg+) - E(Hg)

Hg
--
We focus only on the mercury atom. 

We can employ only quantum mechanical based methods because they are able to describe
electronic structure on neutral and ionized species.

Computational apparatus
~~~~~~~~~~~~~~~~~~~~~~~
conda create -n molmatmodel
conda activate  molmatmodel
conda install -n molmatmodel -c conda-forge nwchem mopac xtb xtb-python pyscf qe ase pymace

quick check of the software
~~~~~~~~~~~~~~~~~~~~~~~~~~~
python -c "import shutil; [print(f'✅ {x}: Ready') if exec(f'try:\n import {x}\nexcept:\n raise') is None else None for x in ['ase', 'xtb', 'pyscf', 'mace']]; [print(f'✅ {bin}: Ready') if shutil.which(bin) else print(f'❌ {bin}: Not found') for bin in ['mopac', 'xtb', 'nwchem', 'pw.x']]"


runs
----
python hg_ie_calculation_01.py > hg_ie_calculation_01.py_logfile

Results and discussion
~~~~~~~~~~~~~~~~~~~~~~

Experiment 10.437 eV


deepseek 
--------
link for script buildup:



