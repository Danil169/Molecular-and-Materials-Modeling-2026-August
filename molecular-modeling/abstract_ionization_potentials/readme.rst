=============================================================================================================================
Ionization potentials of Zn, Cd, Hg and dipole polarizabilities of Zn+, Cd+, Hg+: correlation and relativistic effects (v0.1)
=============================================================================================================================

https://doi.org/10.1016/S0009-2614(99)00665-X

Recalculating published data using computationally cheaper methods.

In this computational exercise we focus on the heaviest atom of mercury.

The first ionization potential represents the energy for removal of the electron from the neutral atom.
It is calculated as the difference between energy of ionized atom (M+) and the neutral atom (M):

IP = E(Hg+) - E(Hg)

We focus only on the mercury atom. 

We can employ only quantum mechanical based methods because they are able to describe
electronic structure on neutral and ionized species.

Computational apparatus
-----------------------
conda create -n molmatmodel
conda activate  molmatmodel
conda install -n molmatmodel -c conda-forge nwchem mopac xtb xtb-python pyscf qe ase pymace

quick check of the software
~~~~~~~~~~~~~~~~~~~~~~~~~~~
python -c "import shutil; [print(f'✅ {x}: Ready') if exec(f'try:\n import {x}\nexcept:\n raise') is None else None for x in ['ase', 'xtb', 'pyscf', 'mace']]; [print(f'✅ {bin}: Ready') if shutil.which(bin) else print(f'❌ {bin}: Not found') for bin in ['mopac', 'xtb', 'nwchem', 'pw.x']]"

which xtb; xtb --version

runs
----
python hg_ie_calculation_01.py > hg_ie_calculation_01.py_logfile

Results and discussion
----------------------

Experimental IP(Hg): **10.437 eV**

.. csv-table::
   :header: "Method / Software", "IP (eV)", "Error vs exp."

   PM7 (MOPAC UHF),            10.518,  +0.08 eV  (+0.8 %)
   GFN1-xTB,                   16.768,  +6.33 eV  (+61 %)  ← fails
   CCSD(T)/def2-SVP (NWChem),  10.018,  −0.42 eV
   CCSD(T)/def2-TZVP (NWChem), 10.139,  −0.30 eV
   CCSD(T)/def2-QZVP (NWChem), 10.325,  −0.11 eV
   ECP-CCSD(T)/def2-TZVP (PySCF), 10.125, −0.31 eV

PM7 gives an excellent result for Hg (0.8 % error) due to well-fitted
semi-empirical parameters for heavy elements.  NWChem CCSD(T) results
improve systematically with basis set size, reaching 10.325 eV with def2-QZVP.
xTB GFN1 is not parameterised for heavy-element ionisation and fails badly.

Full write-up: ``abstract.rst``.
Working files: ``mopac/``, ``xtb/``, ``nwchem/``, ``pyscf/``.

deepseek 
--------
link for AI thread:  https://chat.deepseek.com/share/vpbbo1f5d5b4gzdxsj
