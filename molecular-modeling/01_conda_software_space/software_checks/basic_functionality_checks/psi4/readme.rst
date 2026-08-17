====
PSI4
====

Test PSI4 code.


 conda install -c conda-forge "libxc<7.0.0" --force-reinstall  ??

Instead of forcing down libxc across your global environment setup, 
upgrade psi4 directly to its development release channel. 
This release channel contains updated source code that removes calls to deprecated LibXC variables:

conda install -c psi4/label/dev psi4 --solver=classic

(mace_env) miroi@MIRO:~/work/projects/mol_mat_modeling_schools/Molecular-and-Materials-Modeling-2026-August/molecular-modeling/01_conda_software_space/software_checks/basic_functionality_checks/psi4/.conda install -c psi4/label/dev psi4 --solver=classic
2 channel Terms of Service accepted
Collecting package metadata (current_repodata.json): done
Solving environment: -
... it takes too long ! omitting psi4 code



