====
PSI4
====

Test PSI4 code.


 conda install -c conda-forge "libxc<7.0.0" --force-reinstall  ??

 conda install -c conda-forge "libxc<7.0.0" --force-reinstall
2 channel Terms of Service accepted
Channels:
 - conda-forge
 - defaults
Platform: linux-64
Collecting package metadata (repodata.json): done
Solving environment: / WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
- WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
| WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
/ WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
- WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
\ WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
/ WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
- WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
| WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
WARNING conda.conda_libmamba_solver.solver:_prepare_problems_message(854): Failed to explain problems
Traceback (most recent call last):
  File "/home/miroi/miniconda3/lib/python3.14/site-packages/conda_libmamba_solver/solver.py", line 851, in _prepare_problems_message
    explained_errors = unsolvable.explain_problems(db, problems_format_auto)
libmambapy.bindings.legacy.MambaNativeException: Error parsing version "1.*^". Version contains invalid characters in 1.*^.
failed

LibMambaUnsatisfiableError: Encountered problems while solving:
  - package libxc-6.1.0-cpu_hb3673ca_2 requires libxc-c 6.1.0 cpu_hd8589cd_2, but none of the providers can be installed



Instead of forcing down libxc across your global environment setup, 
upgrade psi4 directly to its development release channel. 
This release channel contains updated source code that removes calls to deprecated LibXC variables:

conda install -c psi4/label/dev psi4 --solver=classic

(mace_env) miroi@MIRO:~/work/projects/mol_mat_modeling_schools/Molecular-and-Materials-Modeling-2026-August/molecular-modeling/01_conda_software_space/software_checks/basic_functionality_checks/psi4/.conda install -c psi4/label/dev psi4 --solver=classic
2 channel Terms of Service accepted
Collecting package metadata (current_repodata.json): done
Solving environment: -
... it takes too long ! omitting psi4 code



