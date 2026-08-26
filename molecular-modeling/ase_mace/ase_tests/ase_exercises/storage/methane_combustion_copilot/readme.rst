Methane combustion
==================

ASE reaxff simulate burning of methane with molecular dynamics
provide a ready-to-run ASE Python script

potential : https://openkim.org/id/Sim_LAMMPS_ReaxFF_ChenowethVanDuinGoddard_2008_CHO__SM_584143153761_001
wget https://openkim.org/files/1759075555/SM_584143153761_001/ffield.reax.cho

mv ffield.reax.cho  CHO.ff

python3 -c "from lammps import lammps; lmp = lammps(); print(lmp.version())"
LAMMPS (7 Feb 2024 - Update 1)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:98)
  using 1 OpenMP thread(s) per MPI task
20240207
Total wall time: 0:00:00

mising reaxff from precompiled package !!!

(myenv) miroi@MIRO:~/work/projects/Molecular-and-Materials-Modeling-2026-May/exercises/storage/methane_combustion_copilot/.lmp -h | grep reaxff

ERROR: Unrecognized pair style 'reaxff' is part of the REAXFF package which is not enabled in this LAMMPS binary. (src/force.cpp:275)
Last command: pair_style reaxff NULL
Total wall time: 0:00:00
