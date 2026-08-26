====================================
Quantum Espresso Software Perfomance 
====================================

checks
------
which pw.x
/home/milias/miniconda3/envs/molmatmodel/bin/pw.x

which mpirun
/home/milias/miniconda3/envs/molmatmodel/bin/mpirun

runs
----

export OMP_NUM_THREADS=1
mpirun -np 4 pw.x -in Tl2.so_scf.in > Tl2.so_scf.in_out_OMP1_MPI4
mpirun -np 8 pw.x -in Tl2.so_scf.in > Tl2.so_scf.in_out_OMP1_MPI8_bltpDesktop

#MPI  wall
4      12.93s
8      13.32s
