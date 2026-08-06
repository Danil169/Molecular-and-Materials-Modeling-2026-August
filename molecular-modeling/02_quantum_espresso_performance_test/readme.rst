==================================
Quantum Espresso benchmark example
==================================

installation of QE
------------------

Ubuntu 24.04
~~~~~~~~~~~~

broken native package 24.04, see https://gitlab.com/QEF/q-e/-/issues/684 !

instead - download the deb-package from  https://github.com/pranabdas/espresso/releases

sudo apt install ./quantum-espresso_7.5-1_amd64.deb

export PATH=/opt/espresso/7.5:$PATH

check on linked libraries: ldd /opt/espresso/7.5/pw.x

which pw.x
/opt/espresso/7.5/pw.x

also provide: 
sudo apt-get install c2x quantum-espresso-data quantum-espresso-data-sssp

check installed QE
~~~~~~~~~~~~~~~~~~

summary: 
dpkg -s quantum-espresso

list files: dpkg -L  quantum-espresso  quantum-espresso-data

installed pseudos: 
ls -l /usr/share/espresso/pseudo/

find QE version
~~~~~~~~~~~~~~~
apt show quantum-espresso

dpkg -l quantum-espresso


QE visualizer
~~~~~~~~~~~~~
sudo apt install xcrysden

displays input, output files

performance test
----------------
We have simple example of 2 molecules, QE inputs is adapted for package's default pseudos.

Consider both OpenMPI processes and OpenMP threads.

launching on 1 MPI process :  wsl2_bash_run.01 1
launching on 2 MPI processes :  wsl2_bash_run.01 2
launching on 4 MPI processes :  wsl2_bash_run.01 4

example on my notebook
~~~~~~~~~~~~~~~~~~~~~~
12th Gen Intel(R) Core(TM) i5-12450H
PWSCF v.7.5 (Ubuntu 24.04)

export OMP_NUM_THREADS=1

N(OpemMPI)   Wall 
1          1m 9.57s
2          35.68s
6          18.18s

challenge
~~~~~~~~~
Ask Google AI on your specific CPU xxx: "how many OpenMP threads and how many OpenMPI processes for xxx"

Increase the cutoffs and run more performance tests with changing both OpenMPI and OpenMP parameters.

