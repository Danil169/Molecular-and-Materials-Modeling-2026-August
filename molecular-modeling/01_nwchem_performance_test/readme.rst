==========================
CH3(.) radical with NWChem
==========================

installation
------------
sudo apt install nwchem

check with
~~~~~~~~~~
lsb_release -a

which nwchem
/usr/bin/nwchem

ldd /usr/bin/nwchem

dpkg -l nwchem
dpkg -S nwchem

dpkg -L nwchem-data

man nwchem

find nwchem version
~~~~~~~~~~~~~~~~~~~
apt show nwchem

dpkg -l nwchem

running NWChem interactively
----------------------------
example from https://nwchemgit.github.io/EPR-pNMR.html

mpirun -np 6 nwchem ch3_zora_b3lyp_prop.nw


or via script with command line:

wsl_bash_run.01

wsl_bash_run.01  4

wsl_bash_run.01  6

some perfomance data
--------------------
12th Gen Intel(R) Core(TM) i5-12450H
NWChem 7.0.2 (Ubuntu 22.04)

Nthr    Wall
2        9.9s
4        6.8s
6        6.6s

NWChem 7.2.2 (Ubuntu 24.04.2)
Nthr    Wall
2       11.5s
4       7.5s
6       7.7s

Task
-----
- find out list of available basis set installed with nwchem package:

dpkg -L nwchem-data

ls /usr/share/nwchem/libraries/*

- use larger basis set than the current 6-311G basis set, make table basis set vs wall time

/usr/share/nwchem/libraries/6-311g  ... 6-311G
6-311gs_polarization  ... 6-311G*


basis    NBF
6-311G   13+3
6-311G*  19+3


