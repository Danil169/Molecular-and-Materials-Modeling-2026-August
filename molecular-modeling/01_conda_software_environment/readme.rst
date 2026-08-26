===============================
Software maintenance with conda
===============================


get miniconda
-------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

add to .bashrc
~~~~~~~~~~~~~~

alias con='source /home/miroi/miniconda3/etc/profile.d/conda.sh; echo -e "miniconda activated :\c"; conda --version'
.con
miniconda activated :conda 26.5.3

working with conda
------------------

installing packages
~~~~~~~~~~~~~~~~~~~
conda activate myenv
conda config --add channels conda-forge
conda install ase nwchem xtb xtb-python mopac nwchem qe pyscf gromacs

which ase
/home/milias/miniconda3/envs/myenv/bin/ase
ase --version
ase-3.29.0

which mopac nwchem pw.x xtb
/home/milias/miniconda3/envs/myenv/bin/mopac
/home/milias/miniconda3/envs/myenv/bin/nwchem
/home/milias/miniconda3/envs/myenv/bin/pw.x
/home/milias/miniconda3/envs/myenv/bin/xtb


Challenge
=========
compare installable sofware versions : conda vs Linux package

mopac
xtb
nwchem
quantum espresso

