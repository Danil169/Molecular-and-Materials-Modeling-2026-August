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
conda install ase nwchem xtb xtb-python mopac nwchem qe pyscf


Challenge
=========
compare installable sofware versions : conda vs Linux package

mopac
xtb
nwchem
quantum espresso

