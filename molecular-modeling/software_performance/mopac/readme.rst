==============
MOPAC for WSL2
==============

installation
------------

for donwload

http://openmopac.net/Download_MOPAC_Executable_Step2.html

https://openmopac.net/download/

wget http://openmopac.net/mopac-23.1.2-linux.tar.gz

wget https://github.com/openmopac/mopac/releases/download/v23.2.3/mopac-23.2.3-linux.run

or
wget https://github.com/openmopac/mopac/releases/download/v23.2.3/mopac-23.2.3-linux.tar.gz


tar xvzf mopac-23.2.3-linux.tar.gz
milias@DESKTOP-7OTLCGO:~/work/software/mopac/mopac-23.2.3-linux/.ls
CITATION.cff  LICENSE  bin/  examples/  include/  lib/

"mini DNA" perfomance test
--------------------------

Notebook:
~~~~~~~~~

Nthr   Job time
1      27.50s
6      18.49s  
12     18.34s


Desktop PC:
~~~~~~~~~~~
12th Gen Intel(R) Core(TM) i7-12700K

Nthr   Job time
3       30.58s
6       10.82s  <-----
10      15.09s
12      15.61s

Nthr    Job tim
4       97.62        
10      12.23
18      12.07

see also http://openmopac.net/Manual/Reducing_computation_time.html
