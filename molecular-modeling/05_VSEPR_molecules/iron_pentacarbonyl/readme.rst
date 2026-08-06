Fe(CO)6 theoretical geometry study
==================================

Avogadro
--------
gives distorted geometry (UFF, LJ)

MOPAC
-----
converged to square pyramidal geometry, which is wrong

it is closed-shell system


NWChem
------
mpirun -np 4 nwchem iron_pentacarbonyl_b3lyp_geopt.nw  > iron_pentacarbonyl_b3lyp_geopt.out

grep -A 2 'Step       Energy      Delta E' *.out

also converged to square pyramidal geometry...

TODO: to investigate it further, force both geometries - square pyramidal and  trigonal bipyramidal - and launch geometry optimizations from them

TODO: apply global minima search tools (CREST ? )



