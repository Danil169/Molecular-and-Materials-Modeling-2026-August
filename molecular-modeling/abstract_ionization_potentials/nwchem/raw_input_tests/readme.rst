=====================
NWChem raw input jobs
=====================

mpirun -np 4 nwchem hg_scf_ecp.nw > hg_scf_ecp.nw_logfile
mpirun -np 4 nwchem hg1plus_scf_ecp.nw > hg1plus_scf_ecp.nw_logfile
mpirun -np 4 nwchem hg1plus_scf-rohf_ecp.nw > hg1plus_scf-rohf_ecp.nw_logfile
