mol delete all
mol load parm7 ./mecp2_MD_traj/Arg133Gly_com.prmtop
mol addfile ./mecp2_MD_traj_fix_all/Arg133Gly_rw_run1.nc first 10000 last 100000 step 10 waitfor all

set com [atomselect top "all"]
set bp [atomselect top "protein"]
set nuc [atomselect top "nucleic"]

file mkdir hbonds

hbonds -sel1 $bp  -dist 3.5 -ang 30  -plot no  -type pair  -writefile yes  -detailout 3c2i_bp-detail.dat -outdir hbonds -log 3c2i_bp.log -outfile 3c2i_bp-hbonds.dat
hbonds -sel1 $bp -sel2 $nuc -dist 3.5 -ang 30 -plot no  -type pair -writefile yes -detailout 3c2i_bp-nuc-detail.dat -outdir hbonds -log 3c2i_bp-nuc.log -outfile 3c2i_bp-nuc-hbonds.dat


exit 
