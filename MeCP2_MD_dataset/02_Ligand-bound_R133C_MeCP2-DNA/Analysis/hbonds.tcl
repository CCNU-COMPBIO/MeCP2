mol delete all
mol load parm7 Arg133Cys_DB06237.com.prmtop
mol addfile Arg133Cys_DB06237_run1.nc first 400 last 4000 step 1 waitfor all

set com [atomselect top "all"]
set bp [atomselect top "protein"]
set nuc [atomselect top "nucleic"]

file mkdir hbonds

hbonds -sel1 $bp  -dist 3.5 -ang 30  -plot no  -type pair  -writefile yes  -detailout 3c2i_bp-detail.dat -outdir hbonds -log 3c2i_bp.log -outfile 3c2i_bp-hbonds.dat
hbonds -sel1 $bp -sel2 $nuc -dist 3.5 -ang 30 -plot no  -type pair -writefile yes -detailout 3c2i_bp-nuc-detail.dat -outdir hbonds -log 3c2i_bp-nuc.log -outfile 3c2i_bp-nuc-hbonds.dat


exit 
