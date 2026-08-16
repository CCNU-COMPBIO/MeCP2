#!/bin/bash
# NaCl_Box_30 NaCl_Box_50 CA_15mM_1264_Box_30 CA_15mM_1264_Box_50 MG_15mM_1264_Box_30 MG_15mM_1264_Box_50
for file in Arg106Gly 
do 
    for run in run1 run2 run3
    do
        mkdir ${file}_${run}
        cd ${file}_${run}
        python ../MMPBSA.py -O -i ../mmpbsa.in -cp ../${file}_com.prmtop  -rp ../${file}_bp.prmtop  -lp ../${file}_dna.prmtop  -y ../../target/${file}_target_${run}.nc
        cd ..
        
    done
done
        



