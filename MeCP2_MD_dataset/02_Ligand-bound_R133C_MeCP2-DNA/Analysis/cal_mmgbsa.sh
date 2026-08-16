#!/bin/bash

for file in Arg133Cys_DB00266 
do
    for run in run1
    do
        echo "========================================"
        echo "Running MMPBSA for ${file}_${run}"
        echo "========================================"

        mkdir -p ${file}_${run}
        cd ${file}_${run} || exit

        python ../MMPBSA.py -O \
            -i ../mmgbsa.in \
            -cp ../../traj/${file}_proDNA.prmtop \
            -rp ../../traj/${file}_protein.prmtop \
            -lp ../../traj/${file}_DNA.prmtop \
            -y ../../traj/${file}_${run}_proDNA.nc \
            -o FINAL_RESULTS_${file}_${run}_protein_DNA.dat \
            -do FINAL_DECOMP_${file}_${run}_protein_DNA.dat

        cd ..

        echo "完成：${file}_${run}"
    done
done


