###########################################################
# 批量RMSF计算脚本（论文最终版）
# 对齐方式：蛋白质骨架 + DNA 骨架
# 计算范围：蛋白质 Cα + DNA 骨架 RMSF
# 适用于：WT + 8个突变体，3次重复 run1/run2/run3
###########################################################

# ===================== 批量体系列表 =====================
foreach {name run} {
    Arg106Gly run1 Arg106Gly run2 Arg106Gly run3
    Arg106Leu run1 Arg106Leu run2 Arg106Leu run3
    Arg106Trp run1 Arg106Trp run2 Arg106Trp run3
    Arg111Gly run1 Arg111Gly run2 Arg111Gly run3
    Arg133Cys run1 Arg133Cys run2 Arg133Cys run3
    Arg133Gly run1 Arg133Gly run2 Arg133Gly run3
    Arg133Leu run1 Arg133Leu run2 Arg133Leu run3
    Glu143Lys run1 Glu143Lys run2 Glu143Lys run3
    WT run1 WT run2 WT run3
} {

    puts "========================================"
    puts "处理体系：$name $run"
    puts "========================================"

    # 加载文件
    mol load parm7 ./${name}_com.prmtop
    mol addfile ./${name}_${run}.nc first 10000 last 100000 step 10 waitfor all

    set total_frames [molinfo top get numframes]
    set nframes [expr $total_frames - 1]

    # ===================== 对齐：蛋白骨架 + DNA骨架 =====================
    set align_sel "(protein and name N CA C O) or (nucleic and backbone)"
    set ref_atoms [atomselect top $align_sel frame 0]
    set mov_atoms [atomselect top $align_sel]
    set all_atoms [atomselect top "all"]

    for {set i 1} {$i <= $nframes} {incr i} {
        $mov_atoms frame $i
        $all_atoms frame $i
        set trans [measure fit $mov_atoms $ref_atoms]
        $all_atoms move $trans
    }

    $ref_atoms delete
    $mov_atoms delete
    $all_atoms delete

    # ===================== 计算 RMSF：蛋白 Cα + DNA 骨架 =====================
    set rmsf_sel "(protein and name CA) or (nucleic and backbone)"
    set rmsf_atoms [atomselect top $rmsf_sel]
    set rmsf_vals [measure rmsf $rmsf_atoms first 1 last $nframes]
    $rmsf_atoms set beta $rmsf_vals

    # ===================== 输出 RMSF 文件 =====================
    set outfile [open ./RMSF/${name}_${run}_complex_RMSF.dat w]
    puts $outfile "# RMSF Calculation Result"
    puts $outfile "# Alignment: Protein Backbone + DNA Backbone"
    puts $outfile "# Calculation: Protein Cα + DNA Backbone"
    puts $outfile "# ResID\tType\tRMSF(Å)"

    set prot_min 1
    set prot_max 69

    # 输出蛋白 Cα RMSF
    for {set i $prot_min} {$i <= $prot_max} {incr i} {
        set sel [atomselect top "protein and name CA and resid $i"]
        if {[$sel num] > 0} {
            set val [lindex [$sel get beta] 0]
            puts $outfile [format "%d\tProtein\t%.4f" $i $val]
        }
        $sel delete
    }

    # 输出 DNA 骨架 RMSF
    set dna_res [lsort -integer -unique [[atomselect top "nucleic and backbone"] get resid]]
    foreach res $dna_res {
        set sel [atomselect top "nucleic and backbone and resid $res"]
        if {[$sel num] > 0} {
            set vals [$sel get beta]
            set sum 0
            foreach v $vals { set sum [expr $sum + $v] }
            set avg [expr $sum / [llength $vals]]
            puts $outfile [format "%d\tDNA\t%.4f" $res $avg]
        }
        $sel delete
    }

    close $outfile
    $rmsf_atoms delete
    mol delete top

    puts "✅ $name $run 完成 → ./RMSF/${name}_${run}_complex_RMSF.dat"
}

exit
