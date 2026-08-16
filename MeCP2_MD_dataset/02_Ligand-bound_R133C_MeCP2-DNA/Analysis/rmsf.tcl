###########################################################
# 批量 RMSF 计算脚本：Arg133Cys 空白体系 + 小分子体系
#
# 逻辑：
# 1. 对齐：蛋白质骨架 + DNA 骨架
# 2. 计算：蛋白质 Cα + DNA 骨架 RMSF
# 3. 输出：
#    - Protein Cα RMSF
#    - DNA backbone per-residue average RMSF
#
# 适用于：
# Arg133Cys_no_ligand run1/run2/run3
# Arg133Cys_DB00694 run1/run2/run3
# Arg133Cys_DB00445 run1/run2/run3
# Arg133Cys_DB06237 run1/run2/run3
#
# 文件格式：
# Arg133Cys_no_ligand.com.prmtop
# Arg133Cys_no_ligand_run1.nc
#
# Arg133Cys_DB00694.com.prmtop
# Arg133Cys_DB00694_run1.nc
#
# 输出目录：
# ./RMSF/
###########################################################


# ===================== 基本参数 =====================

# 图片中对应的所有体系前缀
set systems {
    Arg133Cys_no_ligand
    Arg133Cys_DB00694
    Arg133Cys_DB00445
    Arg133Cys_DB06237
}

set runs {
    run1
    run2
    run3
}

# 参考帧
# 399 对应约 20 ns，因为 4000 帧 = 200 ns，每帧 0.05 ns
set start_frame 399
set ref_frame $start_frame

# 所有体系均为 4000 帧，200 ns
set max_frame 4000

# 每 1 帧计算一次
set step 1


# ===================== 对齐和 RMSF 选择 =====================
# 对齐：蛋白 backbone + DNA backbone
# 蛋白 backbone: N CA C O
# DNA backbone: nucleic and backbone
set align_sel "(protein and name N CA C O) or (nucleic and backbone)"

# RMSF 计算：蛋白 Cα + DNA backbone
set rmsf_sel "(protein and name CA) or (nucleic and backbone)"


# ===================== 输出目录 =====================

if {![file exists "./RMSF"]} {
    file mkdir "./RMSF"
    puts "已创建输出目录：./RMSF"
}


# ===================== 主循环 =====================

foreach sysname $systems {

    set topfile "./${sysname}.com.prmtop"

    if {![file exists $topfile]} {
        puts "跳过 ${sysname}：找不到拓扑文件 $topfile"
        continue
    }

    foreach run $runs {

        set trajfile "./${sysname}_${run}.nc"
        set outfile  "./RMSF/${sysname}_${run}_complex_RMSF.dat"

        if {![file exists $trajfile]} {
            puts "跳过 ${sysname}_${run}：找不到轨迹文件 $trajfile"
            continue
        }

        puts ""
        puts "========================================"
        puts "处理体系：${sysname}_${run}"
        puts "拓扑文件：$topfile"
        puts "轨迹文件：$trajfile"
        puts "输出文件：$outfile"
        puts "========================================"


        # ===================== 加载轨迹 =====================

        mol new $topfile type parm7 waitfor all
        mol addfile $trajfile type netcdf waitfor all

        set molid [molinfo top]
        set total_frames [molinfo $molid get numframes]

        set nf $max_frame
        if {$total_frames < $nf} {
            set nf $total_frames
        }

        set last_frame [expr {$nf - 1}]

        puts "轨迹总帧数：$total_frames"
        puts "实际用于 RMSF 的帧数：$nf"
        puts "使用帧范围：0 到 $last_frame"
        puts "参考帧：$ref_frame"


        # ===================== 检查 ref_frame 是否合理 =====================

        if {$ref_frame >= $nf} {
            puts "错误：参考帧 ref_frame = $ref_frame 超出当前轨迹帧数 $nf"
            mol delete $molid
            continue
        }


        # ===================== 检查原子选择 =====================

        set align_check [atomselect $molid $align_sel frame $ref_frame]
        set rmsf_check  [atomselect $molid $rmsf_sel frame $ref_frame]

        set n_align [$align_check num]
        set n_rmsf  [$rmsf_check num]

        puts "对齐原子数：$n_align"
        puts "RMSF 计算原子数：$n_rmsf"

        if {$n_align == 0} {
            puts "错误：对齐选择为空，请检查 align_sel："
            puts "$align_sel"
            $align_check delete
            $rmsf_check delete
            mol delete $molid
            continue
        }

        if {$n_rmsf == 0} {
            puts "错误：RMSF 选择为空，请检查 rmsf_sel："
            puts "$rmsf_sel"
            $align_check delete
            $rmsf_check delete
            mol delete $molid
            continue
        }

        $align_check delete
        $rmsf_check delete


        # ===================== 对齐：蛋白骨架 + DNA 骨架 =====================
        # 以 ref_frame 为参考帧，将每一帧的 all atoms 按照 align_sel 对齐过去。
        # 移动 all atoms 是为了保证后续 RMSF 基于已对齐轨迹。

        set ref_atoms [atomselect $molid $align_sel frame $ref_frame]
        set mov_atoms [atomselect $molid $align_sel]
        set all_atoms [atomselect $molid "all"]

        for {set i 0} {$i <= $last_frame} {incr i $step} {

            $mov_atoms frame $i
            $all_atoms frame $i

            set trans [measure fit $mov_atoms $ref_atoms]
            $all_atoms move $trans
        }

        $ref_atoms delete
        $mov_atoms delete
        $all_atoms delete


        # ===================== 计算 RMSF：蛋白 Cα + DNA backbone =====================
        # measure rmsf 返回 rmsf_sel 中每个原子的 RMSF 值。
        # 将 RMSF 存入 beta 字段，后续按 resid 提取。

        set rmsf_atoms [atomselect $molid $rmsf_sel]
        set rmsf_vals [measure rmsf $rmsf_atoms first 0 last $last_frame step $step]
        $rmsf_atoms set beta $rmsf_vals


        # ===================== 输出 RMSF 文件 =====================

        set fout [open $outfile w]

        puts $fout "# RMSF Calculation Result"
        puts $fout "# System: ${sysname}_${run}"
        puts $fout "# Alignment: Protein backbone + DNA backbone"
        puts $fout "# Alignment selection: $align_sel"
        puts $fout "# RMSF calculation: Protein C-alpha + DNA backbone"
        puts $fout "# RMSF selection: $rmsf_sel"
        puts $fout "# Frames used: 0 to $last_frame, step $step"
        puts $fout "# Reference frame: $ref_frame"
        puts $fout "# ResID\tType\tRMSF_Angstrom"


        # ===================== 输出蛋白 Cα RMSF =====================
        # 自动获取蛋白 Cα 的 resid，避免固定 1-69 后因为编号不同出错。

        set prot_ca_sel [atomselect $molid "protein and name CA"]
        set prot_resids [lsort -integer -unique [$prot_ca_sel get resid]]
        $prot_ca_sel delete

        foreach resid $prot_resids {

            set sel [atomselect $molid "protein and name CA and resid $resid"]

            if {[$sel num] > 0} {
                set val [lindex [$sel get beta] 0]
                puts $fout [format "%d\tProtein\t%.4f" $resid $val]
            }

            $sel delete
        }


        # ===================== 输出 DNA backbone RMSF =====================
        # DNA 一个残基包含多个 backbone 原子。
        # 对同一个 DNA resid 的 backbone 原子 RMSF 取平均。

        set dna_bb_sel [atomselect $molid "nucleic and backbone"]
        set dna_resids [lsort -integer -unique [$dna_bb_sel get resid]]
        $dna_bb_sel delete

        foreach resid $dna_resids {

            set sel [atomselect $molid "nucleic and backbone and resid $resid"]

            if {[$sel num] > 0} {

                set vals [$sel get beta]

                set sum 0.0
                set n 0

                foreach v $vals {
                    set sum [expr {$sum + $v}]
                    incr n
                }

                if {$n > 0} {
                    set avg [expr {$sum / $n}]
                    puts $fout [format "%d\tDNA\t%.4f" $resid $avg]
                }
            }

            $sel delete
        }


        # ===================== 清理 =====================

        close $fout
        $rmsf_atoms delete
        mol delete $molid

        puts "✅ 完成：${sysname}_${run}"
        puts "输出：$outfile"
    }
}


puts ""
puts "============================================="
puts "  所有 RMSF 计算完成！"
puts "  结果保存在 ./RMSF/ 文件夹"
puts "============================================="

exit
