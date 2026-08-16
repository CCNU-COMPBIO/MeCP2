# ============ VMD 批量RMSD计算脚本 - 论文标准版 ============
# 研究场景：WT + 点突变蛋白-DNA复合物分子动力学模拟
# 功能：
# 1. 对齐：蛋白质骨架 + DNA骨架（刚性区，最稳定）
# 2. 计算：蛋白-DNA 复合物整体骨架 RMSD
# 3. 用于：对比WT与各突变体的复合物构象稳定性
# ==================================================================

# 1. 定义所有体系名称（WT + 点突变）
set systems {
    Arg106Gly
    Arg106Leu
    Arg106Trp
    Arg111Gly
    Arg133Cys
    Arg133Gly
    Arg133Leu
    Glu143Lys
    WT
}

# 2. 定义重复 run
set runs { run1 run2 run3 }

# 3. 参考帧（初始结构 = 第0帧）
set ref_frame 0

# ===================== 【论文标准】对齐选择集：蛋白骨架 + DNA骨架 =====================
set align_sel "(protein and name N CA C O) or (nucleic and backbone)"

# ===================== 【论文标准】RMSD计算集：复合物整体骨架 =====================
set rmsd_sel "(protein and name N CA C O) or (nucleic and backbone)"

set max_frame 100000
set step 10
set max_output 10000

# ===================== 输出目录 =====================
if {![file exists "./RMSD"]} {
    file mkdir "./RMSD"
    puts "创建输出目录：./RMSD"
}

# ===================== 自动批量计算 =====================
foreach sys $systems {
    foreach r $runs {

        set topfile "./${sys}_com.prmtop"
        set trajfile "./${sys}_${r}.nc"
        set outfile "./RMSD/${sys}_${r}_complex_backbone_rmsd.dat"

        if {![file exists $topfile] || ![file exists $trajfile]} {
            puts "跳过：$sys $r （文件缺失）"
            continue
        }

        mol new $topfile waitfor all
        mol addfile $trajfile waitfor all
        set molid [molinfo top]
        set total_frames [molinfo $molid get numframes]

        set nf $max_frame
        if {$total_frames < $nf} { set nf $total_frames }

        # 对齐原子
        set align_sel_current [atomselect $molid $align_sel]
        set align_sel_ref [atomselect $molid $align_sel frame $ref_frame]

        # RMSD原子
        set rmsd_sel_current [atomselect $molid $rmsd_sel]
        set rmsd_sel_ref [atomselect $molid $rmsd_sel frame $ref_frame]

        # 原子检查
        set natoms_align [$align_sel_current num]
        set natoms_rmsd [$rmsd_sel_current num]
        if {$natoms_align == 0 || $natoms_rmsd == 0} {
            puts "错误：$sys $r 原子选择为空"
            mol delete $molid
            continue
        }

        # 输出
        set out [open $outfile w]
        puts $out "# Frame RMSD (Å)  |  Protein-DNA complex backbone RMSD"
        puts $out "# System: $sys | Run: $r | Aligned by protein+DNA backbone"

        set out_index 1
        for {set f 0} {$f < $nf && $out_index <= $max_output} {incr f $step} {

            $align_sel_current frame $f
            $rmsd_sel_current frame $f

            # 对齐：蛋白+DNA骨架
            $align_sel_current move [measure fit $align_sel_current $align_sel_ref]

            # 计算：复合物骨架RMSD
            set rmsd [measure rmsd $rmsd_sel_current $rmsd_sel_ref]

            puts $out "$out_index $rmsd"
            incr out_index
        }

        close $out
        $align_sel_current delete
        $align_sel_ref delete
        $rmsd_sel_current delete
        $rmsd_sel_ref delete
        mol delete $molid

        puts "✅ 完成：$sys $r | 总帧数：$total_frames | 输出点数：[expr {$out_index-1}]"
    }
}

puts "\n============================================="
puts "  蛋白-DNA复合物骨架RMSD 批量计算完成！"
puts "  适用于：WT + 点突变体系对比分析"
puts "  结果文件：./RMSD/*.dat"
puts "=============================================\n"
