# ============================================================
# VMD RMSD 计算脚本：4 个 Arg133Cys-配体体系 × 3 个 run
# ============================================================
# 计算逻辑：
# 1. 对齐：蛋白质骨架 + DNA 骨架
# 2. 计算：蛋白 + DNA 骨架合并后的 total backbone RMSD
# 3. 每个体系的 run1、run2、run3 分别输出
#
# 对应文件命名规则：
#
# 拓扑文件：
#   Arg133Cys_DB00694.com.prmtop
#   Arg133Cys_DB00445.com.prmtop
#   Arg133Cys_DB06213.com.prmtop
#   Arg133Cys_DB06237.com.prmtop
#
# 轨迹文件：
#   Arg133Cys_<DrugBankID>_run1.nc
#   Arg133Cys_<DrugBankID>_run2.nc
#   Arg133Cys_<DrugBankID>_run3.nc
#
# 若某个体系缺少拓扑文件，则跳过该体系；
# 若某个 run 缺少轨迹文件，则跳过该 run，继续计算下一个。
#
# 输出示例：
# ./RMSD/Arg133Cys_DB00694_run1_total_bb_rmsd.dat
# ============================================================


# ===================== 基本参数 =====================

# 需要计算的体系
set systems {
    Arg133Cys_DB00694
    Arg133Cys_DB00445
    Arg133Cys_DB06213
    Arg133Cys_DB06237
    Arg133Cys_no_ligand
}

# 每个体系检查三个独立重复
set runs {
    run1
    run2
    run3
}

# 参考帧
set ref_frame 0

# 最大计算帧数
set max_frame 4000

# 每隔多少帧计算一次
set step 1


# ===================== 对齐和 RMSD 选择 =====================

# 蛋白质 backbone + DNA backbone
#
# 蛋白质：
# N CA C O
#
# DNA：
# P OP1 OP2 C3' C4' C5' O3' O5'
set align_sel \
"(protein and name N CA C O) or (nucleic and name P OP1 OP2 C3' C4' C5' O3' O5')"


# ===================== 输出目录 =====================

if {![file exists "./RMSD"]} {
    file mkdir "./RMSD"
    puts "已创建输出目录：./RMSD"
}


# ===================== 主循环 =====================

foreach sysname $systems {

    set topfile "${sysname}.com.prmtop"

    puts ""
    puts "============================================================"
    puts "检查体系：$sysname"
    puts "============================================================"

    # 如果拓扑不存在，跳过整个体系
    if {![file exists $topfile]} {
        puts "跳过体系：$sysname"
        puts "原因：找不到拓扑文件 $topfile"
        continue
    }

    foreach r $runs {

        set trajfile "${sysname}_${r}.nc"

        set outfile \
        "./RMSD/${sysname}_${r}_total_bb_rmsd.dat"

        # 如果某个 run 的轨迹不存在，跳过该 run
        if {![file exists $trajfile]} {
            puts ""
            puts "跳过：${sysname}_${r}"
            puts "原因：找不到轨迹文件 $trajfile"
            continue
        }

        puts ""
        puts "============================================================"
        puts "开始计算：${sysname}_${r}"
        puts "拓扑文件：$topfile"
        puts "轨迹文件：$trajfile"
        puts "输出文件：$outfile"
        puts "============================================================"


        # ===================== 加载轨迹 =====================

        mol new $topfile type parm7 waitfor all

        mol addfile $trajfile \
        type netcdf \
        waitfor all

        set molid [molinfo top]

        set total_frames \
        [molinfo $molid get numframes]

        # 实际计算帧数取轨迹总帧数和 max_frame 中较小值
        set nf $max_frame

        if {$total_frames < $nf} {
            set nf $total_frames
        }

        puts "轨迹总帧数：$total_frames"
        puts "实际计算帧数：$nf"


        # ===================== 原子选择 =====================

        set ref \
        [atomselect $molid $align_sel frame $ref_frame]

        set cur \
        [atomselect $molid $align_sel]

        set n_ref [$ref num]
        set n_cur [$cur num]

        puts "RMSD 选择原子数：$n_ref"

        # 若原子选择为空，清理后继续下一个 run
        if {$n_ref == 0} {

            puts "错误：RMSD 原子选择为空"
            puts "请检查 align_sel："
            puts "$align_sel"

            $ref delete
            $cur delete

            mol delete $molid

            continue
        }

        if {$n_ref != $n_cur} {

            puts "警告：参考帧与当前选择原子数不一致"

            puts "ref atoms = $n_ref"
            puts "cur atoms = $n_cur"
        }


        # ===================== 打开输出文件 =====================

        set fout [open $outfile w]

        puts $fout \
        "# Frame Protein_DNA_total_backbone_RMSD_Angstrom"


        # ===================== 逐帧计算 RMSD =====================

        for {set f 0} {$f < $nf} {incr f $step} {

            $cur frame $f
            $cur update

            # 将当前帧拟合到本条轨迹的参考帧
            set trans_mat \
            [measure fit $cur $ref]

            $cur move $trans_mat

            # 计算蛋白质和 DNA 合并骨架 RMSD
            set rmsd \
            [measure rmsd $cur $ref]

            puts $fout "$f $rmsd"
        }


        # ===================== 清理 =====================

        close $fout

        $ref delete
        $cur delete

        mol delete $molid

        puts "完成：${sysname}_${r}"
        puts "输出：$outfile"
    }
}


puts ""
puts "============================================================"
puts "所有可用轨迹的 Protein-DNA total backbone RMSD 计算完成"
puts "缺失的拓扑或轨迹已自动跳过"
puts "结果保存在 ./RMSD/ 文件夹"
puts "============================================================"

exit
