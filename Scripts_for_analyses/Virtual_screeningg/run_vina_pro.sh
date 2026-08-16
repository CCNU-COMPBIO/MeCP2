#!/bin/bash

# 创建统一文件夹
mkdir -p output_all
mkdir -p log_file
mkdir -p failed_list

# 批量处理所有 DB 小分子
# 注意：这里会递归搜索 ../03_ligands_pdbqt_new/ 下面所有子文件夹中的 pdbqt
find ../02_ligands_pdbqt_new -type f -name "*.pdbqt" | sort | while read -r f
do
    drugid=$(basename "$(dirname "$f")")
    b=$(basename "$f" .pdbqt)

    echo "========================================"
    echo "正在处理：${drugid}/${b}"

    # 如果已经跑过，就跳过，方便中断后继续
    if [ -s "output_all/${b}.pdbqt" ] && [ -s "log_file/${b}.txt" ]; then
        echo "已存在，跳过：${b}"
        continue
    fi

    # 运行 vina：参数文件名称仍然使用 config.txt
    /home/lenovo/bin/vina \
    --config config.txt \
    --ligand "$f" \
    --out "output_all/${b}.pdbqt" \
    --spacing 0.375 \
    | tee "log_file/${b}.txt"

    # 检查 vina 是否运行失败
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "FAILED: $f" >> failed_list/failed_docking.txt
        echo "对接失败：$f"
    fi

done

echo -e "\n全部完成！"
echo "结果结构 → output_all/"
echo "运行日志 → log_file/"
echo "失败列表 → failed_list/failed_docking.txt"
