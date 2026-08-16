#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ===================== 1. 数据 =====================
data = {
    "Drug": [
        "DB00266", "DB00445", "DB00694", "DB06237", "DB08974"
    ],
    "ΔG_complex": [
        -10.69, -11.60, -10.34, -10.34, -10.48
    ],
    "ΔG_protein": [
        -6.755, -6.673, -6.377, -5.968, -6.345
    ]
}

df = pd.DataFrame(data)

# ===================== 2. 全局参数 =====================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# ===================== 3. 配色 =====================
complex_color = "#4DAF4A"      # 绿色
protein_color = "#F0E442"      # 黄色
edge_color = "black"

# ===================== 4. 创建画布 =====================
fig, ax = plt.subplots(figsize=(7,5), dpi=600)

# ===================== 5. 柱位置 =====================
x = np.arange(len(df))
width = 0.36

ax.bar(
    x - width/2,
    df["ΔG_complex"],
    width=width,
    color=complex_color,
    edgecolor=edge_color,
    linewidth=1,
    label=r'$\Delta G_{complex}$'
)

ax.bar(
    x + width/2,
    df["ΔG_protein"],
    width=width,
    color=protein_color,
    edgecolor=edge_color,
    linewidth=1,
    label=r'$\Delta G_{protein}$'
)

# ===================== 6. 坐标轴 =====================
ax.set_xticks(x)
ax.set_xticklabels(df["Drug"], rotation=0, ha="center", fontweight="bold")

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)


ax.set_ylabel("AutoDock Vina binding affinity (kcal/mol)",
              fontsize=13,
              fontweight="bold")

ax.set_xlabel("Drug molecules",
              fontsize=13,
              fontweight="bold")

# 因为都是负值
ax.set_ylim(-13, 0)

# 关闭图内网格线，仅保留外部黑色边框
ax.grid(False)

# 刻度
ax.tick_params(axis='both',
               labelsize=11,
               width=1.2,
               length=5)

for label in ax.get_xticklabels()+ax.get_yticklabels():
    label.set_fontweight("bold")

# ===================== 7. 图例 =====================
legend = ax.legend(
    frameon=False,
    fontsize=11,
    loc="lower right",
    bbox_to_anchor=(0.98, 0.01)   # 向下移动
)

for t in legend.get_texts():
    t.set_fontweight("bold")

# ===================== 8. 去掉顶部右侧边框 =====================
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

# ===================== 9. 保存 =====================
plt.tight_layout()

output="AutoDock_Vina_binding_affinity_5_drugs_black_frame_no_internal_lines"

plt.savefig(output+".tiff",
            dpi=600,
            bbox_inches="tight")

plt.savefig(output+".pdf",
            bbox_inches="tight")

plt.savefig(output+".png",
            dpi=600,
            bbox_inches="tight")

plt.show()
plt.close()

print("Finished!")
