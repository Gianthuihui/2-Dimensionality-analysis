# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import pairwise_distances
from scipy.stats import pearsonr
from scipy.linalg import eigh


# %%
# Step 1: Load frequency matrix
freq = pd.read_csv("/mnt/f/OneDrive/科研/4_代码/2-Dimensionality-analysis/6-PCoA可视化/data/4级单倍群.csv", index_col=0)


# %%
# Step 2: Calculate Euclidean distance matrix
dist_matrix = pairwise_distances(freq, metric='euclidean')

# Step 3: Classical MDS (PCoA) with eigen decomposition to get variance explained
n = dist_matrix.shape[0]
H = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * H @ (dist_matrix ** 2) @ H

eigvals, eigvecs = eigh(B)
eigvals = eigvals[::-1]
eigvecs = eigvecs[:, ::-1]
coords = eigvecs[:, :2]

# Step 4: Format coordinate data
sample_coords = pd.DataFrame(coords, columns=["PCoA1", "PCoA2"])
sample_coords["Population"] = freq.index

# Variance explained
variance_explained = eigvals / eigvals.sum()
pcoa1_var = round(variance_explained[0] * 100, 2)
pcoa2_var = round(variance_explained[1] * 100, 2)

# 输出前10个主轴解释度（可调）
print("Top 10 axes variance explained (%):")
for i in range(10):
    print(f"PCoA{i+1}: {round(variance_explained[i]*100, 2)}%")

# %%

# Step 5: Merge group info
group_info = pd.read_csv("/mnt/f/OneDrive/科研/4_代码/2-Dimensionality-analysis/6-PCoA可视化/conf/group.csv")
merged_coords = pd.merge(sample_coords, group_info, on="Population", how="inner")

# Step 6: Vector fitting (envfit-like)
haplo_scores = []
for col in freq.columns:
    r1, _ = pearsonr(coords[:, 0], freq[col])
    r2, _ = pearsonr(coords[:, 1], freq[col])
    haplo_scores.append({'Haplogroup': col, 'PCoA1': r1, 'PCoA2': r2})

haplo_scores = pd.DataFrame(haplo_scores)
haplo_scores["VectorLength"] = np.sqrt(haplo_scores["PCoA1"]**2 + haplo_scores["PCoA2"]**2)

# 输出对结构贡献最大的前5个单倍群
top_contributors = haplo_scores.sort_values("VectorLength", ascending=False).head(5)
print("\nTop 5 haplogroup contributors (vector length):")
print(top_contributors[["Haplogroup", "PCoA1", "PCoA2", "VectorLength"]])

# %%
from aquarel import load_theme

# %%
plt.rcParams['font.family'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# Step 7: 绘图（群体-形状，语系-颜色）
plt.figure(figsize=(10, 8))

# 自定义颜色和 marker
marker_styles = ['o', 's', 'D', '^', 'v', '<', '>', 'P', '*', 'X', 'H', 'd', 'p']
color_palette = sns.color_palette("tab10")

# 按语系分组绘制
for i, (lang_group, group_df) in enumerate(merged_coords.groupby('LanguageGroup')):
    lang_color = color_palette[i % len(color_palette)]
    pops = group_df['Population'].unique()

    for j, pop in enumerate(pops):
        plt.rcParams['font.family'] = ['Arial']
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        sub_df = group_df[group_df['Population'] == pop]
        marker = marker_styles[j % len(marker_styles)]
        plt.scatter(
            sub_df["PCoA1"], sub_df["PCoA2"],
            color=lang_color,
            marker=marker,
            s=100,
            label=f"{lang_group}-{pop}"
        )

# 绘制箭头向量
for _, row in haplo_scores.iterrows():
    plt.arrow(0, 0, row["PCoA1"], row["PCoA2"], color='#E57373', alpha=0.6, head_width=0.02)
    plt.text(row["PCoA1"] * 1.1, row["PCoA2"] * 1.1, row["Haplogroup"], color='red', fontsize=8)

# 美化
plt.xlabel(f"PCoA1 ({pcoa1_var}%)")
plt.ylabel(f"PCoA2 ({pcoa2_var}%)")
plt.title("PCoA Biplot: Haplogroup Contribution to Population Structure")
plt.grid(False)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="LanguageGroup-Population", fontsize=8)
plt.tight_layout()
plt.savefig("./1.pdf")
plt.show()

top_contributors.to_csv("/mnt/f/OneDrive/科研/4_代码/2-Dimensionality-analysis/6-PCoA可视化/output/top5_haplogroup_contributors-1.csv", index=False)


