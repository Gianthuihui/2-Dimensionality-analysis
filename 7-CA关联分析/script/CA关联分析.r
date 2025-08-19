#!/usr/bin/env Rscript

# ==================== 1. 加载所需 R 包 ====================
# — 若未安装，请先在 R 控制台运行 install.packages(c("FactoMineR","factoextra")) —
library(FactoMineR)   # CA 核心计算包
library(factoextra)   # CA 可视化辅助包

# ==================== 2. 【请在此处修改】输入输出路径设置 ====================
# —— 将下面三个变量改为你自己的文件与目录路径 —— #
input_file  <- "/mnt/f/OneDrive/科研/4_代码/2-Dimensionality-analysis/7-CA关联分析/data/CA.csv"
output_dir  <- "/mnt/f/OneDrive/科研/4_代码/2-Dimensionality-analysis/7-CA关联分析/output"
plot_prefix <- "CA_result"   # 输出图文件名前缀，可按需修改

# —— 自动创建输出目录，无需手动建文件夹 —— #
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==================== 3. 读取并检查数据 ====================
# 要求：CSV 第一列为行名，其余各列为计数数据
data_raw <- read.csv(input_file, header = TRUE, row.names = 1, check.names = FALSE)
cat(">>> 读取数据：", input_file, "\n")
cat(">>> 数据维度（行 × 列）：", dim(data_raw), "\n")
# 如需查看前几行，可取消下一行注释：
# print(head(data_raw))

# ==================== 4. 卡方检验（Chi-square test） ====================
chi <- chisq.test(data_raw)
cat(">>> Chi-square statistic =", chi$statistic, "\n")
cat(">>> df =", chi$parameter, "\n")
cat(">>> p-value =", chi$p.value, "\n")
chi_df <- data.frame(
  statistic = unname(chi$statistic),
  df        = unname(chi$parameter),
  p_value   = chi$p.value
)
write.csv(chi_df,
          file = file.path(output_dir, paste0(plot_prefix, "_chi_square.csv")),
          row.names = FALSE,
          quote = FALSE)

# ==================== 5. 运行 Correspondence Analysis ====================
# graph = FALSE 表示不自动绘图，后续统一用 factoextra
res.ca <- CA(data_raw, graph = FALSE)

# ==================== 6. 导出特征值 & 方差贡献（Eigenvalues & Variance） ====================
eig_df <- data.frame(
  Dimension        = seq_len(nrow(res.ca$eig)),
  Eigenvalue       = res.ca$eig[, 1],
  VariancePercent  = res.ca$eig[, 2],
  CumVariancePercent = res.ca$eig[, 3]
)
write.csv(eig_df,
          file = file.path(output_dir, paste0(plot_prefix, "_eigenvalues.csv")),
          row.names = FALSE,
          quote = FALSE)

# ==================== 7. 导出 CA 数值结果 ====================
# 7.1 行、列坐标
write.csv(res.ca$row$coord,
          file = file.path(output_dir, paste0(plot_prefix, "_row_coords.csv")),
          quote = FALSE)
write.csv(res.ca$col$coord,
          file = file.path(output_dir, paste0(plot_prefix, "_col_coords.csv")),
          quote = FALSE)

# 7.2 贡献度 (contrib) 与 质量 (cos2)
write.csv(res.ca$row$contrib,
          file = file.path(output_dir, paste0(plot_prefix, "_row_contrib.csv")),
          quote = FALSE)
write.csv(res.ca$col$contrib,
          file = file.path(output_dir, paste0(plot_prefix, "_col_contrib.csv")),
          quote = FALSE)
write.csv(res.ca$row$cos2,
          file = file.path(output_dir, paste0(plot_prefix, "_row_cos2.csv")),
          quote = FALSE)
write.csv(res.ca$col$cos2,
          file = file.path(output_dir, paste0(plot_prefix, "_col_cos2.csv")),
          quote = FALSE)

# ==================== 8. 绘制并保存图形 ====================
# （1）CA 双标图（行 + 列）
pdf(file = file.path(output_dir, paste0(plot_prefix, "_biplot.pdf")),
    width = 8, height = 6)
print(fviz_ca_biplot(res.ca,
                    repel = TRUE,
                    title = "Correspondence Analysis Biplot",
                    axes  = c(1, 2)))
dev.off()

# （2）仅行点图
png(filename = file.path(output_dir, paste0(plot_prefix, "_rows.png")),
    width = 800, height = 600)
print(fviz_ca_row(res.ca,
                 repel = TRUE,
                 title = "CA Row Plot",
                 axes  = c(1, 2)))
dev.off()

# （3）仅列点图
png(filename = file.path(output_dir, paste0(plot_prefix, "_cols.png")),
    width = 800, height = 600)
print(fviz_ca_col(res.ca,
                 repel = TRUE,
                 title = "CA Column Plot",
                 axes  = c(1, 2)))
dev.off()

# ==================== 9. 完成提示 ====================
cat("✅ Correspondence Analysis 完成，所有结果保存在：", output_dir, "\n")
