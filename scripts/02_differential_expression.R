# 02_differential_expression.R
# DESeq2 дифференциальный анализ WT vs KO
# Входные данные: data/processed/dds_kidney.rds
# Выходные данные: results/tables/DEG_kidney_WT_vs_KO.csv
source('config/color_palettes.R')
library(DESeq2); library(apeglm)
