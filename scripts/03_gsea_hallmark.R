# 03_gsea_hallmark.R
# GSEA по Hallmark MSigDB генсетам (fgsea)
# Входные данные: results/tables/DEG_kidney_WT_vs_KO.csv
# Выходные данные: results/tables/GSEA_hallmark_kidney.csv
source('config/color_palettes.R')
library(fgsea); library(dplyr)
