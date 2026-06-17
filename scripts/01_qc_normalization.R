# 01_qc_normalization.R
# QC, фильтрация, DESeq2-нормализация, VST
# Входные данные: data/raw/counts/, data/raw/metadata/
# Выходные данные: data/processed/dds_kidney.rds, vst_kidney.rds
source('config/color_palettes.R')
library(DESeq2)
