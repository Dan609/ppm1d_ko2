# 07_neutrophil_analysis.R
# Анализ нейтрофилов Ppm1d-KO:
# PCA, volcano, heatmap, MA-plot,
# костимуляторные лиганды (Tnfsf9/4-1BBL, Tnfsf4/OX-40L)
# N1/N2 поляризация
source('config/color_palettes.R')
library(DESeq2); library(ggplot2)
