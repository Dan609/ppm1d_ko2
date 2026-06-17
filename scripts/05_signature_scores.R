# 05_signature_scores.R
# Комплексный скоринг сигнатур:
# DDRSenescenceScore, InflammagingScore, RhoActivityScore,
# ProxTubuleScore, MacrophageScore, ToxicityPredispositionScore,
# cGAS-STING, PINK1-mitostress
source('config/color_palettes.R')
library(dplyr); library(tidyr); library(ggplot2)
