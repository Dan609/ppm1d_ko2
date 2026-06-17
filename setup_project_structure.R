# =============================================================
# setup_project_structure.R
#
# ЗАПУСК: из корневой папки проекта (ppm1d_ko_transcriptomics/)
#   source("setup_project_structure.R")
# или из терминала:
#   Rscript setup_project_structure.R
#
# Создаёт всю структуру папок, .gitignore, README.md,
# шаблоны конфигов и Rproj-файл.
# =============================================================

cat("=== Развёртывание структуры проекта PPM1D-KO ===\n\n")

# ----------------------------------------------------------
# 1. СПИСОК ВСЕХ ПАПОК
# ----------------------------------------------------------
dirs <- c(
  # Data
  "data/raw/counts",
  "data/raw/metadata",
  "data/raw/genesets",
  "data/processed",
  
  # Scripts
  "scripts/modeling",
  
  # Results
  "results/tables",
  "results/figures/main",
  "results/figures/supplementary",
  "results/figures/modeling",
  "results/modeling/calibration_output",
  "results/modeling/bifurcation_output",
  
  # Manuscript
  "manuscript/main",
  "manuscript/sections",
  "manuscript/figures_legends",
  "manuscript/references",
  "manuscript/submission/v1_submission",
  "manuscript/submission/v2_revision",
  
  # Notebooks
  "notebooks",
  
  # Config
  "config",
  
  # Tmp
  "tmp/scratch",
  "tmp/old_versions"
)

# Создать все папки
created <- 0
skipped <- 0
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  [+] создана: %s\n", d))
    created <- created + 1
  } else {
    cat(sprintf("  [=] уже есть: %s\n", d))
    skipped <- skipped + 1
  }
}
cat(sprintf("\nПапок создано: %d  |  уже существовало: %d\n\n",
            created, skipped))

# ----------------------------------------------------------
# 2. СОЗДАТЬ PLACEHOLDER-ФАЙЛЫ .gitkeep В ПУСТЫХ ПАПКАХ
#    (чтобы git отслеживал пустые директории)
# ----------------------------------------------------------
gitkeep_dirs <- c(
  "data/raw/counts",
  "data/raw/metadata",
  "data/raw/genesets",
  "results/tables",
  "results/figures/main",
  "results/figures/supplementary",
  "results/figures/modeling",
  "results/modeling/calibration_output",
  "results/modeling/bifurcation_output",
  "manuscript/submission/v1_submission",
  "manuscript/submission/v2_revision",
  "tmp/scratch",
  "tmp/old_versions"
)
for (d in gitkeep_dirs) {
  gk <- file.path(d, ".gitkeep")
  if (!file.exists(gk)) {
    file.create(gk)
  }
}

# ----------------------------------------------------------
# 3. .gitignore
# ----------------------------------------------------------
gitignore_content <- c(
  "# ===== PPM1D-KO Project .gitignore =====",
  "",
  "# R-специфичное",
  ".Rhistory",
  ".Rapp.history",
  ".RData",
  ".Ruserdata",
  "*.Rproj.user/",
  "",
  "# Обработанные данные (генерируются скриптами)",
  "data/processed/",
  "",
  "# Результаты (генерируются скриптами, тяжёлые)",
  "results/",
  "",
  "# Временные файлы",
  "tmp/",
  "",
  "# Большие бинарные файлы данных (хранить в data/raw локально)",
  "*.rds",
  "*.RData",
  "data/raw/counts/*.csv",
  "data/raw/metadata/*.csv",
  "",
  "# НО отслеживать genesets и конфиги (маленькие)",
  "!data/raw/genesets/",
  "!data/raw/genesets/*",
  "",
  "# Системные",
  ".DS_Store",
  "Thumbs.db",
  "desktop.ini",
  "",
  "# Билды рукописи",
  "manuscript/**/*.pdf",
  "manuscript/**/*.log",
  "manuscript/**/*.aux",
  "",
  "# Notebooks outputs",
  "notebooks/*_files/",
  "notebooks/*.html",
  "",
  "# Temp",
  "tmp/",
  "*.tmp",
  "*.bak",
  "~$*"
)
writeLines(gitignore_content, ".gitignore")
cat("[+] Создан .gitignore\n")

# ----------------------------------------------------------
# 4. README.md
# ----------------------------------------------------------
readme_content <- c(
  "# PPM1D/Wip1 Knockout Transcriptomics Project",
  "",
  "## Описание",
  "Bulk RNA-seq анализ нокаутных по PPM1D/Wip1 почек мышей",
  "(C57BL/6, BD27 WT vs BD28/BD29 KO).",
  "Включает биоинформатический пайплайн, математическое моделирование",
  "Rho-GTPase переключателя и черновик рукописи.",
  "",
  "## Структура",
  "```",
  "data/        — исходные данные (counts, metadata, genesets)",
  "scripts/     — R-скрипты анализа и моделирования",
  "results/     — таблицы и рисунки (генерируются, не в git)",
  "manuscript/  — черновик статьи (md/docx разделы)",
  "notebooks/   — Rmd для exploratory анализа",
  "config/      — параметры запуска",
  "tmp/         — временные файлы (не в git)",
  "```",
  "",
  "## Быстрый старт",
  "```r",
  "# 1. Установить зависимости",
  "source('scripts/00_setup.R')",
  "",
  "# 2. Нормализация",
  "source('scripts/01_qc_normalization.R')",
  "",
  "# 3. DESeq2",
  "source('scripts/02_differential_expression.R')",
  "",
  "# 4. Моделирование Rho-переключателя",
  "source('scripts/modeling/run_calibration_pipeline.R')",
  "source('scripts/modeling/rho_bifurcation.R')",
  "results <- run_full_bifurcation_analysis()",
  "```",
  "",
  "## Ключевые выводы",
  "- Wip1KO-почка: «DDR-primed, pre-injured» транскриптомное состояние",
  "- Повышены: Trp53/p21/p16, SASP-цитокины, Rho-сеть",
  "- Математическая модель: KO -> RhoA-dominant аттрактор",
  "- Механизм цисплатин-летальности: наложение стресса на pre-injured фон",
  "",
  "## Данные",
  "- rawcount272829.csv — матрица сырых счётов (12 образцов)",
  "- 4 WT (BD27) + 8 KO (BD28, BD29)",
  "",
  "## Авторы",
  "Bobkov, Bogdanova et al.",
  "",
  "## Статус",
  "🟡 In preparation — результаты получены, рукопись в работе"
)
writeLines(readme_content, "README.md")
cat("[+] Создан README.md\n")

# ----------------------------------------------------------
# 5. Rproj файл
# ----------------------------------------------------------
rproj_name <- basename(getwd())
rproj_content <- c(
  "Version: 1.0",
  "",
  "RestoreWorkspace: No",
  "SaveWorkspace: No",
  "AlwaysSaveHistory: Default",
  "",
  "EnableCodeIndexing: Yes",
  "UseSpacesForTab: Yes",
  "NumSpacesForTab: 2",
  "Encoding: UTF-8",
  "",
  "RnwWeave: Sweave",
  "LaTeX: pdfLaTeX",
  "",
  "AutoAppendNewline: Yes",
  "StripTrailingWhitespace: Yes",
  "",
  "BuildType: Package"
)
writeLines(rproj_content,
           paste0(rproj_name, ".Rproj"))
cat(sprintf("[+] Создан %s.Rproj\n", rproj_name))

# ----------------------------------------------------------
# 6. config/params_kidney.yaml
# ----------------------------------------------------------
yaml_kidney <- c(
  "# Параметры анализа: почки Wip1KO",
  "model_type: kidney_WipKO",
  "",
  "data:",
  "  counts:   data/raw/counts/rawcount272829.csv",
  "  metadata: data/raw/metadata/metadata_kidney.csv",
  "",
  "groups:",
  "  control:   WT",
  "  treatment: KO",
  "",
  "deseq2:",
  "  fdr_threshold:    0.05",
  "  log2fc_threshold: 1.0",
  "  shrinkage:        apeglm",
  "",
  "gsea:",
  "  gmt_hallmark: data/raw/genesets/hallmark_msigdb.gmt",
  "  nperm:        1000",
  "  min_gs_size:  15",
  "  max_gs_size:  500",
  "",
  "modeling:",
  "  n_ic_calibration: 20",
  "  n_ic_full:        50",
  "  seed:             2025",
  "  out_dir:          results/modeling/calibration_output",
  "",
  "output:",
  "  tables:  results/tables",
  "  figures: results/figures"
)
writeLines(yaml_kidney, "config/params_kidney.yaml")
cat("[+] Создан config/params_kidney.yaml\n")

# ----------------------------------------------------------
# 7. config/params_MC38.yaml
# ----------------------------------------------------------
yaml_mc38 <- c(
  "# Параметры анализа: MC38 PPM1D-KO",
  "model_type: MC38_PPM1DKO",
  "",
  "data:",
  "  counts:   data/raw/counts/MC38_rawcounts.csv",
  "  metadata: data/raw/metadata/metadata_MC38.csv",
  "",
  "groups:",
  "  control:   WT",
  "  treatment: KO",
  "",
  "deseq2:",
  "  fdr_threshold:    0.05",
  "  log2fc_threshold: 1.0",
  "  shrinkage:        apeglm",
  "",
  "gsea:",
  "  gmt_hallmark: data/raw/genesets/hallmark_msigdb.gmt",
  "  nperm:        1000",
  "  min_gs_size:  15",
  "  max_gs_size:  500",
  "",
  "modeling:",
  "  n_ic_calibration: 20",
  "  n_ic_full:        50",
  "  seed:             2025",
  "  out_dir:          results/modeling/calibration_output_MC38",
  "",
  "output:",
  "  tables:  results/tables",
  "  figures: results/figures"
)
writeLines(yaml_mc38, "config/params_MC38.yaml")
cat("[+] Создан config/params_MC38.yaml\n")

# ----------------------------------------------------------
# 8. config/params_HT29.yaml
# ----------------------------------------------------------
yaml_ht29 <- c(
  "# Параметры анализа: HT29 PPM1D-KO",
  "model_type: HT29_PPM1DKO",
  "",
  "data:",
  "  counts:   data/raw/counts/HT29_rawcounts.csv",
  "  metadata: data/raw/metadata/metadata_HT29.csv",
  "",
  "groups:",
  "  control:   WT",
  "  treatment: KO",
  "",
  "deseq2:",
  "  fdr_threshold:    0.05",
  "  log2fc_threshold: 1.0",
  "  shrinkage:        apeglm",
  "",
  "gsea:",
  "  gmt_hallmark: data/raw/genesets/hallmark_msigdb.gmt",
  "  nperm:        1000",
  "  min_gs_size:  15",
  "  max_gs_size:  500",
  "",
  "modeling:",
  "  n_ic_calibration: 20",
  "  n_ic_full:        50",
  "  seed:             2025",
  "  out_dir:          results/modeling/calibration_output_HT29",
  "",
  "output:",
  "  tables:  results/tables",
  "  figures: results/figures"
)
writeLines(yaml_ht29, "config/params_HT29.yaml")
cat("[+] Создан config/params_HT29.yaml\n")

# ----------------------------------------------------------
# 9. config/color_palettes.R
# ----------------------------------------------------------
color_palettes <- c(
  "# =================================================",
  "# color_palettes.R",
  "# Единые цветовые палитры для всего проекта",
  "# source('config/color_palettes.R') в начале каждого скрипта",
  "# =================================================",
  "",
  "COLORS <- list(",
  "",
  "  # Генотипы",
  "  genotype = c(",
  "    WT = '#2471a3',",
  "    KO = '#c0392b'",
  "  ),",
  "",
  "  # GTPases",
  "  gtpase = c(",
  "    RhoA  = '#e74c3c',",
  "    Rac1  = '#2980b9',",
  "    Cdc42 = '#27ae60'",
  "  ),",
  "",
  "  # Пути (GSEA / кастомные генсеты)",
  "  pathways = c(",
  "    DDR_Senescence  = '#8e44ad',",
  "    Inflammation    = '#e67e22',",
  "    Rho_GTPases     = '#e74c3c',",
  "    ProxTubule      = '#16a085',",
  "    Macrophage      = '#f39c12',",
  "    Fibroblast      = '#7f8c8d',",
  "    cGAS_STING      = '#2c3e50',",
  "    PINK1_Mitostress = '#d35400'",
  "  ),",
  "",
  "  # Volcano plot",
  "  volcano = c(",
  "    up       = '#c0392b',",
  "    down     = '#2471a3',",
  "    ns       = '#bdc3c7'",
  "  ),",
  "",
  "  # Heatmap (diverging)",
  "  heatmap_low  = '#2471a3',",
  "  heatmap_mid  = '#ffffff',",
  "  heatmap_high = '#c0392b',",
  "",
  "  # Нейтральные",
  "  grey_light = '#ecf0f1',",
  "  grey_mid   = '#95a5a6',",
  "  grey_dark  = '#2c3e50'",
  ")",
  "",
  "# Тема ggplot2 для всего проекта",
  "theme_ppm1d <- function(base_size = 12) {",
  "  theme_bw(base_size = base_size) %+replace%",
  "    theme(",
  "      panel.grid.minor  = element_blank(),",
  "      strip.background  = element_rect(fill = '#f4f4f4',",
  "                                        color = 'grey80'),",
  "      legend.background = element_blank(),",
  "      legend.key        = element_blank(),",
  "      plot.title        = element_text(face = 'bold',",
  "                                        size = base_size + 1),",
  "      plot.subtitle     = element_text(color = 'grey40',",
  "                                        size = base_size - 1)",
  "    )",
  "}",
  "",
  "cat('  [config] Палитры и тема загружены\\n')"
)
writeLines(color_palettes, "config/color_palettes.R")
cat("[+] Создан config/color_palettes.R\n")

# ----------------------------------------------------------
# 10. scripts/00_setup.R
# ----------------------------------------------------------
setup_script <- c(
  "# =================================================",
  "# 00_setup.R",
  "# Установка всех необходимых пакетов",
  "# Запускать один раз при первом развёртывании",
  "# =================================================",
  "",
  "# CRAN пакеты",
  "cran_pkgs <- c(",
  "  'dplyr', 'tidyr', 'ggplot2', 'patchwork',",
  "  'ggrepel', 'ggbeeswarm', 'ggalluvial',",
  "  'pheatmap', 'RColorBrewer', 'viridis',",
  "  'deSolve', 'yaml', 'readr', 'stringr',",
  "  'scales', 'gridExtra', 'cowplot'",
  ")",
  "",
  "# Bioconductor пакеты",
  "bioc_pkgs <- c(",
  "  'DESeq2', 'clusterProfiler', 'fgsea',",
  "  'org.Mm.eg.db', 'org.Hs.eg.db',",
  "  'AnnotationDbi', 'ComplexHeatmap',",
  "  'EnhancedVolcano', 'apeglm', 'BiocParallel'",
  ")",
  "",
  "# Установка CRAN",
  "missing_cran <- cran_pkgs[!cran_pkgs %in%",
  "                            installed.packages()[,'Package']]",
  "if (length(missing_cran) > 0) {",
  "  install.packages(missing_cran, dependencies = TRUE)",
  "}",
  "",
  "# Установка Bioconductor",
  "if (!requireNamespace('BiocManager', quietly = TRUE))",
  "  install.packages('BiocManager')",
  "missing_bioc <- bioc_pkgs[!bioc_pkgs %in%",
  "                            installed.packages()[,'Package']]",
  "if (length(missing_bioc) > 0) {",
  "  BiocManager::install(missing_bioc, update = FALSE)",
  "}",
  "",
  "cat('Все пакеты установлены.\\n')"
)
writeLines(setup_script, "scripts/00_setup.R")
cat("[+] Создан scripts/00_setup.R\n")

# ----------------------------------------------------------
# 11. manuscript/sections/ — шаблоны разделов
# ----------------------------------------------------------
sections <- list(
  
  "manuscript/sections/01_introduction.md" = c(
    "# Introduction",
    "",
    "<!-- STATUS: draft -->",
    "",
    "PPM1D (Wip1) is a serine/threonine phosphatase...",
    "",
    "## Key points to cover",
    "- PPM1D as negative regulator of p53/ATM/CHK1/p38MAPK",
    "- Role in DDR, senescence, aging",
    "- Cisplatin nephrotoxicity — macrophage-mediated mechanism",
    "- Gap: transcriptomic reprogramming of Wip1KO kidney",
    "- Hypothesis: 'DDR-primed pre-injured' state"
  ),
  
  "manuscript/sections/02_methods.md" = c(
    "# Materials and Methods",
    "",
    "<!-- STATUS: complete (2.1-2.4) -->",
    "",
    "## 2.1 Animals",
    "",
    "## 2.2 Cisplatin treatment protocol",
    "",
    "## 2.3 RNA isolation and sequencing",
    "",
    "## 2.4 Bioinformatics pipeline",
    "### 2.4.1 Quality control and normalization",
    "### 2.4.2 Differential expression analysis",
    "### 2.4.3 Gene set enrichment analysis",
    "### 2.4.4 Signature scoring",
    "### 2.4.5 Mathematical modeling of Rho-GTPase switch"
  ),
  
  "manuscript/sections/03_results.md" = c(
    "# Results",
    "",
    "<!-- STATUS: in progress (3.1-3.4 drafted) -->",
    "",
    "## 3.1 Wip1KO kidney exhibits a DDR-primed transcriptome",
    "",
    "## 3.2 Elevated senescence and SASP signatures in Wip1KO",
    "",
    "## 3.3 Rho GTPase network is dysregulated in Wip1KO kidney",
    "",
    "## 3.4 Metabolic reprogramming of proximal tubules",
    "",
    "## 3.5 Mathematical modeling predicts RhoA-dominant attractor",
    "<!-- TODO: add modeling results -->",
    "",
    "## 3.6 Transcriptomic basis of cisplatin hypersensitivity",
    "<!-- TODO -->"
  ),
  
  "manuscript/sections/04_discussion.md" = c(
    "# Discussion",
    "",
    "<!-- STATUS: outline only -->",
    "",
    "## Key points",
    "- DDR-Rho-inflammaging axis as unified mechanism",
    "- Tubuloinflammaging phenotype resembles accelerated aging",
    "- Macrophage polarization as mediator of cisplatin toxicity",
    "- Mathematical model: bistability and KO attractor",
    "- Implications for clinical PPM1D inhibition (oncology)"
  ),
  
  "manuscript/sections/05_abstract.md" = c(
    "# Abstract",
    "",
    "<!-- STATUS: draft -->",
    "",
    "**Background:**",
    "",
    "**Methods:**",
    "",
    "**Results:**",
    "",
    "**Conclusions:**",
    "",
    "**Keywords:** PPM1D, Wip1, kidney injury, cisplatin,",
    "transcriptomics, Rho GTPase, senescence, DDR"
  ),
  
  "manuscript/figures_legends/figure_legends.md" = c(
    "# Figure Legends",
    "",
    "<!-- Подписи к рисункам для рукописи -->",
    "",
    "## Figure 1. Global transcriptomic changes in Wip1KO kidney",
    "**A.** PCA plot... **B.** Volcano plot...",
    "",
    "## Figure 2. Senescence and SASP signatures",
    "",
    "## Figure 3. Rho GTPase network dysregulation",
    "",
    "## Figure 4. Mathematical modeling of RhoA/Rac1 switch",
    "",
    "## Figure 5. Cisplatin hypersensitivity mechanism",
    "",
    "---",
    "## Supplementary Figures",
    "",
    "## Figure S1. Quality control metrics"
  ),
  
  "manuscript/references/ppm1d_papers_notes.md" = c(
    "# Literature Notes — PPM1D Project",
    "",
    "## Key papers",
    "",
    "### PPM1D/Wip1 biology",
    "- Fiscella et al. 1997 — PPM1D discovery",
    "- Bulavin et al. 2002 — Wip1 as p53 phosphatase",
    "- Lindqvist et al. 2009 — Wip1 in DDR recovery",
    "",
    "### Cisplatin nephrotoxicity",
    "- Linkermann et al. 2014 — necroptosis in AKI",
    "- Sharp et al. 2019 — macrophage role",
    "",
    "### Rho GTPase in kidney",
    "- TODO: add refs",
    "",
    "### Senescence & SASP",
    "- Coppé et al. 2008 — SASP definition",
    "- TODO: add refs"
  )
)

for (fpath in names(sections)) {
  writeLines(sections[[fpath]], fpath)
  cat(sprintf("[+] Создан %s\n", fpath))
}

# ----------------------------------------------------------
# 12. notebooks/exploratory_kidney.Rmd
# ----------------------------------------------------------
rmd_content <- c(
  "---",
  "title: 'Exploratory Analysis: Wip1KO Kidney RNA-seq'",
  "author: 'PPM1D Project'",
  "date: '`r Sys.Date()`'",
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_float: true",
  "    code_folding: hide",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = TRUE, warning = FALSE,",
  "                       message = FALSE, fig.width = 8)",
  "source('config/color_palettes.R')",
  "```",
  "",
  "## Данные",
  "```{r load}",
  "# counts <- read.csv('data/raw/counts/rawcount272829.csv',",
  "#                     row.names = 1)",
  "# metadata <- read.csv('data/raw/metadata/metadata_kidney.csv')",
  "```",
  "",
  "## PCA",
  "```{r pca}",
  "# dds <- readRDS('data/processed/dds_kidney.rds')",
  "# plotPCA(vst(dds), intgroup = 'genotype')",
  "```",
  "",
  "## Volcano",
  "```{r volcano}",
  "# res <- readRDS('data/processed/deseq2_results_kidney.rds')",
  "```"
)
writeLines(rmd_content, "notebooks/exploratory_kidney.Rmd")
cat("[+] Создан notebooks/exploratory_kidney.Rmd\n")

# ----------------------------------------------------------
# 13. Скелеты R-скриптов в scripts/
# ----------------------------------------------------------
script_headers <- list(
  
  "scripts/01_qc_normalization.R" = c(
    "# 01_qc_normalization.R",
    "# QC, фильтрация, DESeq2-нормализация, VST",
    "# Входные данные: data/raw/counts/, data/raw/metadata/",
    "# Выходные данные: data/processed/dds_kidney.rds, vst_kidney.rds",
    "source('config/color_palettes.R')",
    "library(DESeq2)"
  ),
  
  "scripts/02_differential_expression.R" = c(
    "# 02_differential_expression.R",
    "# DESeq2 дифференциальный анализ WT vs KO",
    "# Входные данные: data/processed/dds_kidney.rds",
    "# Выходные данные: results/tables/DEG_kidney_WT_vs_KO.csv",
    "source('config/color_palettes.R')",
    "library(DESeq2); library(apeglm)"
  ),
  
  "scripts/03_gsea_hallmark.R" = c(
    "# 03_gsea_hallmark.R",
    "# GSEA по Hallmark MSigDB генсетам (fgsea)",
    "# Входные данные: results/tables/DEG_kidney_WT_vs_KO.csv",
    "# Выходные данные: results/tables/GSEA_hallmark_kidney.csv",
    "source('config/color_palettes.R')",
    "library(fgsea); library(dplyr)"
  ),
  
  "scripts/04_custom_geneset_scoring.R" = c(
    "# 04_custom_geneset_scoring.R",
    "# Скоринг кастомных генсетов:",
    "# Inflammation, Senescence, Small GTPases (Rho/Rac/Cdc42)",
    "source('config/color_palettes.R')",
    "library(dplyr); library(ggplot2)"
  ),
  
  "scripts/05_signature_scores.R" = c(
    "# 05_signature_scores.R",
    "# Комплексный скоринг сигнатур:",
    "# DDRSenescenceScore, InflammagingScore, RhoActivityScore,",
    "# ProxTubuleScore, MacrophageScore, ToxicityPredispositionScore,",
    "# cGAS-STING, PINK1-mitostress",
    "source('config/color_palettes.R')",
    "library(dplyr); library(tidyr); library(ggplot2)"
  ),
  
  "scripts/06_visualization_main.R" = c(
    "# 06_visualization_main.R",
    "# Все основные визуализации для статьи:",
    "# volcano, heatmap top50, circular heatmap,",
    "# barplot генсетов, корреляционная матрица,",
    "# alluvial, chord diagram ген-путь",
    "source('config/color_palettes.R')",
    "library(ggplot2); library(ComplexHeatmap)",
    "library(circlize); library(ggalluvial)"
  ),
  
  "scripts/07_neutrophil_analysis.R" = c(
    "# 07_neutrophil_analysis.R",
    "# Анализ нейтрофилов Ppm1d-KO:",
    "# PCA, volcano, heatmap, MA-plot,",
    "# костимуляторные лиганды (Tnfsf9/4-1BBL, Tnfsf4/OX-40L)",
    "# N1/N2 поляризация",
    "source('config/color_palettes.R')",
    "library(DESeq2); library(ggplot2)"
  ),
  
  "scripts/08_macrophage_analysis.R" = c(
    "# 08_macrophage_analysis.R",
    "# Анализ поляризации макрофагов:",
    "# M1/M2 маркёры, связь с цисплатин-нефротоксичностью",
    "source('config/color_palettes.R')",
    "library(dplyr); library(ggplot2)"
  )
)

for (fpath in names(script_headers)) {
  if (!file.exists(fpath)) {
    writeLines(script_headers[[fpath]], fpath)
    cat(sprintf("[+] Создан %s\n", fpath))
  }
}

# ----------------------------------------------------------
# ИТОГОВЫЙ ОТЧЁТ
# ----------------------------------------------------------
cat("\n")
cat("==============================================\n")
cat("  СТРУКТУРА ПРОЕКТА УСПЕШНО РАЗВЁРНУТА\n")
cat("==============================================\n")
cat(sprintf("  Рабочая директория: %s\n", getwd()))
cat("\n")
cat("  Следующие шаги:\n")
cat("  1. Поместите rawcount272829.csv в data/raw/counts/\n")
cat("  2. Поместите metadata_kidney.csv в data/raw/metadata/\n")
cat("  3. Инициализируйте git:\n")
cat("       git init\n")
cat("       git add .\n")
cat("       git commit -m 'Initial project structure'\n")
cat("  4. Подключите GitHub remote:\n")
cat("       git remote add origin <URL>\n")
cat("       git push -u origin main\n")
cat("  5. Запустите: source('scripts/00_setup.R')\n")
cat("==============================================\n")

