# =================================================
# 00_setup.R
# Установка всех необходимых пакетов
# Запускать один раз при первом развёртывании
# =================================================

# CRAN пакеты
cran_pkgs <- c(
  'dplyr', 'tidyr', 'ggplot2', 'patchwork',
  'ggrepel', 'ggbeeswarm', 'ggalluvial',
  'pheatmap', 'RColorBrewer', 'viridis',
  'deSolve', 'yaml', 'readr', 'stringr',
  'scales', 'gridExtra', 'cowplot'
)

# Bioconductor пакеты
bioc_pkgs <- c(
  'DESeq2', 'clusterProfiler', 'fgsea',
  'org.Mm.eg.db', 'org.Hs.eg.db',
  'AnnotationDbi', 'ComplexHeatmap',
  'EnhancedVolcano', 'apeglm', 'BiocParallel'
)

# Установка CRAN
missing_cran <- cran_pkgs[!cran_pkgs %in%
                            installed.packages()[,'Package']]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

# Установка Bioconductor
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
missing_bioc <- bioc_pkgs[!bioc_pkgs %in%
                            installed.packages()[,'Package']]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE)
}

cat('Все пакеты установлены.\n')
