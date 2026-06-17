# =============================================================================
# 02_heatmaps.R
# PROJECT: ppm1d_ko2 — Signature scoring + heatmaps
#
# ЗАПУСКАЕТСЯ ОТДЕЛЬНО ДЛЯ КАЖДОГО ДАТАСЕТА:
#   Rscript 02_heatmaps.R mouse
#   Rscript 02_heatmaps.R human
#
# Если аргумент не передан — запускает оба последовательно.
#
# Что делает:
#   1. Загружает pipeline_mouse.rds или pipeline_human.rds
#   2. Строит VST-матрицу из norm_counts
#   3. Рассчитывает z-score сигнатурных скоров
#   4. HM1: обзорный хитмап скоров (сигнатуры × образцы)
#   5. HM2: ген-уровневые хитмапы для 4 приоритетных сигнатур
#   6. HM3: корреляционная матрица сигнатур
#   7. Beeswarm-панели ключевых скоров
#   8. Сохраняет vst_mat + sig_scores в RDS
#
# НОВОЕ vs старый скрипт:
#   - PPM1D_KO вместо PPM1A_KO
#   - Два независимых генсета: mouse (lowercase) и human (UPPERCASE)
#   - Единая функция run_for_dataset() — один код для обоих видов
#   - Фигуры в output/figures/mouse/ и output/figures/human/
#   - pipeline_mouse.rds / pipeline_human.rds обновляются по отдельности
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
  library(ggbeeswarm)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

# ── ГЛОБАЛЬНЫЕ ПАРАМЕТРЫ ─────────────────────────────────────────────────────
FDR        <- 0.05
LFC_THRESH <- 1

# ── Цветовые палитры (общие для обоих видов) ──────────────────────────────────
GENO_COLORS <- c(
  "WT"       = "#4393C3",
  "PPM1D_KO" = "#D6604D",   # NEW: было PPM1A_KO
  "PPM1B_KO" = "#74C476",
  "DKO"      = "#9970AB"
)

# Для меток с TNF (только MC-38)
GENO_COLORS_TNF <- c(
  GENO_COLORS,
  "WT\n+TNF"       = "#92C5DE",
  "PPM1D_KO\n+TNF" = "#F4A582",   # NEW: было PPM1A_KO
  "PPM1B_KO\n+TNF" = "#A1D99B"
)

# Категории сигнатур (общие — имена одинаковые, гены разные по виду)
SIG_CATEGORY <- c(
  "Aging_DNAmAge"     = "Aging clock",
  "Aging_LongevityUp" = "Aging clock",
  "Aging_LongevityDn" = "Aging clock",
  "Senescence_SASP"   = "Senescence",
  "Senescence_DDR"    = "Senescence",
  "Inflammaging"      = "Inflammation",
  "cGAS_STING"        = "Inflammation",
  "Rho_Activity"      = "Cytoskeleton",
  "Mito_Stress"       = "Mitochondria"
)

CAT_COLORS <- c(
  "Aging clock"  = "#8DA0CB",
  "Senescence"   = "#E78AC3",
  "Inflammation" = "#FC8D62",
  "Cytoskeleton" = "#66C2A5",
  "Mitochondria" = "#FFD92F"
)

# =============================================================================
# ГЕНСЕТЫ — МЫШЬ (символы нижнего регистра, Mus musculus)
# =============================================================================
SIG_LIST_MOUSE <- list(

  Aging_DNAmAge = c(
    "Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3",
    "Hdac1","Hdac2","Sirt1","Sirt6","Cbx5","Cbx7",
    "Bmi1","Ezh2","Kdm6a","Kdm5c"
  ),

  Aging_LongevityUp = c(
    "Hspa5","Hsp90ab1","Hspa8","Hspa1a","Dnajb1",
    "Foxo3","Sirt1","Sirt3","Sirt6","Ppargc1a",
    "Nfe2l2","Sqstm1","Map1lc3b","Becn1","Atg5",
    "Igfbp1","Igfbp3","Eif4ebp1","Rptor"
  ),

  Aging_LongevityDn = c(
    "Igf1","Igf1r","Ins1","Ins2","Insr","Irs1","Irs2",
    "Pik3ca","Akt1","Mtor","Rps6kb1","Eif4e",
    "Ghr","Gh"
  ),

  Senescence_SASP = c(
    "Il6","Il1b","Il1a","Tnf","Cxcl1","Cxcl2","Cxcl10",
    "Ccl2","Ccl5","Ccl7","Mmp3","Mmp9","Mmp13",
    "Vegfa","Igfbp2","Igfbp4","Igfbp5","Igfbp7",
    "Serpine1","Serpinb2","Cdkn1a","Cdkn2a","Trp53",
    "Gdf15","Hmgb1","Lmnb1","Ilk"
  ),

  Senescence_DDR = c(
    "Atm","Atr","Chek1","Chek2","H2afx","Brca1","Brca2",
    "Rad51","Trp53","Cdkn1a","Cdkn2a","Rb1","E2f1",
    "Ppm1d","Mdm2","Mdm4","Bbc3","Pmaip1"
  ),

  Inflammaging = c(
    "Nfkb1","Nfkb2","Rela","Relb","Ikbkb","Ikbkg",
    "Tnf","Il6","Il1b","Il18","Nlrp3","Pycard","Casp1",
    "Irf3","Irf7","Sting1","Cgas","Mavs","Tbk1",
    "Ccl2","Cxcl10","Ifnb1","Ifng","Stat1","Stat3"
  ),

  cGAS_STING = c(
    "Cgas","Sting1","Tbk1","Irf3","Irf7",
    "Ifnb1","Ifna1","Isg15","Mx1","Mx2","Oas1a",
    "Rnasel","Trex1","Dnase2","Enpp1"
  ),

  Rho_Activity = c(
    "Rhoa","Rhob","Rhoc","Rac1","Rac2","Cdc42",
    "Rock1","Rock2","Pak1","Pak2","Pak4",
    "Limk1","Limk2","Cofilin1","Arpc2","Arpc3",
    "Wasf1","Wasf2","Nwasp",
    "Arhgef1","Arhgef2","Arhgef7","Arhgef12",
    "Arhgap1","Arhgap5","Arhgap35"
  ),

  Mito_Stress = c(
    "Pink1","Prkn","Bnip3","Bnip3l","Fundc1",
    "Mfn1","Mfn2","Opa1","Drp1","Fis1",
    "Ndufs1","Ndufs2","Sdha","Uqcrc1","Cox5a",
    "Tfam","Nrf1","Ppargc1a","Ppargc1b"
  )
)

# =============================================================================
# ГЕНСЕТЫ — ЧЕЛОВЕК (символы UPPERCASE, Homo sapiens)
# NEW: отдельный список — гены ортологичные, регистр официальный HGNC
# =============================================================================
SIG_LIST_HUMAN <- list(

  Aging_DNAmAge = c(
    "DNMT1","DNMT3A","DNMT3B","TET1","TET2","TET3",
    "HDAC1","HDAC2","SIRT1","SIRT6","CBX5","CBX7",
    "BMI1","EZH2","KDM6A","KDM5C"
  ),

  Aging_LongevityUp = c(
    "HSPA5","HSP90AB1","HSPA8","HSPA1A","DNAJB1",
    "FOXO3","SIRT1","SIRT3","SIRT6","PPARGC1A",
    "NFE2L2","SQSTM1","MAP1LC3B","BECN1","ATG5",
    "IGFBP1","IGFBP3","EIF4EBP1","RPTOR"
  ),

  Aging_LongevityDn = c(
    "IGF1","IGF1R","INS","INSR","IRS1","IRS2",
    "PIK3CA","AKT1","MTOR","RPS6KB1","EIF4E",
    "GHR","GH1"
  ),

  Senescence_SASP = c(
    "IL6","IL1B","IL1A","TNF","CXCL1","CXCL2","CXCL10",
    "CCL2","CCL5","CCL7","MMP3","MMP9","MMP13",
    "VEGFA","IGFBP2","IGFBP4","IGFBP5","IGFBP7",
    "SERPINE1","SERPINB2","CDKN1A","CDKN2A","TP53",
    "GDF15","HMGB1","LMNB1","ILK"
  ),

  Senescence_DDR = c(
    "ATM","ATR","CHEK1","CHEK2","H2AX","BRCA1","BRCA2",
    "RAD51","TP53","CDKN1A","CDKN2A","RB1","E2F1",
    "PPM1D","MDM2","MDM4","BBC3","PMAIP1"
  ),

  Inflammaging = c(
    "NFKB1","NFKB2","RELA","RELB","IKBKB","IKBKG",
    "TNF","IL6","IL1B","IL18","NLRP3","PYCARD","CASP1",
    "IRF3","IRF7","STING1","CGAS","MAVS","TBK1",
    "CCL2","CXCL10","IFNB1","IFNG","STAT1","STAT3"
  ),

  cGAS_STING = c(
    "CGAS","STING1","TBK1","IRF3","IRF7",
    "IFNB1","IFNA1","ISG15","MX1","MX2","OAS1",
    "RNASEL","TREX1","DNASE2","ENPP1"
  ),

  Rho_Activity = c(
    "RHOA","RHOB","RHOC","RAC1","RAC2","CDC42",
    "ROCK1","ROCK2","PAK1","PAK2","PAK4",
    "LIMK1","LIMK2","CFL1","ARPC2","ARPC3",
    "WASF1","WASF2","WASL",
    "ARHGEF1","ARHGEF2","ARHGEF7","ARHGEF12",
    "ARHGAP1","ARHGAP5","ARHGAP35"
  ),

  Mito_Stress = c(
    "PINK1","PRKN","BNIP3","BNIP3L","FUNDC1",
    "MFN1","MFN2","OPA1","DNM1L","FIS1",
    "NDUFS1","NDUFS2","SDHA","UQCRC1","COX5A",
    "TFAM","NRF1","PPARGC1A","PPARGC1B"
  )
)

# =============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# =============================================================================

# ── VST из norm_counts (без исходного dds) ────────────────────────────────────
build_vst <- function(norm_counts, coldata, design_formula = ~ genotype) {

  raw_pseudo <- round(as.matrix(norm_counts))
  storage.mode(raw_pseudo) <- "integer"

  cd_deseq <- as.data.frame(coldata) %>%
    tibble::column_to_rownames("sample_id")
  cd_deseq <- cd_deseq[colnames(raw_pseudo), , drop = FALSE]

  stopifnot(identical(colnames(raw_pseudo), rownames(cd_deseq)))

  dds <- DESeqDataSetFromMatrix(
    countData = raw_pseudo,
    colData   = cd_deseq,
    design    = design_formula
  )
  sizeFactors(dds) <- rep(1, ncol(dds))
  vsd <- vst(dds, blind = FALSE)
  assay(vsd)
}

# ── Z-score сигнатурного скора ────────────────────────────────────────────────
score_signature <- function(mat, genes, sig_name) {
  genes_found <- intersect(genes, rownames(mat))
  if (length(genes_found) == 0) {
    warning("Сигнатура '", sig_name, "': ни один ген не найден")
    return(NULL)
  }
  message(sprintf("  %-25s %d / %d генов",
                  sig_name, length(genes_found), length(genes)))

  mean_expr <- colMeans(mat[genes_found, , drop = FALSE])
  z_score   <- scale(mean_expr)[, 1]

  tibble(
    sample_id = names(z_score),
    score     = z_score,
    signature = sig_name,
    n_genes   = length(genes_found)
  )
}

# ── Beeswarm plot одной сигнатуры ─────────────────────────────────────────────
plot_signature <- function(sig_name, scores_df, has_tnf = TRUE) {

  df <- scores_df %>% dplyr::filter(signature == sig_name)
  n_genes <- unique(df$n_genes)

  # Уровни group_label в зависимости от датасета
  if (has_tnf) {
    lvls <- c("WT","PPM1D_KO","PPM1B_KO","DKO",
              "WT\n+TNF","PPM1D_KO\n+TNF","PPM1B_KO\n+TNF")
    pal  <- GENO_COLORS_TNF
  } else {
    lvls <- c("WT","PPM1D_KO","PPM1B_KO")
    pal  <- GENO_COLORS
  }

  df <- df %>%
    dplyr::mutate(
      group_label = factor(
        ifelse(has_tnf & tnf == "TNF12h",
               paste0(as.character(genotype), "\n+TNF"),
               as.character(genotype)),
        levels = lvls
      )
    )

  p <- ggplot(df, aes(x = group_label, y = score, fill = group_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5,
                 color = "grey30", linewidth = 0.4) +
    ggbeeswarm::geom_beeswarm(
      aes(shape = if (has_tnf) tnf else NULL),
      size = 4, cex = 2.5, alpha = 0.9,
      stroke = 0.8, color = "grey20"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title    = sig_name,
      subtitle = paste0("n_genes = ", n_genes,
                        " | z-score of mean VST expression"),
      x = NULL, y = "Z-score"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size = 9),
          plot.title         = element_text(face = "bold", size = 12))

  # Разделитель basal | TNF только для MC-38
  if (has_tnf) {
    p <- p +
      geom_vline(xintercept = 4.5, linetype = "dotted",
                 color = "grey60", linewidth = 0.6) +
      scale_shape_manual(
        values = c("basal" = 21, "TNF12h" = 24),
        labels = c("Basal", "TNF 12h"), name = "Treatment"
      )
  }
  p
}

# ── Фиксированный порядок образцов для хитмапов ───────────────────────────────
# Мышь: WT → PPM1D_KO → PPM1B_KO → DKO → (те же + TNF)
# Человек: WT → PPM1D_KO → PPM1B_KO (TNF нет)

order_samples <- function(coldata, has_tnf = TRUE) {

  geno_levels <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")

  if (has_tnf) {
    # Сначала все basal в нужном порядке, затем TNF в том же порядке
    coldata %>%
      dplyr::mutate(
        geno_ord = factor(genotype, levels = geno_levels),
        tnf_ord  = factor(tnf, levels = c("basal", "TNF12h"))
      ) %>%
      dplyr::arrange(tnf_ord, geno_ord) %>%        # basal блок, потом TNF блок
      dplyr::select(-geno_ord, -tnf_ord)
  } else {
    coldata %>%
      dplyr::mutate(geno_ord = factor(genotype, levels = geno_levels)) %>%
      dplyr::arrange(geno_ord) %>%
      dplyr::select(-geno_ord)
  }
}


# ── HM1: скоры сигнатур × образцы ────────────────────────────────────────────
make_hm1 <- function(scores_all, coldata, has_tnf = TRUE,
                     dataset_label = "") {

  # ── НОВОЕ: фиксируем порядок образцов ──────────────────────────────────────
  coldata <- order_samples(coldata, has_tnf)

  hm_mat <- scores_all %>%
    dplyr::select(signature, sample_id, score) %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = score) %>%
    tibble::column_to_rownames("signature") %>%
    as.matrix()

  hm_mat <- hm_mat[, coldata$sample_id, drop = FALSE]

  # Аннотация колонок
  ann_list <- list(
    Genotype = coldata$genotype,
    col      = list(Genotype = GENO_COLORS[levels(coldata$genotype)])
  )

  if (has_tnf) {
    col_ann <- HeatmapAnnotation(
      Genotype = coldata$genotype,
      TNF      = coldata$tnf,
      col = list(
        Genotype = GENO_COLORS[levels(coldata$genotype)],
        TNF      = c("basal" = "#F7F7F7", "TNF12h" = "#FC8D59")
      ),
      annotation_name_gp = gpar(fontsize = 9),
      annotation_height  = unit(c(4, 4), "mm")
    )
    split_col <- factor(ifelse(coldata$tnf == "basal", "Basal", "TNF 12h"),
                        levels = c("Basal","TNF 12h"))
  } else {
    col_ann <- HeatmapAnnotation(
      Genotype = coldata$genotype,
      col = list(Genotype = GENO_COLORS[levels(coldata$genotype)]),
      annotation_name_gp = gpar(fontsize = 9),
      annotation_height  = unit(4, "mm")
    )
    split_col <- NULL
  }

  # Аннотация строк (категория)
  sig_cat_vec <- SIG_CATEGORY[rownames(hm_mat)]
  row_ann <- rowAnnotation(
    Category = sig_cat_vec,
    col = list(Category = CAT_COLORS),
    annotation_name_gp = gpar(fontsize = 9),
    width = unit(4, "mm")
  )

  col_fun <- colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#2166AC","#92C5DE","white","#F4A582","#D6604D")
  )

  Heatmap(
    hm_mat,
    name               = "Z-score",
    col                = col_fun,
    top_annotation     = col_ann,
    right_annotation   = row_ann,
    column_split       = split_col,
    column_gap         = unit(3, "mm"),
    cluster_columns    = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows       = TRUE,
    show_column_names  = TRUE,
    column_names_gp    = gpar(fontsize = 8),
    row_names_gp       = gpar(fontsize = 9),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", hm_mat[i, j]), x, y,
                gp = gpar(fontsize = 6.5,
                          col = ifelse(abs(hm_mat[i,j]) > 1.2,
                                       "white","grey20")))
    },
    column_title = dataset_label,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Z-score",
      at    = c(-2,-1,0,1,2),
      legend_height = unit(40, "mm")
    )
  )
}

# ── HM2: гены × образцы для одной сигнатуры ──────────────────────────────────
make_hm2_single <- function(sig_name, gene_vec, vst_mat, coldata,
                            has_tnf = TRUE) {

  # ── НОВОЕ: фиксируем порядок образцов ──────────────────────────────────────
  coldata <- order_samples(coldata, has_tnf)

  genes_in <- intersect(gene_vec, rownames(vst_mat))
  if (length(genes_in) < 2) {
    message("  Пропускаем ", sig_name, ": < 2 генов найдено")
    return(NULL)
  }

  mat_z <- t(scale(t(vst_mat[genes_in, coldata$sample_id, drop = FALSE])))

  if (has_tnf) {
    split_col <- factor(ifelse(coldata$tnf == "basal","Basal","TNF 12h"),
                        levels = c("Basal","TNF 12h"))
    tnf_ann   <- coldata$tnf
  } else {
    split_col <- NULL
    tnf_ann   <- NULL
  }

  geno_pal <- GENO_COLORS[levels(coldata$genotype)]

  if (!is.null(tnf_ann)) {
    ca <- HeatmapAnnotation(
      Genotype = coldata$genotype,
      TNF      = tnf_ann,
      col = list(Genotype = geno_pal,
                 TNF = c("basal"="#F7F7F7","TNF12h"="#FC8D59")),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_height  = unit(3, "mm"),
      show_legend        = FALSE
    )
  } else {
    ca <- HeatmapAnnotation(
      Genotype = coldata$genotype,
      col = list(Genotype = geno_pal),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_height  = unit(3, "mm"),
      show_legend        = FALSE
    )
  }

  col_fun_gene <- colorRamp2(c(-2,0,2),
                             c("#2166AC","white","#D6604D"))

  Heatmap(
    mat_z,
    name              = "Row Z-score",
    col               = col_fun_gene,
    top_annotation    = ca,
    column_split      = split_col,
    column_gap        = unit(2, "mm"),
    cluster_columns   = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows      = TRUE,
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 7),
    row_names_gp      = gpar(fontsize = 8, fontface = "italic"),
    column_title      = sig_name,
    column_title_gp   = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(legend_height = unit(30, "mm"))
  )
}

# ── HM3: корреляционная матрица сигнатур ──────────────────────────────────────
make_hm3 <- function(scores_all, dataset_label = "") {

  cor_mat <- scores_all %>%
    dplyr::select(signature, sample_id, score) %>%
    tidyr::pivot_wider(names_from = signature, values_from = score) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix() %>%
    cor(method = "pearson")

  sig_cat_vec <- SIG_CATEGORY[rownames(cor_mat)]
  ann_cor <- rowAnnotation(
    Category = sig_cat_vec,
    col = list(Category = CAT_COLORS),
    annotation_name_gp = gpar(fontsize = 8),
    width = unit(4, "mm")
  )

  col_fun_cor <- colorRamp2(
    c(-1,-0.5,0,0.5,1),
    c("#2166AC","#92C5DE","white","#F4A582","#D6604D")
  )

  Heatmap(
    cor_mat,
    name              = "Pearson r",
    col               = col_fun_cor,
    left_annotation   = ann_cor,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 8),
    row_names_gp      = gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i != j)
        grid.text(sprintf("%.2f", cor_mat[i,j]), x, y,
                  gp = gpar(fontsize = 7,
                            col = ifelse(abs(cor_mat[i,j]) > 0.6,
                                         "white","grey20")))
    },
    column_title    = paste0("Inter-signature correlation — ", dataset_label),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-1,-0.5,0,0.5,1),
      legend_height = unit(35, "mm")
    )
  )
}

# =============================================================================
# ГЛАВНАЯ ФУНКЦИЯ — запускает весь блок для одного датасета
# =============================================================================
run_for_dataset <- function(species) {

  stopifnot(species %in% c("mouse","human"))

  message("\n", paste(rep("═", 56), collapse=""))
  message("  ДАТАСЕТ: ", toupper(species))
  message(paste(rep("═", 56), collapse=""))

  # ── 1. Загрузка ──────────────────────────────────────────────────────────────
  rds_path  <- paste0("output/rds/pipeline_", species, ".rds")
  fig_dir   <- paste0("output/figures/", species, "/")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  pipeline  <- readRDS(rds_path)
  norm_counts <- pipeline$norm_counts
  coldata     <- pipeline$coldata

  # Выбираем генсеты по виду
  sig_list  <- if (species == "mouse") SIG_LIST_MOUSE else SIG_LIST_HUMAN
  has_tnf   <- any(coldata$tnf == "TNF12h")
  dataset_label <- if (species == "mouse") "MC38 (Mus musculus)"
  else                    "HT29 (Homo sapiens)"

  # ── 2. VST ───────────────────────────────────────────────────────────────────
  # Для mouse: дизайн с genotype (TNF — exploratory, не в design)
  design_f <- if (species == "mouse") ~ genotype else ~ genotype
  vst_mat  <- build_vst(norm_counts, coldata, design_f)
  message("✓ VST: ", nrow(vst_mat), " × ", ncol(vst_mat))

  # ── 3. Скоры ─────────────────────────────────────────────────────────────────
  message("\nРасчёт сигнатурных скоров:")
  scores_all <- purrr::map_dfr(
    names(sig_list),
    ~ score_signature(vst_mat, sig_list[[.x]], .x)
  ) %>%
    dplyr::left_join(
      coldata %>% dplyr::select(sample_id, genotype, tnf,
                                n_replicates, exploratory),
      by = "sample_id"
    )

  message("✓ Скоры: ", nrow(scores_all), " строк | ",
          length(unique(scores_all$signature)), " сигнатур")

  # ── 4. Beeswarm панели ───────────────────────────────────────────────────────
  key_sigs <- c("Senescence_SASP","Inflammaging",
                "Rho_Activity","Aging_LongevityUp")

  for (sig in unique(scores_all$signature)) {
    p <- plot_signature(sig, scores_all, has_tnf = has_tnf)
    ggsave(paste0(fig_dir, "Score_", sig, ".png"),
           p, width = 9, height = 5, dpi = 300)
  }
  message("✓ Beeswarm: ", length(unique(scores_all$signature)), " файлов")

  # Сводная панель 4 ключевых скора
  panel_plots <- lapply(key_sigs, plot_signature,
                        scores_df = scores_all, has_tnf = has_tnf)
  p_panel <- patchwork::wrap_plots(panel_plots, ncol = 2) +
    patchwork::plot_annotation(
      title    = paste0("Signature scores — ", dataset_label),
      subtitle = "Z-score of mean VST expression",
      theme    = theme(plot.title    = element_text(size=13, face="bold"),
                       plot.subtitle = element_text(size=10, color="grey40"))
    )
  ggsave(paste0(fig_dir, "Score_panel_4x.pdf"),
         p_panel, width = 16, height = 10)
  ggsave(paste0(fig_dir, "Score_panel_4x.png"),
         p_panel, width = 16, height = 10, dpi = 300)
  message("✓ Панель 4 скора сохранена")

  # ── 5. HM1: обзорный хитмап скоров ──────────────────────────────────────────
  hm1 <- make_hm1(scores_all, coldata, has_tnf = has_tnf,
                  dataset_label = dataset_label)

  for (ext in c("pdf","png")) {
    fname <- paste0(fig_dir, "HM1_signature_scores_overview.", ext)
    if (ext == "pdf") pdf(fname, width = 12, height = 5.5)
    else              png(fname, width = 3600, height = 1650, res = 300)
    draw(hm1)
    dev.off()
  }
  message("✓ HM1 сохранён")

  # ── 6. HM2: ген-уровневые хитмапы (4 приоритетных) ──────────────────────────
  priority_sigs <- sig_list[c("Senescence_SASP","Inflammaging",
                              "cGAS_STING","Rho_Activity")]

  hm2_list <- lapply(names(priority_sigs), function(nm) {
    make_hm2_single(nm, priority_sigs[[nm]], vst_mat, coldata,
                    has_tnf = has_tnf)
  })
  names(hm2_list) <- names(priority_sigs)
  hm2_list <- Filter(Negate(is.null), hm2_list)

  # PNG по одному на сигнатуру
  for (nm in names(hm2_list)) {
    png(paste0(fig_dir, "HM2_genes_", nm, ".png"),
        width = 2800, height = 2200, res = 300)
    draw(hm2_list[[nm]])
    dev.off()
  }
  message("✓ HM2: ", length(hm2_list), " файлов сохранено")

  # ── 7. HM3: корреляционная матрица ───────────────────────────────────────────
  hm3 <- make_hm3(scores_all, dataset_label)

  for (ext in c("pdf","png")) {
    fname <- paste0(fig_dir, "HM3_signature_correlation.", ext)
    if (ext == "pdf") pdf(fname, width = 7, height = 6.5)
    else              png(fname, width = 2100, height = 1950, res = 300)
    draw(hm3)
    dev.off()
  }
  message("✓ HM3 сохранён")

  # ── 8. Таблица скоров ────────────────────────────────────────────────────────
  dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
  write_csv(scores_all,
            paste0("output/tables/signature_scores_", species, ".csv"))
  message("✓ Таблица скоров: output/tables/signature_scores_", species, ".csv")

  # ── 9. Обновляем RDS ─────────────────────────────────────────────────────────
  pipeline$vst_mat    <- vst_mat
  pipeline$sig_scores <- scores_all
  pipeline$sig_list   <- sig_list
  saveRDS(pipeline, rds_path)
  message("✓ RDS обновлён: +vst_mat, +sig_scores, +sig_list")

  # ── Итог ─────────────────────────────────────────────────────────────────────
  message("\n── Файлы (", species, ") ──")
  message("  ", fig_dir, "Score_*.png             [",
          length(unique(scores_all$signature)), " файлов]")
  message("  ", fig_dir, "Score_panel_4x.png")
  message("  ", fig_dir, "HM1_signature_scores_overview.png")
  message("  ", fig_dir, "HM2_genes_*.png              [",
          length(hm2_list), " файлов]")
  message("  ", fig_dir, "HM3_signature_correlation.png")

  invisible(pipeline)
}

# =============================================================================
# ТОЧКА ВХОДА
# =============================================================================

# ── Для запуска из RStudio вручную — раскомментируй нужную строку: ──────────
args <- "mouse"
# args <- "human"
# args <- c("mouse", "human")   # оба сразу


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # Без аргументов — запускаем оба датасета
  message("Аргумент не передан — запускаем mouse + human")
  run_for_dataset("mouse")
  if (file.exists("output/rds/pipeline_human.rds")) {
    run_for_dataset("human")
  } else {
    message("⚠ pipeline_human.rds не найден — пропускаем HT29")
  }
} else if (args[1] %in% c("mouse","human")) {
  run_for_dataset(args[1])
} else {
  stop("Неизвестный аргумент: '", args[1],
       "'. Используй: Rscript 02_heatmaps.R mouse  или  human")
}


# ── Вариант 1: загрузить из pipeline RDS (если 02 уже запускался ранее) ───────
pipeline    <- readRDS("output/rds/pipeline_mouse.rds")
scores_all  <- pipeline$sig_scores

# ── Вариант 2: если pipeline ещё не содержит sig_scores — запусти 02 целиком ──
# source("scripts/02_heatmaps.R")   # раскомментируй и выполни

# Быстрая проверка в скрипте 02_heatmaps.R
scores_wide <- scores_all %>%
  dplyr::select(sample_id, signature, score) %>%
  tidyr::pivot_wider(names_from = signature, values_from = score)

cor(scores_wide$Aging_LongevityUp,
    scores_wide$Senescence_DDR,
    use = "complete.obs")
# Если r > 0.5 → это стресс-ответ, не longevity

scores_wide <- scores_all %>%
  dplyr::select(sample_id, signature, score) %>%
  tidyr::pivot_wider(names_from = signature, values_from = score)

# Точечная проверка LongevityUp vs DDR
r <- cor(scores_wide$Aging_LongevityUp,
         scores_wide$Senescence_DDR,
         use = "complete.obs")

cat(sprintf(
  "LongevityUp ~ Senescence_DDR: r = %.3f\n%s\n",
  r,
  dplyr::case_when(
    r >  0.5 ~ "→ СТРЕСС-ОТВЕТ: высокий LongevityUp коактивируется с DDR",
    r < -0.5 ~ "→ НАСТОЯЩИЙ anti-aging: LongevityUp антикоррелирует с DDR",
    TRUE     ~ "→ Слабая связь, интерпретация неоднозначна"
  )
))

# Полная матрица корреляций всех сигнатур (для HM3 и Discussion)
sig_cols <- names(scores_wide)[-1]   # убираем sample_id
cor_mat  <- cor(scores_wide[, sig_cols], use = "complete.obs")

# Вывод только строки LongevityUp — отсортировано по убыванию r
sort(cor_mat["Aging_LongevityUp", ], decreasing = TRUE) %>%
  round(3) %>%
  print()


message("\n══════════════════════════════════════════════════")
message("✓ Скрипт 02_heatmaps.R завершён")
message("══════════════════════════════════════════════════")
