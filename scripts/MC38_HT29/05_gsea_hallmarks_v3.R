# =============================================================================
# 05_gsea_hallmarks_v3.R
#
# ИЗМЕНЕНИЯ ОТНОСИТЕЛЬНО v2:
#   - ИСПРАВЛЕНО: WipKO = Wip1/PPM1D KO (не PPM1A!)
#     PPM1A — отдельная фосфатаза, не связанная с Wip1/ATM-осью
#   - НОВОЕ: Двухвидовой анализ: мышь (MC38) + человек (HT29)
#     Каждый вид анализируется своими генсетами msigdbr, результаты
#     затем приводятся к общей таблице для кросс-видовых сравнений
#   - НОВОЕ: ортологовый мэппинг мышь → человек для NES-сравнения
#   - СТРУКТУРА КОНТРАСТОВ:
#
#   МЫШЬ (MC38, ПЕРВЫЙ ПРИОРИТЕТ):
#     Baseline:
#       Wip1KO_vs_WT         (BD2/5/20 vs BD1/3/4)
#       Ppm1bKO_vs_WT        (BD6/7/8  vs BD1/3/4)
#       DoubleKO_vs_WT       (BD9/57   vs BD1/3/4)
#     TNF 12h (exploratory):
#       WT_TNF_vs_WT         (BD12 vs BD1/3/4)
#       Wip1KO_TNF_vs_Wip1KO (BD13 vs BD2/5/20)
#       Ppm1bKO_TNF_vs_Ppm1bKO (BD14 vs BD6/7/8)
#     Interaction:
#       Interaction_Wip1KO_x_TNF
#       Interaction_Ppm1bKO_x_TNF
#
#   ЧЕЛОВЕК (HT29, ВТОРОЙ ПРИОРИТЕТ):
#     Baseline:
#       HT29_Wip1KO_vs_WT    (BD11/22/23 vs BD21/58/10)
#       HT29_Ppm1bKO_vs_WT   (BD24/25/26 vs BD21/58/10)
#
#   МЕТРИКА РАНЖИРОВАНИЯ:
#     rank = sign(log2FC) × (−log10(pvalue))
#     Объединяет направление и статистическую достоверность.
#
#   NOTA BENE: msigdbr для мыши и человека — разные объекты!
#   Для мыши: msigdbr(species = "Mus musculus", category = "H")
#   Для человека: msigdbr(species = "Homo sapiens", category = "H")
#   Названия путей (HALLMARK_*) идентичны → возможен прямой NES-сравнение.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(circlize)
})

# ── 0. КОНФИГУРАЦИЯ ──────────────────────────────────────────────────────────
OUTPUT_RDS    <- "output/rds/"
OUTPUT_TABLES <- "output/tables/"
OUTPUT_FIG    <- "output/figures/"
dir.create(OUTPUT_RDS,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_FIG,    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_FIG, "gsea_enrichment_plots"), showWarnings = FALSE)

# ── 1. ЗАГРУЗКА РЕЗУЛЬТАТОВ ──────────────────────────────────────────────────
# Ожидаем, что 01_deseq2_pipeline_v3.R сохранил:
#   output/rds/deseq2_all_results.rds
# Структура объекта:
#   $mouse$results_base      — именованный список tibble (мышь, baseline)
#   $mouse$results_tnf       — именованный список tibble (мышь, TNF)
#   $mouse$results_int       — именованный список tibble (мышь, interaction)
#   $human$results_base      — именованный список tibble (человек, baseline)

pipeline <- readRDS(file.path(OUTPUT_RDS, "pipeline_mouse.rds"))

# Словарь всех контрастов с аннотацией вида и типа
contrasts_meta <- bind_rows(
  tibble(
    name    = names(pipeline$mouse$results_base),
    species = "mouse",
    type    = "baseline"
  ),
  tibble(
    name    = names(pipeline$mouse$results_tnf),
    species = "mouse",
    type    = "tnf_response"
  ),
  tibble(
    name    = names(pipeline$mouse$results_int),
    species = "mouse",
    type    = "interaction"
  ),
  tibble(
    name    = names(pipeline$human$results_base),
    species = "human",
    type    = "baseline"
  )
)
message("Всего контрастов: ", nrow(contrasts_meta))
print(contrasts_meta)

# ── 2. HALLMARK ГЕНСЕТЫ ───────────────────────────────────────────────────────
# Создаём отдельные объекты для мыши и человека.
# Названия путей идентичны → NES напрямую сопоставимы между видами.

make_hallmark_objects <- function(species_name) {
  raw <- msigdbr(species = species_name, category = "H")

  # fgsea: named list pathway → gene_symbol
  gs_list <- raw %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gs_name) %>%
    summarise(genes = list(gene_symbol), .groups = "drop") %>%
    deframe()

  # clusterProfiler: TERM2GENE data.frame
  term2gene <- raw %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    as.data.frame()

  list(gs_list = gs_list, term2gene = term2gene)
}

hallmark_mouse <- make_hallmark_objects("Mus musculus")
hallmark_human <- make_hallmark_objects("Homo sapiens")
message("Hallmark путей мышь:   ", length(hallmark_mouse$gs_list))
message("Hallmark путей человек: ", length(hallmark_human$gs_list))

# ── 3. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ───────────────────────────────────────────────

# 3A. Построение ранжированного вектора
# Входные данные: tibble с колонками gene_id, log2FoldChange, pvalue
build_ranked_vec <- function(de_tbl) {
  de_tbl %>%
    dplyr::select(gene_id, log2FoldChange, pvalue) %>%
    dplyr::filter(!is.na(log2FoldChange)) %>%
    dplyr::mutate(
      pvalue_clean = dplyr::case_when(
        is.na(pvalue) | pvalue <= 0 ~ 1,
        pvalue > 1                  ~ 1,
        TRUE                        ~ pvalue
      ),
      rank_metric = sign(log2FoldChange) * (-log10(pvalue_clean))
    ) %>%
    dplyr::arrange(desc(rank_metric)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    { setNames(.$rank_metric, .$gene_id) }
}

# 3B. fgsea для одного контраста
run_fgsea_contrast <- function(de_tbl, contrast_name, pathways,
                               min_size = 10, max_size = 500, n_perm = 1000) {

  message("\n>>> fgsea: ", contrast_name)

  ranked_vec <- build_ranked_vec(de_tbl)
  message("  Генов в ранге: ", length(ranked_vec),
          " | Диапазон: [", round(min(ranked_vec),2),
          " ; ", round(max(ranked_vec),2), "]")

  set.seed(42)
  res <- tryCatch(
    fgsea(
      pathways    = pathways,
      stats       = ranked_vec,
      minSize     = min_size,
      maxSize     = max_size,
      eps         = 0,
      nPermSimple = n_perm
    ),
    error = function(e) {
      warning("fgsea ошибка (", contrast_name, "): ", e$message)
      return(NULL)
    }
  )
  if (is.null(res)) return(NULL)

  res %>%
    as.data.frame() %>%
    dplyr::mutate(
      contrast      = contrast_name,
      NES_direction = ifelse(NES > 0, "enriched_up", "enriched_down"),
      sig           = padj < 0.05,
      pathway_short = gsub("HALLMARK_", "", pathway)
    ) %>%
    dplyr::arrange(padj)
}

# ── 4. ЗАПУСК GSEA ПО ВСЕМ КОНТРАСТАМ ────────────────────────────────────────

run_all_gsea <- function(results_list, pathways, species_tag, type_tag) {
  lapply(names(results_list), function(nm) {
    res <- run_fgsea_contrast(
      de_tbl        = results_list[[nm]],
      contrast_name = nm,
      pathways      = pathways
    )
    if (!is.null(res)) res$species <- species_tag
    if (!is.null(res)) res$contrast_type <- type_tag
    res
  }) %>%
    setNames(names(results_list)) %>%
    Filter(Negate(is.null), .)
}

gsea_results_list <- c(
  run_all_gsea(pipeline$mouse$results_base, hallmark_mouse$gs_list, "mouse", "baseline"),
  run_all_gsea(pipeline$mouse$results_tnf,  hallmark_mouse$gs_list, "mouse", "tnf_response"),
  run_all_gsea(pipeline$mouse$results_int,  hallmark_mouse$gs_list, "mouse", "interaction")
 # run_all_gsea(pipeline$human$results_base, hallmark_human$gs_list, "human", "baseline")
)

gsea_combined <- bind_rows(gsea_results_list) %>%
  dplyr::mutate(
   # species       = factor(species,       levels = c("mouse", "human")), #
    contrast_type = factor(contrast_type, levels = c("baseline","tnf_response","interaction"))
  )

write.csv(gsea_combined,
          file.path(OUTPUT_TABLES, "gsea_hallmarks_all_contrasts_v3.csv"),
          row.names = FALSE)
message("\n✓ GSEA завершена. Строк: ", nrow(gsea_combined),
        " | Значимых (padj<0.05): ", sum(gsea_combined$sig, na.rm = TRUE))

# ── 5. ВИЗУАЛИЗАЦИИ ───────────────────────────────────────────────────────────

# ── 5A. NES HEATMAP: ВСЕ ЗНАЧИМЫЕ ПУТИ × ВСЕ КОНТРАСТЫ ─────────────────────
# Отдельно для мыши и человека (разные гены, разный масштаб NES)

make_nes_heatmap <- function(gsea_df, species_label) {

  sig_paths <- gsea_df %>%
    dplyr::filter(sig) %>%
    dplyr::pull(pathway_short) %>%
    unique()

  if (length(sig_paths) == 0) {
    message("! Нет значимых путей для ", species_label)
    return(invisible(NULL))
  }

  nes_mat <- gsea_df %>%
    dplyr::filter(pathway_short %in% sig_paths) %>%
    dplyr::select(pathway_short, contrast, NES) %>%
    pivot_wider(names_from = contrast, values_from = NES) %>%
    column_to_rownames("pathway_short") %>%
    as.matrix()
  nes_mat[is.na(nes_mat)] <- 0

  # Аннотация столбцов — тип контраста и вид
  meta_cols <- gsea_df %>%
    dplyr::distinct(contrast, contrast_type) %>%
    dplyr::filter(contrast %in% colnames(nes_mat))

  col_type <- setNames(
    meta_cols$contrast_type,
    meta_cols$contrast
  )
  col_type_vec <- as.character(col_type[colnames(nes_mat)])

  type_colors <- c(
    "baseline"     = "#AEC6E8",
    "tnf_response" = "#FDAE6B",
    "interaction"  = "#C994C7"
  )

  col_anno <- HeatmapAnnotation(
    Type = col_type_vec,
    col  = list(Type = type_colors),
    annotation_name_side = "left",
    annotation_legend_param = list(Type = list(title = "Contrast type"))
  )

  nes_col_fun <- colorRamp2(
    breaks = c(-3, -1.5, 0, 1.5, 3),
    colors = c("#053061","#92C5DE","#F7F7F7","#F4A582","#67001F")
  )

  Heatmap(
    nes_mat,
    name             = "NES",
    col              = nes_col_fun,
    top_annotation   = col_anno,
    cluster_rows     = TRUE,
    cluster_columns  = FALSE,
    show_row_dend    = TRUE,
    row_names_gp     = gpar(fontsize = 8),
    column_names_gp  = gpar(fontsize = 8),
    column_names_rot = 45,

    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- nes_mat[i, j]
      if (!is.na(v) && abs(v) >= 1.5)
        grid.text(round(v, 1), x, y,
                  gp = gpar(fontsize = 7,
                            col = ifelse(abs(v) > 2.5, "white", "black")))
    },

    column_title     = paste0("Hallmark GSEA NES — ", species_label),
    column_title_gp  = gpar(fontsize = 11, fontface = "bold"),
    width  = unit(max(12, ncol(nes_mat) * 1.5), "cm"),
    height = unit(max(8,  nrow(nes_mat) * 0.45), "cm"),

    heatmap_legend_param = list(
      title         = "NES",
      at            = c(-3,-2,-1,0,1,2,3),
      legend_height = unit(4, "cm")
    )
  )
}

# Мышь
ht_mouse <- make_nes_heatmap(
  dplyr::filter(gsea_combined, species == "mouse"),
  "Mouse MC38"
)
pdf(file.path(OUTPUT_FIG, "BIOINF_gsea_NES_heatmap_mouse.pdf"),
    width = 18, height = 12)
draw(ht_mouse)
dev.off()

# Человек
ht_human <- make_nes_heatmap(
  dplyr::filter(gsea_combined, species == "human"),
  "Human HT29"
)
pdf(file.path(OUTPUT_FIG, "BIOINF_gsea_NES_heatmap_human.pdf"),
    width = 14, height = 10)
draw(ht_human)
dev.off()
message("✓ NES heatmaps сохранены")

# ── 5B. КРОСС-ВИДОВОЙ DOTPLOT: мышь vs человек для BASELINE ─────────────────
# Ключевой рисунок: какие пути значимы у обоих видов при Wip1KO?
# Это усиливает трансляционную значимость результатов.

cross_species <- gsea_combined %>%
  dplyr::filter(
    contrast_type == "baseline",
    grepl("Wip1KO_vs_WT|HT29_Wip1KO_vs_WT", contrast)
  )

# Пути, значимые хотя бы в одном виде
paths_cross <- cross_species %>%
  dplyr::filter(sig) %>%
  dplyr::pull(pathway_short) %>%
  unique()

# Все пути в обоих видах (включая незначимые — для полной матрицы)
plot_cross <- cross_species %>%
  dplyr::filter(pathway_short %in% paths_cross) %>%
  dplyr::mutate(
    species_label = ifelse(species == "mouse",
                           "Mouse MC38 (Wip1KO vs WT)",
                           "Human HT29 (Wip1KO vs WT)"),
    # Консерватизм: помечаем пути, значимые у обоих
    conserved = pathway_short %in% (
      cross_species %>%
        dplyr::filter(sig) %>%
        dplyr::group_by(pathway_short) %>%
        dplyr::filter(n_distinct(species) == 2) %>%
        dplyr::pull(pathway_short) %>%
        unique()
    )
  )

p_cross <- ggplot(
  plot_cross,
  aes(x     = species_label,
      y     = reorder(pathway_short, NES),
      size  = -log10(padj + 1e-10),
      color = NES,
      shape = conserved)
) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#D6604D",
    midpoint = 0,
    name     = "NES"
  ) +
  scale_size_continuous(range = c(2, 9), name = "-log10(FDR)") +
  scale_shape_manual(
    values = c(`TRUE` = 18, `FALSE` = 16),
    labels = c("Species-specific", "Conserved (both species)"),
    name   = "Conservation"
  ) +
  geom_vline(xintercept = 1.5, color = "grey70", linetype = "dashed") +
  labs(
    title    = "Cross-species GSEA: Wip1KO vs WT — Mouse MC38 vs Human HT29",
    subtitle = "Diamond = pathway significant in BOTH species | FDR < 0.05 in ≥ 1 species",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_crossspecies_Wip1KO.pdf"),
       p_cross, width = 13, height = 10)
ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_crossspecies_Wip1KO.png"),
       p_cross, width = 13, height = 10, dpi = 300)
message("✓ Cross-species dotplot сохранён")

# ── 5C. BARPLOT: Double_KO vs Single KOs (мышь) ──────────────────────────────
# Аддитивность vs синергия: усиливает ли двойной KO эффекты одиночных?

target_contrasts_mouse <- c("Wip1KO_vs_WT", "Ppm1bKO_vs_WT", "DoubleKO_vs_WT")

baseline_compare <- gsea_combined %>%
  dplyr::filter(
    species  == "mouse",
    contrast %in% target_contrasts_mouse
  )

paths_show <- baseline_compare %>%
  dplyr::filter(sig) %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_max(abs(NES), n = 20) %>%
  dplyr::ungroup() %>%
  dplyr::pull(pathway_short) %>%
  unique()

plot_comp <- baseline_compare %>%
  dplyr::filter(pathway_short %in% paths_show) %>%
  dplyr::mutate(
    NES_plot = ifelse(sig, NES, 0),
    contrast = factor(contrast,
                      levels = c("Wip1KO_vs_WT",
                                 "Ppm1bKO_vs_WT",
                                 "DoubleKO_vs_WT"))
  )

p_comp <- ggplot(
  plot_comp,
  aes(x = NES_plot,
      y = reorder(pathway_short, NES),
      fill = contrast)
) +
  geom_col(
    position = position_dodge(0.75),
    width    = 0.65,
    color    = "white",
    linewidth = 0.25
  ) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "Wip1KO_vs_WT"    = "#C51B7D",
      "Ppm1bKO_vs_WT"   = "#4393C3",
      "DoubleKO_vs_WT"  = "#E08214"
    ),
    labels = c(
      "Wip1KO_vs_WT"    = "Wip1KO vs WT (PPM1D/Wip1 KO)",
      "Ppm1bKO_vs_WT"   = "Ppm1bKO vs WT (PPM1B KO)",
      "DoubleKO_vs_WT"  = "Double KO vs WT"
    )
  ) +
  labs(
    title    = "Hallmark GSEA: NES comparison — Single KO vs Double KO (Mouse MC38)",
    subtitle = "NES = 0 for non-significant paths | FDR < 0.05",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL,
    fill     = "Genotype"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8)
  )

ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_double_vs_single_KO_NES_v3.pdf"),
       p_comp, width = 13, height = 9)
ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_double_vs_single_KO_NES_v3.png"),
       p_comp, width = 13, height = 9, dpi = 300)
message("✓ Double vs Single KO barplot сохранён")

# ── 5D. TNF INTERACTION DOTPLOT (мышь) ───────────────────────────────────────

interaction_gsea <- gsea_combined %>%
  dplyr::filter(
    species       == "mouse",
    contrast_type == "interaction",
    sig
  ) %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_max(abs(NES), n = 15) %>%
  dplyr::ungroup()

if (nrow(interaction_gsea) > 0) {
  p_inter <- ggplot(
    interaction_gsea,
    aes(x    = NES,
        y    = reorder(pathway_short, NES),
        fill = NES_direction,
        size = -log10(padj + 1e-10))
  ) +
    geom_point(shape = 21, color = "grey40", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    facet_wrap(~ contrast, ncol = 2) +
    scale_fill_manual(
      values = c(
        "enriched_up"   = "#D6604D",
        "enriched_down" = "#2166AC"
      ),
      labels = c(
        "enriched_up"   = "Hyper-response to TNF in KO",
        "enriched_down" = "Hypo-response to TNF in KO"
      ),
      name = "TNF × KO interaction"
    ) +
    scale_size_continuous(range = c(3, 10), name = "-log10(FDR)") +
    labs(
      title    = "Hallmark GSEA — Interaction: TNF response altered in KO (Mouse MC38)",
      subtitle = "Up = pathway MORE activated by TNF in KO than WT | FDR < 0.05",
      x = "NES (interaction term)", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "top"
    )

  ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_TNF_interaction_mouse_v3.pdf"),
         p_inter, width = 14, height = 9)
  ggsave(file.path(OUTPUT_FIG, "BIOINF_gsea_TNF_interaction_mouse_v3.png"),
         p_inter, width = 14, height = 9, dpi = 300)
  message("✓ Interaction plot сохранён")
} else {
  message("! Interaction: нет sig. путей padj<0.05. ",
          "Попробуйте padj<0.1 для exploratory.")
}

# ── 5E. MOUNTAIN PLOTS (clusterProfiler) для BASELINE ────────────────────────
# Классические enrichment plots для топ-6 путей в каждом baseline-контрасте

run_cp_gsea <- function(de_tbl, term2gene) {
  gene_list <- build_ranked_vec(de_tbl)
  GSEA(
    geneList      = gene_list,
    TERM2GENE     = term2gene,
    pvalueCutoff  = 0.25,
    pAdjustMethod = "BH",
    verbose       = FALSE,
    minGSSize     = 10,
    maxGSSize     = 500,
    eps           = 1e-10
  )
}

save_mountain_plots <- function(results_list, term2gene, species_tag) {
  for (nm in names(results_list)) {
    message("Mountain plot: ", nm)
    cp_obj <- tryCatch(
      run_cp_gsea(results_list[[nm]], term2gene),
      error = function(e) { warning(nm, ": ", e$message); NULL }
    )
    if (is.null(cp_obj)) next

    top6 <- cp_obj@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::slice_max(abs(NES), n = 6) %>%
      dplyr::pull(ID)

    if (length(top6) == 0) {
      message("  Нет sig. путей для mountain plot: ", nm)
      next
    }

    p_e <- gseaplot2(
      x          = cp_obj,
      geneSetID  = top6,
      pvalue_by  = "p.adjust",
      title      = paste0("GSEA Hallmark — ", nm, " [", species_tag, "]")
    )
    fn <- file.path(OUTPUT_FIG, "gsea_enrichment_plots",
                    paste0("GSEA_mountain_", species_tag, "_", nm, ".pdf"))
    ggsave(fn, p_e, width = 10, height = 8)
    message("  ✓ Saved: ", basename(fn))
  }
}

save_mountain_plots(pipeline$mouse$results_base, hallmark_mouse$term2gene, "mouse")
save_mountain_plots(pipeline$human$results_base, hallmark_human$term2gene, "human")

# ── 6. КРОСС-ВИДОВАЯ ТАБЛИЦА КОНСЕРВИРОВАННЫХ ПУТЕЙ ─────────────────────────
# Для раздела Results: какие пути значимы у ОБОИХ видов при Wip1KO?

conserved_wip1ko <- gsea_combined %>%
  dplyr::filter(
    sig,
    grepl("Wip1KO_vs_WT|HT29_Wip1KO_vs_WT", contrast)
  ) %>%
  dplyr::group_by(pathway_short) %>%
  dplyr::filter(n_distinct(species) == 2) %>%   # должны быть оба вида
  dplyr::summarise(
    NES_mouse    = NES[species == "mouse"],
    padj_mouse   = padj[species == "mouse"],
    NES_human    = NES[species == "human"],
    padj_human   = padj[species == "human"],
    direction    = ifelse(mean(NES) > 0, "UP", "DOWN"),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(desc(abs(NES_mouse)))

write.csv(conserved_wip1ko,
          file.path(OUTPUT_TABLES, "gsea_conserved_Wip1KO_mouse_human.csv"),
          row.names = FALSE)
message("\n✓ Консервированных путей Wip1KO (мышь ∩ человек): ",
        nrow(conserved_wip1ko))

# ── 7. ФИНАЛЬНОЕ СОХРАНЕНИЕ ──────────────────────────────────────────────────
saveRDS(gsea_combined, file.path(OUTPUT_RDS, "gsea_hallmarks_all_v3.rds"))

message("\n══════════════════════════════════════════════")
message("  GSEA PIPELINE v3 COMPLETE")
message("  Видов:              ", n_distinct(gsea_combined$species))
message("  Контрастов:         ", n_distinct(gsea_combined$contrast))
message("  Путей всего:        ", n_distinct(gsea_combined$pathway))
message("  Значимых (FDR<.05): ",
        sum(gsea_combined$sig, na.rm = TRUE),
        " записей в ",
        n_distinct(gsea_combined$contrast[gsea_combined$sig]),
        " контрастах")
message("══════════════════════════════════════════════")
