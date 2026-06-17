# =============================================================================
# 04b_zscore_visualization_v3.R
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Z-score визуализация экспрессии генов сигнатурных панелей.
#   Позволяет видеть ОТНОСИТЕЛЬНЫЕ изменения независимо от абсолютного уровня.
#
#   ЗАДАЧИ:
#   1. Heatmap топ-DEG (все 4 генотипа × TNF)
#   2. Панельные heatmaps (DDR, SASP, Rho, Inflammaging)
#   3. TNF-interaction heatmap
#   4. DKO-специфические гены
#   5. Dotplot всех модульных скоров (из 04_module_clocks)
#   6. TNF-delta barplot
#   7. DKO additivity scatter
#
# ИСПРАВЛЕНИЯ v3:
#   - реальные имена генотипов: PPM1D_KO / PPM1B_KO / DKO
#   - vst_mat — матрица (не SE), берём напрямую из pipeline
#   - coldata — tibble с sample_id, genotype, tnf, group
#   - scores_all загружается из 04 (или пересчитывается здесь)
#   - убраны все ссылки на WipKO / PPM1A_KO / cond_colors старые
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
})

# =============================================================================
# 0. Директории
# =============================================================================
dir.create("output/figures/zscore",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/module_clocks", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",           recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Загрузка данных
# =============================================================================
pipeline  <- readRDS("output/rds/pipeline_mouse.rds")
vst_mat   <- as.matrix(pipeline$norm_counts)   # гены × образцы
coldata   <- pipeline$coldata                   # sample_id, genotype, tnf, group
sig_list  <- pipeline$sig_list

# scores_all из скрипта 04 (если уже есть в памяти — пропустить)
if (!exists("scores_all")) {
  message("scores_all не найден в памяти — пересчитываем из модулей 04")
  source("scripts/Double_KO/04_module_clocks_v3.R")   # подгружает scores_all
}

# =============================================================================
# 2. Цвета и порядки (единые для всего скрипта)
# =============================================================================
geno_order <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")

geno_colors <- c(
  "WT"       = "#2166AC",
  "PPM1D_KO" = "#D73027",
  "PPM1B_KO" = "#F46D43",
  "DKO"      = "#762A83"
)

tnf_colors <- c(
  "basal"  = "#F0F0F0",
  "TNF12h" = "#FC8D59"
)

# Логический порядок образцов
sample_order <- coldata %>%
  mutate(
    grp = paste0(genotype, "_", tnf),
    grp = factor(grp, levels = c(
      "WT_basal",       "PPM1D_KO_basal",  "PPM1B_KO_basal", "DKO_basal",
      "WT_TNF12h",      "PPM1D_KO_TNF12h", "PPM1B_KO_TNF12h"
    ))
  ) %>%
  arrange(grp, sample_id) %>%
  filter(sample_id %in% colnames(vst_mat)) %>%
  pull(sample_id)

# Аннотация столбцов (для ComplexHeatmap)
ann_df <- coldata %>%
  filter(sample_id %in% sample_order) %>%
  arrange(match(sample_id, sample_order)) %>%
  select(sample_id, genotype, tnf) %>%
  mutate(
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf,      levels = c("basal", "TNF12h"))
  )

ha_top <- HeatmapAnnotation(
  Genotype = ann_df$genotype,
  TNF      = ann_df$tnf,
  col = list(
    Genotype = geno_colors,
    TNF      = tnf_colors
  ),
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_side = "left",
  simple_anno_size     = unit(4, "mm")
)

# Разбивка колонок basal | TNF12h
col_split <- factor(
  ifelse(ann_df$tnf == "TNF12h", "TNF 12h", "Basal"),
  levels = c("Basal", "TNF 12h")
)

# Z-score палитра
zscore_col <- colorRamp2(
  c(-2.5, -1, 0, 1, 2.5),
  c("#053061", "#92C5DE", "#F7F7F7", "#F4A582", "#67001F")
)

# =============================================================================
# 3. Вспомогательные функции
# =============================================================================

# Z-score матрица по генам
make_zscore_matrix <- function(genes) {
  g   <- intersect(genes, rownames(vst_mat))
  if (length(g) < 2) return(NULL)
  sub <- vst_mat[g, sample_order, drop = FALSE]
  t(scale(t(sub)))   # z-score по строкам (генам)
}

# Универсальная функция рисования heatmap
draw_zscore_heatmap <- function(genes,
                                title,
                                filename,
                                row_split_vec = NULL,
                                w = 14, h = NULL) {

  z_mat <- make_zscore_matrix(genes)

  if (is.null(z_mat) || nrow(z_mat) == 0) {
    message("⚠ Нет генов для: ", title, " — пропускаем")
    return(invisible(NULL))
  }

  n <- nrow(z_mat)
  if (is.null(h)) h <- max(8, n * 0.22)

  ht <- Heatmap(
    z_mat,
    name              = "Z-score",
    col               = zscore_col,
    top_annotation    = ha_top,
    column_split      = col_split,
    column_title_gp   = gpar(fontsize = 10, fontface = "bold"),
    column_gap        = unit(4, "mm"),
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    row_split         = row_split_vec,
    show_row_dend     = (n <= 80),
    show_column_dend  = FALSE,
    row_names_gp      = gpar(fontsize = ifelse(n > 60, 6, 8)),
    column_names_gp   = gpar(fontsize = 8),
    column_names_rot  = 45,
    border            = TRUE,
    heatmap_legend_param = list(
      title         = "Z-score",
      at            = c(-2.5, -1, 0, 1, 2.5),
      labels        = c("−2.5", "−1", "0", "1", "2.5"),
      legend_height = unit(3.5, "cm")
    )
  )

  out_base <- file.path("output/figures/zscore", filename)
  pdf(paste0(out_base, ".pdf"), width = w, height = h)
  draw(ht,
       column_title    = title,
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side   = "bottom",
       annotation_legend_side = "bottom")
  dev.off()

  png(paste0(out_base, ".png"), width = w, height = h,
      units = "in", res = 300, bg = "white")
  draw(ht,
       column_title    = title,
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side   = "bottom",
       annotation_legend_side = "bottom")
  dev.off()

  message("✓ Heatmap: ", filename, "  (", n, " генов)")
  invisible(ht)
}

# =============================================================================
# 4. Генсеты
# =============================================================================
# Ожидаем: ddr_panel, sasp_panel, rho_panel, inflammaging_panel
# (из sig_list или source)
ddr_panel         <- sig_list[["Senescence_DDR"]]
sasp_panel        <- sig_list[["Senescence_SASP"]]
rho_panel         <- sig_list[["Rho_Activity"]]
inflammaging_panel <- sig_list[["Inflammaging"]]
cgassting_panel   <- sig_list[["cGAS_STING"]]
mito_panel        <- sig_list[["Mito_Stress"]]

# =============================================================================
# 5. HEATMAP A: Топ-50 DEG (все генотипы, basal)
# =============================================================================
message("\n── Heatmap A: топ DEG ──")

# Собираем все baseline контрасты
all_base_results <- pipeline$results_base   # именованный список tibble

top50_genes <- map_dfr(all_base_results, ~ as_tibble(.x), .id = "contrast") %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) >= 1,
    !contrast %in% c("DKO_vs_PPM1D_KO", "DKO_vs_PPM1B_KO")
  ) %>%
  group_by(gene_id) %>%
  summarise(max_lfc = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop") %>%
  slice_max(max_lfc, n = 50, with_ties = FALSE) %>%
  pull(gene_id)

message("  Топ DEG генов: ", length(top50_genes))

draw_zscore_heatmap(
  genes    = top50_genes,
  title    = "Top-50 DEG (basal) — 4 genotypes + TNF context",
  filename = "BIOINF_HM_A_top50_DEG",
  h        = 16
)

# =============================================================================
# 6. HEATMAP B: DDR панель
# =============================================================================
message("\n── Heatmap B: DDR ──")
draw_zscore_heatmap(
  genes    = ddr_panel,
  title    = "DDR signature genes — all conditions",
  filename = "BIOINF_HM_B_DDR_panel",
  h        = 14
)

# =============================================================================
# 7. HEATMAP C: SASP панель
# =============================================================================
message("\n── Heatmap C: SASP ──")
draw_zscore_heatmap(
  genes    = sasp_panel,
  title    = "SASP cytokines — KO genotypes + TNF stimulation",
  filename = "BIOINF_HM_C_SASP_panel"
)

# =============================================================================
# 8. HEATMAP D: Rho GTPase панель
# =============================================================================
message("\n── Heatmap D: Rho ──")
draw_zscore_heatmap(
  genes    = rho_panel,
  title    = "Rho GTPase network — 4 genotypes + TNF",
  filename = "BIOINF_HM_D_Rho_panel",
  h        = 18
)

# =============================================================================
# 9. HEATMAP E: Inflammaging панель
# =============================================================================
message("\n── Heatmap E: Inflammaging ──")
draw_zscore_heatmap(
  genes    = inflammaging_panel,
  title    = "Inflammaging signature — all conditions",
  filename = "BIOINF_HM_E_Inflammaging_panel"
)

# =============================================================================
# 10. HEATMAP F: TNF-interaction DEGs
# =============================================================================
message("\n── Heatmap F: TNF interaction DEGs ──")

results_int <- pipeline$results_tnf_interaction   # именованный список

if (!is.null(results_int) && length(results_int) > 0) {

  int_genes_list <- map(results_int, function(df) {
    df %>% filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(gene_id)
  })

  int_genes_list <- keep(int_genes_list, ~ length(.x) > 0)
  int_all        <- unique(unlist(int_genes_list))

  message("  Interaction DEGs всего: ", length(int_all))

  if (length(int_all) >= 5) {
    # row_split: к какому контрасту относится ген
    row_split_vec <- sapply(int_all, function(g) {
      hits <- names(int_genes_list)[sapply(int_genes_list, function(gs) g %in% gs)]
      if (length(hits) == 0) "shared" else hits[1]
    })

    draw_zscore_heatmap(
      genes         = int_all,
      title         = "TNF-interaction DEGs: altered TNF response in KO vs WT",
      filename      = "BIOINF_HM_F_TNF_interaction",
      row_split_vec = row_split_vec,
      h             = max(10, length(int_all) * 0.25)
    )
  } else {
    message("  ⚠ Менее 5 interaction DEGs при FDR<0.05, |LFC|≥1")
    message("  → Попробуй: padj < 0.1 или |LFC| ≥ 0.5 для эксплораторного heatmap")
  }

} else {
  message("  ⚠ results_tnf_interaction не найден в pipeline — пропускаем")
}

# =============================================================================
# 11. HEATMAP G: DKO-специфические гены
# =============================================================================
message("\n── Heatmap G: DKO-специфические DEG ──")

get_sig_genes <- function(contrast_name, lfc_thr = 1) {
  df <- pipeline$results_base[[contrast_name]]
  if (is.null(df)) return(character(0))
  df %>% filter(padj < 0.05, abs(log2FoldChange) >= lfc_thr) %>% pull(gene_id)
}

dko_sig  <- get_sig_genes("DKO_vs_WT")
ppm1d_sig <- get_sig_genes("PPM1D_KO_vs_WT")
ppm1b_sig <- get_sig_genes("PPM1B_KO_vs_WT")

# DKO-специфические = в DKO, но НЕ одновременно в обоих single KOs
dko_specific <- setdiff(dko_sig, intersect(ppm1d_sig, ppm1b_sig))
message("  DKO-специфических DEG: ", length(dko_specific))

if (length(dko_specific) >= 5) {
  draw_zscore_heatmap(
    genes    = dko_specific,
    title    = "DKO-specific DEGs (absent in both PPM1D_KO and PPM1B_KO)",
    filename = "BIOINF_HM_G_DKO_specific",
    h        = max(8, length(dko_specific) * 0.22)
  )
  write.csv(
    data.frame(gene_id = dko_specific),
    "output/tables/DKO_specific_DEGs.csv",
    row.names = FALSE
  )
} else {
  message("  ⚠ Меньше 5 DKO-специфических генов")
}

# =============================================================================
# 12. VIZ-P1: Dotplot модульных скоров (все модули × условия)
# =============================================================================
message("\n── VIZ-P1: Dotplot модульных скоров ──")

scores_summary_plot <- scores_all %>%
  group_by(module, category, genotype, tnf) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sem        = sd(score, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf, levels = c("basal", "TNF12h"))
  )

p_dot <- ggplot(
  scores_summary_plot,
  aes(x = genotype, y = mean_score,
      color = genotype, shape = tnf)
) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_score - sem, ymax = mean_score + sem),
    width = 0.25, linewidth = 0.5, alpha = 0.7
  ) +
  geom_point(size = 3.5, stroke = 0.8) +
  scale_color_manual(values = geno_colors, name = "Genotype") +
  scale_shape_manual(
    values = c("basal" = 16, "TNF12h" = 17),
    name   = "Stimulation"
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  labs(
    title    = "Module scores — mean ± SEM by genotype",
    subtitle = "Circle = basal; triangle = TNF12h",
    x        = NULL,
    y        = "Mean z-score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 7),
    strip.text       = element_text(size = 8, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("output/figures/module_clocks/VIZ_P1_dotplot_modules.png",
       p_dot, width = 16, height = 13, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_P1_dotplot_modules.pdf",
       p_dot, width = 16, height = 13)
message("✓ VIZ-P1 dotplot сохранён")

# =============================================================================
# 13. VIZ-P2: TNF-delta barplot
# =============================================================================
message("\n── VIZ-P2: TNF-delta barplot ──")

tnf_delta <- scores_all %>%
  filter(genotype != "DKO") %>%          # у DKO нет TNF
  group_by(module, genotype, tnf) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = tnf, values_from = mean_score) %>%
  mutate(
    delta    = TNF12h - basal,
    genotype = factor(genotype, levels = setdiff(geno_order, "DKO"))
  )

p_delta <- ggplot(
  tnf_delta,
  aes(x = module, y = delta, fill = genotype)
) +
  geom_col(position = position_dodge(0.75), width = 0.65, color = "white") +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  scale_fill_manual(values = geno_colors, name = "Genotype") +
  labs(
    title    = "TNF-driven module score shift (Δ) by genotype",
    subtitle = "Δ = mean(TNF12h) − mean(basal)  |  DKO excluded (no TNF data)",
    x        = NULL,
    y        = "Δ z-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  )

ggsave("output/figures/module_clocks/VIZ_P2_TNF_delta_barplot.png",
       p_delta, width = 14, height = 6, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_P2_TNF_delta_barplot.pdf",
       p_delta, width = 14, height = 6)
message("✓ VIZ-P2 TNF-delta barplot сохранён")

write.csv(tnf_delta, "output/tables/TNF_delta_by_module_genotype.csv",
          row.names = FALSE)

# =============================================================================
# 14. VIZ-P3: DKO additivity scatter
# =============================================================================
message("\n── VIZ-P3: DKO additivity scatter ──")

dko_additivity <- scores_all %>%
  filter(tnf == "basal") %>%
  group_by(module, genotype) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = mean_score) %>%
  mutate(
    additive_pred = PPM1D_KO + PPM1B_KO - WT,     # аддитивная модель
    synergy_score = DKO - additive_pred,           # отклонение
    module        = as.character(module)
  )

# Сохраняем таблицу
write.csv(dko_additivity,
          "output/tables/DKO_additivity_test.csv",
          row.names = FALSE)

p_add <- ggplot(
  dko_additivity,
  aes(x = additive_pred, y = DKO, label = module, color = synergy_score)
) +
  # Линия идеальной аддитивности
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.8) +
  # Доверительная полоса ±0.2 вокруг линии аддитивности
  geom_ribbon(
    aes(ymin = additive_pred - 0.2, ymax = additive_pred + 0.2),
    fill = "grey85", alpha = 0.4, color = NA,
    inherit.aes = FALSE,
    data = data.frame(
      additive_pred = range(dko_additivity$additive_pred, na.rm = TRUE)
    )
  ) +
  geom_point(size = 4.5, alpha = 0.9) +
  ggrepel::geom_label_repel(
    size = 3, box.padding = 0.4,
    max.overlaps = 20, seed = 42,
    label.size = 0.2
  ) +
  scale_color_gradient2(
    low      = "#4575B4",
    mid      = "grey80",
    high     = "#D73027",
    midpoint = 0,
    name     = "Synergy\n(DKO − additive)",
    limits   = c(-1, 1),
    oob      = scales::squish
  ) +
  labs(
    title    = "Double KO: observed vs. additive prediction",
    subtitle = "Additive = PPM1D_KO + PPM1B_KO − WT  |  Grey band = ±0.2 additivity zone",
    x        = "Additive prediction (z-score)",
    y        = "DKO observed (z-score)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave("output/figures/module_clocks/VIZ_P3_DKO_additivity.png",
       p_add, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_P3_DKO_additivity.pdf",
       p_add, width = 8, height = 7)
message("✓ VIZ-P3 DKO additivity scatter сохранён")

# =============================================================================
# 15. Итог
# =============================================================================
message("\n", strrep("=", 70))
message("✓✓ 04b COMPLETE")
message(strrep("=", 70))
message("\nHeatmaps (output/figures/zscore/):")
message("  BIOINF_HM_A — топ-50 DEG")
message("  BIOINF_HM_B — DDR панель")
message("  BIOINF_HM_C — SASP панель")
message("  BIOINF_HM_D — Rho GTPase панель")
message("  BIOINF_HM_E — Inflammaging панель")
message("  BIOINF_HM_F — TNF-interaction DEGs")
message("  BIOINF_HM_G — DKO-специфические гены")
message("\nModule clock plots (output/figures/module_clocks/):")
message("  VIZ_P1 — dotplot модульных скоров")
message("  VIZ_P2 — TNF-delta barplot")
message("  VIZ_P3 — DKO additivity scatter")
message("\nТаблицы (output/tables/):")
message("  DKO_specific_DEGs.csv")
message("  TNF_delta_by_module_genotype.csv")
message("  DKO_additivity_test.csv")
