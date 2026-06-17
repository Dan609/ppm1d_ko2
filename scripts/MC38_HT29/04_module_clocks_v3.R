# =============================================================================
# # =============================================================================
# 04_module_clocks_v3.R
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   "Module clocks" — многомерные биологические часы:
#   DDR, Inflammaging, Rho, SASP, ProxTubule и др. как спицы радара.
#   Позволяет увидеть ПРОФИЛЬ состояния клетки целиком.
#
#   ВОПРОСЫ:
#   1. DOUBLE KO (DKO): синергия vs. аддитивность vs. компенсация?
#   2. TNF-ТРАЕКТОРИЯ: длина и направление вектора basal→TNF12h
#      (гипотеза: у KO вектор длиннее, направлен в сторону SASP)
#   3. ФАЗОВЫЕ ПОРТРЕТЫ: DDR×Inflammaging и Rho×SASP
#
# ИСПРАВЛЕНИЯ v3:
#   - PPM1D_KO / DKO / PPM1B_KO — реальные имена из coldata
#   - group_order использует реальные уровни: DKO_basal, PPM1D_KO_basal, etc.
#   - убран mod_scores (не определён), везде scores_all
#   - phase_data: правильные имена колонок после pivot_wider
#   - убраны GENO_COLORS / WipKO / Double_KO
#   - score_signature определена внутри скрипта (не зависит от 03)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fmsb)
  library(ggrepel)
  library(patchwork)
  library(ggridges)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# =============================================================================
# 0. Создаём папки
# =============================================================================
dir.create("output/figures/module_clocks", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",               recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Загружаем pipeline
# =============================================================================
pipeline <- readRDS("output/rds/pipeline_mouse.rds")

vst_mat  <- as.matrix(pipeline$norm_counts)   # матрица VST: гены × образцы
coldata  <- pipeline$coldata                   # tibble: sample_id, genotype, tnf, group, ...
sig_list <- pipeline$sig_list                  # генсеты из 02/03 скриптов

message("✓ vst_mat: ", nrow(vst_mat), " генов × ", ncol(vst_mat), " образцов")
message("✓ sig_list генсеты: ", paste(names(sig_list), collapse = ", "))

# Упорядочиваем столбцы матрицы по coldata
stopifnot(all(coldata$sample_id %in% colnames(vst_mat)))
vst_mat <- vst_mat[, coldata$sample_id]

# =============================================================================
# 2. Реальные уровни факторов (из coldata)
# =============================================================================
# genotype: "WT" "PPM1D_KO" "PPM1B_KO" "DKO"
# group:    "WT_basal" "PPM1D_KO_basal" "PPM1B_KO_basal" "DKO_basal"
#           "WT_TNF12h" "PPM1D_KO_TNF12h" "PPM1B_KO_TNF12h"
# tnf:      "basal" "TNF12h"

geno_order  <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")
group_order <- c(
  "WT_basal", "PPM1D_KO_basal", "PPM1B_KO_basal", "DKO_basal",
  "WT_TNF12h", "PPM1D_KO_TNF12h", "PPM1B_KO_TNF12h"
)

geno_colors <- c(
  "WT"       = "#2166AC",
  "PPM1D_KO" = "#D73027",
  "PPM1B_KO" = "#F46D43",
  "DKO"      = "#762A83"
)

cat_colors <- c(
  "Senescence"   = "#E41A1C",
  "Inflammation" = "#FF7F00",
  "Cytoskeleton" = "#984EA3",
  "Organ health" = "#4DAF4A",
  "Mitochondria" = "#377EB8",
  "Aging clock"  = "#A65628",
  "Other"        = "#CCCCCC"
)

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# =============================================================================
# 3. Модули: берём из sig_list + добавляем кастомные сенесцентные
# =============================================================================

# Загружаем кастомные генсеты из 03 скрипта если они сохранены
sen_gs_path <- "data/genesets/senescence_genesets.rds"
if (file.exists(sen_gs_path)) {
  sen_gs <- readRDS(sen_gs_path)
} else {
  sen_gs <- list()
  message("⚠ senescence_genesets.rds не найден — используем только sig_list")
}

# Объединяем: sig_list имеет приоритет
all_gs <- c(sen_gs, sig_list)
all_gs <- all_gs[!duplicated(names(all_gs))]   # убираем дубликаты (приоритет sig_list)

message("✓ Всего генсетов для скоринга: ", length(all_gs))

# Определяем рабочие модули — берём что есть
module_candidates <- c(
  "Senescence_DDR", "sen_ddr",
  "Senescence_SASP", "sen_sasp",
  "sen_core",
  "sen_ois",
  "sen_mitochondria",
  "Inflammaging",
  "cGAS_STING",
  "Rho_Activity",
  "ProxTubule",
  "Mito_Stress",
  "Tyshkovskiy_DNAmAge",   "Aging_DNAmAge",
  "Tyshkovskiy_LongevityUp",   "Aging_LongevityUp",
  "Tyshkovskiy_LongevityDown", "Aging_LongevityDown"
)

# Берём первый найденный синоним для каждого модуля
module_alias <- list(
  DDR            = c("Senescence_DDR",  "sen_ddr"),
  SASP           = c("Senescence_SASP", "sen_sasp"),
  Sen_Core       = "sen_core",
  Sen_OIS        = "sen_ois",
  Sen_Mito       = "sen_mitochondria",
  Inflammaging   = "Inflammaging",
  cGAS_STING     = "cGAS_STING",
  RhoActivity    = "Rho_Activity",
  ProxTubule     = "ProxTubule",
  Mito_Stress    = "Mito_Stress",
  Tysh_DNAmAge       = c("Tyshkovskiy_DNAmAge",      "Aging_DNAmAge"),
  Tysh_LongevityUp   = c("Tyshkovskiy_LongevityUp",  "Aging_LongevityUp"),
  Tysh_LongevityDown = c("Tyshkovskiy_LongevityDown","Aging_LongevityDown")
)

modules <- list()
for (mod_name in names(module_alias)) {
  for (alias in module_alias[[mod_name]]) {
    if (alias %in% names(all_gs)) {
      modules[[mod_name]] <- all_gs[[alias]]
      break
    }
  }
}

# Убираем не найденные
modules <- Filter(Negate(is.null), modules)
message("✓ Найдено рабочих модулей: ", length(modules), ": ",
        paste(names(modules), collapse = ", "))

# Категории для каждого модуля
module_categories <- c(
  DDR            = "Senescence",
  SASP           = "Senescence",
  Sen_Core       = "Senescence",
  Sen_OIS        = "Senescence",
  Sen_Mito       = "Senescence",
  Inflammaging   = "Inflammation",
  cGAS_STING     = "Inflammation",
  RhoActivity    = "Cytoskeleton",
  ProxTubule     = "Organ health",
  Mito_Stress    = "Mitochondria",
  Tysh_DNAmAge       = "Aging clock",
  Tysh_LongevityUp   = "Aging clock",
  Tysh_LongevityDown = "Aging clock"
)

# =============================================================================
# 4. Z-score нормализация матрицы и скоринг модулей
# =============================================================================
vst_z <- t(scale(t(vst_mat)))

# Функция скоринга (средний z-score по генсету)
score_signature_local <- function(mat_z, geneset, sig_name) {
  genes_found <- intersect(geneset, rownames(mat_z))
  if (length(genes_found) < 3) {
    message("⚠ ", sig_name, ": только ", length(genes_found),
            " генов в матрице — пропускаем")
    return(NULL)
  }
  scores <- colMeans(mat_z[genes_found, , drop = FALSE], na.rm = TRUE)
  tibble(
    sample_id = names(scores),
    score     = as.numeric(scores),
    module    = sig_name,
    n_genes   = length(genes_found)
  )
}

# Считаем скоры для всех модулей
scores_raw <- map_dfr(names(modules),
                      ~ score_signature_local(vst_z, modules[[.x]], .x))

message("✓ Скоры посчитаны: ", nrow(scores_raw), " строк")

# =============================================================================
# 5. Join с coldata
# =============================================================================
scores_all <- scores_raw %>%
  left_join(
    coldata %>%
      select(sample_id, genotype, tnf, group,
             any_of(c("n_replicates", "exploratory"))),
    by = "sample_id"
  ) %>%
  mutate(
    # СНАЧАЛА module как character — unlist на случай если list
    module   = as.character(unlist(module)),
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf,      levels = c("basal", "TNF12h")),
    # category до factor(module) — пока module ещё character
    category = module_categories[module],
    category = replace_na(category, "Other"),
    category = factor(category, levels = names(cat_colors)),
    # module — в фактор последним
    module   = factor(module, levels = names(modules))
  )

# Диагностика
message("  Строк: ",          nrow(scores_all))
message("  NA в genotype: ",  sum(is.na(scores_all$genotype)))
message("  NA в tnf: ",       sum(is.na(scores_all$tnf)))
message("  NA в module: ",    sum(is.na(scores_all$module)))
message("  NA в category: ",  sum(is.na(scores_all$category)))
print(dplyr::count(scores_all, module, category))

stopifnot(
  "NA в genotype" = sum(is.na(scores_all$genotype)) == 0,
  "NA в tnf"      = sum(is.na(scores_all$tnf))      == 0
)

# Экспорт
write_csv(scores_all, "output/tables/module_scores_all_samples.csv")

# =============================================================================
# 6. Сводная таблица mean ± SD
# =============================================================================
scores_summary <- scores_all %>%
  group_by(module, category, genotype, tnf) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score   = sd(score,  na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

write_csv(scores_summary, "output/tables/module_scores_summary.csv")
message("✓ Таблицы экспортированы")

# =============================================================================
# 7. VIZ-A: Радарный график (basal, средние по генотипу)
# =============================================================================

# Оставляем только модули с достаточным числом генов (уже отфильтровано)
radar_data <- scores_summary %>%
  filter(tnf == "basal") %>%
  select(genotype, module, mean_score) %>%
  pivot_wider(names_from = module, values_from = mean_score) %>%
  column_to_rownames("genotype")

# DKO может отсутствовать в geno_order если нет basal — проверяем
geno_in_radar <- intersect(geno_order, rownames(radar_data))

radar_max <- apply(radar_data, 2, max, na.rm = TRUE) * 1.3
radar_min <- apply(radar_data, 2, min, na.rm = TRUE) * 1.3
radar_fmsb <- rbind(radar_max, radar_min, radar_data[geno_in_radar, ])

radar_cols_border <- geno_colors[geno_in_radar]
radar_cols_fill   <- adjustcolor(radar_cols_border, alpha.f = 0.15)

draw_radar <- function() {
  par(mar = c(2, 2, 3, 2))
  fmsb::radarchart(
    radar_fmsb,
    axistype   = 1,
    pcol       = radar_cols_border,
    pfcol      = radar_cols_fill,
    plwd       = 2.5,
    plty       = 1,
    cglcol     = "grey70",
    cglty      = 1,
    cglwd      = 0.5,
    axislabcol = "grey40",
    vlcex      = 0.78,
    title      = "Module scores — basal (mean z-score)"
  )
  legend(
    x      = 1.15, y = 1.25,
    legend = geno_in_radar,
    col    = radar_cols_border,
    lwd    = 2.5, lty = 1,
    bty    = "n", cex = 0.85
  )
}

png("output/figures/module_clocks/VIZ_A_radar_basal.png",
    width = 9, height = 8, units = "in", res = 300, bg = "white")
draw_radar()
dev.off()

pdf("output/figures/module_clocks/VIZ_A_radar_basal.pdf",
    width = 9, height = 8)
draw_radar()
dev.off()

message("✓ VIZ-A Radar сохранён")

# =============================================================================
# 8. VIZ-B: Ridge plot — распределение скоров, только basal
# =============================================================================
p_ridge <- scores_all %>%
  filter(tnf == "basal") %>%
  mutate(
    genotype = factor(genotype, levels = rev(geno_order)),
    module   = factor(module,   levels = rev(levels(scores_all$module)))
  ) %>%
  ggplot(aes(x = score, y = module,
             fill = genotype, color = genotype)) +
  ggridges::geom_density_ridges(
    alpha           = 0.55,
    scale           = 0.85,
    jittered_points = TRUE,
    point_size      = 2,
    point_alpha     = 0.9,
    position        = ggridges::position_raincloud(
      width  = 0.05,
      height = 0.15
    )
  ) +
  scale_fill_manual(values  = geno_colors, name = "Genotype") +
  scale_color_manual(values = geno_colors, name = "Genotype") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Module score distributions — basal",
    subtitle = "Mean z-score across signature genes | MC38 cells",
    x = "Module score (z-score)", y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(size = 10)
  )

ggsave("output/figures/module_clocks/VIZ_B_ridge_basal.png",
       p_ridge, width = 11, height = 9, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_B_ridge_basal.pdf",
       p_ridge, width = 11, height = 9)

message("✓ VIZ-B Ridge сохранён")

# =============================================================================
# 9. VIZ-C: Slope plot — basal → TNF12h по генотипам
# =============================================================================
# DKO не имеет TNF-образца, поэтому для него только basal-точка
slope_data <- scores_summary %>%
  mutate(
    tnf      = factor(tnf, levels = c("basal", "TNF12h")),
    genotype = factor(genotype, levels = geno_order)
  )

p_slope <- ggplot(slope_data,
                  aes(x     = tnf,
                      y     = mean_score,
                      group = interaction(genotype, module),
                      color = genotype)) +
  geom_line(aes(linetype = genotype), linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_score - sd_score,
        ymax = mean_score + sd_score),
    width = 0.1, alpha = 0.5
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  scale_color_manual(values = geno_colors, name = "Genotype") +
  scale_linetype_manual(
    values = c("WT"       = "solid",
               "PPM1D_KO" = "dashed",
               "PPM1B_KO" = "dotdash",
               "DKO"      = "dotted"),
    name = "Genotype"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  labs(
    title    = "TNF effect on module scores by genotype",
    subtitle = "Mean ± SD; basal → TNF12h | DKO: basal only",
    x = NULL, y = "Module score (z-score)"
  ) +
  theme_pub +
  theme(
    strip.text      = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "bottom",
    axis.text.x     = element_text(size = 9)
  )

ggsave("output/figures/module_clocks/VIZ_C_slope_tnf_effect.png",
       p_slope, width = 14, height = 10, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_C_slope_tnf_effect.pdf",
       p_slope, width = 14, height = 10)

message("✓ VIZ-C Slope сохранён")

# =============================================================================
# 10. VIZ-D: Heatmap всех скоров — фиксированный порядок колонок
# =============================================================================
score_mat <- scores_all %>%
  select(sample_id, module, score) %>%
  pivot_wider(names_from = module, values_from = score) %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()   # модули × образцы

# Порядок образцов по группам
sample_order_df <- coldata %>%
  filter(sample_id %in% colnames(score_mat)) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group, sample_id)

sample_order <- sample_order_df$sample_id
score_mat_ord <- score_mat[, sample_order]

# Аннотация столбцов
col_ann_df <- sample_order_df %>%
  select(genotype, tnf) %>%
  mutate(
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf,      levels = c("basal", "TNF12h"))
  ) %>%
  as.data.frame()
rownames(col_ann_df) <- sample_order

# Аннотация строк
row_ann_df <- data.frame(
  category = module_categories[rownames(score_mat_ord)],
  row.names = rownames(score_mat_ord)
)
row_ann_df$category[is.na(row_ann_df$category)] <- "Other"

# Разделители между группами — надёжный вариант
gap_positions <- vapply(
  group_order,
  function(g) sum(as.character(sample_order_df$group) == g),
  integer(1)
)
gap_positions <- head(cumsum(gap_positions), -1)
gap_positions <- gap_positions[gap_positions > 0]   # убираем нули (группы без образцов)

message("gap_positions: ", paste(gap_positions, collapse = ", "))

# Цвета аннотаций — именованные векторы
used_cat   <- unique(row_ann_df$category)
ann_colors <- list(
  genotype = geno_colors[geno_order],
  tnf      = c("basal" = "#F0F0F0", "TNF12h" = "#FC8D59"),
  category = cat_colors[used_cat]
)

# Палитра heatmap
score_range <- max(abs(range(score_mat_ord, na.rm = TRUE)), na.rm = TRUE)
score_range <- ceiling(score_range * 10) / 10
pal100      <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
my_breaks   <- seq(-score_range, score_range, length.out = 101)

pheatmap::pheatmap(
  score_mat_ord,
  annotation_col    = col_ann_df,
  annotation_row    = row_ann_df,
  annotation_colors = ann_colors,
  color             = pal100,
  breaks            = my_breaks,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  gaps_col          = gap_positions,
  show_colnames     = TRUE,
  cellwidth         = 22,
  cellheight        = 18,
  fontsize          = 9,
  fontsize_row      = 9,
  fontsize_col      = 7.5,
  border_color      = NA,
  main              = "Module scores — all samples (ordered by genotype × treatment)",
  filename          = "output/figures/module_clocks/VIZ_D_heatmap_scores.png",
  width = 14, height = 7
)

message("✓ VIZ-D Heatmap сохранён")

# =============================================================================
# 11. VIZ-E: Фазовый портрет DDR × Inflammaging
# =============================================================================

# Проверяем наличие обоих модулей
if (all(c("DDR", "Inflammaging") %in% names(modules))) {

  phase_ddr <- scores_all %>%
    filter(module %in% c("DDR", "Inflammaging")) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    group_by(genotype, tnf) %>%
    summarise(
      DDR_mean = mean(DDR,          na.rm = TRUE),
      Inf_mean = mean(Inflammaging, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    pivot_wider(
      names_from  = tnf,
      values_from = c(DDR_mean, Inf_mean)
    )

  # Реальные имена колонок: DDR_mean_basal, DDR_mean_TNF12h, Inf_mean_basal, Inf_mean_TNF12h
  message("Колонки phase_ddr: ", paste(names(phase_ddr), collapse = ", "))

  # Генотипы с TNF-данными (WT, PPM1D_KO, PPM1B_KO)
  arrows_ddr <- phase_ddr %>%
    filter(!is.na(DDR_mean_TNF12h)) %>%
    mutate(genotype = factor(genotype, levels = geno_order))

  # DKO — только basal
  dko_ddr <- phase_ddr %>%
    filter(genotype == "DKO")

  p_phase_ddr <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # Стрелки basal → TNF12h
    geom_segment(
      data = arrows_ddr,
      aes(x    = DDR_mean_basal, y    = Inf_mean_basal,
          xend = DDR_mean_TNF12h, yend = Inf_mean_TNF12h,
          color = genotype),
      arrow     = arrow(length = unit(0.28, "cm"), type = "closed"),
      linewidth = 1.3, alpha = 0.85
    ) +
    # Точки basal (генотипы со стрелками)
    geom_point(
      data = arrows_ddr,
      aes(x = DDR_mean_basal, y = Inf_mean_basal, color = genotype),
      size = 5, shape = 16
    ) +
    # Точки TNF12h (конец стрелки)
    geom_point(
      data = arrows_ddr,
      aes(x = DDR_mean_TNF12h, y = Inf_mean_TNF12h, color = genotype),
      size = 3.5, shape = 1, stroke = 1.5
    ) +
    # DKO — треугольник, только basal
    geom_point(
      data  = dko_ddr,
      aes(x = DDR_mean_basal, y = Inf_mean_basal),
      color = geno_colors["DKO"],
      size  = 5, shape = 17
    ) +
    # Подписи
    ggrepel::geom_label_repel(
      data = phase_ddr,
      aes(x     = DDR_mean_basal,
          y     = Inf_mean_basal,
          label = genotype,
          color = genotype),
      size          = 3.2,
      box.padding   = 0.5,
      show.legend   = FALSE,
      seed          = 42
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = "● basal   ○ TNF12h   ▲ DKO basal only",
             hjust = 1.05, vjust = -0.5,
             size = 2.6, color = "grey40") +
    scale_color_manual(values = geno_colors, name = "Genotype") +
    coord_cartesian(clip = "off") +
    labs(
      title    = "Phase portrait: DDR × Inflammaging",
      subtitle = "● basal  ○ TNF12h  ▲ DKO (no TNF data) | MC38",
      x = "DDR module score (z-score)",
      y = "Inflammaging module score (z-score)"
    ) +
    theme_pub +
    theme(plot.margin = margin(10, 25, 10, 10))

  ggsave("output/figures/module_clocks/VIZ_E_phase_DDR_Inflammaging.png",
         p_phase_ddr, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_E_phase_DDR_Inflammaging.pdf",
         p_phase_ddr, width = 7, height = 6)

  write_csv(phase_ddr, "output/tables/phase_DDR_Inflammaging.csv")
  message("✓ VIZ-E Phase DDR×Inflammaging сохранён")

} else {
  message("⚠ VIZ-E пропущен: DDR или Inflammaging не найдены в modules")
}

# =============================================================================
# 12. VIZ-F: Фазовый портрет RhoActivity × SASP
# =============================================================================

if (all(c("RhoActivity", "SASP") %in% names(modules))) {

  phase_rho <- scores_all %>%
    filter(module %in% c("RhoActivity", "SASP")) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    group_by(genotype, tnf) %>%
    summarise(
      Rho_mean  = mean(RhoActivity, na.rm = TRUE),
      SASP_mean = mean(SASP,        na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    pivot_wider(
      names_from  = tnf,
      values_from = c(Rho_mean, SASP_mean)
    )

  message("Колонки phase_rho: ", paste(names(phase_rho), collapse = ", "))

  arrows_rho <- phase_rho %>%
    filter(!is.na(Rho_mean_TNF12h)) %>%
    mutate(genotype = factor(genotype, levels = geno_order))

  dko_rho <- phase_rho %>%
    filter(genotype == "DKO")

  p_phase_rho <- ggplot() +
    # Квадранты — биологические метки
    annotate("text", x = -Inf, y =  Inf,
             label = "Low Rho\nHigh SASP",
             hjust = -0.1, vjust = 1.3,
             size = 2.8, color = "grey60", fontface = "italic") +
    annotate("text", x =  Inf, y =  Inf,
             label = "High Rho\nHigh SASP",
             hjust = 1.1, vjust = 1.3,
             size = 2.8, color = "#B2182B", fontface = "bold.italic") +
    annotate("text", x = -Inf, y = -Inf,
             label = "Low Rho\nLow SASP",
             hjust = -0.1, vjust = -0.3,
             size = 2.8, color = "grey60", fontface = "italic") +
    annotate("text", x =  Inf, y = -Inf,
             label = "High Rho\nLow SASP",
             hjust = 1.1, vjust = -0.3,
             size = 2.8, color = "grey50", fontface = "italic") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # Стрелки
    geom_segment(
      data = arrows_rho,
      aes(x    = Rho_mean_basal,  y    = SASP_mean_basal,
          xend = Rho_mean_TNF12h, yend = SASP_mean_TNF12h,
          color = genotype),
      arrow     = arrow(length = unit(0.28, "cm"), type = "closed"),
      linewidth = 1.3, alpha = 0.85
    ) +
    # Точки basal
    geom_point(
      data = arrows_rho,
      aes(x = Rho_mean_basal, y = SASP_mean_basal, color = genotype),
      size = 5, shape = 16
    ) +
    # Точки TNF12h
    geom_point(
      data = arrows_rho,
      aes(x = Rho_mean_TNF12h, y = SASP_mean_TNF12h, color = genotype),
      size = 3.5, shape = 1, stroke = 1.5
    ) +
    # DKO — только basal
    geom_point(
      data  = dko_rho,
      aes(x = Rho_mean_basal, y = SASP_mean_basal),
      color = geno_colors["DKO"],
      size  = 5, shape = 17
    ) +
    # Подписи
    ggrepel::geom_label_repel(
      data = phase_rho,
      aes(x     = Rho_mean_basal,
          y     = SASP_mean_basal,
          label = genotype,
          color = genotype),
      size          = 3.2,
      box.padding   = 0.5,
      show.legend   = FALSE,
      seed          = 7
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = "● basal   ○ TNF12h   ▲ DKO basal only",
             hjust = 1.05, vjust = -0.5,
             size = 2.6, color = "grey40") +
    scale_color_manual(values = geno_colors, name = "Genotype") +
    coord_cartesian(clip = "off") +
    labs(
      title    = "Phase portrait: RhoActivity × SASP",
      subtitle = "● basal  ○ TNF12h  ▲ DKO (no TNF data) | MC38",
      x = "Rho Activity module score (z-score)",
      y = "SASP module score (z-score)"
    ) +
    theme_pub +
    theme(plot.margin = margin(10, 25, 10, 10))

  ggsave("output/figures/module_clocks/VIZ_F_phase_Rho_SASP.png",
         p_phase_rho, width = 7.5, height = 6.5, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_F_phase_Rho_SASP.pdf",
         p_phase_rho, width = 7.5, height = 6.5)

  write_csv(phase_rho, "output/tables/phase_Rho_SASP.csv")
  message("✓ VIZ-F Phase Rho×SASP сохранён")

} else {
  message("⚠ VIZ-F пропущен: RhoActivity или SASP не найдены в modules")
}

# =============================================================================
# 13. VIZ-G: DKO additivity test — панель по модулям
# =============================================================================
additivity <- scores_summary %>%
  filter(tnf == "basal") %>%
  select(module, genotype, mean_score) %>%
  pivot_wider(names_from = genotype, values_from = mean_score)

# Проверяем наличие всех колонок
req_geno <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")
miss     <- setdiff(req_geno, names(additivity))

if (length(miss) > 0) {
  message("⚠ VIZ-G: отсутствуют генотипы: ", paste(miss, collapse = ", "))
} else {

  additivity <- additivity %>%
    mutate(
      expected  = (PPM1D_KO - WT) + (PPM1B_KO - WT) + WT,
      observed  = DKO,
      ratio     = observed / expected,
      interp    = case_when(
        ratio > 1.2  ~ "Synergy (>120%)",
        ratio < 0.8  ~ "Compensation (<80%)",
        TRUE         ~ "Additive (80–120%)"
      )
    )

  write_csv(additivity, "output/tables/dko_additivity_test.csv")

  p_additivity <- additivity %>%
    select(module, observed, expected, wt = WT) %>%
    pivot_longer(c(observed, expected),
                 names_to = "type", values_to = "score") %>%
    ggplot(aes(x = module, y = score, fill = type)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_hline(
      data = additivity %>% select(module, wt = WT),
      aes(yintercept = wt),
      linetype = "dashed", color = geno_colors["WT"],
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = c("observed" = "#762A83", "expected" = "#C994C7"),
      labels = c("observed" = "Observed DKO",
                 "expected" = "Expected additive"),
      name = NULL
    ) +
    labs(
      title    = "DKO: observed vs. expected additive module score",
      subtitle = "Purple = observed DKO | Light = PPM1D_KO + PPM1B_KO additive model | Dashed = WT",
      x = NULL, y = "Mean z-score"
    ) +
    theme_pub +
    theme(
      axis.text.x     = element_text(angle = 40, hjust = 1),
      legend.position = "top"
    )

  ggsave("output/figures/module_clocks/VIZ_G_dko_additivity.png",
         p_additivity, width = 12, height = 6, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_G_dko_additivity.pdf",
         p_additivity, width = 12, height = 6)

  message("✓ VIZ-G DKO additivity сохранён")
}

# =============================================================================
# 14. Итог
# =============================================================================
message("\n══════════════════════════════════════════════════════════")
message("✓ 04_module_clocks_v3.R ЗАВЕРШЁН")
message("  Таблицы:")
message("    output/tables/module_scores_all_samples.csv")
message("    output/tables/module_scores_summary.csv")
message("    output/tables/phase_DDR_Inflammaging.csv")
message("    output/tables/phase_Rho_SASP.csv")
message("    output/tables/dko_additivity_test.csv")
message("  Рисунки:")
message("    VIZ_A — радар (basal, mean по генотипу)")
message("    VIZ_B — ridge plot (basal)")
message("    VIZ_C — slope plot (basal→TNF12h)")
message("    VIZ_D — heatmap всех образцов")
message("    VIZ_E — фазовый портрет DDR×Inflammaging")
message("    VIZ_F — фазовый портрет Rho×SASP")
message("    VIZ_G — DKO additivity test")
message("══════════════════════════════════════════════════════════")
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   "Module clocks" — многомерные биологические часы:
#   DDR, Inflammaging, Rho, SASP, ProxTubule и др. как спицы радара.
#   Позволяет увидеть ПРОФИЛЬ состояния клетки целиком.
#
#   ВОПРОСЫ:
#   1. DOUBLE KO (DKO): синергия vs. аддитивность vs. компенсация?
#   2. TNF-ТРАЕКТОРИЯ: длина и направление вектора basal→TNF12h
#      (гипотеза: у KO вектор длиннее, направлен в сторону SASP)
#   3. ФАЗОВЫЕ ПОРТРЕТЫ: DDR×Inflammaging и Rho×SASP
#
# ИСПРАВЛЕНИЯ v3:
#   - PPM1D_KO / DKO / PPM1B_KO — реальные имена из coldata
#   - group_order использует реальные уровни: DKO_basal, PPM1D_KO_basal, etc.
#   - убран mod_scores (не определён), везде scores_all
#   - phase_data: правильные имена колонок после pivot_wider
#   - убраны GENO_COLORS / WipKO / Double_KO
#   - score_signature определена внутри скрипта (не зависит от 03)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fmsb)
  library(ggrepel)
  library(patchwork)
  library(ggridges)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# =============================================================================
# 0. Создаём папки
# =============================================================================
dir.create("output/figures/module_clocks", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",               recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Загружаем pipeline
# =============================================================================
pipeline <- readRDS("output/rds/deseq2_all_results.rds")

vst_mat  <- as.matrix(pipeline$norm_counts)   # матрица VST: гены × образцы
coldata  <- pipeline$coldata                   # tibble: sample_id, genotype, tnf, group, ...
sig_list <- pipeline$sig_list                  # генсеты из 02/03 скриптов

message("✓ vst_mat: ", nrow(vst_mat), " генов × ", ncol(vst_mat), " образцов")
message("✓ sig_list генсеты: ", paste(names(sig_list), collapse = ", "))

# Упорядочиваем столбцы матрицы по coldata
stopifnot(all(coldata$sample_id %in% colnames(vst_mat)))
vst_mat <- vst_mat[, coldata$sample_id]

# =============================================================================
# 2. Реальные уровни факторов (из coldata)
# =============================================================================
# genotype: "WT" "PPM1D_KO" "PPM1B_KO" "DKO"
# group:    "WT_basal" "PPM1D_KO_basal" "PPM1B_KO_basal" "DKO_basal"
#           "WT_TNF12h" "PPM1D_KO_TNF12h" "PPM1B_KO_TNF12h"
# tnf:      "basal" "TNF12h"

geno_order  <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")
group_order <- c(
  "WT_basal", "PPM1D_KO_basal", "PPM1B_KO_basal", "DKO_basal",
  "WT_TNF12h", "PPM1D_KO_TNF12h", "PPM1B_KO_TNF12h"
)

geno_colors <- c(
  "WT"       = "#2166AC",
  "PPM1D_KO" = "#D73027",
  "PPM1B_KO" = "#F46D43",
  "DKO"      = "#762A83"
)

cat_colors <- c(
  "Senescence"   = "#E41A1C",
  "Inflammation" = "#FF7F00",
  "Cytoskeleton" = "#984EA3",
  "Organ health" = "#4DAF4A",
  "Mitochondria" = "#377EB8",
  "Aging clock"  = "#A65628",
  "Other"        = "#CCCCCC"
)

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# =============================================================================
# 3. Модули: берём из sig_list + добавляем кастомные сенесцентные
# =============================================================================

# Загружаем кастомные генсеты из 03 скрипта если они сохранены
sen_gs_path <- "data/genesets/senescence_genesets.rds"
if (file.exists(sen_gs_path)) {
  sen_gs <- readRDS(sen_gs_path)
} else {
  sen_gs <- list()
  message("⚠ senescence_genesets.rds не найден — используем только sig_list")
}

# Объединяем: sig_list имеет приоритет
all_gs <- c(sen_gs, sig_list)
all_gs <- all_gs[!duplicated(names(all_gs))]   # убираем дубликаты (приоритет sig_list)

message("✓ Всего генсетов для скоринга: ", length(all_gs))

# Определяем рабочие модули — берём что есть
module_candidates <- c(
  "Senescence_DDR", "sen_ddr",
  "Senescence_SASP", "sen_sasp",
  "sen_core",
  "sen_ois",
  "sen_mitochondria",
  "Inflammaging",
  "cGAS_STING",
  "Rho_Activity",
  "ProxTubule",
  "Mito_Stress",
  "Tyshkovskiy_DNAmAge",   "Aging_DNAmAge",
  "Tyshkovskiy_LongevityUp",   "Aging_LongevityUp",
  "Tyshkovskiy_LongevityDown", "Aging_LongevityDown"
)

# Берём первый найденный синоним для каждого модуля
module_alias <- list(
  DDR            = c("Senescence_DDR",  "sen_ddr"),
  SASP           = c("Senescence_SASP", "sen_sasp"),
  Sen_Core       = "sen_core",
  Sen_OIS        = "sen_ois",
  Sen_Mito       = "sen_mitochondria",
  Inflammaging   = "Inflammaging",
  cGAS_STING     = "cGAS_STING",
  RhoActivity    = "Rho_Activity",
  ProxTubule     = "ProxTubule",
  Mito_Stress    = "Mito_Stress",
  Tysh_DNAmAge       = c("Tyshkovskiy_DNAmAge",      "Aging_DNAmAge"),
  Tysh_LongevityUp   = c("Tyshkovskiy_LongevityUp",  "Aging_LongevityUp"),
  Tysh_LongevityDown = c("Tyshkovskiy_LongevityDown","Aging_LongevityDown")
)

modules <- list()
for (mod_name in names(module_alias)) {
  for (alias in module_alias[[mod_name]]) {
    if (alias %in% names(all_gs)) {
      modules[[mod_name]] <- all_gs[[alias]]
      break
    }
  }
}

# Убираем не найденные
modules <- Filter(Negate(is.null), modules)
message("✓ Найдено рабочих модулей: ", length(modules), ": ",
        paste(names(modules), collapse = ", "))

# Категории для каждого модуля
module_categories <- c(
  DDR            = "Senescence",
  SASP           = "Senescence",
  Sen_Core       = "Senescence",
  Sen_OIS        = "Senescence",
  Sen_Mito       = "Senescence",
  Inflammaging   = "Inflammation",
  cGAS_STING     = "Inflammation",
  RhoActivity    = "Cytoskeleton",
  ProxTubule     = "Organ health",
  Mito_Stress    = "Mitochondria",
  Tysh_DNAmAge       = "Aging clock",
  Tysh_LongevityUp   = "Aging clock",
  Tysh_LongevityDown = "Aging clock"
)

# =============================================================================
# 4. Z-score нормализация матрицы и скоринг модулей
# =============================================================================
vst_z <- t(scale(t(vst_mat)))

# Функция скоринга (средний z-score по генсету)
score_signature_local <- function(mat_z, geneset, sig_name) {
  genes_found <- intersect(geneset, rownames(mat_z))
  if (length(genes_found) < 3) {
    message("⚠ ", sig_name, ": только ", length(genes_found),
            " генов в матрице — пропускаем")
    return(NULL)
  }
  scores <- colMeans(mat_z[genes_found, , drop = FALSE], na.rm = TRUE)
  tibble(
    sample_id = names(scores),
    score     = as.numeric(scores),
    module    = sig_name,
    n_genes   = length(genes_found)
  )
}

# Считаем скоры для всех модулей
scores_raw <- map_dfr(names(modules),
                      ~ score_signature_local(vst_z, modules[[.x]], .x))

message("✓ Скоры посчитаны: ", nrow(scores_raw), " строк")

# =============================================================================
# 5. Join с coldata
# =============================================================================
scores_all <- scores_raw %>%
  left_join(
    coldata %>%
      select(sample_id, genotype, tnf, group,
             any_of(c("n_replicates", "exploratory"))),
    by = "sample_id"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf,      levels = c("basal", "TNF12h")),
    module   = factor(module,   levels = names(modules)),
    category = factor(
      module_categories[as.character(module)],
      levels = names(cat_colors)
    )
  )

# Диагностика
message("\n── Проверка scores_all ──────────────────────────────────")
message("  Строк: ", nrow(scores_all))
message("  NA в genotype: ", sum(is.na(scores_all$genotype)))
message("  NA в tnf: ",      sum(is.na(scores_all$tnf)))
print(count(scores_all, module))

stopifnot(
  "NA в genotype" = sum(is.na(scores_all$genotype)) == 0,
  "NA в tnf"      = sum(is.na(scores_all$tnf))      == 0
)

# Экспорт
write_csv(scores_all, "output/tables/module_scores_all_samples.csv")

# =============================================================================
# 6. Сводная таблица mean ± SD
# =============================================================================
scores_summary <- scores_all %>%
  group_by(module, category, genotype, tnf) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score   = sd(score,  na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

write_csv(scores_summary, "output/tables/module_scores_summary.csv")
message("✓ Таблицы экспортированы")

# =============================================================================
# 7. VIZ-A: Радарный график (basal, средние по генотипу)
# =============================================================================

# Оставляем только модули с достаточным числом генов (уже отфильтровано)
radar_data <- scores_summary %>%
  filter(tnf == "basal") %>%
  select(genotype, module, mean_score) %>%
  pivot_wider(names_from = module, values_from = mean_score) %>%
  column_to_rownames("genotype")

# DKO может отсутствовать в geno_order если нет basal — проверяем
geno_in_radar <- intersect(geno_order, rownames(radar_data))

radar_max <- apply(radar_data, 2, max, na.rm = TRUE) * 1.3
radar_min <- apply(radar_data, 2, min, na.rm = TRUE) * 1.3
radar_fmsb <- rbind(radar_max, radar_min, radar_data[geno_in_radar, ])

radar_cols_border <- geno_colors[geno_in_radar]
radar_cols_fill   <- adjustcolor(radar_cols_border, alpha.f = 0.15)

draw_radar <- function() {
  par(mar = c(2, 2, 3, 2))
  fmsb::radarchart(
    radar_fmsb,
    axistype   = 1,
    pcol       = radar_cols_border,
    pfcol      = radar_cols_fill,
    plwd       = 2.5,
    plty       = 1,
    cglcol     = "grey70",
    cglty      = 1,
    cglwd      = 0.5,
    axislabcol = "grey40",
    vlcex      = 0.78,
    title      = "Module scores — basal (mean z-score)"
  )
  legend(
    x      = 1.15, y = 1.25,
    legend = geno_in_radar,
    col    = radar_cols_border,
    lwd    = 2.5, lty = 1,
    bty    = "n", cex = 0.85
  )
}

png("output/figures/module_clocks/VIZ_A_radar_basal.png",
    width = 9, height = 8, units = "in", res = 300, bg = "white")
draw_radar()
dev.off()

pdf("output/figures/module_clocks/VIZ_A_radar_basal.pdf",
    width = 9, height = 8)
draw_radar()
dev.off()

message("✓ VIZ-A Radar сохранён")

# =============================================================================
# 8. VIZ-B: Ridge plot — распределение скоров, только basal
# =============================================================================
p_ridge <- scores_all %>%
  filter(tnf == "basal") %>%
  mutate(
    genotype = factor(genotype, levels = rev(geno_order)),
    module   = factor(module,   levels = rev(levels(scores_all$module)))
  ) %>%
  ggplot(aes(x = score, y = module,
             fill = genotype, color = genotype)) +
  ggridges::geom_density_ridges(
    alpha           = 0.55,
    scale           = 0.85,
    jittered_points = TRUE,
    point_size      = 2,
    point_alpha     = 0.9,
    position        = ggridges::position_raincloud(
      width  = 0.05,
      height = 0.15
    )
  ) +
  scale_fill_manual(values  = geno_colors, name = "Genotype") +
  scale_color_manual(values = geno_colors, name = "Genotype") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Module score distributions — basal",
    subtitle = "Mean z-score across signature genes | MC38 cells",
    x = "Module score (z-score)", y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(size = 10)
  )

ggsave("output/figures/module_clocks/VIZ_B_ridge_basal.png",
       p_ridge, width = 11, height = 9, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_B_ridge_basal.pdf",
       p_ridge, width = 11, height = 9)

message("✓ VIZ-B Ridge сохранён")

# =============================================================================
# 9. VIZ-C: Slope plot — basal → TNF12h по генотипам
# =============================================================================
# DKO не имеет TNF-образца, поэтому для него только basal-точка
slope_data <- scores_summary %>%
  mutate(
    tnf      = factor(tnf, levels = c("basal", "TNF12h")),
    genotype = factor(genotype, levels = geno_order)
  )

p_slope <- ggplot(slope_data,
                  aes(x     = tnf,
                      y     = mean_score,
                      group = interaction(genotype, module),
                      color = genotype)) +
  geom_line(aes(linetype = genotype), linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_score - sd_score,
        ymax = mean_score + sd_score),
    width = 0.1, alpha = 0.5
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  scale_color_manual(values = geno_colors, name = "Genotype") +
  scale_linetype_manual(
    values = c("WT"       = "solid",
               "PPM1D_KO" = "dashed",
               "PPM1B_KO" = "dotdash",
               "DKO"      = "dotted"),
    name = "Genotype"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  labs(
    title    = "TNF effect on module scores by genotype",
    subtitle = "Mean ± SD; basal → TNF12h | DKO: basal only",
    x = NULL, y = "Module score (z-score)"
  ) +
  theme_pub +
  theme(
    strip.text      = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "bottom",
    axis.text.x     = element_text(size = 9)
  )

ggsave("output/figures/module_clocks/VIZ_C_slope_tnf_effect.png",
       p_slope, width = 14, height = 10, dpi = 300, bg = "white")
ggsave("output/figures/module_clocks/VIZ_C_slope_tnf_effect.pdf",
       p_slope, width = 14, height = 10)

message("✓ VIZ-C Slope сохранён")

# =============================================================================
# 10. VIZ-D: Heatmap всех скоров — фиксированный порядок колонок
# =============================================================================
score_mat <- scores_all %>%
  select(sample_id, module, score) %>%
  pivot_wider(names_from = module, values_from = score) %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()   # модули × образцы

# Порядок образцов по группам
sample_order_df <- coldata %>%
  filter(sample_id %in% colnames(score_mat)) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group, sample_id)

sample_order <- sample_order_df$sample_id
score_mat_ord <- score_mat[, sample_order]

# Аннотация столбцов
col_ann_df <- sample_order_df %>%
  select(genotype, tnf) %>%
  mutate(
    genotype = factor(genotype, levels = geno_order),
    tnf      = factor(tnf,      levels = c("basal", "TNF12h"))
  ) %>%
  as.data.frame()
rownames(col_ann_df) <- sample_order

# Аннотация строк
row_ann_df <- data.frame(
  category = module_categories[rownames(score_mat_ord)],
  row.names = rownames(score_mat_ord)
)
row_ann_df$category[is.na(row_ann_df$category)] <- "Other"

# Разделители между группами
gap_positions <- sample_order_df %>%
  mutate(grp = as.character(group)) %>%
  group_by(grp) %>%
  summarise(n = n(), .groups = "drop") %>%
  # сохраняем порядок group_order
  slice(match(intersect(group_order, unique(grp)), grp)) %>%
  mutate(gap = cumsum(n)) %>%
  pull(gap) %>%
  head(-1)

# Цвета аннотаций — именованные векторы
used_cat   <- unique(row_ann_df$category)
ann_colors <- list(
  genotype = geno_colors[geno_order],
  tnf      = c("basal" = "#F0F0F0", "TNF12h" = "#FC8D59"),
  category = cat_colors[used_cat]
)

# Палитра heatmap
score_range <- max(abs(range(score_mat_ord, na.rm = TRUE)), na.rm = TRUE)
score_range <- ceiling(score_range * 10) / 10
pal100      <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
my_breaks   <- seq(-score_range, score_range, length.out = 101)

pheatmap::pheatmap(
  score_mat_ord,
  annotation_col    = col_ann_df,
  annotation_row    = row_ann_df,
  annotation_colors = ann_colors,
  color             = pal100,
  breaks            = my_breaks,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  gaps_col          = gap_positions,
  show_colnames     = TRUE,
  cellwidth         = 22,
  cellheight        = 18,
  fontsize          = 9,
  fontsize_row      = 9,
  fontsize_col      = 7.5,
  border_color      = NA,
  main              = "Module scores — all samples (ordered by genotype × treatment)",
  filename          = "output/figures/module_clocks/VIZ_D_heatmap_scores.png",
  width = 14, height = 7
)

message("✓ VIZ-D Heatmap сохранён")

# =============================================================================
# 11. VIZ-E: Фазовый портрет DDR × Inflammaging
# =============================================================================

# Проверяем наличие обоих модулей
if (all(c("DDR", "Inflammaging") %in% names(modules))) {

  phase_ddr <- scores_all %>%
    filter(module %in% c("DDR", "Inflammaging")) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    group_by(genotype, tnf) %>%
    summarise(
      DDR_mean = mean(DDR,          na.rm = TRUE),
      Inf_mean = mean(Inflammaging, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    pivot_wider(
      names_from  = tnf,
      values_from = c(DDR_mean, Inf_mean)
    )

  # Реальные имена колонок: DDR_mean_basal, DDR_mean_TNF12h, Inf_mean_basal, Inf_mean_TNF12h
  message("Колонки phase_ddr: ", paste(names(phase_ddr), collapse = ", "))

  # Генотипы с TNF-данными (WT, PPM1D_KO, PPM1B_KO)
  arrows_ddr <- phase_ddr %>%
    filter(!is.na(DDR_mean_TNF12h)) %>%
    mutate(genotype = factor(genotype, levels = geno_order))

  # DKO — только basal
  dko_ddr <- phase_ddr %>%
    filter(genotype == "DKO")

  p_phase_ddr <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # Стрелки basal → TNF12h
    geom_segment(
      data = arrows_ddr,
      aes(x    = DDR_mean_basal, y    = Inf_mean_basal,
          xend = DDR_mean_TNF12h, yend = Inf_mean_TNF12h,
          color = genotype),
      arrow     = arrow(length = unit(0.28, "cm"), type = "closed"),
      linewidth = 1.3, alpha = 0.85
    ) +
    # Точки basal (генотипы со стрелками)
    geom_point(
      data = arrows_ddr,
      aes(x = DDR_mean_basal, y = Inf_mean_basal, color = genotype),
      size = 5, shape = 16
    ) +
    # Точки TNF12h (конец стрелки)
    geom_point(
      data = arrows_ddr,
      aes(x = DDR_mean_TNF12h, y = Inf_mean_TNF12h, color = genotype),
      size = 3.5, shape = 1, stroke = 1.5
    ) +
    # DKO — треугольник, только basal
    geom_point(
      data  = dko_ddr,
      aes(x = DDR_mean_basal, y = Inf_mean_basal),
      color = geno_colors["DKO"],
      size  = 5, shape = 17
    ) +
    # Подписи
    ggrepel::geom_label_repel(
      data = phase_ddr,
      aes(x     = DDR_mean_basal,
          y     = Inf_mean_basal,
          label = genotype,
          color = genotype),
      size          = 3.2,
      box.padding   = 0.5,
      show.legend   = FALSE,
      seed          = 42
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = "● basal   ○ TNF12h   ▲ DKO basal only",
             hjust = 1.05, vjust = -0.5,
             size = 2.6, color = "grey40") +
    scale_color_manual(values = geno_colors, name = "Genotype") +
    coord_cartesian(clip = "off") +
    labs(
      title    = "Phase portrait: DDR × Inflammaging",
      subtitle = "● basal  ○ TNF12h  ▲ DKO (no TNF data) | MC38",
      x = "DDR module score (z-score)",
      y = "Inflammaging module score (z-score)"
    ) +
    theme_pub +
    theme(plot.margin = margin(10, 25, 10, 10))

  ggsave("output/figures/module_clocks/VIZ_E_phase_DDR_Inflammaging.png",
         p_phase_ddr, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_E_phase_DDR_Inflammaging.pdf",
         p_phase_ddr, width = 7, height = 6)

  write_csv(phase_ddr, "output/tables/phase_DDR_Inflammaging.csv")
  message("✓ VIZ-E Phase DDR×Inflammaging сохранён")

} else {
  message("⚠ VIZ-E пропущен: DDR или Inflammaging не найдены в modules")
}

# =============================================================================
# 12. VIZ-F: Фазовый портрет RhoActivity × SASP
# =============================================================================

if (all(c("RhoActivity", "SASP") %in% names(modules))) {

  phase_rho <- scores_all %>%
    filter(module %in% c("RhoActivity", "SASP")) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    group_by(genotype, tnf) %>%
    summarise(
      Rho_mean  = mean(RhoActivity, na.rm = TRUE),
      SASP_mean = mean(SASP,        na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    pivot_wider(
      names_from  = tnf,
      values_from = c(Rho_mean, SASP_mean)
    )

  message("Колонки phase_rho: ", paste(names(phase_rho), collapse = ", "))

  arrows_rho <- phase_rho %>%
    filter(!is.na(Rho_mean_TNF12h)) %>%
    mutate(genotype = factor(genotype, levels = geno_order))

  dko_rho <- phase_rho %>%
    filter(genotype == "DKO")

  p_phase_rho <- ggplot() +
    # Квадранты — биологические метки
    annotate("text", x = -Inf, y =  Inf,
             label = "Low Rho\nHigh SASP",
             hjust = -0.1, vjust = 1.3,
             size = 2.8, color = "grey60", fontface = "italic") +
    annotate("text", x =  Inf, y =  Inf,
             label = "High Rho\nHigh SASP",
             hjust = 1.1, vjust = 1.3,
             size = 2.8, color = "#B2182B", fontface = "bold.italic") +
    annotate("text", x = -Inf, y = -Inf,
             label = "Low Rho\nLow SASP",
             hjust = -0.1, vjust = -0.3,
             size = 2.8, color = "grey60", fontface = "italic") +
    annotate("text", x =  Inf, y = -Inf,
             label = "High Rho\nLow SASP",
             hjust = 1.1, vjust = -0.3,
             size = 2.8, color = "grey50", fontface = "italic") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # Стрелки
    geom_segment(
      data = arrows_rho,
      aes(x    = Rho_mean_basal,  y    = SASP_mean_basal,
          xend = Rho_mean_TNF12h, yend = SASP_mean_TNF12h,
          color = genotype),
      arrow     = arrow(length = unit(0.28, "cm"), type = "closed"),
      linewidth = 1.3, alpha = 0.85
    ) +
    # Точки basal
    geom_point(
      data = arrows_rho,
      aes(x = Rho_mean_basal, y = SASP_mean_basal, color = genotype),
      size = 5, shape = 16
    ) +
    # Точки TNF12h
    geom_point(
      data = arrows_rho,
      aes(x = Rho_mean_TNF12h, y = SASP_mean_TNF12h, color = genotype),
      size = 3.5, shape = 1, stroke = 1.5
    ) +
    # DKO — только basal
    geom_point(
      data  = dko_rho,
      aes(x = Rho_mean_basal, y = SASP_mean_basal),
      color = geno_colors["DKO"],
      size  = 5, shape = 17
    ) +
    # Подписи
    ggrepel::geom_label_repel(
      data = phase_rho,
      aes(x     = Rho_mean_basal,
          y     = SASP_mean_basal,
          label = genotype,
          color = genotype),
      size          = 3.2,
      box.padding   = 0.5,
      show.legend   = FALSE,
      seed          = 7
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = "● basal   ○ TNF12h   ▲ DKO basal only",
             hjust = 1.05, vjust = -0.5,
             size = 2.6, color = "grey40") +
    scale_color_manual(values = geno_colors, name = "Genotype") +
    coord_cartesian(clip = "off") +
    labs(
      title    = "Phase portrait: RhoActivity × SASP",
      subtitle = "● basal  ○ TNF12h  ▲ DKO (no TNF data) | MC38",
      x = "Rho Activity module score (z-score)",
      y = "SASP module score (z-score)"
    ) +
    theme_pub +
    theme(plot.margin = margin(10, 25, 10, 10))

  ggsave("output/figures/module_clocks/VIZ_F_phase_Rho_SASP.png",
         p_phase_rho, width = 7.5, height = 6.5, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_F_phase_Rho_SASP.pdf",
         p_phase_rho, width = 7.5, height = 6.5)

  write_csv(phase_rho, "output/tables/phase_Rho_SASP.csv")
  message("✓ VIZ-F Phase Rho×SASP сохранён")

} else {
  message("⚠ VIZ-F пропущен: RhoActivity или SASP не найдены в modules")
}

# =============================================================================
# 13. VIZ-G: DKO additivity test — панель по модулям
# =============================================================================
additivity <- scores_summary %>%
  filter(tnf == "basal") %>%
  select(module, genotype, mean_score) %>%
  pivot_wider(names_from = genotype, values_from = mean_score)

# Проверяем наличие всех колонок
req_geno <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")
miss     <- setdiff(req_geno, names(additivity))

if (length(miss) > 0) {
  message("⚠ VIZ-G: отсутствуют генотипы: ", paste(miss, collapse = ", "))
} else {

  additivity <- additivity %>%
    mutate(
      expected  = (PPM1D_KO - WT) + (PPM1B_KO - WT) + WT,
      observed  = DKO,
      ratio     = observed / expected,
      interp    = case_when(
        ratio > 1.2  ~ "Synergy (>120%)",
        ratio < 0.8  ~ "Compensation (<80%)",
        TRUE         ~ "Additive (80–120%)"
      )
    )

  write_csv(additivity, "output/tables/dko_additivity_test.csv")

  p_additivity <- additivity %>%
    select(module, observed, expected, wt = WT) %>%
    pivot_longer(c(observed, expected),
                 names_to = "type", values_to = "score") %>%
    ggplot(aes(x = module, y = score, fill = type)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_hline(
      data = additivity %>% select(module, wt = WT),
      aes(yintercept = wt),
      linetype = "dashed", color = geno_colors["WT"],
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = c("observed" = "#762A83", "expected" = "#C994C7"),
      labels = c("observed" = "Observed DKO",
                 "expected" = "Expected additive"),
      name = NULL
    ) +
    labs(
      title    = "DKO: observed vs. expected additive module score",
      subtitle = "Purple = observed DKO | Light = PPM1D_KO + PPM1B_KO additive model | Dashed = WT",
      x = NULL, y = "Mean z-score"
    ) +
    theme_pub +
    theme(
      axis.text.x     = element_text(angle = 40, hjust = 1),
      legend.position = "top"
    )

  ggsave("output/figures/module_clocks/VIZ_G_dko_additivity.png",
         p_additivity, width = 12, height = 6, dpi = 300, bg = "white")
  ggsave("output/figures/module_clocks/VIZ_G_dko_additivity.pdf",
         p_additivity, width = 12, height = 6)

  message("✓ VIZ-G DKO additivity сохранён")
}

# =============================================================================
# 14. Итог
# =============================================================================
message("\n══════════════════════════════════════════════════════════")
message("✓ 04_module_clocks_v3.R ЗАВЕРШЁН")
message("  Таблицы:")
message("    output/tables/module_scores_all_samples.csv")
message("    output/tables/module_scores_summary.csv")
message("    output/tables/phase_DDR_Inflammaging.csv")
message("    output/tables/phase_Rho_SASP.csv")
message("    output/tables/dko_additivity_test.csv")
message("  Рисунки:")
message("    VIZ_A — радар (basal, mean по генотипу)")
message("    VIZ_B — ridge plot (basal)")
message("    VIZ_C — slope plot (basal→TNF12h)")
message("    VIZ_D — heatmap всех образцов")
message("    VIZ_E — фазовый портрет DDR×Inflammaging")
message("    VIZ_F — фазовый портрет Rho×SASP")
message("    VIZ_G — DKO additivity test")
message("══════════════════════════════════════════════════════════")
