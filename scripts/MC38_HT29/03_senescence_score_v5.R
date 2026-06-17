# =============================================================================
# 03_senescence_score_v5.R
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   PPM1D/Wip1 (PPM1A_KO в этом эксперименте) — серин/треониновая фосфатаза,
#   негативный регулятор сенесценции: дефосфорилирует p53 (Ser15), CHK1/2,
#   p38MAPK, "выключая" DDR-сигналинг.
#
#   ВОПРОСЫ:
#   1. DOUBLE KO ADDITIVITY: аддитивный или синергетический эффект PPM1A+PPM1B?
#   2. TNF КАК "ВТОРОЙ УДАР": KO-фон + TNF = форсированная сенесценция?
#   3. SASP DELTA: амплификация SASP-ответа в KO vs WT при TNF-стимуляции.
#
# ИСПРАВЛЕНИЯ v5:
#   - PPM1A_KO вместо WipKO во всех метках и палитрах
#   - join через coldata$sample_id + coldata$group (структура из скрипта 01)
#   - функция score_signature определена ОДИН РАЗ, ДО первого вызова
#   - all_scores собирается в один проход без пересборки
#   - condition выводится из coldata$group (не пересчитывается)
#   - генсеты — мышиные символы (MC38, Mus musculus)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(patchwork)
})

# =============================================================================
# 0. Создаём папки
# =============================================================================
dir.create("output/figures/senescence", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",             recursive = TRUE, showWarnings = FALSE)
dir.create("data/genesets",             recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Загрузка pipeline из скрипта 01
# =============================================================================
pipeline <- readRDS("output/rds/pipeline_mouse.rds")

# norm_counts — агрегированная матрица 19330 × 14
vst_mat      <- as.matrix(pipeline$norm_counts)
coldata      <- pipeline$coldata        # tibble: sample_id, genotype, tnf, group, ...
results_base <- pipeline$results_base
results_tnf  <- pipeline$results_tnf

message("✓ vst_mat: ", nrow(vst_mat), " генов × ", ncol(vst_mat), " образцов")
message("✓ coldata образцы: ", paste(coldata$sample_id, collapse = ", "))
message("✓ coldata группы:  ", paste(levels(coldata$group), collapse = ", "))

# Проверяем совпадение имён
stopifnot(
  "colnames(vst_mat) не совпадают с coldata$sample_id" =
    all(sort(colnames(vst_mat)) == sort(coldata$sample_id))
)

# Упорядочиваем столбцы матрицы по coldata
vst_mat <- vst_mat[, coldata$sample_id]

# =============================================================================
# 2. Определяем порядок образцов для визуализаций
# =============================================================================
sample_order <- c(
  "BD01", "BD03", "BD04",   # WT basal
  "BD02", "BD05", "BD20",   # PPM1A_KO basal
  "BD06", "BD07", "BD08",   # PPM1B_KO basal
  "BD09", "BD57",            # Double_KO basal
  "BD12", "BD13", "BD14"    # TNF series
)
sample_order <- intersect(sample_order, colnames(vst_mat))

# =============================================================================
# 3. Z-score нормализация по генам (строкам)
# =============================================================================
vst_z <- t(scale(t(vst_mat)))
vst_z <- vst_z[, sample_order]

message("✓ vst_z: ", ncol(vst_z), " образцов в порядке: ",
        paste(colnames(vst_z), collapse = ", "))

# =============================================================================
# 4. Генсеты сенесценции (мышиные символы, Mus musculus)
# =============================================================================

# ── 4.1 Ядро сенесценции ─────────────────────────────────────────────────────
sen_core <- c(
  "Cdkn1a", "Cdkn2a", "Cdkn2b", "Cdkn1b",
  "Trp53",  "Rb1",    "Rbl1",   "Rbl2",
  "Mdm2",   "Mdm4",   "Pml",
  "Hmga1",  "Hmga2",
  "Bhlhe40",          # Dec1 — официальный символ
  "Lmnb1",  "Lmnb2"
)

# ── 4.2 SASP-цитокины ─────────────────────────────────────────────────────────
sen_sasp <- c(
  "Il6",    "Il1a",   "Il1b",
  "Il10",   "Il13",   "Il15",   "Il18",
  "Ccl2",   "Ccl3",   "Ccl5",   "Ccl7",    "Ccl8",
  "Cxcl1",  "Cxcl2",  "Cxcl5",  "Cxcl10",  "Cxcl12", "Cx3cl1",
  "Tnf",    "Tnfsf10","Tnfsf4", "Tnfsf9",
  "Mmp1",   "Mmp3",   "Mmp9",   "Mmp10",
  "Mmp12",  "Mmp13",  "Mmp14",
  "Vegfa",  "Hgf",
  "Igfbp2", "Igfbp3", "Igfbp5", "Igfbp7",
  "Serpine1",
  "Timp1",  "Timp2",
  "Il6ra",  "Il1r1",  "Tnfrsf1a"
)

# ── 4.3 DDR-маркёры ───────────────────────────────────────────────────────────
sen_ddr <- c(
  "H2ax",   "Atm",    "Atr",
  "Dna2",   "Rad17",  "Rpa1",   "Rpa2",
  "Chek1",  "Chek2",
  "Brca1",  "Brca2",  "Rad51",
  "Trp53bp1",         # 53BP1 — официальный символ у мыши
  "Mdc1",   "Rnf8",   "Rnf168",
  "Bbc3",             # PUMA
  "Pmaip1",           # NOXA
  "Bax",    "Bcl2",
  "Parp1",  "Xrcc1",  "Lig3",
  "Mlh1",   "Msh2",   "Msh6",
  "Nbn",    "Mre11a", "Rad50"
)

# ── 4.4 Митохондриальный стресс / cGAS-STING ──────────────────────────────────
sen_mitochondria <- c(
  "Cgas",   "Sting1", "Tbk1",
  "Irf3",   "Irf7",
  "Mx1",    "Mx2",
  "Ifit1",  "Ifit2",  "Ifit3",
  "Pink1",  "Prkn",
  "Bnip3",  "Bnip3l", "Fundc1",
  "Dnm1l",            # Drp1 — официальный символ
  "Fis1",
  "Mfn1",   "Mfn2",   "Opa1",
  "Ndufs1", "Ndufs2", "Ndufs3",
  "Sdha",   "Sdhb",
  "Cycs",
  "Cox4i1", "Cox5a",
  "Atp5f1a",          # Atp5a1 — официальный символ
  "Sod2",   "Cat",    "Gpx4",
  "Nfe2l2"
)

# ── 4.5 OIS (онкоген-индуцированная сенесценция) ─────────────────────────────
sen_ois <- c(
  "Hras",   "Kras",   "Nras",
  "Braf",   "Raf1",
  "E2f1",   "E2f3",
  "Ccnd1",  "Ccnd2",  "Ccnd3",
  "Cdk4",   "Cdk6",
  "Ccne1",  "Ccne2",
  "Cdk2",
  "Pcna",   "Mcm2",   "Mcm7",
  "Mki67"             # Ki67 — официальный символ
)

# Объединяем все генсеты в именованный список
senescence_genesets <- list(
  sen_core         = sen_core,
  sen_sasp         = sen_sasp,
  sen_ddr          = sen_ddr,
  sen_mitochondria = sen_mitochondria,
  sen_ois          = sen_ois
)

# Сохраняем для переиспользования
saveRDS(senescence_genesets, "data/genesets/senescence_genesets.rds")

# =============================================================================
# 5. Проверка покрытия генсетов в матрице
# =============================================================================
check_coverage <- function(geneset, name, mat) {
  found   <- intersect(geneset, rownames(mat))
  missing <- setdiff(geneset, rownames(mat))
  pct     <- round(100 * length(found) / length(geneset))
  message(sprintf("  %-20s: %d/%d генов (%.0f%%)",
                  name, length(found), length(geneset), pct))
  if (length(missing) > 0 && length(missing) <= 8)
    message("    Не найдены: ", paste(missing, collapse = ", "))
  invisible(found)
}

message("\n── Покрытие генсетов в vst_mat ──────────────────────────")
walk2(senescence_genesets, names(senescence_genesets),
      ~ check_coverage(.x, .y, vst_mat))

# Фильтруем каждый генсет до генов, присутствующих в матрице
senescence_genesets_filt <- map(senescence_genesets,
                                ~ intersect(.x, rownames(vst_mat)))

message("\n── Финальный размер генсетов (filtered) ─────────────────")
walk2(senescence_genesets_filt, names(senescence_genesets_filt),
      ~ message(sprintf("  %-20s: %d генов", .y, length(.x))))

# =============================================================================
# 6. Если в pipeline есть sig_list — добавляем туда наши генсеты
#    (для совместимости с downstream скриптами 04–06)
# =============================================================================
if (!is.null(pipeline$sig_list)) {
  sig_list <- pipeline$sig_list
} else {
  sig_list <- list()
}

# Добавляем / перезаписываем сенесцентные генсеты
sig_list <- modifyList(sig_list, senescence_genesets_filt)
message("\n✓ sig_list содержит ", length(sig_list), " генсетов: ",
        paste(names(sig_list), collapse = ", "))

# =============================================================================
# 7. Функция скоринга — средний z-score по генсету
#    (определяется ОДИН РАЗ здесь)
# =============================================================================
score_signature <- function(mat_z, geneset, sig_name) {
  genes_found <- intersect(geneset, rownames(mat_z))
  if (length(genes_found) < 3) {
    warning("Генсет '", sig_name, "': найдено < 3 генов, пропускаем")
    return(NULL)
  }
  scores <- colMeans(mat_z[genes_found, , drop = FALSE], na.rm = TRUE)
  tibble(
    sample    = names(scores),
    score     = as.numeric(scores),
    signature = sig_name,
    n_genes   = length(genes_found)
  )
}

# =============================================================================
# 8. Считаем скоры для всех генсетов в sig_list
# =============================================================================
all_scores_raw <- map_dfr(names(sig_list),
                          ~ score_signature(vst_z, sig_list[[.x]], .x))

message("✓ Скоры посчитаны: ", nrow(all_scores_raw), " строк, ",
        n_distinct(all_scores_raw$sample), " образцов, ",
        n_distinct(all_scores_raw$signature), " сигнатур")

# =============================================================================
# 9. Присоединяем метаданные из coldata
#    coldata содержит: sample_id, genotype, tnf, group, n_replicates, exploratory
# =============================================================================

# Маппинг group → короткий label для визуализаций
# Смотрим реальные уровни group в coldata
print(levels(coldata$group))
print(levels(coldata$genotype))

group_to_condition <- c(
  "WT_basal"          = "WT",
  "PPM1D_KO_basal"    = "PPM1D_KO",
  "PPM1B_KO_basal"    = "PPM1B_KO",
  "DKO_basal"         = "DKO",          # <-- было Double_KO_basal
  "WT_TNF12h"         = "WT_TNF",
  "PPM1D_KO_TNF12h"   = "PPM1D_KO_TNF",
  "PPM1B_KO_TNF12h"   = "PPM1B_KO_TNF"
)

condition_levels <- c(
  "WT", "PPM1D_KO", "PPM1B_KO", "DKO",
  "WT_TNF", "PPM1D_KO_TNF", "PPM1B_KO_TNF"
)

# Палитры под реальные имена
condition_colors <- c(
  "WT"            = "#4DAF4A",
  "PPM1D_KO"      = "#E41A1C",
  "PPM1B_KO"      = "#FF7F00",
  "DKO"           = "#984EA3",
  "WT_TNF"        = "#A6D96A",
  "PPM1D_KO_TNF"  = "#FB9A99",
  "PPM1B_KO_TNF"  = "#FDBF6F"
)

genotype_colors <- c(
  "WT"       = "#4DAF4A",
  "PPM1D_KO" = "#E41A1C",
  "PPM1B_KO" = "#FF7F00",
  "DKO"      = "#984EA3"
)


all_scores <- all_scores_raw %>%
  left_join(
    coldata %>%
      select(any_of(c("sample", "sample_id",
                      "genotype", "tnf",
                      "group", "n_replicates", "exploratory"))),
    by = setNames(
      intersect(c("sample", "sample_id"), names(coldata)),
      "sample"
    )
  ) %>%
  mutate(
    condition = factor(
      recode(as.character(group), !!!group_to_condition),
      levels = condition_levels
    )
  )

# Проверка
message("NA в condition: ", sum(is.na(all_scores$condition)))
print(table(all_scores$condition, useNA = "always"))

# Проверка
n_na_cond    <- sum(is.na(all_scores$condition))
n_na_geno    <- sum(is.na(all_scores$genotype))
n_samples    <- n_distinct(all_scores$sample)

message("\n── Проверка join ────────────────────────────────────────")
message("  NA в condition: ", n_na_cond)
message("  NA в genotype:  ", n_na_geno)
message("  Уникальных образцов: ", n_samples)
print(table(all_scores$condition, useNA = "always"))

stopifnot(
  "NA в condition после join" = n_na_cond == 0,
  "NA в genotype после join"  = n_na_geno == 0,
  "Не все 14 образцов"        = n_samples == 14
)

# Сохраняем
write_csv(all_scores, "output/tables/senescence_scores_all.csv")
message("✓ Таблица сохранена: output/tables/senescence_scores_all.csv")

# =============================================================================
# 10. Цветовые палитры
# =============================================================================
condition_colors <- c(
  "WT"            = "#4DAF4A",
  "PPM1A_KO"      = "#E41A1C",
  "PPM1B_KO"      = "#FF7F00",
  "Double_KO"     = "#984EA3",
  "WT_TNF"        = "#A6D96A",
  "PPM1A_KO_TNF"  = "#FB9A99",
  "PPM1B_KO_TNF"  = "#FDBF6F"
)

genotype_colors <- c(
  "WT"         = "#4DAF4A",
  "PPM1A_KO"   = "#E41A1C",
  "PPM1B_KO"   = "#FF7F00",
  "Double_KO"  = "#984EA3"
)

# =============================================================================
# 11. Сводная матрица для heatmap (средние по condition-группам)
# =============================================================================
score_matrix <- all_scores %>%
  group_by(signature, condition) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_score) %>%
  column_to_rownames("signature") %>%
  as.matrix()

# Биологически осмысленный порядок строк
sig_order_preferred <- c(
  "sen_core", "sen_ois",
  "Senescence_DDR", "Senescence_SASP",
  "sen_sasp", "Inflammaging",
  "sen_ddr",  "cGAS_STING",
  "sen_mitochondria", "Mito_Stress",
  "Rho_Activity",
  "ProxTubule",
  "Tyshkovskiy_DNAmAge",
  "Tyshkovskiy_LongevityUp",
  "Tyshkovskiy_LongevityDown"
)
sig_order <- intersect(sig_order_preferred, rownames(score_matrix))
# Добавляем оставшиеся генсеты (не вошедшие в preferred)
sig_order <- c(sig_order,
               setdiff(rownames(score_matrix), sig_order))

col_order <- intersect(condition_levels, colnames(score_matrix))
score_matrix <- score_matrix[sig_order, col_order]

message("\n── score_matrix: ", nrow(score_matrix), " × ", ncol(score_matrix))

# =============================================================================
# 12. Fig_SEN_01: Heatmap средних скоров по группам
# =============================================================================
col_fun <- colorRamp2(
  c(-1.5, -0.75, 0, 0.75, 1.5),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

# Реальные имена условий в матрице (после всех исправлений)
cond_in_matrix <- colnames(score_matrix)

# Именованный вектор цветов — имена ДОЛЖНЫ совпадать со значениями аннотации
col_vec <- condition_colors[cond_in_matrix]
names(col_vec) <- cond_in_matrix   # <-- ЭТО БЫЛО ПОТЕРЯНО

# Аннотация столбцов — передаём вектор значений, col = список с именованным вектором
col_ann <- HeatmapAnnotation(
  Condition = cond_in_matrix,                      # вектор значений
  col       = list(Condition = col_vec),           # именованный вектор цветов
  show_legend          = TRUE,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9)
)

# Аннотация строк
module_map <- c(
  "sen_core"         = "Senescence Core",
  "sen_ois"          = "Senescence Core",
  "Senescence_DDR"   = "Senescence Core",
  "Senescence_SASP"  = "SASP / Inflammation",
  "sen_sasp"         = "SASP / Inflammation",
  "Inflammaging"     = "SASP / Inflammation",
  "sen_ddr"          = "DDR",
  "cGAS_STING"       = "DDR",
  "sen_mitochondria" = "Mito / Stress",
  "Mito_Stress"      = "Mito / Stress",
  "Rho_Activity"     = "Cytoskeleton",
  "ProxTubule"       = "Kidney",
  "Tyshkovskiy_DNAmAge"       = "Aging Clocks",
  "Tyshkovskiy_LongevityUp"   = "Aging Clocks",
  "Tyshkovskiy_LongevityDown" = "Aging Clocks",
  "Aging_DNAmAge"    = "Aging Clocks",
  "Aging_LongevityUp"   = "Aging Clocks",
  "Aging_LongevityDown" = "Aging Clocks"
)

module_colors <- c(
  "Senescence Core"     = "#D9534F",
  "SASP / Inflammation" = "#F0AD4E",
  "DDR"                 = "#5B9BD5",
  "Mito / Stress"       = "#70AD47",
  "Cytoskeleton"        = "#9B59B6",
  "Kidney"              = "#1ABC9C",
  "Aging Clocks"        = "#95A5A6",
  "Other"               = "#CCCCCC"
)

row_modules <- module_map[rownames(score_matrix)]
row_modules[is.na(row_modules)] <- "Other"

# Именованный вектор для rowAnnotation — имена = значения аннотации
used_modules      <- unique(row_modules)
row_module_colors <- module_colors[used_modules]
names(row_module_colors) <- used_modules   # <-- то же самое правило

row_ann <- rowAnnotation(
  Module = row_modules,
  col    = list(Module = row_module_colors),
  show_legend          = TRUE,
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 9)
)

png("output/figures/senescence/Fig_SEN_01_heatmap_group_means.png",
    width = 2600, height = 2000, res = 200)

ht <- Heatmap(
  score_matrix,
  name             = "Mean\nZ-score",
  col              = col_fun,
  top_annotation   = col_ann,
  right_annotation = row_ann,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_names_side   = "left",
  row_names_gp     = gpar(fontsize = 10),
  column_names_gp  = gpar(fontsize = 10),
  column_names_rot = 35,
  rect_gp          = gpar(col = "white", lwd = 0.8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", score_matrix[i, j]),
      x, y,
      gp = gpar(
        fontsize = 8,
        col = ifelse(abs(score_matrix[i, j]) > 0.8, "white", "black")
      )
    )
  },
  heatmap_legend_param = list(
    direction = "vertical",
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht, padding = unit(c(5, 15, 5, 5), "mm"))
dev.off()
message("✓ Fig_SEN_01 сохранён")

# =============================================================================
# 13. Fig_SEN_02: Barplot — скоры по отдельным образцам
# =============================================================================
key_sigs <- c(
  "sen_core", "sen_sasp", "sen_ddr", "sen_mitochondria", "sen_ois",
  # Добавляем из sig_list если есть
  intersect(c("Senescence_SASP", "Senescence_DDR",
              "Inflammaging", "cGAS_STING",
              "Rho_Activity", "Mito_Stress"),
            names(sig_list))
)
key_sigs <- unique(key_sigs)

p_bar <- all_scores %>%
  filter(signature %in% key_sigs) %>%
  mutate(
    signature = factor(signature, levels = key_sigs),
    sample    = factor(sample, levels = sample_order)
  ) %>%
  ggplot(aes(x = sample, y = score, fill = condition)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed",
             color = "grey40") +
  scale_fill_manual(values = condition_colors, name = "Condition",
                    drop = FALSE) +
  facet_wrap(~ signature, scales = "free_y",
             ncol = min(4, ceiling(sqrt(length(key_sigs))))) +
  labs(
    title    = "Senescence & related signature scores per sample",
    subtitle = "Mean z-score of normalized expression | MC38 cells",
    x = NULL, y = "Mean Z-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    strip.text      = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

ggsave("output/figures/senescence/Fig_SEN_02_barplot_per_sample.png",
       p_bar, width = 18, height = 10, dpi = 200)
message("✓ Fig_SEN_02 сохранён")

# =============================================================================
# 14. Fig_SEN_03: TNF delta-score
#     delta = TNF_score − basal_score (per genotype, per signature)
# =============================================================================

# Образцы с TNF: BD12 (WT), BD13 (PPM1A_KO), BD14 (PPM1B_KO)
# Baseline mean для каждого генотипа:
#   WT baseline       = mean(BD01, BD03, BD04)
#   PPM1A_KO baseline = mean(BD02, BD05, BD20)
#   PPM1B_KO baseline = mean(BD06, BD07, BD08)

# Baseline means по генотипам
baseline_means <- all_scores %>%
  filter(tnf == "basal",
         genotype %in% c("WT", "PPM1A_KO", "PPM1B_KO"),
         signature %in% key_sigs) %>%
  group_by(genotype, signature) %>%
  summarise(baseline_mean = mean(score, na.rm = TRUE), .groups = "drop")

# TNF-образцы
tnf_scores <- all_scores %>%
  filter(tnf == "TNF12h",
         signature %in% key_sigs) %>%
  left_join(baseline_means, by = c("genotype", "signature")) %>%
  mutate(
    delta    = score - baseline_mean,
    genotype = factor(genotype,
                      levels = c("WT", "PPM1A_KO", "PPM1B_KO"))
  )

p_delta <- tnf_scores %>%
  ggplot(aes(x = genotype, y = delta, fill = genotype)) +
  geom_col(width = 0.65, color = "grey20", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label = sprintf("%.2f", delta),
                y = delta + sign(delta) * 0.03),
            size = 3, vjust = 0) +
  scale_fill_manual(
    values = genotype_colors[c("WT", "PPM1A_KO", "PPM1B_KO")],
    name   = "Genotype"
  ) +
  facet_wrap(~ signature, scales = "free_y",
             ncol = min(4, ceiling(sqrt(length(key_sigs))))) +
  labs(
    title    = "TNF-induced score change (Δ = TNF − basal mean)",
    subtitle = "Positive = TNF activates signature | 'SASP storm' hypothesis",
    x = NULL, y = "Δ Mean Z-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    strip.text      = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

ggsave("output/figures/senescence/Fig_SEN_03_TNF_delta_score.png",
       p_delta, width = 16, height = 8, dpi = 200)
message("✓ Fig_SEN_03 сохранён")

# =============================================================================
# 15. Fig_SEN_04: Double KO additivity test
#     Ожидаемое аддитивное = (PPM1A_KO − WT) + (PPM1B_KO − WT) + WT
# =============================================================================
baseline_geno <- all_scores %>%
  filter(tnf == "basal",
         signature %in% names(senescence_genesets_filt)) %>%
  group_by(genotype, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

# Сводим в wide-формат для расчёта
additivity_wide <- baseline_geno %>%
  pivot_wider(names_from = genotype, values_from = mean_score)

# Проверяем наличие нужных колонок
required_geno <- c("WT", "PPM1A_KO", "PPM1B_KO", "Double_KO")
missing_geno  <- setdiff(required_geno, colnames(additivity_wide))
if (length(missing_geno) > 0) {
  warning("Отсутствуют генотипы для additivity test: ",
          paste(missing_geno, collapse = ", "))
} else {

  additivity_test <- additivity_wide %>%
    mutate(
      expected_additive = (PPM1A_KO - WT) + (PPM1B_KO - WT) + WT,
      observed          = Double_KO,
      additivity_ratio  = observed / expected_additive,
      interpretation    = case_when(
        additivity_ratio > 1.2  ~ "Synergy (>120%)",
        additivity_ratio < 0.8  ~ "Compensation (<80%)",
        TRUE                    ~ "Additive (80–120%)"
      )
    )

  print(additivity_test %>%
          select(signature, WT, PPM1A_KO, PPM1B_KO,
                 Double_KO, expected_additive,
                 additivity_ratio, interpretation))

  write_csv(additivity_test,
            "output/tables/double_ko_additivity_test.csv")
  message("✓ Additivity test сохранён: output/tables/double_ko_additivity_test.csv")

  # Визуализация
  p_additivity <- additivity_test %>%
    select(signature, observed, expected = expected_additive, wt = WT) %>%
    pivot_longer(c(observed, expected),
                 names_to = "type", values_to = "score") %>%
    ggplot(aes(x = signature, y = score, fill = type)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_hline(data = . %>% distinct(signature, wt),
               aes(yintercept = wt),
               linetype = "dashed", color = "#4DAF4A",
               linewidth = 0.5) +
    scale_fill_manual(
      values = c("observed" = "#984EA3", "expected" = "#D4A0D8"),
      labels = c("observed" = "Observed Double KO",
                 "expected" = "Expected additive"),
      name   = NULL
    ) +
    labs(
      title    = "Double KO: observed vs expected additive senescence score",
      subtitle = "Purple = observed | Light = additive model | Dashed = WT mean",
      x = NULL, y = "Mean Z-score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 40, hjust = 1),
      strip.text      = element_text(face = "bold"),
      legend.position = "top",
      plot.title      = element_text(face = "bold")
    )

  ggsave("output/figures/senescence/Fig_SEN_04_dKO_additivity.png",
         p_additivity, width = 12, height = 6, dpi = 200)
  message("✓ Fig_SEN_04 сохранён")
}

# =============================================================================
# 16. Fig_SEN_05: Сводная панель (Panel A + B + C)
# =============================================================================

# Panel A — boxplot скоров по генотипу + TNF
p_panel_a <- all_scores %>%
  filter(signature %in% names(senescence_genesets_filt)) %>%
  mutate(
    genotype  = factor(genotype, levels = c("WT","PPM1A_KO","PPM1B_KO","Double_KO")),
    tnf_label = factor(ifelse(tnf == "TNF12h", "TNF 12h", "Baseline"),
                       levels = c("Baseline", "TNF 12h"))
  ) %>%
  ggplot(aes(x = genotype, y = score,
             fill = genotype, shape = tnf_label)) +
  geom_boxplot(
    aes(group = interaction(genotype, tnf_label)),
    width = 0.45, outlier.shape = NA, alpha = 0.45,
    position = position_dodge(0.75), color = "grey30"
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 3, alpha = 0.9, aes(color = genotype)
  ) +
  facet_wrap(~ signature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values  = genotype_colors, name = "Genotype") +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  scale_shape_manual(
    values = c("Baseline" = 21, "TNF 12h" = 24),
    name   = "Treatment"
  ) +
  labs(
    title    = "A. Senescence module scores: genotype × TNF",
    subtitle = "▲ = +TNF 12h | MC38 cells",
    x = NULL, y = "Mean Z-score (module)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    strip.text      = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

# Panel B — TNF delta для ключевых сигнатур
p_panel_b <- tnf_scores %>%
  ggplot(aes(x = signature, y = delta, fill = genotype)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_fill_manual(
    values = genotype_colors[c("WT", "PPM1A_KO", "PPM1B_KO")],
    name   = "Genotype"
  ) +
  labs(
    title    = "B. TNF-induced SASP amplification (Δ = TNF − baseline)",
    subtitle = "Larger bar = stronger TNF response",
    x = NULL, y = "Δ Z-score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

# Сборка панели
p_full <- (p_panel_a / p_panel_b) +
  plot_layout(heights = c(2.5, 1)) +
  plot_annotation(
    title    = "Senescence landscape: PPM1A KO, PPM1B KO, Double KO ± TNF",
    subtitle = "MC38 murine colon carcinoma | PPM1A = Wip1/PPM1D",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("output/figures/senescence/Fig_SEN_05_full_panel.png",
       p_full, width = 16, height = 16, dpi = 200)
ggsave("output/figures/senescence/Fig_SEN_05_full_panel.pdf",
       p_full, width = 16, height = 16)
message("✓ Fig_SEN_05 (сводная панель) сохранён")

# =============================================================================
# 17. Итог
# =============================================================================
message("\n══════════════════════════════════════════════════════════")
message("✓ 03_senescence_score_v5.R ЗАВЕРШЁН")
message("  Таблицы:")
message("    output/tables/senescence_scores_all.csv")
message("    output/tables/double_ko_additivity_test.csv")
message("  Рисунки:")
message("    Fig_SEN_01 — heatmap группы (средние z-score)")
message("    Fig_SEN_02 — barplot по отдельным образцам")
message("    Fig_SEN_03 — TNF delta-score по генотипам")
message("    Fig_SEN_04 — Double KO additivity test")
message("    Fig_SEN_05 — сводная панель A+B")
message("══════════════════════════════════════════════════════════")
