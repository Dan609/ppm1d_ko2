# =============================================================================
# 01_deseq2_pipeline.R
# PROJECT: ppm1d_ko2 — Unified phosphatase KO transcriptomics pipeline
#
# DATASETS:
#   A) MC-38  (Mus musculus)   — ПЕРВЫЙ приоритет
#      BD01,BD03,BD04  = WT         (n=3, basal)
#      BD02,BD05,BD20  = PPM1D_KO   (n=3, basal)   ← ИСПРАВЛЕНО: было PPM1A_KO
#      BD06,BD07,BD08  = PPM1B_KO   (n=3, basal)
#      BD09,BD57       = DKO        (n=2, basal)
#      BD12            = WT + TNF12h        (n=1, exploratory)
#      BD13            = PPM1D_KO + TNF12h  (n=1, exploratory)
#      BD14            = PPM1B_KO + TNF12h  (n=1, exploratory)
#
#   B) HT29  (Homo sapiens)    — ВТОРОЙ приоритет
#      BD21,BD58,BD10  = WT         (n=3)
#      BD11,BD22,BD23  = PPM1D_KO   (n=3)
#      BD24,BD25,BD26  = PPM1B_KO   (n=3)
#
# ЛОГИКА: DESeq2 уже был запущен, результаты лежат в xlsx/tsv.
#   1. Загружаем нормализованные счёты
#   2. Загружаем DE-результаты, стандартизуем колонки
#   3. PPM1D_KO везде (не PPM1A!)
#   4. Сохраняем два RDS-объекта: pipeline_mouse и pipeline_human
#      + общий объект pipeline_all для кросс-видового анализа
#
# ВАЖНО: WIP1 = PPM1D (официальный символ HGNC:9277).
#         PPM1A — ДРУГАЯ фосфатаза. Ошибка исправлена во всём пайплайне.
# =============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ── ГЛОБАЛЬНЫЕ ПАРАМЕТРЫ ─────────────────────────────────────────────────────
FDR        <- 0.05
LFC_THRESH <- 1          # |log2FC| ≥ 1

# NEW: раздельные директории для двух видов
DATA_DIR_MOUSE  <- "data/mouse"    # MC-38: xlsx/tsv DE-результаты + norm_counts
DATA_DIR_HUMAN  <- "data/human"    # HT29:  xlsx/tsv DE-результаты + norm_counts

DATA_DIR_MOUSE  <- "data/new"    # MC-38: xlsx/tsv DE-результаты + norm_counts
DATA_DIR_HUMAN  <- "data/new"    # HT29:  xlsx/tsv DE-результаты + norm_counts

# NEW: цветовая палитра генотипов — единая для обоих датасетов
GENO_COLORS <- c(
  "WT"        = "#4393C3",   # синий
  "PPM1D_KO"  = "#D6604D",   # красный   ← ИСПРАВЛЕНО с PPM1A_KO
  "PPM1B_KO"  = "#74C476",   # зелёный
  "DKO"       = "#9970AB"    # фиолетовый (только MC-38)
)

# ── СОЗДАНИЕ ДИРЕКТОРИЙ ──────────────────────────────────────────────────────
for (d in c("output/rds", "output/tables",
            "output/figures/mouse", "output/figures/human",
            "data/processed")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
message("✓ Директории созданы")

# =============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# =============================================================================

# ── Загрузка одного DE-файла (xlsx или tsv) ──────────────────────────────────
read_de_file <- function(filename,
                         contrast_name,
                         is_exploratory = FALSE,
                         data_dir) {

  path <- file.path(data_dir, filename)
  if (!file.exists(path)) {
    warning("Файл не найден: ", path)
    return(NULL)
  }

  ext <- tools::file_ext(filename)

  df <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path)
  } else {
    read.csv(path,
             sep         = "\t",
             row.names   = 1,
             check.names = FALSE) %>%
      tibble::rownames_to_column("gene_id")
  }

  df <- df %>% as_tibble(.name_repair = "unique")

  # Стандартизуем имя колонки с ID генов
  if (!"gene_id" %in% colnames(df)) {
    df <- df %>% dplyr::rename(gene_id = 1)
  }

  # Стандартизуем имена колонок DESeq2
  df <- df %>%
    dplyr::rename_with(~ dplyr::recode(.,
                                       "log2FC"  = "log2FoldChange",
                                       "Log2FC"  = "log2FoldChange",
                                       "p_value" = "pvalue",
                                       "Pvalue"  = "pvalue",
                                       "FDR"     = "padj"
    ))

  # Проверка обязательных колонок
  required <- c("gene_id", "log2FoldChange", "pvalue", "padj")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0)
    warning(contrast_name, " — отсутствуют колонки: ",
            paste(missing, collapse = ", "))

  df %>%
    dplyr::mutate(
      contrast    = contrast_name,
      exploratory = is_exploratory,
      sig         = !is.na(padj) & padj < FDR & abs(log2FoldChange) >= LFC_THRESH,
      direction   = dplyr::case_when(
        sig & log2FoldChange > 0 ~ "up",
        sig & log2FoldChange < 0 ~ "down",
        TRUE                      ~ "ns"
      )
    )
}

# ── Агрегация технических реплик (lanes) → биологические образцы ─────────────
# BD01_L001_001 → BD01
aggregate_lanes <- function(norm_counts_raw) {
  extract_sample_id <- function(x) sub("_L\\d+.*$", "", x)

  lane_ids   <- colnames(norm_counts_raw)
  sample_ids <- sapply(lane_ids, extract_sample_id)

  unique_samples <- unique(sample_ids)
  message("  Lanes: ", ncol(norm_counts_raw),
          " → образцов: ", length(unique_samples))

  agg <- sapply(unique_samples, function(sid) {
    cols <- which(sample_ids == sid)
    if (length(cols) == 1) norm_counts_raw[, cols]
    else                   rowMeans(norm_counts_raw[, cols])
  })

  agg <- as.data.frame(agg)
  rownames(agg) <- rownames(norm_counts_raw)
  agg
}

# ── Volcano plot (универсальная функция) ──────────────────────────────────────
make_volcano <- function(de_df,
                         title_str,
                         top_n      = 15,
                         fc_thresh  = LFC_THRESH,
                         fdr_thresh = FDR) {

  de_plot <- de_df %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj + 1e-300),
      color_group = dplyr::case_when(
        sig & log2FoldChange > 0 ~ "up",
        sig & log2FoldChange < 0 ~ "down",
        TRUE                      ~ "ns"
      )
    )

  # Подписываем топ-N по значимости
  top_genes <- de_plot %>%
    dplyr::filter(sig) %>%
    dplyr::slice_max(neg_log10_padj, n = top_n) %>%
    dplyr::pull(gene_id)

  de_plot <- de_plot %>%
    dplyr::mutate(label = ifelse(gene_id %in% top_genes, gene_id, ""))

  n_up   <- sum(de_plot$color_group == "up")
  n_down <- sum(de_plot$color_group == "down")

  ggplot(de_plot, aes(x = log2FoldChange, y = neg_log10_padj,
                      color = color_group)) +
    geom_point(data = . %>% dplyr::filter(color_group == "ns"),
               size = 0.5, alpha = 0.25, color = "grey75") +
    geom_point(data = . %>% dplyr::filter(color_group != "ns"),
               size = 1.6, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 2.8, max.overlaps = 20,
                             segment.color = "grey50", segment.size = 0.3,
                             show.legend = FALSE) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_thresh),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(
      values = c("up" = "#D6604D", "down" = "#4393C3", "ns" = "grey75"),
      labels = c("up"   = paste0("Up (", n_up, ")"),
                 "down" = paste0("Down (", n_down, ")"),
                 "ns"   = "NS"),
      name = NULL
    ) +
    labs(title    = title_str,
         subtitle = paste0("FDR<", fdr_thresh, ", |log2FC|≥", fc_thresh,
                           "  |  ↑ ", n_up, "  ↓ ", n_down),
         x = "log2 Fold Change",
         y = "-log10(FDR)") +
    theme_bw(base_size = 11) +
    theme(legend.position  = "top",
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, color = "grey40"))
}

# ── DEG summary barplot (универсальная функция) ───────────────────────────────
make_deg_barplot <- function(results_list, title_str) {

  deg_summary <- bind_rows(results_list) %>%
    dplyr::filter(direction != "ns") %>%
    dplyr::group_by(contrast, exploratory, direction) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::mutate(
      n_plot   = ifelse(direction == "down", -n, n),
      contrast = factor(contrast, levels = rev(unique(contrast)))
    )

  ggplot(deg_summary, aes(x = n_plot, y = contrast, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
    geom_text(aes(label = abs(n),
                  x = n_plot + ifelse(direction == "up", 20, -20)),
              size = 3,
              hjust = ifelse(deg_summary$direction == "up", 0, 1)) +
    # Рамка для exploratory
    {
      expl <- deg_summary %>%
        dplyr::filter(exploratory) %>%
        dplyr::distinct(contrast)
      if (nrow(expl) > 0)
        geom_tile(data = expl,
                  aes(x = 0, y = contrast, width = Inf, height = 1),
                  fill = NA, color = "orange", linewidth = 0.8,
                  inherit.aes = FALSE)
    } +
    scale_fill_manual(values = c("up" = "#D6604D", "down" = "#4393C3"),
                      name = "Direction") +
    scale_x_continuous(labels = abs,
                       name   = "Number of DEGs (|log2FC|≥1, FDR<0.05)") +
    labs(title    = title_str,
         subtitle = "Orange frame = exploratory (n=1)",
         y        = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position    = "top")
}

# =============================================================================
# БЛОК A: MC-38 (МЫШb)
# =============================================================================
message("\n══ БЛОК A: MC-38 (Mus musculus) ══")

# ── A1. Нормализованные счёты ─────────────────────────────────────────────────
norm_counts_raw_mouse <- read.csv(
  file.path(DATA_DIR_MOUSE, "normalized_counts_first.tsv"),
  sep = "\t", row.names = 1, check.names = FALSE
)
message("MC-38 сырая матрица: ",
        nrow(norm_counts_raw_mouse), " × ", ncol(norm_counts_raw_mouse))

norm_counts_mouse <- aggregate_lanes(norm_counts_raw_mouse)
message("MC-38 агрегировано: ",
        nrow(norm_counts_mouse), " × ", ncol(norm_counts_mouse))

write.csv(norm_counts_mouse, "data/processed/norm_counts_mouse.csv")

# ── A2. Coldata MC-38 ─────────────────────────────────────────────────────────
# NEW: PPM1D_KO (не PPM1A_KO!) + добавлена колонка organism
coldata_mouse <- tibble(
  sample_id = c("BD01","BD03","BD04",
                "BD02","BD05","BD20",
                "BD06","BD07","BD08",
                "BD09","BD57",
                "BD12","BD13","BD14"),

  organism  = "mouse",             # NEW

  cell_line = "MC38",              # NEW

  # NEW: WIP1 = PPM1D — исправлено везде
  genotype  = factor(
    c("WT","WT","WT",
      "PPM1D_KO","PPM1D_KO","PPM1D_KO",
      "PPM1B_KO","PPM1B_KO","PPM1B_KO",
      "DKO","DKO",
      "WT","PPM1D_KO","PPM1B_KO"),
    levels = c("WT","PPM1D_KO","PPM1B_KO","DKO")
  ),

  tnf = factor(
    c(rep("basal", 11), "TNF12h","TNF12h","TNF12h"),
    levels = c("basal","TNF12h")
  ),

  group = factor(paste(
    c("WT","WT","WT",
      "PPM1D_KO","PPM1D_KO","PPM1D_KO",
      "PPM1B_KO","PPM1B_KO","PPM1B_KO",
      "DKO","DKO",
      "WT","PPM1D_KO","PPM1B_KO"),
    c(rep("basal",11), rep("TNF12h",3)),
    sep = "_"
  )),

  n_replicates = c(3,3,3, 3,3,3, 3,3,3, 2,2, 1,1,1),
  exploratory  = n_replicates < 2
)

message("✓ coldata_mouse: ", nrow(coldata_mouse), " образцов")
print(coldata_mouse %>% dplyr::select(sample_id, genotype, tnf,
                                      n_replicates, exploratory))

# Синхронизация порядка столбцов матрицы
common_samples_mouse <- intersect(coldata_mouse$sample_id,
                                  colnames(norm_counts_mouse))
norm_counts_mouse    <- norm_counts_mouse[, common_samples_mouse]
coldata_mouse        <- coldata_mouse %>%
  dplyr::filter(sample_id %in% common_samples_mouse) %>%
  dplyr::arrange(match(sample_id, common_samples_mouse))

# ── A3. DE-результаты MC-38 ───────────────────────────────────────────────────
# BASELINE (confirmatory, n≥2)
# NEW: имена контрастов — PPM1D_KO (не PPM1A_KO)

results_mouse_base <- list()

results_mouse_base[["PPM1D_KO_vs_WT"]] <- read_de_file(
  "wt_vs_WipKo_results.xlsx", "PPM1D_KO_vs_WT",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

results_mouse_base[["PPM1B_KO_vs_WT"]] <- read_de_file(
  "wt_vs_ppm1b_KO_results.xlsx", "PPM1B_KO_vs_WT",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

results_mouse_base[["DKO_vs_WT"]] <- read_de_file(
  "wt_vs_double_KO_results.xlsx", "DKO_vs_WT",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

# NEW: прямые сравнения фосфатаз между собой
# (если файлы уже есть — подключаем; если нет — NULL, будут добавлены позже)
results_mouse_base[["PPM1D_KO_vs_PPM1B_KO"]] <- read_de_file(
  "PPM1D_KO_vs_PPM1B_KO_results.xlsx", "PPM1D_KO_vs_PPM1B_KO",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

results_mouse_base[["DKO_vs_PPM1D_KO"]] <- read_de_file(
  "DKO_vs_PPM1D_KO_results.xlsx", "DKO_vs_PPM1D_KO",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

results_mouse_base[["DKO_vs_PPM1B_KO"]] <- read_de_file(
  "DKO_vs_PPM1B_KO_results.xlsx", "DKO_vs_PPM1B_KO",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

# TNF basal reference
results_mouse_base[["WT_TNF_vs_WT_basal"]] <- read_de_file(
  "wt_vs_wt_TNF_12h.tsv", "WT_TNF_vs_WT_basal",
  is_exploratory = FALSE, data_dir = DATA_DIR_MOUSE)

# EXPLORATORY (n=1 TNF-образцы)
results_mouse_tnf <- list()

results_mouse_tnf[["PPM1D_KO_TNF_vs_WT_basal"]] <- read_de_file(
  "wt_vs_WipKo__TNF_12h_results.xlsx", "PPM1D_KO_TNF_vs_WT_basal",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1B_KO_TNF_vs_WT_basal"]] <- read_de_file(
  "wt_vs_ppm1b_KO__TNF_12h_results.xlsx", "PPM1B_KO_TNF_vs_WT_basal",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1D_KO_TNF_vs_PPM1D_KO_basal"]] <- read_de_file(
  "WipKo_vs_WipKo__TNF_12h_results.xlsx", "PPM1D_KO_TNF_vs_PPM1D_KO_basal",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1D_KO_TNF_vs_WT_TNF"]] <- read_de_file(
  "wt_TNF_12h_vs_WipKo__TNF_12h_results.xlsx", "PPM1D_KO_TNF_vs_WT_TNF",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1B_KO_TNF_vs_WT_TNF"]] <- read_de_file(
  "wt_TNF_12h_vs_ppm1b_KO__TNF_12h_results.xlsx", "PPM1B_KO_TNF_vs_WT_TNF",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1B_KO_TNF_vs_PPM1B_KO_basal"]] <- read_de_file(
  "ppm1b_KO_vs_ppm1b_KO__TNF_12h_results.xlsx", "PPM1B_KO_TNF_vs_PPM1B_KO_basal",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

results_mouse_tnf[["PPM1D_KO_TNF_vs_PPM1B_KO_TNF"]] <- read_de_file(
  "WipKo__TNF_12h_vs_ppm1b_KO__TNF_12h_results.xlsx", "PPM1D_KO_TNF_vs_PPM1B_KO_TNF",
  is_exploratory = TRUE, data_dir = DATA_DIR_MOUSE)

# Убираем NULL (файлы не найдены)
results_mouse_base <- Filter(Negate(is.null), results_mouse_base)
results_mouse_tnf  <- Filter(Negate(is.null), results_mouse_tnf)

message("✓ MC-38 baseline-контрастов: ", length(results_mouse_base))
message("✓ MC-38 exploratory-контрастов: ", length(results_mouse_tnf))

# =============================================================================
# БЛОК B: HT29 (ЧЕЛОВЕК)   — NEW блок
# =============================================================================
message("\n══ БЛОК B: HT29 (Homo sapiens) ══")

# ── B1. Нормализованные счёты ─────────────────────────────────────────────────
norm_counts_raw_human <- read.csv(
  file.path(DATA_DIR_HUMAN, "normalized_counts_first.tsv"),
  sep = "\t", row.names = 1, check.names = FALSE
)
message("HT29 сырая матрица: ",
        nrow(norm_counts_raw_human), " × ", ncol(norm_counts_raw_human))

norm_counts_human <- aggregate_lanes(norm_counts_raw_human)
message("HT29 агрегировано: ",
        nrow(norm_counts_human), " × ", ncol(norm_counts_human))

write.csv(norm_counts_human, "data/processed/norm_counts_human.csv")

# ── B2. Coldata HT29 ──────────────────────────────────────────────────────────
coldata_human <- tibble(
  sample_id = c("BD21","BD58","BD10",
                "BD11","BD22","BD23",
                "BD24","BD25","BD26"),

  organism  = "human",
  cell_line = "HT29",

  genotype  = factor(
    c("WT","WT","WT",
      "PPM1D_KO","PPM1D_KO","PPM1D_KO",
      "PPM1B_KO","PPM1B_KO","PPM1B_KO"),
    levels = c("WT","PPM1D_KO","PPM1B_KO")
  ),

  # HT29 — без TNF
  tnf          = factor("basal", levels = c("basal","TNF12h")),
  group        = genotype,
  n_replicates = 3L,
  exploratory  = FALSE
)

message("✓ coldata_human: ", nrow(coldata_human), " образцов")
print(coldata_human %>% dplyr::select(sample_id, genotype, n_replicates))

# Синхронизация
common_samples_human <- intersect(coldata_human$sample_id,
                                  colnames(norm_counts_human))
norm_counts_human    <- norm_counts_human[, common_samples_human]
coldata_human        <- coldata_human %>%
  dplyr::filter(sample_id %in% common_samples_human) %>%
  dplyr::arrange(match(sample_id, common_samples_human))

# ── B3. DE-результаты HT29 ───────────────────────────────────────────────────
results_human_base <- list()

results_human_base[["PPM1D_KO_vs_WT"]] <- read_de_file(
  "HT29_PPM1D_KO_vs_WT_results.xlsx", "PPM1D_KO_vs_WT",
  is_exploratory = FALSE, data_dir = DATA_DIR_HUMAN)

results_human_base[["PPM1B_KO_vs_WT"]] <- read_de_file(
  "HT29_PPM1B_KO_vs_WT_results.xlsx", "PPM1B_KO_vs_WT",
  is_exploratory = FALSE, data_dir = DATA_DIR_HUMAN)

results_human_base[["PPM1D_KO_vs_PPM1B_KO"]] <- read_de_file(
  "HT29_PPM1D_KO_vs_PPM1B_KO_results.xlsx", "PPM1D_KO_vs_PPM1B_KO",
  is_exploratory = FALSE, data_dir = DATA_DIR_HUMAN)

results_human_base <- Filter(Negate(is.null), results_human_base)

message("✓ HT29 baseline-контрастов: ", length(results_human_base))

# =============================================================================
# СВОДНЫЕ ТАБЛИЦЫ DEG
# =============================================================================

summarise_de <- function(results_list, dataset_label) {
  bind_rows(results_list) %>%
    dplyr::group_by(contrast, exploratory) %>%
    dplyr::summarise(
      total_genes = n(),
      sig_up      = sum(direction == "up",   na.rm = TRUE),
      sig_down    = sum(direction == "down",  na.rm = TRUE),
      sig_total   = sum(sig,                  na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(dataset = dataset_label) %>%
    dplyr::arrange(exploratory, contrast)
}

summary_mouse <- summarise_de(c(results_mouse_base, results_mouse_tnf), "MC38_mouse")

# ── Диагностика: что нашлось в results_human_base ────────────────────────────
message("Найдено HT29 контрастов: ", length(results_human_base))
message("Имена: ", paste(names(results_human_base), collapse = ", "))

# Проверяем файлы в папке human
message("Файлы в data/human/:")
print(list.files(DATA_DIR_HUMAN))


summary_human <- summarise_de(results_human_base,                        "HT29_human")

message("\n── DEG summary: MC-38 ──")
print(summary_mouse, n = 20)

message("\n── DEG summary: HT29 ──")
print(summary_human, n = 10)

write.csv(bind_rows(summary_mouse, summary_human),
          "output/tables/DEG_summary_all.csv", row.names = FALSE)

write.csv(bind_rows(summary_mouse),
          "output/tables/DEG_summary_all.csv", row.names = FALSE)

message("✓ Сводная таблица: output/tables/DEG_summary_all.csv")

# =============================================================================
# ВИЗУАЛИЗАЦИИ — MC-38
# =============================================================================
message("\n══ Визуализации: MC-38 ══")

# ── PCA MC-38 ─────────────────────────────────────────────────────────────────
mat_m <- t(as.matrix(norm_counts_mouse))
mat_m <- mat_m[, apply(mat_m, 2, var) > 0]
pca_m <- prcomp(mat_m, scale. = TRUE, center = TRUE)
pca_m_var <- summary(pca_m)$importance[2, ] * 100

pca_df_mouse <- as.data.frame(pca_m$x[, 1:3]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(coldata_mouse, by = "sample_id")

p_pca_mouse <- ggplot(pca_df_mouse,
                      aes(x = PC1, y = PC2,
                          fill  = genotype,
                          shape = tnf,
                          label = sample_id)) +
  geom_point(size = 5, color = "white", stroke = 1.2, alpha = 0.95) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20,
                           segment.color = "grey60") +
  scale_fill_manual(values = GENO_COLORS, name = "Genotype") +
  scale_shape_manual(values = c("basal" = 21, "TNF12h" = 24),
                     labels = c("Basal", "TNF 12h"),
                     name   = "Treatment") +
  labs(title    = "PCA — MC38 (n=14)",
       subtitle = paste0("PC1: ", round(pca_m_var[1], 1),
                         "% | PC2: ", round(pca_m_var[2], 1), "%"),
       x = paste0("PC1 (", round(pca_m_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_m_var[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank()) +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 5)),
         shape = guide_legend(override.aes = list(fill = "grey50", size = 5)))

ggsave("output/figures/mouse/QC_PCA_MC38.pdf",
       p_pca_mouse, width = 9, height = 7)
ggsave("output/figures/mouse/QC_PCA_MC38.png",
       p_pca_mouse, width = 9, height = 7, dpi = 300)
message("✓ PCA MC-38 сохранён")

# ── Volcano MC-38 (baseline, панель) ─────────────────────────────────────────
volcano_mouse_base <- lapply(names(results_mouse_base), function(nm) {
  make_volcano(results_mouse_base[[nm]], nm)
})

p_vol_mouse <- patchwork::wrap_plots(volcano_mouse_base, ncol = 2) +
  patchwork::plot_annotation(
    title = "Volcano plots — MC38 baseline contrasts",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("output/figures/mouse/Volcano_baseline_MC38.pdf",
       p_vol_mouse, width = 14, height = 6 * ceiling(length(results_mouse_base) / 2))
ggsave("output/figures/mouse/Volcano_baseline_MC38.png",
       p_vol_mouse, width = 14, height = 6 * ceiling(length(results_mouse_base) / 2),
       dpi = 300)
message("✓ Volcano MC-38 baseline сохранён")

# ── DEG barplot MC-38 ─────────────────────────────────────────────────────────
p_bar_mouse <- make_deg_barplot(
  c(results_mouse_base, results_mouse_tnf),
  "DEG counts — MC38 (Mus musculus)"
)
ggsave("output/figures/mouse/DEG_counts_MC38.png",
       p_bar_mouse, width = 10, height = 8, dpi = 300)
message("✓ DEG barplot MC-38 сохранён")

# =============================================================================
# ВИЗУАЛИЗАЦИИ — HT29
# =============================================================================
message("\n══ Визуализации: HT29 ══")

# ── PCA HT29 ──────────────────────────────────────────────────────────────────
mat_h <- t(as.matrix(norm_counts_human))
mat_h <- mat_h[, apply(mat_h, 2, var) > 0]
pca_h <- prcomp(mat_h, scale. = TRUE, center = TRUE)
pca_h_var <- summary(pca_h)$importance[2, ] * 100

pca_df_human <- as.data.frame(pca_h$x[, 1:3]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(coldata_human, by = "sample_id")

p_pca_human <- ggplot(pca_df_human,
                      aes(x = PC1, y = PC2,
                          fill  = genotype,
                          label = sample_id)) +
  geom_point(size = 5, shape = 21, color = "white",
             stroke = 1.2, alpha = 0.95) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 15,
                           segment.color = "grey60") +
  scale_fill_manual(values = GENO_COLORS, name = "Genotype") +
  labs(title    = "PCA — HT29 (n=9)",
       subtitle = paste0("PC1: ", round(pca_h_var[1], 1),
                         "% | PC2: ", round(pca_h_var[2], 1), "%"),
       x = paste0("PC1 (", round(pca_h_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_h_var[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("output/figures/human/QC_PCA_HT29.pdf",
       p_pca_human, width = 8, height = 6)
ggsave("output/figures/human/QC_PCA_HT29.png",
       p_pca_human, width = 8, height = 6, dpi = 300)
message("✓ PCA HT29 сохранён")

# ── Volcano HT29 ─────────────────────────────────────────────────────────────
volcano_human <- lapply(names(results_human_base), function(nm) {
  make_volcano(results_human_base[[nm]], paste0("HT29 — ", nm))
})

p_vol_human <- patchwork::wrap_plots(volcano_human, ncol = 2) +
  patchwork::plot_annotation(
    title = "Volcano plots — HT29 baseline contrasts",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("output/figures/human/Volcano_baseline_HT29.pdf",
       p_vol_human, width = 14, height = 6 * ceiling(length(results_human_base) / 2))
ggsave("output/figures/human/Volcano_baseline_HT29.png",
       p_vol_human, width = 14, height = 6 * ceiling(length(results_human_base) / 2),
       dpi = 300)
message("✓ Volcano HT29 сохранён")

# ── DEG barplot HT29 ─────────────────────────────────────────────────────────
p_bar_human <- make_deg_barplot(
  results_human_base,
  "DEG counts — HT29 (Homo sapiens)"
)
ggsave("output/figures/human/DEG_counts_HT29.png",
       p_bar_human, width = 8, height = 5, dpi = 300)
message("✓ DEG barplot HT29 сохранён")

# =============================================================================
# СОХРАНЕНИЕ RDS-ОБЪЕКТОВ
# =============================================================================
message("\n══ Сохранение RDS ══")

# Объект для MC-38
pipeline_mouse <- list(
  organism     = "mouse",
  cell_line    = "MC38",
  norm_counts  = norm_counts_mouse,
  coldata      = coldata_mouse,
  results_base = results_mouse_base,
  results_tnf  = results_mouse_tnf,
  all_de       = bind_rows(c(results_mouse_base, results_mouse_tnf)),
  params       = list(
    FDR              = FDR,
    LFC_THRESH       = LFC_THRESH,
    date             = Sys.Date(),
    note_ppm1d       = "WIP1 = PPM1D (HGNC:9277); NOT PPM1A",
    note_exploratory = "TNF contrasts n=1: trend direction only",
    contrasts_base   = names(results_mouse_base),
    contrasts_tnf    = names(results_mouse_tnf)
  )
)

# Объект для HT29
pipeline_human <- list(
  organism     = "human",
  cell_line    = "HT29",
  norm_counts  = norm_counts_human,
  coldata      = coldata_human,
  results_base = results_human_base,
  results_tnf  = list(),           # TNF не предусмотрен для HT29
  all_de       = bind_rows(results_human_base),
  params       = list(
    FDR            = FDR,
    LFC_THRESH     = LFC_THRESH,
    date           = Sys.Date(),
    note_ppm1d     = "WIP1 = PPM1D (HGNC:9277); NOT PPM1A",
    contrasts_base = names(results_human_base)
  )
)

# NEW: объединённый объект для кросс-видового анализа (скрипт 05)
pipeline_all <- list(
  mouse        = pipeline_mouse,
  human        = pipeline_human,
  params       = list(
    FDR        = FDR,
    LFC_THRESH = LFC_THRESH,
    date       = Sys.Date(),
    # Контрасты, присутствующие в ОБОИХ датасетах — основа X1–X4
    shared_contrasts = intersect(names(results_mouse_base),
                                 names(results_human_base))
  )
)

saveRDS(pipeline_mouse, "output/rds/pipeline_mouse.rds")
saveRDS(pipeline_human, "output/rds/pipeline_human.rds")
saveRDS(pipeline_all,   "output/rds/pipeline_all.rds")

message("✓ pipeline_mouse.rds сохранён")
message("✓ pipeline_human.rds сохранён")
message("✓ pipeline_all.rds   сохранён")

# ── Финальная проверка ────────────────────────────────────────────────────────
test_m <- readRDS("output/rds/pipeline_mouse.rds")
test_h <- readRDS("output/rds/pipeline_human.rds")
test_a <- readRDS("output/rds/pipeline_all.rds")

message("\n══ ФИНАЛЬНАЯ ПРОВЕРКА ══")
message("MC-38  norm_counts : ",
        nrow(test_m$norm_counts), " × ", ncol(test_m$norm_counts))
message("MC-38  baseline    : ", length(test_m$results_base),
        " | exploratory: ",       length(test_m$results_tnf))
message("MC-38  total DEG rows: ", nrow(test_m$all_de))
message("HT29   norm_counts : ",
        nrow(test_h$norm_counts), " × ", ncol(test_h$norm_counts))
message("HT29   baseline    : ", length(test_h$results_base))
message("HT29   total DEG rows: ", nrow(test_h$all_de))
message("Shared contrasts   : ",
        paste(test_a$params$shared_contrasts, collapse = ", "))

message("\n══════════════════════════════════════════════════")
message("✓ Скрипт 01 завершён. Файлы:")
message("  RDS:     output/rds/pipeline_mouse.rds")
message("           output/rds/pipeline_human.rds")
message("           output/rds/pipeline_all.rds")
message("  Figures: output/figures/mouse/QC_PCA_MC38.png")
message("           output/figures/mouse/Volcano_baseline_MC38.png")
message("           output/figures/mouse/DEG_counts_MC38.png")
message("           output/figures/human/QC_PCA_HT29.png")
message("           output/figures/human/Volcano_baseline_HT29.png")
message("           output/figures/human/DEG_counts_HT29.png")
message("  Tables:  output/tables/DEG_summary_all.csv")
message("══════════════════════════════════════════════════")
