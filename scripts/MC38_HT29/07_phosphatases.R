# =============================================================================
# SCRIPT 07: PHOSPHATASE LANDSCAPE ANALYSIS
# Версия: v2 — исправлен нейминг (PPM1D/Wip1 ≠ PPM1A),
#              адаптирован под новый coldata (мышь MC38 + человек HT29)
#
# Цель: экспрессия фосфатаз в KO-моделях — валидация нокаута,
#       компенсаторная регуляция, ремоделинг фосфатазного ландшафта
#
# ИСПРАВЛЕНИЕ нейминга (обязательно к прочтению):
#   - BD02/BD05/BD20 = MC38 Wip1 KO = нокаут PPM1D (ген: Ppm1d у мыши / PPM1D у человека)
#   - В старом скрипте эта группа ошибочно называлась PPM1A_KO
#   - PPM1A (ген Ppm1a / PPM1A) — это другая фосфатаза (PP2Cα), не связанная с DDR
#   - В данном скрипте используется корректный нейминг: Wip1_KO / PPM1D_KO
#
# Входные данные (из скрипта 01):
#   pipeline$vst_mat   — VST-матрица [гены × образцы]
#   pipeline$coldata   — метаданные с колонками:
#                        sample_id, genotype, tnf, group, organism
#
# Генотипы в coldata (исправленные):
#   WT, Wip1_KO, PPM1B_KO, Double_KO  — у мыши (organism == "mouse")
#   WT, Wip1_KO, PPM1B_KO             — у человека (organism == "human")
#
# Выходные данные:
#   output/figures_phosphatase/PHOS_A_heatmap_all_{mouse|human}.png
#   output/figures_phosphatase/PHOS_B_heatmap_PPM_{mouse|human}.png
#   output/figures_phosphatase/PHOS_C_KO_validation_{mouse|human}.png
#   output/figures_phosphatase/PHOS_D_volcano_DE_{mouse|human}.png
#   output/figures_phosphatase/PHOS_E_dotplot_landscape_{mouse|human}.png
#   output/figures_phosphatase/PHOS_F_slope_TNF_mouse.png
#   output/figures_phosphatase/PHOS_G_crossspecies_comparison.png
#   output/tables_phosphatase/00_phosphatase_detected_{mouse|human}.csv
#   output/tables_phosphatase/01_phosphatase_DE_allcontrasts_{mouse|human}.csv
#   output/tables_phosphatase/02_phosphatase_summary_wide_{mouse|human}.csv
#   output/tables_phosphatase/03_phosphatase_significant_DE_{mouse|human}.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(pheatmap); library(RColorBrewer); library(tibble)
  library(stringr); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(scales); library(purrr)
  library(DESeq2)
})

dir.create("output/figures_phosphatase", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables_phosphatase",  showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0. ЗАГРУЗКА ДАННЫХ ИЗ ПАЙПЛАЙНА
# =============================================================================

pipeline <- readRDS("output/rds/pipeline_mouse.rds")

vst_mat <- if (!is.null(pipeline$vst_mat)) {
  pipeline$vst_mat
} else if (!is.null(pipeline$norm_counts)) {
  as.matrix(pipeline$norm_counts)
} else {
  stop("Не найдена VST/norm матрица в pipeline (ожидаем vst_mat или norm_counts)")
}
if (!is.matrix(vst_mat)) vst_mat <- as.matrix(vst_mat)

coldata <- pipeline$coldata

# ── Восстановление dds ───────────────────────────────────────────────────────
# Приоритет: pipeline$dds → отдельный rds → пересборка из raw counts
dds <- NULL

if (!is.null(pipeline$dds) && is(pipeline$dds, "DESeqDataSet")) {
  dds <- pipeline$dds
  message("✓ dds загружен из pipeline$dds")
} else if (file.exists("output/rds/dds_mouse.rds")) {
  dds <- readRDS("output/rds/dds_mouse.rds")
  message("✓ dds загружен из output/rds/dds_mouse.rds")
} else if (!is.null(pipeline$raw_counts) || !is.null(pipeline$counts)) {
  # Пересборка из raw counts если они есть в pipeline
  raw <- if (!is.null(pipeline$raw_counts)) pipeline$raw_counts
  else pipeline$counts
  message("⚙ dds отсутствует — пересборка из raw counts...")
  tryCatch({
    cd_dds <- coldata %>%
      dplyr::filter(sample_id %in% colnames(raw)) %>%
      dplyr::mutate(
        genotype = factor(as.character(genotype)),
        genotype = relevel(genotype, ref = "WT")
      ) %>%
      as.data.frame()
    rownames(cd_dds) <- cd_dds$sample_id
    raw_sub <- raw[, cd_dds$sample_id]
    dds <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(raw_sub)),
      colData   = cd_dds,
      design    = ~ genotype
    )
    dds <- DESeq(dds, quiet = TRUE)
    saveRDS(dds, "output/rds/dds_mouse.rds")
    message("✓ dds пересобран и сохранён → output/rds/dds_mouse.rds")
  }, error = function(e) {
    message("⚠ Не удалось пересобрать dds: ", e$message)
    message("  PHOS-D (volcano) будет пропущен")
    dds <<- NULL
  })
} else {
  message("⚠ dds недоступен и raw counts не найдены — PHOS-D будет пропущен")
}

# ── Проверка обязательных колонок ────────────────────────────────────────────
required_cols <- c("sample_id", "genotype", "tnf", "group")
missing_cols  <- setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0)
  stop("coldata не содержит колонок: ", paste(missing_cols, collapse = ", "))

has_organism <- "organism" %in% colnames(coldata)
if (!has_organism) {
  message("⚠ Колонка 'organism' не найдена в coldata — предполагаем только мышь (MC38)")
  coldata$organism <- "mouse"
}

message("✓ vst_mat: ",  nrow(vst_mat), " генов × ", ncol(vst_mat), " образцов")
message("✓ coldata: ",  nrow(coldata), " образцов")
message("  Организмы: ", paste(sort(unique(coldata$organism)), collapse = ", "))
message("  Генотипы:  ", paste(sort(unique(coldata$genotype)),  collapse = ", "))
message("  dds статус: ", ifelse(is.null(dds), "НЕДОСТУПЕН (PHOS-D пропущен)", "OK"))


# =============================================================================
# 1. ПАНЕЛЬ ФОСФАТАЗ — ДВЕ ВЕРСИИ (МЫШЬ / ЧЕЛОВЕК)
# =============================================================================
# Символы генов: мышиные (строчные) vs. человеческие (прописные первая буква)
# Соответствие: Ppm1d (мышь) <-> PPM1D (человек), и т.д.

phosphatase_genes_mouse <- c(
  # PPM/PP2C — серин/треониновые, Mg2+-зависимые
  # ВАЖНО: Ppm1d = Wip1, основной субъект нашего проекта (НЕ Ppm1a!)
  "Ppm1d",  # Wip1 — главный KO в проекте
  "Ppm1a",  # PP2Cα — другая фосфатаза (НЕ нокаутирована у мыши)
  "Ppm1b",  # PP2Cβ — второй KO
  "Ppm1e", "Ppm1f", "Ppm1g", "Ppm1h", "Ppm1j",
  # PPP — PP1, PP2A, PP2B (кальциневрин)
  "Ppp1ca", "Ppp1cb", "Ppp1cc",
  "Ppp2ca", "Ppp2cb",
  "Ppp3ca", "Ppp3cb",
  # PTP — тирозиновые фосфатазы
  "Ptpn1", "Ptpn2", "Ptpn3", "Ptpn4", "Ptpn5",
  "Ptpn6", "Ptpn7", "Ptpn9",
  "Ptpn11", "Ptpn12", "Ptpn13", "Ptpn14",
  "Ptpn18", "Ptpn20", "Ptpn21", "Ptpn22", "Ptpn23",
  "Ptpra", "Ptprb", "Ptprc", "Ptprd", "Ptpre", "Ptprf",
  "Ptprg", "Ptprh", "Ptprj", "Ptprk", "Ptprm",
  "Ptprn", "Ptprn2", "Ptpro", "Ptprq", "Ptprs",
  "Ptprt", "Ptpru", "Ptprz1",
  # ACP/ALP — кислые/щелочные фосфатазы
  "Acp1", "Acp2", "Acp3", "Acp5",
  "Alpl", "Alpi",
  # DUSP — MAP-киназные фосфатазы (прямые мишени Wip1/p38MAPK-оси)
  "Dusp1", "Dusp4", "Dusp6", "Dusp10", "Dusp16",
  # Slingshot — Rho/Actin ось (связь Wip1 → Rho → Cofilin)
  "Ssh1", "Ssh2", "Ssh3",
  # PTEN
  "Pten"
)

phosphatase_genes_human <- c(
  # PPM/PP2C
  "PPM1D",  # Wip1 — KO у человека (HT29 WipKo)
  "PPM1A",  # PP2Cα — другая фосфатаза (НЕ нокаутирована)
  "PPM1B",  # PP2Cβ — KO у человека (HT29 PPM1BKO)
  "PPM1E", "PPM1F", "PPM1G", "PPM1H", "PPM1J",
  # PPP
  "PPP1CA", "PPP1CB", "PPP1CC",
  "PPP2CA", "PPP2CB",
  "PPP3CA", "PPP3CB",
  # PTP non-receptor
  "PTPN1", "PTPN2", "PTPN3", "PTPN4", "PTPN5",
  "PTPN6", "PTPN7", "PTPN9",
  "PTPN11", "PTPN12", "PTPN13", "PTPN14",
  "PTPN18", "PTPN20", "PTPN21", "PTPN22", "PTPN23",
  # PTP receptor
  "PTPRA", "PTPRB", "PTPRC", "PTPRD", "PTPRE", "PTPRF",
  "PTPRG", "PTPRH", "PTPRJ", "PTPRK", "PTPRM",
  "PTPRN", "PTPRN2", "PTPRO", "PTPRQ", "PTPRS",
  "PTPRT", "PTPRU", "PTPRZ1",
  # ACP/ALP
  "ACP1", "ACP2", "ACP5",
  "ALPL", "ALPI",
  # DUSP
  "DUSP1", "DUSP4", "DUSP6", "DUSP10", "DUSP16",
  # Slingshot
  "SSH1", "SSH2", "SSH3",
  # PTEN
  "PTEN"
)

# Функция классификации (работает для обоих видов)
classify_phosphatase <- function(genes) {
  dplyr::case_when(
    str_detect(genes, regex("^Ppm|^PPM",   ignore_case = FALSE)) ~ "PPM/PP2C",
    str_detect(genes, regex("^Ppp|^PPP",   ignore_case = FALSE)) ~ "PPP",
    str_detect(genes, regex("^Ptpn|^PTPN", ignore_case = FALSE)) ~ "PTPN (non-receptor)",
    str_detect(genes, regex("^Ptpr|^PTPR", ignore_case = FALSE)) ~ "PTPR (receptor)",
    str_detect(genes, regex("^Acp|^ACP|^Alp|^ALP", ignore_case = FALSE)) ~ "ACP/ALP",
    str_detect(genes, regex("^Dusp|^DUSP", ignore_case = FALSE)) ~ "DUSP",
    str_detect(genes, regex("^Ssh|^SSH",   ignore_case = FALSE)) ~ "Slingshot",
    str_detect(genes, regex("^Pten$|^PTEN$", ignore_case = FALSE)) ~ "PTEN",
    TRUE ~ "Other"
  )
}

class_colors <- c(
  "PPM/PP2C"            = "#E41A1C",
  "PPP"                 = "#377EB8",
  "PTPN (non-receptor)" = "#4DAF4A",
  "PTPR (receptor)"     = "#984EA3",
  "ACP/ALP"             = "#FF7F00",
  "DUSP"                = "#A65628",
  "Slingshot"           = "#F781BF",
  "PTEN"                = "#999999",
  "Other"               = "#CCCCCC"
)

# =============================================================================
# 2. ПАРАМЕТРЫ ВИЗУАЛИЗАЦИИ — УНИФИЦИРОВАННЫЕ ДЛЯ ОБОИХ ВИДОВ
# =============================================================================

# ИСПРАВЛЕНО: PPM1D_KO вместо PPM1A_KO
geno_order_mouse <- c("WT", "PPM1D_KO", "PPM1B_KO", "DKO")
geno_order_human <- c("WT", "PPM1D_KO", "PPM1B_KO")

geno_colors <- c(
  "WT"        = "#4DAF4A",
  "PPM1D_KO"   = "#E41A1C",   # PPM1D KO — красный (был PPM1A_KO)
  "PPM1B_KO"  = "#FF7F00",
  "DKO" = "#984EA3"
)

condition_colors_mouse <- c(
  "WT"          = "#4DAF4A",
  "PPM1D_KO"     = "#E41A1C",
  "PPM1B_KO"    = "#FF7F00",
  "DKO"   = "#984EA3",
  "WT_TNF"      = "#A6D96A",
  "PPM1D_KO_TNF" = "#FB9A99",
  "PPM1B_KO_TNF"= "#FDBF6F"
)

condition_colors_human <- c(
  "WT"       = "#4DAF4A",
  "PPM1D_KO_KO"  = "#E41A1C",
  "PPM1B_KO" = "#FF7F00"
)

# Тема для публикации
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey93"),
    strip.text       = element_text(face = "bold", size = 9),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    legend.position  = "right"
  )

# =============================================================================
# 3. ЯДРО АНАЛИЗА — ФУНКЦИЯ run_phosphatase_analysis(organism)
#    Запускается отдельно для "mouse" и "human"
# =============================================================================

run_phosphatase_analysis <- function(organism) {

  stopifnot(organism %in% c("mouse", "human"))
  org_label <- ifelse(organism == "mouse", "MC38 (mouse)", "HT29 (human)")
  message("\n══════════════════════════════════════════════════")
  message("  Organism: ", org_label)
  message("══════════════════════════════════════════════════")

  # --- 3.1 Параметры, специфичные для вида ---
  phos_genes_all <- if (organism == "mouse") phosphatase_genes_mouse
  else                     phosphatase_genes_human

  # Гены KO-мишеней (ИСПРАВЛЕНО: PPM1D = Wip1, не PPM1A!)
  ko_targets <- if (organism == "mouse") c("Ppm1d", "Ppm1b")
  else                     c("PPM1D", "PPM1B")

  geno_order <- if (organism == "mouse") geno_order_mouse
  else                     geno_order_human

  cond_colors <- if (organism == "mouse") condition_colors_mouse
  else                     condition_colors_human

  # --- 3.2 Подматрица образцов текущего вида ---
  coldata_org <- coldata %>%
    dplyr::filter(organism == !!organism) %>%
    dplyr::mutate(
      genotype  = factor(genotype, levels = geno_order),
      tnf       = factor(tnf, levels = c("basal", "TNF12h"))
    )

  available_samples <- intersect(coldata_org$sample_id, colnames(vst_mat))
  if (length(available_samples) == 0)
    stop("Нет образцов организма '", organism, "' в vst_mat!")
  message("  Образцы (", length(available_samples), "): ",
          paste(available_samples, collapse = ", "))

  # Канонический порядок образцов
  sample_order <- coldata_org %>%
    dplyr::filter(sample_id %in% available_samples) %>%
    dplyr::arrange(tnf, genotype) %>%
    dplyr::pull(sample_id)

  vst_org <- vst_mat[, sample_order]

  # --- 3.3 Детектированные фосфатазы ---
  phos_detected <- phos_genes_all[phos_genes_all %in% rownames(vst_org)]
  phos_missing  <- setdiff(phos_genes_all, rownames(vst_org))

  message(sprintf("  Фосфатаз обнаружено: %d / %d",
                  length(phos_detected), length(phos_genes_all)))
  if (length(phos_missing) > 0 && length(phos_missing) <= 15)
    message("  Отсутствуют: ", paste(phos_missing, collapse = ", "))

  phos_class <- setNames(classify_phosphatase(phos_detected), phos_detected)

  write.csv(
    data.frame(gene  = phos_detected,
               class = phos_class,
               ko_target = phos_detected %in% ko_targets),
    sprintf("output/tables_phosphatase/00_phosphatase_detected_%s.csv", organism),
    row.names = FALSE
  )

  # --- 3.4 Подматрица фосфатаз + z-score ---
  phos_mat   <- vst_org[phos_detected, ]
  phos_mat_z <- t(scale(t(phos_mat)))

  # Длинный формат
  phos_long <- phos_mat_z %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "zscore") %>%
    dplyr::left_join(coldata_org, by = "sample_id") %>%
    dplyr::mutate(
      class     = phos_class[gene],
      genotype  = factor(genotype, levels = geno_order),
      is_target = gene %in% ko_targets
    )

  # --- 3.5 Аннотации для heatmap ---
  # Аннотация столбцов
  col_ann_df <- coldata_org %>%
    dplyr::filter(sample_id %in% sample_order) %>%
    dplyr::select(sample_id, genotype, tnf) %>%
    dplyr::mutate(genotype = as.character(genotype),
                  tnf      = as.character(tnf)) %>%
    tibble::column_to_rownames("sample_id") %>%
    .[sample_order, ]

  # Вычисляем позиции для gaps_col (между генотипами)
  gap_positions <- which(diff(as.integer(factor(col_ann_df$genotype,
                                                levels = geno_order))) != 0)

  ann_colors <- list(
    genotype  = geno_colors[intersect(names(geno_colors), unique(col_ann_df$genotype))],
    tnf       = c("basal" = "#F0F0F0", "TNF12h" = "#FC8D59"),
    Class     = class_colors,
    KO_target = c("KO gene" = "#E41A1C", "Other" = "#CCCCCC")
  )

  # =============================================================================
  # VIZ PHOS-A: HEATMAP — все фосфатазы × все образцы
  # =============================================================================

  row_ann_all <- data.frame(
    Class     = phos_class[rownames(phos_mat_z)],
    KO_target = ifelse(rownames(phos_mat_z) %in% ko_targets, "KO gene", "Other"),
    row.names = rownames(phos_mat_z)
  )

  # ── ПАТЧ: пересобрать col_ann_df_ord прямо из coldata ────────────────────────
  # Не зависит от внешнего окружения; гарантирует все 4 генотипа в аннотации

  # 1. Нормализуем geno_colors — добавляем DKO если его нет
  if (!"DKO" %in% names(geno_colors) && "DKO" %in% unique(coldata$genotype)) {
    geno_colors["DKO"] <- "#8856A7"
  }
  # Аналогично для Double_KO (старый алиас)
  if (!"Double_KO" %in% names(geno_colors) && "Double_KO" %in% unique(coldata$genotype)) {
    geno_colors["Double_KO"] <- "#8856A7"
  }

  # 2. Пересобираем col_ann_df_ord из coldata напрямую
  col_ann_df_ord <- coldata %>%
    dplyr::filter(sample_id %in% colnames(vst_mat)) %>%
    dplyr::select(sample_id, genotype, tnf) %>%
    dplyr::mutate(
      genotype = as.character(genotype),   # снимаем старый factor с кривыми уровнями
      genotype = factor(genotype, levels = geno_order)
    ) %>%
    tibble::column_to_rownames("sample_id")

  # 3. Порядок образцов: WT → PPM1D_KO → PPM1B_KO → DKO, внутри — по tnf
  sample_order <- coldata %>%
    dplyr::filter(sample_id %in% colnames(vst_mat)) %>%
    dplyr::mutate(
      genotype = factor(as.character(genotype), levels = geno_order),
      tnf      = factor(tnf, levels = c("basal", "TNF12h"))
    ) %>%
    dplyr::arrange(genotype, tnf) %>%
    dplyr::pull(sample_id)

  # 4. Позиции gaps_col: между группами генотипов
  gap_positions <- coldata %>%
    dplyr::filter(sample_id %in% sample_order) %>%
    dplyr::mutate(genotype = factor(as.character(genotype), levels = geno_order)) %>%
    dplyr::arrange(factor(sample_id, levels = sample_order)) %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::mutate(cumsum = cumsum(n)) %>%
    dplyr::pull(cumsum) %>%
    head(-1)   # убираем последний (после последней группы gap не нужен)

  # 5. Цвета аннотации — только реально присутствующие уровни
  present_genotypes <- levels(droplevels(col_ann_df_ord$genotype))
  present_genotypes <- present_genotypes[!is.na(present_genotypes)]

  ann_phos_colors <- list(
    Class     = class_colors,
    KO_target = c("KO gene" = "#D73027", "Other" = "#CCCCCC"),
    genotype  = geno_colors[present_genotypes],   # <── только присутствующие
    tnf       = c("basal" = "#F0F0F0", "TNF12h" = "#FC8D59")
  )

  # Диагностика перед построением:
  cat("\n── col_ann_df_ord genotype levels:\n")
  print(table(col_ann_df_ord$genotype, useNA = "always"))
  cat("── geno_colors покрывает:", paste(names(ann_phos_colors$genotype), collapse=", "), "\n")
  cat("── sample_order (", length(sample_order), "образцов):", paste(sample_order, collapse=" "), "\n\n")
  # ─────────────────────────────────────────────────────────────────────────────

  pheatmap::pheatmap(
    phos_mat_z[, sample_order],
    annotation_col    = col_ann_df,
    annotation_row    = row_ann_all,
    annotation_colors = ann_colors,
    color    = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks   = seq(-2.5, 2.5, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col     = gap_positions,
    show_rownames = TRUE,
    fontsize_row  = 7,
    fontsize_col  = 7,
    border_color  = NA,
    cellwidth  = 25,
   # cellheight = 10,
    main = sprintf("Phosphatase expression landscape — %s (z-score)", org_label),
    filename = sprintf("output/figures_phosphatase/PHOS_A_heatmap_all_%s.png", organism),
    width = 16, height = 14
  )
  message("  ✓ PHOS-A heatmap all")

  # =============================================================================
  # VIZ PHOS-B: ZOOM — только PPM-семейство
  # =============================================================================

  ppm_prefix <- if (organism == "mouse") "^Ppm" else "^PPM"
  ppm_genes  <- phos_detected[str_detect(phos_detected, ppm_prefix)]

  if (length(ppm_genes) >= 2) {
    ppm_mat_z  <- phos_mat_z[ppm_genes, sample_order]
    row_ann_ppm <- data.frame(
      KO_target = ifelse(ppm_genes %in% ko_targets, "KO gene", "Other"),
      row.names  = ppm_genes
    )
    pheatmap::pheatmap(
      ppm_mat_z,
      annotation_col    = col_ann_df,
      annotation_row    = row_ann_ppm,
      annotation_colors = ann_colors[c("genotype","tnf","KO_target")],
      color    = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks   = seq(-2.5, 2.5, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      gaps_col     = gap_positions,
      fontsize_row  = 10,
      fontsize_col  = 9,
      border_color  = NA,
      cellwidth  = 28,
      cellheight = 22,
      # ИСПРАВЛЕНО: заголовок отражает PPM1D (Wip1), а не PPM1A
      main = sprintf("PPM/PP2C family [incl. PPM1D/Wip1 & PPM1B] — %s", org_label),
      filename = sprintf("output/figures_phosphatase/PHOS_B_heatmap_PPM_%s.png", organism),
      width = 13, height = max(5, length(ppm_genes) * 0.7 + 3)
    )
    message("  ✓ PHOS-B heatmap PPM")
  } else {
    message("  ⚠ PHOS-B пропущен — < 2 PPM-генов")
  }

  # =============================================================================
  # VIZ PHOS-C: BARPLOT валидации нокаута — PPM1D и PPM1B
  # =============================================================================

  target_expr <- phos_long %>%
    dplyr::filter(is_target, tnf == "basal") %>%
    dplyr::group_by(gene, genotype) %>%
    dplyr::summarise(
      mean_z = mean(zscore, na.rm = TRUE),
      sd_z   = sd(zscore,   na.rm = TRUE),
      n      = n(),
      se_z   = ifelse(n > 1, sd_z / sqrt(n), 0),
      .groups = "drop"
    )

  # Метки для заголовка: корректное название генов
  ko_labels <- if (organism == "mouse") {
    c("Ppm1d" = "Ppm1d (Wip1)", "Ppm1b" = "Ppm1b")
  } else {
    c("PPM1D" = "PPM1D (Wip1)", "PPM1B" = "PPM1B")
  }
  target_expr$gene_label <- ko_labels[target_expr$gene]
  phos_long_targets <- phos_long %>%
    dplyr::filter(is_target, tnf == "basal") %>%
    dplyr::mutate(gene_label = ko_labels[gene])

  p_ko_val <- ggplot(target_expr,
                     aes(x = genotype, y = mean_z, fill = genotype)) +
    geom_col(width = 0.65, alpha = 0.9, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                  width = 0.25, linewidth = 0.7, color = "grey30") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_jitter(
      data = phos_long_targets,
      aes(x = genotype, y = zscore, color = genotype),
      width = 0.12, size = 2.8, alpha = 0.85, show.legend = FALSE
    ) +
    geom_hline(
      data = target_expr %>% dplyr::filter(genotype == "WT"),
      aes(yintercept = mean_z),
      linetype = "dotted", color = "grey30", linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values  = geno_colors) +
    scale_color_manual(values = geno_colors) +
    facet_wrap(~ gene_label, scales = "free_y", ncol = 2) +
    labs(
      # ИСПРАВЛЕНО: PPM1D (Wip1), не PPM1A
      title    = sprintf("KO validation: PPM1D (Wip1) and PPM1B expression — %s", org_label),
      subtitle = "Basal conditions | VST z-score | dotted = WT mean | bars = mean ± SE",
      x = NULL, y = "VST z-score",
      fill = "Genotype"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    sprintf("output/figures_phosphatase/PHOS_C_KO_validation_%s.png", organism),
    p_ko_val, width = 9, height = 5, dpi = 300, bg = "white"
  )
  ggsave(
    sprintf("output/figures_phosphatase/PHOS_C_KO_validation_%s.pdf", organism),
    p_ko_val, width = 9, height = 5
  )
  message("  ✓ PHOS-C KO validation")

  # =============================================================================
  # VIZ PHOS-D: VOLCANO — DE фосфатаз по всем контрастам
  # =============================================================================

  # Извлекаем результаты DESeq2 для контрастов данного вида
  # Ожидаем, что в dds есть поле metadata с organism, или используем отдельный dds
  # Если dds единый для всех образцов — фильтруем результаты по генам вида

  FACTOR_NAME <- "genotype"
  available_contrasts <- tryCatch(
    DESeq2::resultsNames(dds),
    error = function(e) {
      message("  ⚠ dds недоступен: ", e$message)
      character(0)
    }
  )

  if (length(available_contrasts) > 0) {
    contrast_names <- available_contrasts[
      str_detect(available_contrasts, FACTOR_NAME) &
        !str_detect(available_contrasts, "Intercept")
    ]
    if (length(contrast_names) == 0)
      contrast_names <- available_contrasts[!str_detect(available_contrasts, "Intercept")]

    message("  Контрасты: ", paste(contrast_names, collapse = ", "))

    phos_de <- purrr::map_dfr(contrast_names, function(cn) {
      tryCatch({
        res_tmp <- DESeq2::results(dds, name = cn,
                                   independentFiltering = TRUE,
                                   cooksCutoff = TRUE)
        as.data.frame(res_tmp) %>%
          rownames_to_column("gene") %>%
          dplyr::filter(gene %in% phos_detected) %>%
          dplyr::mutate(
            contrast  = cn,
            class     = phos_class[gene],
            is_target = gene %in% ko_targets,
            sig = dplyr::case_when(
              padj < 0.05 & log2FoldChange >  1 ~ "Up",
              padj < 0.05 & log2FoldChange < -1 ~ "Down",
              TRUE ~ "NS"
            )
          )
      }, error = function(e) {
        message("  Пропускаем контраст ", cn, ": ", e$message)
        NULL
      })
    })

    write.csv(phos_de,
              sprintf("output/tables_phosphatase/01_phosphatase_DE_allcontrasts_%s.csv",
                      organism),
              row.names = FALSE)

    phos_de_plot <- phos_de %>%
      dplyr::filter(!is.na(padj), !is.na(log2FoldChange))

    if (nrow(phos_de_plot) > 0) {
      de_colors <- c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")
      label_data <- phos_de_plot %>% dplyr::filter(sig != "NS" | is_target)

      p_volcano <- ggplot(
        phos_de_plot,
        aes(x = log2FoldChange, y = -log10(padj),
            color = sig, size = is_target)
      ) +
        geom_point(alpha = 0.75) +
        geom_vline(xintercept = c(-1, 1),
                   linetype = "dashed", color = "grey40", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", color = "grey40", linewidth = 0.5) +
        scale_color_manual(values = de_colors, name = "DE status") +
        scale_size_manual(
          # ИСПРАВЛЕНО: PPM1D/PPM1B, не PPM1A/PPM1B
          values = c("TRUE" = 4.5, "FALSE" = 2.5),
          guide  = "none"
        ) +
        labs(
          title    = sprintf("DE of phosphatase genes — %s", org_label),
          subtitle = "FDR < 0.05, |log2FC| ≥ 1 | large dots = PPM1D (Wip1) / PPM1B",
          x = "log2 Fold Change (vs. WT)", y = "-log10(FDR)"
        ) +
        theme_pub

      if (nrow(label_data) > 0) {
        p_volcano <- p_volcano +
          ggrepel::geom_label_repel(
            data = label_data, aes(label = gene),
            size = 2.8, max.overlaps = 30,
            box.padding = 0.35, seed = 42, show.legend = FALSE
          )
      }

      use_facet <- length(unique(phos_de_plot$contrast)) > 1
      if (use_facet)
        p_volcano <- p_volcano +
        facet_wrap(~ contrast, scales = "free", ncol = 2)

      n_cont <- length(unique(phos_de_plot$contrast))
      fig_h  <- if (use_facet) ceiling(n_cont / 2) * 4.5 + 1.5 else 6
      fig_w  <- if (use_facet) 14 else 8

      ggsave(
        sprintf("output/figures_phosphatase/PHOS_D_volcano_DE_%s.png", organism),
        p_volcano, width = fig_w, height = fig_h, dpi = 300, bg = "white"
      )
      ggsave(
        sprintf("output/figures_phosphatase/PHOS_D_volcano_DE_%s.pdf", organism),
        p_volcano, width = fig_w, height = fig_h
      )
      message("  ✓ PHOS-D volcano")

      # Сводные таблицы
      summary_wide <- phos_de %>%
        dplyr::filter(!is.na(padj)) %>%
        dplyr::select(gene, class, contrast, log2FoldChange, padj, sig) %>%
        tidyr::pivot_wider(
          names_from  = contrast,
          values_from = c(log2FoldChange, padj, sig)
        ) %>%
        dplyr::arrange(class, gene)

      sig_table <- phos_de %>%
        dplyr::filter(sig != "NS", !is.na(padj)) %>%
        dplyr::select(gene, class, contrast, log2FoldChange, padj, sig) %>%
        dplyr::arrange(class, contrast, padj)

      write.csv(summary_wide,
                sprintf("output/tables_phosphatase/02_phosphatase_summary_wide_%s.csv",
                        organism),
                row.names = FALSE)
      write.csv(sig_table,
                sprintf("output/tables_phosphatase/03_phosphatase_significant_DE_%s.csv",
                        organism),
                row.names = FALSE)
      message(sprintf("  ✓ Таблицы DE | значимых: %d", nrow(sig_table)))

    } else {
      message("  ⚠ phos_de_plot пуст — PHOS-D пропущен")
    }
  } else {
    message("  ⚠ dds недоступен — DE-анализ пропущен")
  }

  # =============================================================================
  # VIZ PHOS-E: DOT PLOT — средний z-score × генотип (только basal)
  # =============================================================================

  dot_data <- phos_long %>%
    dplyr::filter(tnf == "basal") %>%
    dplyr::group_by(gene, genotype, class) %>%
    dplyr::summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

  gene_order_dot <- dot_data %>%
    dplyr::filter(genotype != "WT") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(max_abs = max(abs(mean_z)), .groups = "drop") %>%
    dplyr::arrange(desc(max_abs)) %>%
    dplyr::pull(gene)
  dot_data$gene <- factor(dot_data$gene, levels = rev(gene_order_dot))

  p_dot <- ggplot(dot_data,
                  aes(x = genotype, y = gene,
                      fill = mean_z, size = abs(mean_z))) +
    geom_point(shape = 21, color = "white", stroke = 0.3) +
    scale_fill_distiller(
      palette = "RdBu", direction = -1,
      limits = c(-2.5, 2.5), oob = scales::squish,
      name = "Mean z-score"
    ) +
    scale_size_continuous(range = c(1, 8), name = "|z-score|") +
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(limits = geno_order) +
    labs(
      title    = sprintf("Phosphatase landscape — basal conditions — %s", org_label),
      subtitle = "Dot size and color = mean VST z-score | PPM1D (Wip1) highlighted",
      x = "Genotype", y = NULL
    ) +
    theme_pub +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      axis.text.y  = element_text(size = 7),
      strip.text.y = element_text(angle = 0, size = 8),
      panel.spacing = unit(0.2, "lines")
    )

  ggsave(
    sprintf("output/figures_phosphatase/PHOS_E_dotplot_landscape_%s.png", organism),
    p_dot, width = 9,
    height = max(10, length(phos_detected) * 0.18 + 4),
    dpi = 300, bg = "white"
  )
  message("  ✓ PHOS-E dot plot landscape")

  # =============================================================================
  # VIZ PHOS-F: SLOPE PLOT — basal → TNF (только мышь: есть TNF-образцы)
  # =============================================================================

  if (organism == "mouse") {
    phos_slope <- phos_long %>%
      dplyr::filter(genotype != "Double_KO") %>%  # нет TNF данных для Double KO
      dplyr::group_by(gene, genotype, tnf, class) %>%
      dplyr::summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

    responsive_genes <- phos_slope %>%
      tidyr::pivot_wider(names_from = tnf, values_from = mean_z) %>%
      dplyr::filter(!is.na(TNF12h)) %>%
      dplyr::mutate(delta = abs(TNF12h - basal)) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(max_delta = max(delta, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(max_delta > 0.5) %>%
      dplyr::pull(gene)

    if (length(responsive_genes) > 0) {
      p_slope <- phos_slope %>%
        dplyr::filter(gene %in% responsive_genes) %>%
        ggplot(aes(x = tnf, y = mean_z, color = genotype, group = genotype)) +
        geom_line(linewidth = 1.1, alpha = 0.85) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(
          data = . %>% dplyr::filter(tnf == "TNF12h"),
          aes(label = genotype),
          size = 2.5, nudge_x = 0.15,
          show.legend = FALSE, seed = 42
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        scale_color_manual(values = geno_colors) +
        scale_x_discrete(limits = c("basal", "TNF12h")) +
        facet_wrap(~ gene, scales = "free_y", ncol = 5) +
        labs(
          title    = "Phosphatase TNF-response trajectories — MC38 (mouse)",
          subtitle = "Genes with Δz > 0.5 in ≥1 genotype | slope = direction of change",
          x = "Condition", y = "Mean VST z-score", color = "Genotype"
        ) +
        theme_pub +
        theme(strip.text = element_text(size = 7.5))

      ggsave(
        "output/figures_phosphatase/PHOS_F_slope_TNF_mouse.png",
        p_slope,
        width  = max(14, ceiling(length(responsive_genes) / 5) * 2 + 4),
        height = ceiling(length(responsive_genes) / 5) * 2.5 + 2,
        dpi = 300, bg = "white"
      )
      message("  ✓ PHOS-F slope TNF (", length(responsive_genes),
              " responsive genes)")
    } else {
      message("  ⚠ PHOS-F пропущен — нет генов с Δz > 0.5 при TNF")
    }
  }

  # Возвращаем объекты для кросс-видового анализа
  invisible(list(
    phos_long    = phos_long,
    phos_mat_z   = phos_mat_z,
    phos_class   = phos_class,
    ko_targets   = ko_targets,
    geno_order   = geno_order,
    organism     = organism
  ))
}

# =============================================================================
# 4. ЗАПУСК АНАЛИЗА ДЛЯ КАЖДОГО ОРГАНИЗМА
# =============================================================================

results_mouse <- run_phosphatase_analysis("mouse")
results_human <- run_phosphatase_analysis("human")

# =============================================================================
# 5. VIZ PHOS-G: КРОСС-ВИДОВОЕ СРАВНЕНИЕ — PPM-семейство
#    Сравниваем z-score PPM1D и PPM1B в Wip1_KO и PPM1B_KO
#    у мыши (MC38) и человека (HT29)
# =============================================================================

# Общий пул генов: мышиный символ → человеческий ортолог
orthologs_ppm <- tibble::tribble(
  ~mouse_gene, ~human_gene, ~protein,
  "Ppm1d",     "PPM1D",     "PPM1D/Wip1",
  "Ppm1b",     "PPM1B",     "PPM1B",
  "Ppm1a",     "PPM1A",     "PPM1A",
  "Ppm1e",     "PPM1E",     "PPM1E",
  "Ppm1f",     "PPM1F",     "PPM1F",
  "Ppm1g",     "PPM1G",     "PPM1G"
)

# Данные мыши — KO для PPM1D и PPM1B, basal
mouse_ppm <- results_mouse$phos_long %>%
  dplyr::filter(
    gene %in% orthologs_ppm$mouse_gene,
    tnf == "basal"
  ) %>%
  dplyr::group_by(gene, genotype) %>%
  dplyr::summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(orthologs_ppm, by = c("gene" = "mouse_gene")) %>%
  dplyr::mutate(organism = "mouse (MC38)")

# Данные человека
human_ppm <- results_human$phos_long %>%
  dplyr::filter(
    gene %in% orthologs_ppm$human_gene,
    tnf == "basal"
  ) %>%
  dplyr::group_by(gene, genotype) %>%
  dplyr::summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(orthologs_ppm, by = c("gene" = "human_gene")) %>%
  dplyr::mutate(organism = "human (HT29)")

crossspecies_ppm <- dplyr::bind_rows(mouse_ppm, human_ppm) %>%
  dplyr::filter(!is.na(protein))

if (nrow(crossspecies_ppm) > 0) {
  p_cross <- ggplot(crossspecies_ppm,
                    aes(x = genotype, y = mean_z,
                        fill = genotype, alpha = organism)) +
    geom_col(position = "dodge", width = 0.7, color = "white", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = geno_colors) +
    scale_alpha_manual(values = c("mouse (MC38)" = 1.0, "human (HT29)" = 0.6)) +
    facet_wrap(~ protein, scales = "free_y", ncol = 3) +
    labs(
      title    = "PPM/PP2C family: cross-species comparison (MC38 mouse vs. HT29 human)",
      subtitle = "Mean VST z-score | basal conditions | PPM1D = Wip1 (NOT PPM1A!)",
      x = "Genotype", y = "Mean VST z-score",
      fill = "Genotype", alpha = "Organism"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    "output/figures_phosphatase/PHOS_G_crossspecies_comparison.png",
    p_cross, width = 14, height = 8, dpi = 300, bg = "white"
  )
  message("✓ PHOS-G cross-species comparison")
} else {
  message("⚠ PHOS-G пропущен — нет общих PPM генов между видами")
}

# =============================================================================
# 6. ИТОГОВЫЙ ОТЧЁТ
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  PHOSPHATASE LANDSCAPE ANALYSIS v2 — SUMMARY\n")
cat("  ИСПРАВЛЕНИЕ: Wip1 KO = PPM1D (не PPM1A!)\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("  Mouse (MC38): %d фосфатаз обнаружено\n",
            length(results_mouse$phos_class)))
cat(sprintf("  Human (HT29): %d фосфатаз обнаружено\n",
            length(results_human$phos_class)))
cat("\n  Figures → output/figures_phosphatase/\n")
cat("    PHOS_A — heatmap все фосфатазы (mouse + human)\n")
cat("    PHOS_B — zoom PPM/PP2C семейство (mouse + human)\n")
cat("    PHOS_C — KO validation: PPM1D (Wip1) + PPM1B (mouse + human)\n")
cat("    PHOS_D — volcano DE фосфатаз (mouse + human)\n")
cat("    PHOS_E — dot plot landscape (mouse + human)\n")
cat("    PHOS_F — slope TNF-response (mouse only)\n")
cat("    PHOS_G — cross-species PPM family comparison\n")
cat("\n  Tables  → output/tables_phosphatase/\n")
cat("    00 — detected phosphatases list\n")
cat("    01 — DE all contrasts\n")
cat("    02 — summary wide\n")
cat("    03 — significant DE only\n")
cat("═══════════════════════════════════════════════════════════\n\n")
