# =============================================================================
# 06_advanced_visualizations_v3.R
#
# ИЗМЕНЕНИЯ ОТНОСИТЕЛЬНО v2:
#   - ИСПРАВЛЕНО: PPM1A_KO → Wip1KO везде в именах, подписях, палитрах
#     (PPM1D/Wip1 — мишень скрипта; PPM1A — несвязанная фосфатаза)
#   - НОВОЕ: поддержка двух видов — мышь (MC38) и человек (HT29)
#   - НОВОЕ: VIZ-7 — кросс-видовое сравнение z-score профилей
#   - Double_KO только у мыши — честно отражено в facet-структуре
#
# БИОЛОГИЧЕСКИЙ КОНТЕКСТ:
#   Скрипт визуализирует результаты scoring-анализа (скрипт 03/04):
#   z-scores транскриптомных сигнатур (DDR, SASP, Inflammaging, RhoGTPase,
#   cGAS–STING, Mito_Stress) по всем генотипам и условиям.
#   Главный вопрос: насколько PPM1D/Wip1 и PPM1B уникально или совместно
#   регулируют эти программы? Обнаруживается ли эффект в обоих видах?
#
# ТРЕБУЕТ в окружении / файловой системе:
#   output/rds/scores_all_v3.rds      — tidy tibble скоров (оба вида)
#   output/rds/vst_matrices_v3.rds    — list: $mouse и $human VST-матрицы
#   output/rds/deseq2_all_results.rds — list: $mouse и $human DE-результаты
#   output/rds/coldata_v3.rds         — метаданные образцов (оба вида)
#   output/rds/sig_list.rds           — именованный список генсетов
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)      # VIZ-1: ridgeline plot
  library(fmsb)          # VIZ-3: radar chart
  library(UpSetR)        # VIZ-4: upset plot
  library(ggalluvial)    # VIZ-5: alluvial/sankey
  library(ggpubr)        # VIZ-6: stat_compare_means
  library(patchwork)     # компоновка панелей
  library(ggrepel)       # подписи без перекрытия
  library(writexl)       # экспорт Excel
  library(ComplexHeatmap)
  library(circlize)
})

# ── 0. ЗАГРУЗКА ДАННЫХ ────────────────────────────────────────────────────────
scores_all <- readRDS("output/rds/scores_all_v3.rds")
vst_list   <- readRDS("output/rds/vst_matrices_v3.rds")
pipeline   <- readRDS("output/rds/deseq2_all_results.rds")
coldata    <- readRDS("output/rds/coldata_v3.rds")
sig_list   <- readRDS("output/rds/sig_list.rds")
gsea_combined <- readRDS("output/rds/gsea_hallmarks_all_v3.rds")

# Разделяем VST-матрицы по виду
vst_mouse <- vst_list$mouse   # гены × образцы (мышь, MC38)
vst_human <- vst_list$human   # гены × образцы (человек, HT29)

# Проверка структуры scores_all:
#   Ожидаем колонки: species, sample_id, signature, score, genotype, tnf
stopifnot(all(c("species","sample_id","signature","score",
                "genotype","tnf") %in% colnames(scores_all)))

cat("Виды:", unique(scores_all$species), "\n")
cat("Генотипы мышь:", unique(scores_all$genotype[scores_all$species=="mouse"]), "\n")
cat("Генотипы человек:", unique(scores_all$genotype[scores_all$species=="human"]), "\n")
cat("Сигнатуры:", unique(scores_all$signature), "\n")

# ── 0A. КОНФИГУРАЦИЯ ПУТЕЙ И ДИРЕКТОРИЙ ──────────────────────────────────────
for (d in c("output/figures/senescence",
            "output/figures/gsea",
            "output/figures/cross_species",
            "output/tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ── 0B. ЕДИНАЯ ЦВЕТОВАЯ ПАЛИТРА ПРОЕКТА ─────────────────────────────────────
# БИОЛОГИЧЕСКИЙ СМЫСЛ ПАЛИТРЫ:
#   Синий = WT (контроль, "здоровье")
#   Красный = Wip1KO (главный генотип, PPM1D/Wip1 — регулятор p53/ATM)
#   Зелёный = Ppm1bKO (PPM1B — TAK1/NF-κB ось, отдельная от Wip1)
#   Фиолетовый = DoubleKO (синергия обеих фосфатаз)
pal_geno_mouse <- c(
  "WT"       = "#4393C3",   # синий — контроль
  "Wip1KO"   = "#D6604D",   # красный — PPM1D/Wip1 KO (главный генотип)
  "Ppm1bKO"  = "#74C476",   # зелёный — PPM1B KO
  "DoubleKO" = "#9970AB"    # фиолетовый — двойной KO
)

# Человек (только WT, Wip1KO, Ppm1bKO — без Double)
pal_geno_human <- c(
  "WT"      = "#4393C3",
  "Wip1KO"  = "#D6604D",
  "Ppm1bKO" = "#74C476"
)

pal_tnf <- c(
  "basal"  = "white",
  "TNF12h" = "#FC8D59"
)

pal_species <- c(
  "mouse" = "#2B6CB0",   # тёмно-синий
  "human" = "#C05621"    # тёмно-оранжевый
)

# Единая тема для всех публикационных рисунков
theme_pub <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# VIZ-1: RIDGELINE PLOT — распределение z-scores по генотипам
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Ridgeline показывает не просто среднее, а ФОРМУ распределения.
#   При малом n (3 образца на группу) violin/boxplot нединформативны,
#   а ridgeline честно показывает разброс и положение каждой точки.
#
#   ОЖИДАНИЯ:
#   - Wip1KO: правое смещение для DDR, SASP, Inflammaging
#     (PPM1D/Wip1 — прямой дефосфорилатор ATM, CHK1/2, p53 и p38MAPK)
#   - Ppm1bKO: умеренный сдвиг, преимущественно NF-κB/Inflammaging
#     (PPM1B дефосфорилирует TAK1 — апстрим IKK-комплекса)
#   - DoubleKO: максимальный сдвиг, особенно Inflammaging + RhoActivity
#     (синергия двух нокаутов при отсутствии двух тормозов одновременно)
#   - TNF (пунктир vs сплошная): ожидаем большее смещение у KO,
#     что означает «гиперответ» на TNF — потеря фосфатазного ограничения
# =============================================================================

# Ключевые сигнатуры для основного рисунка
ridge_sigs_main <- c(
  "DDR",          # DNA damage response — p53/ATM ось
  "SASP",         # Secretory phenotype — IL-6, IL-1β, MMP
  "Inflammaging", # NF-κB / хроническое воспаление
  "RhoActivity",  # Rho GTPase — актиновый цитоскелет, миграция
  "cGAS_STING",   # Цитозольная ДНК → интерфероновый ответ
  "Mito_Stress"   # Митохондриальный стресс / митофагия
)

# Проверяем наличие сигнатур в данных
ridge_sigs_avail <- intersect(ridge_sigs_main, unique(scores_all$signature))
cat("Доступно для ridgeline:", ridge_sigs_avail, "\n")

# ── VIZ-1A: Мышь — ridgeline по генотипам ─────────────────────────────────────
p_ridge_mouse <- scores_all %>%
  dplyr::filter(
    species   == "mouse",
    signature %in% ridge_sigs_avail
  ) %>%
  dplyr::mutate(
    genotype  = factor(genotype,
                       levels = c("WT","Wip1KO","Ppm1bKO","DoubleKO")),
    tnf       = factor(tnf, levels = c("basal","TNF12h")),
    signature = factor(signature, levels = ridge_sigs_avail)
  ) %>%
  ggplot(aes(
    x     = score,
    y     = genotype,
    fill  = genotype,
    color = genotype
  )) +
  # Ridgeline с linetype = TNF-условие
  # Сплошная = basal, пунктирная = TNF12h
  ggridges::geom_density_ridges(
    aes(linetype = tnf),
    alpha           = 0.50,
    scale           = 1.15,
    rel_min_height  = 0.01,
    linewidth       = 0.65,
    jittered_points = TRUE,   # точки = отдельные образцы (n=2-3)
    point_size      = 2.5,
    point_alpha     = 0.95
  ) +
  # Вертикальная линия z=0 — граница "активации"
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "grey40",
    linewidth  = 0.45
  ) +
  scale_fill_manual(values = pal_geno_mouse, name = "Genotype") +
  scale_color_manual(values = pal_geno_mouse, guide = "none") +
  scale_linetype_manual(
    values = c("basal" = 1, "TNF12h" = 2),
    name   = "Stimulation"
  ) +
  facet_wrap(~signature, ncol = 3, scales = "free_x") +
  labs(
    title    = "Senescence & inflammatory signature scores — Mouse MC38",
    subtitle = "Solid = Basal | Dashed = TNF 12h | Points = individual samples",
    x        = "Signature z-score",
    y        = NULL
  ) +
  theme_pub

ggsave("output/figures/senescence/viz1_ridgeline_mouse.pdf",
       p_ridge_mouse, width = 14, height = 10)
ggsave("output/figures/senescence/viz1_ridgeline_mouse.png",
       p_ridge_mouse, width = 14, height = 10, dpi = 300)
message("✓ VIZ-1A ridgeline (мышь) сохранён")

# ── VIZ-1B: Человек — ridgeline (только basal, нет TNF-образцов) ─────────────
# ПРИМЕЧАНИЕ: у HT29 нет TNF-стимуляции → только basal
# DoubleKO у человека также отсутствует → 3 генотипа

p_ridge_human <- scores_all %>%
  dplyr::filter(
    species   == "human",
    tnf       == "basal",
    signature %in% ridge_sigs_avail
  ) %>%
  dplyr::mutate(
    genotype  = factor(genotype,
                       levels = c("WT","Wip1KO","Ppm1bKO")),
    signature = factor(signature, levels = ridge_sigs_avail)
  ) %>%
  ggplot(aes(
    x     = score,
    y     = genotype,
    fill  = genotype,
    color = genotype
  )) +
  ggridges::geom_density_ridges(
    alpha           = 0.50,
    scale           = 1.1,
    rel_min_height  = 0.01,
    linewidth       = 0.65,
    jittered_points = TRUE,
    point_size      = 2.5,
    point_alpha     = 0.95
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.45) +
  scale_fill_manual(values = pal_geno_human, name = "Genotype") +
  scale_color_manual(values = pal_geno_human, guide = "none") +
  facet_wrap(~signature, ncol = 3, scales = "free_x") +
  labs(
    title    = "Senescence & inflammatory signature scores — Human HT29",
    subtitle = "Basal condition only | Points = individual samples (n=3 per genotype)",
    x        = "Signature z-score",
    y        = NULL
  ) +
  theme_pub

ggsave("output/figures/senescence/viz1_ridgeline_human.pdf",
       p_ridge_human, width = 14, height = 8)
ggsave("output/figures/senescence/viz1_ridgeline_human.png",
       p_ridge_human, width = 14, height = 8, dpi = 300)
message("✓ VIZ-1B ridgeline (человек) сохранён")

# ── VIZ-1C: Violin + boxplot — МЫШЬ (basal), публикационный вариант ─────────
# КОГДА ИСПОЛЬЗОВАТЬ VIOLIN ВМЕСТО RIDGELINE:
#   Violin лучше при n≥5, ridgeline — при n=2-3.
#   Здесь violin + jitter (при n=3 он честнее violin-формы),
#   facet: signature (строки) × TNF (столбцы)

p_violin_mouse <- scores_all %>%
  dplyr::filter(
    species   == "mouse",
    signature %in% ridge_sigs_avail
  ) %>%
  dplyr::mutate(
    genotype  = factor(genotype,
                       levels = c("WT","Wip1KO","Ppm1bKO","DoubleKO")),
    tnf       = factor(tnf,
                       levels = c("basal","TNF12h"),
                       labels = c("Basal","TNF 12h")),
    signature = factor(signature, levels = ridge_sigs_avail)
  ) %>%
  ggplot(aes(x = genotype, y = score, fill = genotype)) +

  # Violin — форма распределения (при n=3 чисто информационная)
  geom_violin(
    alpha     = 0.40,
    scale     = "width",
    trim      = FALSE,
    linewidth = 0.45,
    color     = "grey35"
  ) +
  # Медиана + IQR
  geom_boxplot(
    width         = 0.14,
    outlier.shape = NA,
    fill          = "white",
    alpha         = 0.85,
    linewidth     = 0.4
  ) +
  # Отдельные точки — при n=3 это главный информационный элемент
  geom_jitter(
    aes(color = genotype),
    width = 0.10,
    size  = 2.8,
    alpha = 0.95
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey45", linewidth = 0.35) +

  scale_fill_manual(values  = pal_geno_mouse, name = "Genotype") +
  scale_color_manual(values = pal_geno_mouse, guide = "none") +

  # facet: signature (строки) × TNF-условие (столбцы)
  # DoubleKO появляется только в Basal — честно
  facet_grid(signature ~ tnf, scales = "free_y") +

  labs(
    title    = "Transcriptional signature scores — Mouse MC38 (all genotypes)",
    subtitle = "DoubleKO: basal only (n=2) | Jitter = individual samples",
    x        = NULL,
    y        = "Signature z-score"
  ) +
  theme_pub +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 8),
    strip.text.y = element_text(size = 8)
  )

ggsave("output/figures/senescence/viz1c_violin_mouse.pdf",
       p_violin_mouse, width = 11, height = 16)
ggsave("output/figures/senescence/viz1c_violin_mouse.png",
       p_violin_mouse, width = 11, height = 16, dpi = 300)
message("✓ VIZ-1C violin (мышь, все контрасты) сохранён")

# =============================================================================
# VIZ-2: BUBBLE CHART — GSEA enrichment (мышь + человек)
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Bubble chart одновременно кодирует три измерения:
#   X = генотип/контраст, Y = биологическая программа,
#   Размер = количество генов в пути (gene set size),
#   Цвет = NES (direction + magnitude), Прозрачность = –log10(FDR).
#
#   Это важнее, чем простой barplot NES, потому что:
#   - Большой размер bubble при слабом NES = путь статистически обоснован,
#     но биологически мал по величине изменений
#   - Маленький размер при высоком NES = специфический, узкий путь
#     (может быть артефактом малого генсета — осторожно!)
# =============================================================================

# Вспомогательная функция bubble chart для произвольного набора контрастов
make_bubble_gsea <- function(gsea_df, contrasts_vec, title_str, n_top = 12) {

  # Форматируем короткие имена путей
  df_plot <- gsea_df %>%
    dplyr::filter(contrast %in% contrasts_vec, padj < 0.05) %>%
    dplyr::mutate(
      pathway_label = gsub("HALLMARK_", "", pathway) %>%
        gsub("_", " ", .) %>%
        stringr::str_to_title() %>%
        stringr::str_trunc(40),
      neg_log_padj  = -log10(padj + 1e-10),
      contrast      = factor(contrast, levels = contrasts_vec)
    ) %>%
    dplyr::group_by(contrast) %>%
    dplyr::slice_max(abs(NES), n = n_top) %>%
    dplyr::ungroup()

  if (nrow(df_plot) == 0) {
    message("  ! Нет значимых путей для: ", title_str)
    return(invisible(NULL))
  }

  # Порядок путей по среднему NES
  path_order <- df_plot %>%
    dplyr::group_by(pathway_label) %>%
    dplyr::summarise(mean_NES = mean(NES)) %>%
    dplyr::arrange(mean_NES) %>%
    dplyr::pull(pathway_label)

  df_plot <- df_plot %>%
    dplyr::mutate(
      pathway_label = factor(pathway_label, levels = path_order)
    )

  ggplot(df_plot,
         aes(x     = contrast,
             y     = pathway_label,
             size  = size,
             color = NES,
             alpha = neg_log_padj)) +
    geom_point() +
    # Значение NES внутри пузыря
    geom_text(
      aes(label = round(NES, 1)),
      size     = 2.4,
      color    = "white",
      fontface = "bold"
    ) +
    scale_size_continuous(
      name   = "Gene set size",
      range  = c(4, 16),
      breaks = c(50, 100, 200, 300)
    ) +
    scale_color_gradient2(
      low      = "#2166AC",    # синий = подавление
      mid      = "grey90",
      high     = "#D6302F",    # красный = активация
      midpoint = 0,
      name     = "NES",
      limits   = c(-3, 3),
      oob      = scales::squish
    ) +
    scale_alpha_continuous(
      name  = "-log10(FDR)",
      range = c(0.4, 1.0)
    ) +
    labs(
      title    = title_str,
      subtitle = "Bubble size = gene set size | Color = NES | Opacity = –log10(FDR)",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1,
                                  face = "bold", size = 10),
      axis.text.y  = element_text(size = 9),
      legend.position = "right",
      plot.title   = element_text(face = "bold")
    )
}

# ── VIZ-2A: Мышь, baseline ────────────────────────────────────────────────────
p_bub_mouse_base <- make_bubble_gsea(
  gsea_df       = dplyr::filter(gsea_combined, species == "mouse"),
  contrasts_vec = c("Wip1KO_vs_WT", "Ppm1bKO_vs_WT", "DoubleKO_vs_WT"),
  title_str     = "Hallmark GSEA — KO vs WT (Mouse MC38, basal)"
)
ggsave("output/figures/gsea/viz2a_bubble_mouse_baseline.pdf",
       p_bub_mouse_base, width = 10, height = 12)
ggsave("output/figures/gsea/viz2a_bubble_mouse_baseline.png",
       p_bub_mouse_base, width = 10, height = 12, dpi = 300)
message("✓ VIZ-2A bubble chart (мышь, baseline) сохранён")

# ── VIZ-2B: Мышь, TNF-контрасты ──────────────────────────────────────────────
p_bub_mouse_tnf <- make_bubble_gsea(
  gsea_df       = dplyr::filter(gsea_combined, species == "mouse"),
  contrasts_vec = c("WT_TNF_vs_WT",
                    "Wip1KO_TNF_vs_Wip1KO",
                    "Ppm1bKO_TNF_vs_Ppm1bKO"),
  title_str     = "Hallmark GSEA — TNF stimulation contrasts (Mouse MC38)"
)
ggsave("output/figures/gsea/viz2b_bubble_mouse_tnf.pdf",
       p_bub_mouse_tnf, width = 10, height = 12)
ggsave("output/figures/gsea/viz2b_bubble_mouse_tnf.png",
       p_bub_mouse_tnf, width = 10, height = 12, dpi = 300)
message("✓ VIZ-2B bubble chart (мышь, TNF) сохранён")

# ── VIZ-2C: Человек, baseline ─────────────────────────────────────────────────
p_bub_human_base <- make_bubble_gsea(
  gsea_df       = dplyr::filter(gsea_combined, species == "human"),
  contrasts_vec = c("HT29_Wip1KO_vs_WT", "HT29_Ppm1bKO_vs_WT"),
  title_str     = "Hallmark GSEA — KO vs WT (Human HT29, basal)"
)
ggsave("output/figures/gsea/viz2c_bubble_human_baseline.pdf",
       p_bub_human_base, width = 9, height = 11)
ggsave("output/figures/gsea/viz2c_bubble_human_baseline.png",
       p_bub_human_base, width = 9, height = 11, dpi = 300)
message("✓ VIZ-2C bubble chart (человек, baseline) сохранён")

# =============================================================================
# VIZ-3: RADAR / SPIDER CHART — транскриптомный «профиль» каждого генотипа
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Радар позволяет одним взглядом оценить, какой генотип сдвинут по каким
#   биологическим осям одновременно. Ключевые паттерны:
#
#   - «Wip1KO shape»: большая площадь по DDR + SASP + Inflammaging
#     → PPM1D/Wip1 одновременно тормозит p53/ATM И NF-κB (через IκB)
#   - «Ppm1bKO shape»: преимущественно Inflammaging + RhoActivity
#     → PPM1B/PP2Cβ тормозит TAK1→IKK и, возможно, PAK-сигналинг
#   - «DoubleKO shape»: максимальная площадь по всем осям
#     → аддитивный или синергетический эффект
#
#   Нормализация в [0,1] по каждой сигнатуре необходима потому, что
#   абсолютные значения z-scores несопоставимы между сигнатурами разного
#   размера и компоновки.
# =============================================================================

make_radar_chart <- function(scores_df, species_tag, pal_geno,
                             geno_levels, outfile_prefix) {

  radar_data <- scores_df %>%
    dplyr::filter(species == species_tag, tnf == "basal") %>%
    dplyr::group_by(genotype, signature) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = signature, values_from = mean_score)

  # Приводим порядок генотипов
  geno_present <- intersect(geno_levels, radar_data$genotype)
  radar_data   <- radar_data %>%
    dplyr::filter(genotype %in% geno_present)

  radar_mat <- radar_data %>%
    tibble::column_to_rownames("genotype") %>%
    as.data.frame()

  # Нормализация [0,1] по каждой сигнатуре
  radar_norm <- as.data.frame(
    apply(radar_mat, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) return(rep(0.5, length(x)))
      (x - rng[1]) / diff(rng)
    })
  )
  rownames(radar_norm) <- rownames(radar_mat)

  # fmsb требует: строка 1 = максимум, строка 2 = минимум
  radar_fmsb <- rbind(
    rep(1, ncol(radar_norm)),
    rep(0, ncol(radar_norm)),
    radar_norm[geno_present[geno_present %in% rownames(radar_norm)], ]
  )
  rownames(radar_fmsb)[1:2] <- c("max","min")

  col_radar <- pal_geno[rownames(radar_fmsb)[-(1:2)]]

  # Функция отрисовки (вызывается для PDF и PNG)
  draw_radar_fn <- function() {
    par(mar = c(2, 2, 3, 2))
    fmsb::radarchart(
      radar_fmsb,
      axistype    = 1,
      pcol        = col_radar,
      pfcol       = adjustcolor(col_radar, alpha.f = 0.15),
      plwd        = 2.5,
      plty        = 1,
      cglcol      = "grey80",
      cglty       = 1,
      cglwd       = 0.5,
      axislabcol  = "grey40",
      vlcex       = 0.88,
      caxislabels = c("0","0.25","0.5","0.75","1"),
      title       = paste0("Signature profile — ", species_tag,
                           " (basal, normalized)")
    )
    legend("topright",
           legend = names(col_radar),
           col    = col_radar,
           lwd    = 2.5,
           bty    = "n",
           cex    = 0.9)
  }

  png(paste0(outfile_prefix, ".png"),
      width = 8, height = 7, units = "in", res = 300, bg = "white")
  draw_radar_fn()
  dev.off()

  pdf(paste0(outfile_prefix, ".pdf"), width = 8, height = 7)
  draw_radar_fn()
  dev.off()

  message("✓ Radar chart (", species_tag, ") сохранён")
}

# Мышь
make_radar_chart(
  scores_df    = scores_all,
  species_tag  = "mouse",
  pal_geno     = pal_geno_mouse,
  geno_levels  = c("WT","Wip1KO","Ppm1bKO","DoubleKO"),
  outfile_prefix = "output/figures/senescence/viz3_radar_mouse"
)

# Человек
make_radar_chart(
  scores_df    = scores_all,
  species_tag  = "human",
  pal_geno     = pal_geno_human,
  geno_levels  = c("WT","Wip1KO","Ppm1bKO"),
  outfile_prefix = "output/figures/senescence/viz3_radar_human"
)

# =============================================================================
# VIZ-4: UPSET PLOT — пересечения DEG между генотипами
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   UpSet — это лучший способ показать пересечения >3 множеств (Venn
#   нечитаем при n>3). Ключевые вопросы:
#
#   а) Какова «общая корзина» генов — значимых у всех KO?
#      Это гены, регулируемые обеими фосфатазами независимо.
#      Кандидаты: p53/NF-κB мишени (оба тормоза отключены).
#
#   б) Сколько генов уникально для DoubleKO?
#      Если много → синергетический эффект (не просто объединение).
#      Это ключевой аргумент для статьи о Double_KO.
#
#   в) Перекрытие мышь–человек для Wip1KO:
#      Консервированные гены = наиболее трансляционно значимые мишени.
# =============================================================================

# Вспомогательная функция извлечения DEG-сетов
extract_deg_set <- function(de_tbl,
                            padj_thr = 0.05,
                            lfc_thr  = 1) {
  de_tbl %>%
    dplyr::filter(
      !is.na(padj),
      padj < padj_thr,
      abs(log2FoldChange) >= lfc_thr
    ) %>%
    dplyr::pull(gene_id)
}

# ── VIZ-4A: Мышь — пересечение трёх baseline KO ──────────────────────────────
deg_sets_mouse <- list(
  "Wip1KO vs WT"   = extract_deg_set(pipeline$mouse$results_base[["Wip1KO_vs_WT"]]),
  "Ppm1bKO vs WT"  = extract_deg_set(pipeline$mouse$results_base[["Ppm1bKO_vs_WT"]]),
  "DoubleKO vs WT" = extract_deg_set(pipeline$mouse$results_base[["DoubleKO_vs_WT"]])
)

# Убираем пустые сеты
deg_sets_mouse <- Filter(function(x) length(x) > 0, deg_sets_mouse)
cat("DEG-сеты мышь:\n")
sapply(deg_sets_mouse, length)

draw_upset_mouse <- function() {
  UpSetR::upset(
    UpSetR::fromList(deg_sets_mouse),
    sets            = rev(names(deg_sets_mouse)),
    order.by        = "freq",
    mb.ratio        = c(0.60, 0.40),
    point.size      = 3.8,
    line.size       = 1.1,
    mainbar.y.label = "DEGs in intersection",
    sets.x.label    = "Total DEGs per contrast",
    text.scale      = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
    # Выделяем два биологически важных пересечения:
    queries = list(
      # Гены специфичные для DoubleKO = синергетические мишени
      list(
        query      = intersects,
        params     = list("DoubleKO vs WT"),
        color      = "#9970AB",
        active     = TRUE,
        query.name = "DoubleKO-unique (synergy)"
      ),
      # Гены во всех трёх KO = «ядро» фосфатаз-регулируемого транскриптома
      list(
        query      = intersects,
        params     = list("Wip1KO vs WT",
                          "Ppm1bKO vs WT",
                          "DoubleKO vs WT"),
        color      = "#D73027",
        active     = TRUE,
        query.name = "Core (all 3 KO)"
      )
    )
  )
}

png("output/figures/senescence/viz4a_upset_mouse_baseline.png",
    width = 10, height = 6, units = "in", res = 300, bg = "white")
draw_upset_mouse()
dev.off()

pdf("output/figures/senescence/viz4a_upset_mouse_baseline.pdf",
    width = 10, height = 6)
draw_upset_mouse()
dev.off()
message("✓ VIZ-4A UpSet (мышь, baseline) сохранён")

# ── VIZ-4B: Человек — пересечение двух baseline KO ───────────────────────────
deg_sets_human <- list(
  "HT29 Wip1KO vs WT"  = extract_deg_set(
    pipeline$human$results_base[["HT29_Wip1KO_vs_WT"]]),
  "HT29 Ppm1bKO vs WT" = extract_deg_set(
    pipeline$human$results_base[["HT29_Ppm1bKO_vs_WT"]])
)
deg_sets_human <- Filter(function(x) length(x) > 0, deg_sets_human)

if (length(deg_sets_human) >= 2) {
  png("output/figures/senescence/viz4b_upset_human_baseline.png",
      width = 8, height = 5, units = "in", res = 300, bg = "white")
  UpSetR::upset(
    UpSetR::fromList(deg_sets_human),
    sets       = rev(names(deg_sets_human)),
    order.by   = "freq",
    point.size = 3.8,
    line.size  = 1.1,
    text.scale = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1)
  )
  dev.off()
  message("✓ VIZ-4B UpSet (человек, baseline) сохранён")
}

# ── VIZ-4C: Кросс-видовое — Wip1KO мышь ∩ человек ───────────────────────────
# БИОЛОГИЧЕСКИЙ СМЫСЛ: пересечение мышь-человек = эволюционно консервированные
# мишени Wip1/PPM1D. Это наиболее важные кандидаты для трансляционной медицины.
#
# ВАЖНО: нужен ортологовый мэппинг мышь → человек (HGNC символы)
# Используем biomaRt или готовую таблицу ортологов

# Вариант без biomaRt: переводим мышиные символы в верхний регистр
# (грубое, но работающее приближение для ~70% генов)
mouse_wip1ko_deg <- toupper(
  extract_deg_set(pipeline$mouse$results_base[["Wip1KO_vs_WT"]])
)
human_wip1ko_deg <- extract_deg_set(
  pipeline$human$results_base[["HT29_Wip1KO_vs_WT"]]
)

cross_species_sets <- list(
  "Mouse Wip1KO (uppercase)" = mouse_wip1ko_deg,
  "Human Wip1KO"             = human_wip1ko_deg
)

if (all(sapply(cross_species_sets, length) > 0)) {
  conserved_genes <- intersect(mouse_wip1ko_deg, human_wip1ko_deg)
  message("Консервированных DEG Wip1KO (мышь ∩ человек, приближение): ",
          length(conserved_genes))
  write.csv(
    data.frame(gene = conserved_genes),
    "output/tables/conserved_DEG_Wip1KO_mouse_human.csv",
    row.names = FALSE
  )

  png("output/figures/cross_species/viz4c_upset_crossspecies_Wip1KO.png",
      width = 7, height = 5, units = "in", res = 300, bg = "white")
  UpSetR::upset(
    UpSetR::fromList(cross_species_sets),
    sets       = rev(names(cross_species_sets)),
    order.by   = "freq",
    point.size = 4,
    line.size  = 1.2,
    text.scale = c(1.5, 1.3, 1.3, 1.1, 1.4, 1.2),
    main.bar.color = "#D6604D",
    sets.bar.color = c("#C05621","#2B6CB0"),
    mainbar.y.label = "Shared DEGs",
    sets.x.label    = "Total Wip1KO DEGs"
  )
  dev.off()
  message("✓ VIZ-4C UpSet кросс-видовой Wip1KO сохранён")
}

# =============================================================================
# VIZ-5: ALLUVIAL / SANKEY — поток генов через генотипы
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Alluvial показывает «судьбу» каждого гена по мере нарастания нокаута:
#   WT → Wip1KO → Ppm1bKO → DoubleKO.
#
#   ОЖИДАЕМЫЕ ПАТТЕРНЫ:
#   1. «Консистентные UP» — гены, которые повышены у всех KO.
#      Это прямые мишени фосфатаз (SASP-цитокины, p53-мишени).
#   2. «DoubleKO-специфические» — гены, появляющиеся только при DoubleKO.
#      Это синергетически регулируемые гены (оба тормоза нужны одновременно).
#   3. «Перекрёстные» — UP у Wip1KO, но DOWN у Ppm1bKO.
#      Сигнализирует об антагонизме двух фосфатаз в отношении конкретного пути.
#
#   ЦВЕТ ПОТОКА = pathway-принадлежность гена.
#   Это позволяет видеть, как именно SASP-гены (красный) vs Rho-гены (фиолетовый)
#   ведут себя при переходе между генотипами.
# =============================================================================

# Классификация: UP / DOWN / NS для каждого гена в каждом контрасте
classify_deg_col <- function(de_tbl,
                             padj_thr = 0.05,
                             lfc_thr  = 1) {
  dplyr::case_when(
    is.na(de_tbl$padj)                           ~ "NS",
    de_tbl$padj >= padj_thr                      ~ "NS",
    de_tbl$log2FoldChange >= lfc_thr             ~ "UP",
    de_tbl$log2FoldChange <= -lfc_thr            ~ "DOWN",
    TRUE                                          ~ "NS"
  )
}

# Берём только мышь, baseline-контрасты с тремя генотипами
df_W  <- pipeline$mouse$results_base[["Wip1KO_vs_WT"]]
df_P  <- pipeline$mouse$results_base[["Ppm1bKO_vs_WT"]]
df_D  <- pipeline$mouse$results_base[["DoubleKO_vs_WT"]]

# Проверяем совпадение порядка генов
stopifnot(identical(df_W$gene_id, df_P$gene_id),
          identical(df_W$gene_id, df_D$gene_id))

# Сборка alluvial-таблицы с annotation к pathway
alluvial_genes <- tibble::tibble(
  gene_sym = df_W$gene_id,
  Wip1KO   = classify_deg_col(df_W),
  Ppm1bKO  = classify_deg_col(df_P),
  DoubleKO = classify_deg_col(df_D)
) %>%
  dplyr::filter(Wip1KO != "NS" | Ppm1bKO != "NS" | DoubleKO != "NS") %>%
  dplyr::mutate(
    # Приоритет: DDR > SASP > Inflammaging > Rho > cGAS_STING > Other
    pathway = dplyr::case_when(
      gene_sym %in% sig_list[["DDR"]]          ~ "DDR",
      gene_sym %in% sig_list[["SASP"]]         ~ "SASP",
      gene_sym %in% sig_list[["Inflammaging"]] ~ "Inflammaging",
      gene_sym %in% sig_list[["RhoActivity"]]  ~ "Rho signaling",
      gene_sym %in% sig_list[["cGAS_STING"]]   ~ "cGAS–STING",
      TRUE                                      ~ "Other DEG"
    )
  )

cat("Alluvial: всего генов =", nrow(alluvial_genes), "\n")
print(table(alluvial_genes$pathway))

# Агрегация по категориям
alluvial_counts <- alluvial_genes %>%
  dplyr::count(Wip1KO, Ppm1bKO, DoubleKO, pathway) %>%
  dplyr::filter(n >= 3) %>%   # убираем единичные случаи (шум)
  dplyr::mutate(
    dplyr::across(
      c(Wip1KO, Ppm1bKO, DoubleKO),
      ~ factor(.x, levels = c("UP", "DOWN", "NS"))
    )
  )

# Цвета pathway-категорий
path_colors <- c(
  "DDR"           = "#F781BF",   # розовый — p53/ATM
  "SASP"          = "#E41A1C",   # красный — секреторный фенотип
  "Inflammaging"  = "#FF7F00",   # оранжевый — NF-κB/хроническое воспаление
  "Rho signaling" = "#984EA3",   # фиолетовый — актиновый цитоскелет
  "cGAS–STING"    = "#4DAF4A",   # зелёный — интерфероновый ответ
  "Other DEG"     = "#BBBBBB"    # серый — прочее
)

p_alluvial <- ggplot(
  alluvial_counts,
  aes(
    axis1 = Wip1KO,
    axis2 = Ppm1bKO,
    axis3 = DoubleKO,
    y     = n,
    fill  = pathway
  )
) +
  ggalluvial::geom_alluvium(
    alpha    = 0.72,
    width    = 1/4,
    knot.pos = 0.4
  ) +
  ggalluvial::geom_stratum(
    width     = 1/4,
    fill      = "grey92",
    color     = "grey35",
    linewidth = 0.4
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size     = 3.8,
    fontface = "bold",
    color    = "grey15"
  ) +
  scale_x_discrete(
    limits = c("Wip1KO", "Ppm1bKO", "DoubleKO"),
    labels = c("Wip1KO\nvs WT\n(PPM1D KO)",
               "Ppm1bKO\nvs WT\n(PPM1B KO)",
               "DoubleKO\nvs WT"),
    expand = c(0.12, 0.12)
  ) +
  scale_fill_manual(values = path_colors, name = "Pathway / signature") +
  labs(
    title    = "Gene fate across phosphatase knockouts — Mouse MC38",
    subtitle = "Flow of DEGs (UP / DOWN / NS) | Colour = pathway membership",
    y = "Number of genes",
    x = NULL
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 11, face = "bold"))

ggsave("output/figures/senescence/viz5_alluvial_mouse.png",
       p_alluvial, width = 11, height = 7, dpi = 300, bg = "white")
ggsave("output/figures/senescence/viz5_alluvial_mouse.pdf",
       p_alluvial, width = 11, height = 7)
message("✓ VIZ-5 Alluvial (мышь) сохранён")

# =============================================================================
# VIZ-6: VIOLIN + JITTER — экспрессия ключевых генов-маркёров
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Показываем нормализованную (VST) экспрессию ключевых «нарративных» генов —
#   тех, через которые рассказывается история статьи.
#
#   ВЫБОР ГЕНОВ ПО КАТЕГОРИЯМ:
#   - Il6, Tnf, Il1b, Ccl2: «ядро» SASP-воспаления; Il6 — главный SASP-цитокин
#     и NF-κB-мишень; оба PPM1D и PPM1B тормозят его активацию
#   - Trp53, Cdkn1a(p21): маркёры p53-активации / остановки клеточного цикла
#   - Cdkn2a(p16): маркёр необратимой сенесценции (INK4a-ось)
#   - Rhoa, Rock1: маркёры Rho-активности; Rock1 фосфорилирует LIMK → стресс-волокна
#   - Cgas, Sting1: маркёры цитоплазматической ДНК-сенсинговой оси
#
#   Wilcoxon test vs WT: при n=3 это единственный корректный непараметрический
#   тест (t-test слишком чувствителен к форме распределения при малом n).
# =============================================================================

# Ключевые гены — одинаковые для обоих видов (в верхнем регистре для человека)
key_genes_mouse <- c(
  "Il6","Tnf","Il1b","Ccl2","Cxcl10",   # SASP / NF-κB
  "Trp53","Cdkn1a","Cdkn2a",             # p53 / сенесценция
  "Rhoa","Rock1","Pak1",                  # Rho-сигналинг
  "Cgas","Sting1"                         # cGAS–STING
)
key_genes_human <- c(
  "IL6","TNF","IL1B","CCL2","CXCL10",
  "TP53","CDKN1A","CDKN2A",
  "RHOA","ROCK1","PAK1",
  "CGAS","STING1"
)

# Функция violin-plot для одного вида
make_violin_genes <- function(vst_mat, key_genes_vec, coldata_df,
                              species_tag, pal_geno, geno_levels,
                              outfile_prefix) {

  genes_avail <- intersect(key_genes_vec, rownames(vst_mat))
  if (length(genes_avail) == 0) {
    message("! Нет ключевых генов в VST-матрице для ", species_tag)
    return(invisible(NULL))
  }

  # tidy-формат
  violin_df <- vst_mat[genes_avail, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      -gene,
      names_to  = "sample_id",
      values_to = "vst"
    ) %>%
    dplyr::left_join(
      coldata_df %>%
        dplyr::select(sample_id, genotype, tnf),
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      genotype = factor(genotype, levels = geno_levels),
      gene     = factor(gene, levels = genes_avail)
    )

  comparisons_list <- lapply(
    setdiff(geno_levels, "WT"),
    function(g) c("WT", g)
  )

  # ── Basal violin ────────────────────────────────────────────────────────────
  p_basal <- violin_df %>%
    dplyr::filter(tnf == "basal") %>%
    ggplot(aes(x = genotype, y = vst, fill = genotype, color = genotype)) +
    geom_violin(alpha = 0.40, trim = FALSE, scale = "width",
                linewidth = 0.45, color = "grey30") +
    geom_boxplot(width = 0.14, outlier.shape = NA,
                 fill = "white", alpha = 0.85, linewidth = 0.4) +
    geom_jitter(width = 0.10, size = 2.8, alpha = 0.95) +
    ggpubr::stat_compare_means(
      comparisons = comparisons_list,
      method      = "wilcox.test",
      label       = "p.signif",
      size        = 3,
      bracket.size = 0.4
    ) +
    scale_fill_manual(values  = pal_geno) +
    scale_color_manual(values = pal_geno) +
    facet_wrap(~gene, scales = "free_y", nrow = 2) +
    labs(
      title    = paste0("VST-normalized expression — ",
                        species_tag, " (basal)"),
      subtitle = "Wilcoxon vs WT: * p<0.05, ** p<0.01, *** p<0.001",
      x = NULL, y = "VST expression",
      fill = "Genotype", color = "Genotype"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          legend.position = "none")

  ggsave(paste0(outfile_prefix, "_basal.pdf"),
         p_basal, width = 16, height = 7)
  ggsave(paste0(outfile_prefix, "_basal.png"),
         p_basal, width = 16, height = 7, dpi = 300)

  message("✓ VIZ-6 violin (", species_tag, ", basal) сохранён")
}

# Мышь
make_violin_genes(
  vst_mat     = vst_mouse,
  key_genes_vec = key_genes_mouse,
  coldata_df  = dplyr::filter(coldata, species == "mouse"),
  species_tag = "Mouse MC38",
  pal_geno    = pal_geno_mouse,
  geno_levels = c("WT","Wip1KO","Ppm1bKO","DoubleKO"),
  outfile_prefix = "output/figures/senescence/viz6_violin_mouse"
)

# Человек
make_violin_genes(
  vst_mat     = vst_human,
  key_genes_vec = key_genes_human,
  coldata_df  = dplyr::filter(coldata, species == "human"),
  species_tag = "Human HT29",
  pal_geno    = pal_geno_human,
  geno_levels = c("WT","Wip1KO","Ppm1bKO"),
  outfile_prefix = "output/figures/senescence/viz6_violin_human"
)

# ── VIZ-6B: TNF-эффект — мышь, цитокины ──────────────────────────────────────
# БИОЛОГИЧЕСКИЙ СМЫСЛ:
#   Показываем, как TNF изменяет экспрессию цитокинов в каждом генотипе.
#   PPM1D/Wip1 дефосфорилирует IκBα → без него NF-κB активируется сильнее.
#   ОЖИДАНИЕ: Wip1KO+TNF > WT+TNF для Il6, Tnf, Ccl2 (hyperactivation).

sasp_genes_tnf <- intersect(
  c("Il6","Tnf","Il1b","Cxcl10","Ccl2"),
  rownames(vst_mouse)
)

if (length(sasp_genes_tnf) > 0) {

  violin_tnf_df <- vst_mouse[sasp_genes_tnf, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene,
                        names_to  = "sample_id",
                        values_to = "vst") %>%
    dplyr::left_join(
      dplyr::filter(coldata, species == "mouse") %>%
        dplyr::select(sample_id, genotype, tnf),
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      genotype = factor(genotype,
                        levels = c("WT","Wip1KO","Ppm1bKO","DoubleKO")),
      tnf      = factor(tnf, levels = c("basal","TNF12h"),
                        labels = c("Basal","TNF 12h")),
      gene     = factor(gene, levels = sasp_genes_tnf)
    )

  p_tnf_effect <- violin_tnf_df %>%
    ggplot(aes(x = tnf, y = vst, fill = tnf)) +
    geom_violin(alpha = 0.50, trim = TRUE) +
    geom_boxplot(width = 0.15, alpha = 0.90,
                 outlier.shape = NA, color = "grey20") +
    geom_jitter(aes(color = genotype), width = 0.10, size = 2.2, alpha = 0.9) +
    # Соединяем basal–TNF для одного образца
    geom_line(aes(group = sample_id), alpha = 0.30, color = "grey55") +
    ggpubr::stat_compare_means(method = "wilcox.test",
                               label = "p.signif", size = 3) +
    scale_fill_manual(values  = pal_tnf) +
    scale_color_manual(values = pal_geno_mouse) +
    facet_grid(gene ~ genotype, scales = "free_y") +
    labs(
      title    = "TNF-induced SASP cytokine expression — Mouse MC38",
      subtitle = "Lines = basal→TNF12h within sample | Color = genotype",
      x = "Treatment", y = "VST expression",
      fill = "Treatment", color = "Genotype"
    ) +
    theme_pub

  ggsave("output/figures/senescence/viz6b_tnf_effect_mouse.pdf",
         p_tnf_effect, width = 14, height = 10)
  ggsave("output/figures/senescence/viz6b_tnf_effect_mouse.png",
         p_tnf_effect, width = 14, height = 10, dpi = 300)
  message("✓ VIZ-6B TNF-effect violin (мышь) сохранён")
}

# =============================================================================
# VIZ-7: КРОСС-ВИДОВОЙ HEATMAP — z-score профили мышь vs человек
#
# БИОЛОГИЧЕСКИЙ СМЫСЛ (НОВЫЙ РИСУНОК):
#   Главный вопрос второго приоритета: воспроизводится ли транскриптомный
#   фенотип Wip1KO у человека (HT29) так же, как у мыши (MC38)?
#
#   Если да → PPM1D/Wip1 регулирует универсальную, эволюционно консервированную
#   программу (p53/NF-κB), а не мышиноспецифический артефакт.
#
#   СТРУКТУРА HEATMAP:
#   Строки = транскриптомные сигнатуры (DDR, SASP, Inflammaging, Rho, cGAS, Mito)
#   Столбцы = образцы, сгруппированные по виду → генотипу
#   Значения = z-score каждого образца по каждой сигнатуре
#   Аннотации = вид (цвет панели) + генотип (цвет полосы)
# =============================================================================

# Собираем матрицу z-scores для всех образцов обоих видов (только basal)
cs_matrix_df <- scores_all %>%
  dplyr::filter(tnf == "basal") %>%
  dplyr::select(species, sample_id, genotype, signature, score) %>%
  tidyr::pivot_wider(names_from = signature, values_from = score) %>%
  dplyr::arrange(species, genotype)

# Матрица значений (образцы × сигнатуры)
sig_cols <- setdiff(names(cs_matrix_df),
                    c("species","sample_id","genotype"))
cs_mat <- cs_matrix_df %>%
  tibble::column_to_rownames("sample_id") %>%
  dplyr::select(all_of(sig_cols)) %>%
  as.matrix()

# Аннотация строк для ComplexHeatmap
row_anno_df <- cs_matrix_df %>%
  dplyr::select(sample_id, species, genotype) %>%
  tibble::column_to_rownames("sample_id")

species_colors <- c("mouse" = "#2B6CB0", "human" = "#C05621")
geno_colors    <- c("WT"      = "#4393C3",
                    "Wip1KO"  = "#D6604D",
                    "Ppm1bKO" = "#74C476",
                    "DoubleKO"= "#9970AB")

row_anno <- ComplexHeatmap::rowAnnotation(
  Species  = row_anno_df$species,
  Genotype = row_anno_df$genotype,
  col      = list(
    Species  = species_colors,
    Genotype = geno_colors
  ),
  annotation_legend_param = list(
    Species  = list(title = "Species"),
    Genotype = list(title = "Genotype")
  )
)

# Цветовая шкала для z-score
zscore_col_fun <- circlize::colorRamp2(
  breaks = c(-3, -1.5, 0, 1.5, 3),
  colors = c("#053061","#92C5DE","#F7F7F7","#F4A582","#67001F")
)

ht_cross <- ComplexHeatmap::Heatmap(
  cs_mat,
  name                   = "z-score",
  col                    = zscore_col_fun,
  right_annotation       = row_anno,

  # Кластеризация образцов внутри видов отдельно
  cluster_rows           = TRUE,
  cluster_columns        = TRUE,
  show_row_dend          = TRUE,
  show_column_dend       = TRUE,

  row_names_gp           = gpar(fontsize = 8),
  column_names_gp        = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot       = 30,

  # Подписи значений внутри ячеек
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- cs_mat[i, j]
    if (!is.na(v) && abs(v) >= 1.5)
      grid.text(round(v, 1), x, y,
                gp = gpar(fontsize = 7,
                          col = ifelse(abs(v) > 2.5, "white", "black")))
  },

  column_title    = "Transcriptional signature z-scores: Mouse MC38 vs Human HT29 (basal)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),

  width  = unit(12, "cm"),
  height = unit(max(8, nrow(cs_mat) * 0.55), "cm"),

  heatmap_legend_param = list(
    title         = "z-score",
    at            = c(-3,-2,-1,0,1,2,3),
    legend_height = unit(4, "cm")
  )
)

pdf("output/figures/cross_species/viz7_crossspecies_zscore_heatmap.pdf",
    width = 14, height = max(7, nrow(cs_mat) * 0.4 + 4))
ComplexHeatmap::draw(ht_cross)
dev.off()

png("output/figures/cross_species/viz7_crossspecies_zscore_heatmap.png",
    width = 14, height = max(7, nrow(cs_mat) * 0.4 + 4),
    units = "in", res = 300, bg = "white")
ComplexHeatmap::draw(ht_cross)
dev.off()
message("✓ VIZ-7 кросс-видовой heatmap z-scores сохранён")

# =============================================================================
# ЭКСПОРТ ТАБЛИЦ
# =============================================================================

# Функция форматирования GSEA-таблицы для публикации
format_gsea_table <- function(gsea_df, contrasts_vec) {
  gsea_df %>%
    dplyr::filter(contrast %in% contrasts_vec, padj < 0.05) %>%
    dplyr::mutate(
      pathway_clean = gsub("HALLMARK_","", pathway) %>%
        gsub("_"," ",.) %>%
        stringr::str_to_title(),
      NES     = round(NES, 3),
      ES      = round(ES, 3),
      pval    = signif(pval, 3),
      padj    = signif(padj, 3),
      direction = ifelse(NES > 0, "UP", "DOWN"),
      contrast  = factor(contrast, levels = contrasts_vec)
    ) %>%
    dplyr::arrange(contrast, desc(abs(NES))) %>%
    dplyr::select(
      Species   = species,
      Contrast  = contrast,
      Pathway   = pathway_clean,
      NES, ES,
      Direction = direction,
      `Gene set size` = size,
      pval, padj
    )
}

tbl_mouse_base <- format_gsea_table(
  gsea_combined,
  c("Wip1KO_vs_WT","Ppm1bKO_vs_WT","DoubleKO_vs_WT")
)
tbl_mouse_tnf  <- format_gsea_table(
  gsea_combined,
  c("WT_TNF_vs_WT","Wip1KO_TNF_vs_Wip1KO","Ppm1bKO_TNF_vs_Ppm1bKO")
)
tbl_human_base <- format_gsea_table(
  gsea_combined,
  c("HT29_Wip1KO_vs_WT","HT29_Ppm1bKO_vs_WT")
)

# Robust hits: пути, значимые в ≥2 контрастах
tbl_robust <- bind_rows(tbl_mouse_base, tbl_mouse_tnf, tbl_human_base) %>%
  dplyr::group_by(Pathway) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::mutate(
    n_contrasts    = dplyr::n(),
    mean_NES       = round(mean(NES), 3),
    consistent_dir = dplyr::n_distinct(Direction) == 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_contrasts), desc(abs(mean_NES)))

writexl::write_xlsx(
  list(
    "Mouse_baseline"    = as.data.frame(tbl_mouse_base),
    "Mouse_TNF"         = as.data.frame(tbl_mouse_tnf),
    "Human_baseline"    = as.data.frame(tbl_human_base),
    "Robust_hits_2plus" = as.data.frame(tbl_robust)
  ),
  path = "output/tables/GSEA_Hallmark_results_v3.xlsx"
)

write.csv(tbl_robust,
          "output/tables/Table_S_GSEA_robust_hits_v3.csv",
          row.names = FALSE)

message("✓ Таблицы экспортированы: GSEA_Hallmark_results_v3.xlsx (4 листа)")

# =============================================================================
# ИТОГ
# =============================================================================
message("\n══════════════════════════════════════════════════════════════════")
message("✓ 06_advanced_visualizations_v3.R завершён")
message("  VIZ-1 — Ridgeline: z-score распределения, мышь + человек")
message("  VIZ-2 — Bubble:    GSEA enrichment, мышь (basal+TNF) + человек")
message("  VIZ-3 — Radar:     транскриптомный профиль, мышь + человек")
message("  VIZ-4 — UpSet:     DEG-пересечения, мышь + человек + кросс-вид")
message("  VIZ-5 — Alluvial:  поток генов WT→Wip1KO→Ppm1bKO→DoubleKO (мышь)")
message("  VIZ-6 — Violin:    маркёрные гены, мышь + человек + TNF-эффект")
message("  VIZ-7 — Heatmap:   кросс-видовое z-score сравнение (НОВОЕ)")
message("══════════════════════════════════════════════════════════════════")
