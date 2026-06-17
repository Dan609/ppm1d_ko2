# =================================================
# color_palettes.R
# Единые цветовые палитры для всего проекта
# source('config/color_palettes.R') в начале каждого скрипта
# =================================================

COLORS <- list(

  # Генотипы
  genotype = c(
    WT = '#2471a3',
    KO = '#c0392b'
  ),

  # GTPases
  gtpase = c(
    RhoA  = '#e74c3c',
    Rac1  = '#2980b9',
    Cdc42 = '#27ae60'
  ),

  # Пути (GSEA / кастомные генсеты)
  pathways = c(
    DDR_Senescence  = '#8e44ad',
    Inflammation    = '#e67e22',
    Rho_GTPases     = '#e74c3c',
    ProxTubule      = '#16a085',
    Macrophage      = '#f39c12',
    Fibroblast      = '#7f8c8d',
    cGAS_STING      = '#2c3e50',
    PINK1_Mitostress = '#d35400'
  ),

  # Volcano plot
  volcano = c(
    up       = '#c0392b',
    down     = '#2471a3',
    ns       = '#bdc3c7'
  ),

  # Heatmap (diverging)
  heatmap_low  = '#2471a3',
  heatmap_mid  = '#ffffff',
  heatmap_high = '#c0392b',

  # Нейтральные
  grey_light = '#ecf0f1',
  grey_mid   = '#95a5a6',
  grey_dark  = '#2c3e50'
)

# Тема ggplot2 для всего проекта
theme_ppm1d <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = '#f4f4f4',
                                        color = 'grey80'),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.title        = element_text(face = 'bold',
                                        size = base_size + 1),
      plot.subtitle     = element_text(color = 'grey40',
                                        size = base_size - 1)
    )
}

cat('  [config] Палитры и тема загружены\n')
