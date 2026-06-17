# PPM1D/Wip1 Knockout Transcriptomics Project

## Описание
Bulk RNA-seq анализ нокаутных по PPM1D/Wip1 почек мышей
(C57BL/6, BD27 WT vs BD28/BD29 KO).
Включает биоинформатический пайплайн, математическое моделирование
Rho-GTPase переключателя и черновик рукописи.

## Структура
```
data/        — исходные данные (counts, metadata, genesets)
scripts/     — R-скрипты анализа и моделирования
results/     — таблицы и рисунки (генерируются, не в git)
manuscript/  — черновик статьи (md/docx разделы)
notebooks/   — Rmd для exploratory анализа
config/      — параметры запуска
tmp/         — временные файлы (не в git)
```

## Быстрый старт
```r
# 1. Установить зависимости
source('scripts/00_setup.R')

# 2. Нормализация
source('scripts/01_qc_normalization.R')

# 3. DESeq2
source('scripts/02_differential_expression.R')

# 4. Моделирование Rho-переключателя
source('scripts/modeling/run_calibration_pipeline.R')
source('scripts/modeling/rho_bifurcation.R')
results <- run_full_bifurcation_analysis()
```

## Ключевые выводы
- Wip1KO-почка: «DDR-primed, pre-injured» транскриптомное состояние
- Повышены: Trp53/p21/p16, SASP-цитокины, Rho-сеть
- Математическая модель: KO -> RhoA-dominant аттрактор
- Механизм цисплатин-летальности: наложение стресса на pre-injured фон

## Данные
- rawcount272829.csv — матрица сырых счётов (12 образцов)
- 4 WT (BD27) + 8 KO (BD28, BD29)

## Авторы
Bobkov, Bogdanova et al.

## Статус
🟡 In preparation — результаты получены, рукопись в работе
