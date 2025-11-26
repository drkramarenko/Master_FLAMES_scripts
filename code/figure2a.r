#!/usr/bin/env Rscript

# Input: ST3 Genetic correlations: DCM, HCM & cardiac MRI traits
# Output: Figure 2a

# Package versions used:
# readxl   1.4.5
# dplyr    1.1.4
# reshape2 1.4.4
# ggplot2  3.5.2

library(readxl)    # read_excel()
library(dplyr)     # %>%, mutate(), across()
library(reshape2)  # melt()
library(ggplot2)   # ggplot2 plotting

# --- Paths and input file ---
PROJECT_DIR <- "/path/to/project"
FIG_INPUT_DIR <- file.path(PROJECT_DIR, "figures/for_figures")
FIG_OUTPUT_DIR <- file.path(PROJECT_DIR, "figures/output")

HEATMAP_XLSX <- file.path(FIG_INPUT_DIR, "heatmap_main.xlsx")
HEATMAP_SHEET <- "DCM_HCM"

# --- Read and prepare data ---
gen_cor_DCM_HCM <- read_excel(HEATMAP_XLSX, sheet = HEATMAP_SHEET) %>%
  mutate(across(-rows, as.numeric))

# Long format for ggplot
long_data <- melt(
  gen_cor_DCM_HCM,
  id.vars = "rows",
  variable.name = "Variable",
  value.name = "Correlation"
)

# Ensure rows and columns use the same order
order_levels <- c("DCM", "Ecc", "LVESVi", "LVconc", "HCM")
long_data$rows     <- factor(long_data$rows,     levels = order_levels)
long_data$Variable <- factor(long_data$Variable, levels = order_levels)

# --- Create heatmap ---
DCM_HCM_mat <- ggplot(long_data, aes(x = Variable, y = rows, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low  = "#078eff",   # blue
    mid  = "#f2f7f3",   # neutral
    high = "#B71C1C",   # red
    midpoint = 0,
    limits   = c(-1, 1)

  ) +
  labs(title = "", x = "", y = "") +
  theme_minimal()

print(DCM_HCM_mat)

# --- Save outputs ---
dir.create(FIG_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "DCM_HCM_mat.png"),
  plot     = DCM_HCM_mat,
  width    = 6,
  height   = 4,
  dpi      = 450,
  bg       = "transparent"
)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "DCM_HCM_mat.pdf"),
  plot     = DCM_HCM_mat,
  width    = 6,
  height   = 4,
  device   = "pdf",
  bg       = "transparent"
)
