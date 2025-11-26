#!/usr/bin/env Rscript

# Input: ST4 Comparison of lead variant effects: DCM vs HCM
# Output: Figure 2b
# ggplot2 3.5.2
# dplyr   1.1.4   (if used for preprocessing)
# data.table / readr – optional depending on how loci_for_mr_selected is loaded

library(ggplot2)

# --- Input object ---
# loci_for_mr_selected must contain:
#   DCM_BETA, HCM_BETA, DCM_SE, HCM_SE, SNP_label
# Load your data here if needed:
# loci_for_mr_selected <- readRDS(file.path(PROJECT_DIR, "data/processed/loci_for_mr_selected.rds"))

# --- Plot ---
effect_concordance_plot <- ggplot(
  loci_for_mr_selected,
  aes(x = DCM_BETA, y = HCM_BETA, color = SNP_label)
) +
  # Regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray40",
    linetype = "dashed",
    linewidth = 0.5,
    fullrange = TRUE
  ) +
  # Points
  geom_point(alpha = 0.7, size = 3) +
  
  # Manual color scale
  scale_color_manual(
    values = c(
      "DCM"  = "#7030A0",
      "HCM"  = "#00B050",
      "Both" = "red4"
    )
  ) +
  
  # Error bars
  geom_errorbar(
    aes(ymin = HCM_BETA - HCM_SE, ymax = HCM_BETA + HCM_SE),
    width = 0,
    alpha  = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = DCM_BETA - DCM_SE, xmax = DCM_BETA + DCM_SE),
    height = 0,
    alpha  = 0.5
  ) +

  # Reference lines
  geom_vline(xintercept = 0, color = "gray60") +
  geom_hline(yintercept = 0, color = "gray60") +
  
  # Labels
  labs(
    x = "DCM Beta",
    y = "HCM Beta",
    title = "",
    color = "Locus label"
  ) +
  
  theme_classic() +
  xlim(-0.4, 0.4) +
  ylim(-0.4, 0.55) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title  = element_text(size = 10),
    axis.text   = element_text(size = 7)
  )

# Display plot
print(effect_concordance_plot)

# --- Save outputs ---
ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "figure2b_effect_concordance.png"),
  plot     = effect_concordance_plot,
  width    = 6,
  height   = 4,
  dpi      = 450,
  bg       = "transparent"
)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "figure2b_effect_concordance.pdf"),
  plot     = effect_concordance_plot,
  width    = 6,
  height   = 4,
  device   = "pdf",
  bg       = "transparent"
)
