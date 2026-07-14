# new_lake_plotting.R
# Load analytics objects and save manuscript figure PNGs with explicit names.

Packages <- c(
  "cowplot", "tidyverse"
)
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) install.packages(missing_packages)
lapply(Packages, library, character.only = TRUE)

load("data/new_lake_analytics.RData")
# Manuscript base theme (match legacy Plotting.R aesthetics)
theme_base <- theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())
ggsave("fig2a.png", plot = fig2a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig2b.png", plot = fig2b, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig2c.png", plot = fig2c, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig2d.png", plot = fig2d + theme_base, width = 12, height = 10, dpi = 600, units = "cm")

ggsave("fig3a.png", plot = fig3a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig3b.png", plot = fig3b, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig3c.png", plot = fig3c, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig3d.png", plot = fig3d + theme_base, width = 12, height = 10, dpi = 600, units = "cm")

ggsave("supp1.png", plot = plot_grid(supp.a, supp.b, nrow = 2), width = 25, height = 15, dpi = 600, units = "cm")

# Save combined figures for convenience
ggsave("fig2_panel.png", plot = plot_grid(fig2a, fig2b, fig2c, fig2d, nrow = 2), width = 25, height = 20, dpi = 600, units = "cm")
ggsave("fig3_panel.png", plot = plot_grid(fig3a, fig3b, fig3c, fig3d, nrow = 2), width = 25, height = 20, dpi = 600, units = "cm")

message("new_lake_plotting.R complete: figures saved with explicit manuscript filenames")
