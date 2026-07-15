# new_lake_plotting.R
# Save all manuscript figures as high-resolution PNGs.

Packages <- c("cowplot", "ggplot2")
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) install.packages(missing_packages)
lapply(Packages, library, character.only = TRUE)

load("data/new_lake_analytics.RData")

dir.create("Newfigs", showWarnings = FALSE)

# Figure 1: Shannon diversity and LCBD (4-panel)
ggsave("Newfigs/Figure1.png",
       plot = plot_grid(fig1a, fig1b, fig1c, fig1d, nrow = 2),
       width = 25, height = 20, dpi = 600, units = "cm")

# Figure 2: Beta diversity components
ggsave("Newfigs/Figure2.png",
       plot = fig2,
       width = 15, height = 12, dpi = 600, units = "cm")

# Figure 3: Relative density change vs body mass
ggsave("Newfigs/Figure3.png",
       plot = fig3,
       width = 15, height = 12, dpi = 600, units = "cm")

# Figure 4: CWM by fish presence
ggsave("Newfigs/Figure4.png",
       plot = fig4,
       width = 12, height = 12, dpi = 600, units = "cm")

# Figure S1: elevation and site counts
ggsave("Newfigs/FigureS1.png",
       plot = plot_grid(fig_s1a, fig_s1b, nrow = 1),
       width = 25, height = 12, dpi = 600, units = "cm")

# Figure S2: species density + occupancy (2-panel)
ggsave("Newfigs/FigureS2.png",
       plot = plot_grid(fig_s2a, fig_s2b, nrow = 2),
       width = 30, height = 30, dpi = 600, units = "cm")

# Figure S3: body mass of absent species
ggsave("Newfigs/FigureS3.png",
       plot = fig_s3,
       width = 12, height = 12, dpi = 600, units = "cm")

# Figure S4: CWM vs elevation by fish presence
ggsave("Newfigs/FigureS4.png",
       plot = fig_s4,
       width = 25, height = 12, dpi = 600, units = "cm")

message("new_lake_plotting.R complete: all manuscript figures saved to Newfigs/")






