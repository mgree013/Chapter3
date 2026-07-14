# Title: The influence of non-native fish on stream macroinvertebrate and lake zooplankton communities along elevational gradients
# Authors: Matthew Green, David Herbst, and Kurt Anderson

# plotting Code

library(cowplot)
library(tidyverse)

# This script assumes plot objects are already defined in the environment.
# Source code/Lake.analysis.R first before running this file.

###########################################################################
# Figure 1
fig1_grid <- plot_grid(new.fig1a, fig1a, new.prop.a, new.prop.a2, nrow = 2)
ggsave("fig1.png", plot = fig1_grid, width = 24, height = 20, dpi = 600, units = "cm")

ggsave("fig1_new_a.png", plot = new.fig1a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig1_a.png", plot = fig1a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig1_prop_a.png", plot = new.prop.a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig1_prop_a2.png", plot = new.prop.a2, width = 12, height = 10, dpi = 600, units = "cm")

###########################################################################
# Figure 2
fig2_grid <- plot_grid(fig2a, fig2b, fig2c, fig2d, nrow = 2)
ggsave("fig2.png", plot = fig2_grid, width = 25, height = 20, dpi = 600, units = "cm")

###########################################################################
# Figure 3
fig3_grid <- plot_grid(fig3a, fig3b, fig3c, fig3d, nrow = 2)
ggsave("fig3.png", plot = fig3_grid, width = 25, height = 20, dpi = 600, units = "cm")

###########################################################################
# Figure 4
fig4_grid <- plot_grid(fig5a, fig6a, nrow = 1)
ggsave("fig4.png", plot = fig5a, width = 12, height = 10, dpi = 600, units = "cm")
ggsave("fig6.png", plot = fig6a, width = 20, height = 15, dpi = 600, units = "cm")

###########################################################################
# Beta/Figure 5
fig5_grid <- plot_grid(beta.a, beta.stream.a, beta.stream.b, nrow = 3)
ggsave("beta_plots.png", plot = fig5_grid, width = 25, height = 30, dpi = 600, units = "cm")
ggsave("betas.png", plot = beta.a, width = 15, height = 10, dpi = 600, units = "cm")

###########################################################################
# Supplemental Figure S1
supp_grid <- plot_grid(supp.a, supp.b, nrow = 1)
ggsave("supp1.png", plot = supp_grid, width = 25, height = 15, dpi = 600, units = "cm")
