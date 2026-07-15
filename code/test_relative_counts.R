# test_relative_counts.R
# Test whether using subsample-corrected COUNTS (instead of absolute density)
# for the species-level relative change analysis restores the body-size effect
# weakened by the corrected density formula (Table 3 / Figure 3).
#
# Rationale:
#   zoop_density = Counts * sample_vol / (tow_depth * tow_number * 33.02)
#   For RELATIVE change (fish/fishless), any per-site scalar cancels IF the
#   scalar is uncorrelated with fish presence. sample_vol is the risky term;
#   tow_depth and tow_number are ecological covariates that could differ
#   between fish and fishless sites. Using raw Counts removes all three scalars
#   and isolates the biological signal.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(bbmle); library(cowplot)
})

load("data/new_lake_processed.RData")

av_zoop_bm <- av_zoop_body_size_new %>% rownames_to_column("Taxon")

theme_ms <- theme(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(), panel.background = element_blank()
)

# Shared taxon exclusions (same as analytics script Fig 3)
exclude_taxa <- c("Collotheca_spp.", "Eurycercus_lamellatus", "Lecane_spp.",
                  "Monostyla_spp.", "Polyarthra_vulgaris",
                  "Simocephalus_serrulatus", "nauplii")

env_abund_filt <- env_abund %>%
  left_join(av_zoop_bm, by = c("Species_Name" = "Taxon")) %>%
  dplyr::rename(Taxon = Species_Name) %>%
  filter(
    lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
    HA >= 0.5, lake_max_depth > 3,
    !Taxon %in% exclude_taxa,
    !is.na(Body_mass_ug)
  )

# ── Helper: compute relative change and fit models ────────────────────────────
fit_rel_change <- function(data, value_col, label) {
  rc <- data %>%
    group_by(Taxon, actual_fish_presence, Body_mass_ug) %>%
    summarise(Mean_val = mean(log(.data[[value_col]] + 1), na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = actual_fish_presence, values_from = Mean_val) %>%
    filter(is.finite(No), is.finite(Yes), No > 0, Yes > 0) %>%
    mutate(relative_change = Yes / No,
           method = label)

  m_bm   <- glm(relative_change ~ log(Body_mass_ug + 1),
                family = gaussian, data = rc)
  m_null <- glm(relative_change ~ 1,
                family = gaussian, data = rc)

  aic_tab <- bbmle::AICtab(m_bm, m_null, weights = TRUE)
  r2_bm   <- round((m_bm$null.deviance - m_bm$deviance) / m_bm$null.deviance, 3)
  coef_bm <- round(coef(m_bm)[2], 4)
  wi_bm   <- round(as.data.frame(aic_tab)$weight[1], 3)
  daic_bm <- round(as.data.frame(aic_tab)$dAIC[1], 2)

  list(rc = rc, m_bm = m_bm, m_null = m_null, aic_tab = aic_tab,
       r2 = r2_bm, coef = coef_bm, wi = wi_bm, dAIC = daic_bm,
       n_spp = nrow(rc))
}

# ── Method 1: subsample-corrected COUNTS ─────────────────────────────────────
res_counts <- fit_rel_change(env_abund_filt, "Counts", "Counts")

# ── Method 2: absolute density (corrected formula) ───────────────────────────
res_density <- fit_rel_change(env_abund_filt, "zoop_density", "Density")

# ── Summary table ────────────────────────────────────────────────────────────
cat("\n=======================================================\n")
cat("  COMPARISON: Relative change ~ body mass\n")
cat("=======================================================\n\n")

cat(sprintf("%-30s %6s %8s %6s %8s %8s\n",
            "Method", "n_spp", "~BM wi", "dAIC", "R2", "slope"))
cat(strrep("-", 70), "\n")
cat(sprintf("%-30s %6d %8.3f %6.2f %8.3f %8.4f\n",
            "Counts (subsample-corrected)",
            res_counts$n_spp, res_counts$wi, res_counts$dAIC,
            res_counts$r2, res_counts$coef))
cat(sprintf("%-30s %6d %8.3f %6.2f %8.3f %8.4f\n",
            "Density (corrected formula)",
            res_density$n_spp, res_density$wi, res_density$dAIC,
            res_density$r2, res_density$coef))
cat(sprintf("%-30s %6s %8s %6s %8s %8s\n",
            "MANUSCRIPT (geometric)", "~28", "0.936", "5.4", "0.130", "-"))
cat(strrep("-", 70), "\n\n")

cat("AIC table — COUNTS method:\n")
print(res_counts$aic_tab)
cat("Slope coefficient:", res_counts$coef, "\n")
cat("Negative slope = larger-bodied species decline more in fish sites\n\n")

cat("AIC table — DENSITY method:\n")
print(res_density$aic_tab)
cat("Slope coefficient:", res_density$coef, "\n\n")

# ── Interpretation ────────────────────────────────────────────────────────────
cat("=======================================================\n")
cat("  INTERPRETATION\n")
cat("=======================================================\n\n")

if (res_counts$wi > 0.7 && res_counts$coef < 0) {
  cat("COUNTS method SUPPORTS the manuscript conclusion:\n")
  cat("  Larger-bodied species have greater relative decline in fish sites.\n")
  cat("  Model weight (wi =", res_counts$wi, ") vs density method (wi =",
      res_density$wi, ").\n")
  cat("  Using counts restores the body-size signal lost by sample_vol noise.\n")
} else if (res_counts$coef < 0 && res_counts$wi > res_density$wi) {
  cat("COUNTS method shows STRONGER (but not conclusive) body-size signal:\n")
  cat("  Negative slope retained, wi improved from", res_density$wi,
      "to", res_counts$wi, "\n")
  cat("  Suggests sample_vol variation was adding noise to the density signal.\n")
} else {
  cat("COUNTS method does NOT clearly restore the body-size conclusion.\n")
  cat("  The weakened body-size effect persists even without sample_vol.\n")
  cat("  Consider: tow_depth/tow_number may also vary with fish presence.\n")
}

# ── Figures ───────────────────────────────────────────────────────────────────
make_fig3 <- function(rc_df, subtitle) {
  rc_df %>%
    ggplot(aes(x = log(Body_mass_ug + 1), y = relative_change)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, colour = "steelblue") +
    geom_hline(yintercept = 1, linetype = "dotted") +
    ylab("Relative Change in Zooplankton\n(fish sites / fishless sites)") +
    xlab(expression("Log Body Mass (" * mu * "g)")) +
    ggtitle(subtitle) +
    theme_ms
}

fig_counts <- make_fig3(res_counts$rc,
  paste0("Using COUNTS  (wi=", res_counts$wi, ", R²=", res_counts$r2, ")"))
fig_density <- make_fig3(res_density$rc,
  paste0("Using DENSITY (wi=", res_density$wi, ", R²=", res_density$r2, ")"))

panel <- plot_grid(fig_counts, fig_density, nrow = 1,
                   labels = c("Counts", "Density"))

dir.create("Newfigs", showWarnings = FALSE)
ggsave("Newfigs/Fig3_counts_vs_density_test.png",
       plot = panel, width = 28, height = 12, dpi = 300, units = "cm")
message("Test figure saved: Newfigs/Fig3_counts_vs_density_test.png")
