# new_lake_analytics.R
# Build all manuscript figures (Figs 1-4, S1-S4) and models (Tables 1-4, S3).

Packages <- c(
  "tidyverse", "vegan", "viridis", "FD", "betareg", "bbmle", "performance", "betapart"
)
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) install.packages(missing_packages)
lapply(Packages, library, character.only = TRUE)

load("data/new_lake_processed.RData")

theme_ms <- theme(
  axis.line        = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border     = element_blank(),
  panel.background = element_blank()
)

# ── env_div: filtered working dataset ────────────────────────────────────────
env_div <- env %>%
  left_join(local_diversity, by = c("lake_id", "survey_date")) %>%
  filter(!lake_id %in% c(11505, 42219, 71257, 71282)) %>%
  mutate(Com.Size = log(Com.Size + 1), Biomass = log(Sum.Biomass + 1)) %>%
  filter(lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
         HA >= 0.5, lake_max_depth > 3)

# ── Figure S1: fish presence vs elevation + site counts ──────────────────────
fig_s1a <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = lake_elevation_nbr,
             fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  xlab("Fish Presence") + ylab("Elevation (m)") + ggtitle("a)") +
  theme_ms + theme(legend.position = "none")

fig_s1b <- env_div %>%
  ggplot(aes(x = actual_fish_presence, fill = actual_fish_presence)) +
  geom_bar(stat = "count") +
  stat_count(geom = "text", colour = "black", size = 3.5,
             aes(label = after_stat(count)), vjust = -0.25) +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  xlab("Fish Presence") + ylab("Number of Lake Sites") + ggtitle("b)") +
  theme_ms + theme(legend.position = "none")

# ── Figure 1: Shannon diversity and LCBD by fish presence (a,b) and elevation (c,d) ──
env_div_model <- env_div %>% filter(!is.na(N1), !is.na(betas.LCBD))

fig1a <- env_div_model %>%
  ggplot(aes(x = actual_fish_presence, y = N1,
             fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  xlab("Fish Presence") + ylab("Species (Shannon) Diversity") + ggtitle("a)") +
  theme_ms + theme(legend.position = "none")

fig1b <- env_div_model %>%
  ggplot(aes(x = actual_fish_presence, y = betas.LCBD,
             fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  xlab("Fish Presence") + ylab("Beta-Diversity (LCBD)") + ggtitle("b)") +
  theme_ms + theme(legend.position = "none")

fig1c <- env_div_model %>%
  ggplot(aes(x = lake_elevation_nbr, y = N1)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") + ylab("Species (Shannon) Diversity") + ggtitle("c)") +
  theme_ms

fig1d <- env_div_model %>%
  ggplot(aes(x = lake_elevation_nbr, y = betas.LCBD)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") + ylab("Beta-diversity (LCBD)") + ggtitle("d)") +
  theme_ms

# ── Table 1 models: Shannon diversity and LCBD ───────────────────────────────
m_N1_elev     <- glm(N1 ~ lake_elevation_nbr,
                     family = gaussian, data = env_div_model)
m_N1_fish     <- glm(N1 ~ actual_fish_presence,
                     family = gaussian, data = env_div_model)
m_N1_interact <- glm(N1 ~ lake_elevation_nbr * actual_fish_presence,
                     family = gaussian, data = env_div_model)
m_N1_null     <- glm(N1 ~ 1,
                     family = gaussian, data = env_div_model)

m_LCBD_elev     <- betareg(betas.LCBD ~ lake_elevation_nbr,     data = env_div_model)
m_LCBD_fish     <- betareg(betas.LCBD ~ actual_fish_presence,   data = env_div_model)
m_LCBD_interact <- betareg(betas.LCBD ~ lake_elevation_nbr * actual_fish_presence,
                            data = env_div_model)
m_LCBD_null     <- betareg(betas.LCBD ~ 1,                      data = env_div_model)

table1_N1   <- bbmle::AICtab(m_N1_elev, m_N1_fish, m_N1_interact, m_N1_null,
                              weights = TRUE, sort = FALSE)
table1_LCBD <- bbmle::AICtab(m_LCBD_elev, m_LCBD_fish, m_LCBD_interact, m_LCBD_null,
                              weights = TRUE, sort = FALSE)

# ── Figure 2: Beta diversity components (betapart) ───────────────────────────
species_cols <- setdiff(colnames(site.sp.quad), c("lake_id", "survey_date"))

sp_filt <- sp_abund_env %>%
  filter(lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
         HA >= 0.5, lake_max_depth > 3) %>%
  unite("site_id", lake_id, survey_date, remove = FALSE)

avail_spp <- intersect(species_cols, colnames(sp_filt))
spp_mat <- sp_filt %>%
  dplyr::select(all_of(avail_spp)) %>%
  mutate(across(everything(), as.numeric)) %>%
  replace(is.na(.), 0)

keep <- rowSums(spp_mat) > 0
spp_mat <- spp_mat[keep, , drop = FALSE]
sp_filt <- sp_filt[keep, ]
rownames(spp_mat) <- sp_filt$site_id

lake.dist <- beta.pair.abund(spp_mat, index.family = "bray")

beta_df <- data.frame(
  dist.bal  = as.vector(lake.dist$beta.bray.bal),
  dist.gra  = as.vector(lake.dist$beta.bray.gra),
  dist.bray = as.vector(lake.dist$beta.bray)
) %>% pivot_longer(everything(), names_to = "comp", values_to = "beta")

fig2 <- beta_df %>%
  ggplot(aes(x = comp, y = beta, fill = comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(
    limits = c("dist.bal", "dist.gra", "dist.bray"),
    labels = c(expression(beta[bal]), expression(beta[gra]), expression(beta[bray]))
  ) +
  xlab("Turnover and Nestedness Components") +
  ylab(expression("Zooplankton " * beta * "-diversity")) +
  theme_ms + theme(legend.position = "none")

# Table 2: betareg comparing bal vs gra components
beta_comp_df <- beta_df %>%
  filter(comp %in% c("dist.bal", "dist.gra")) %>%
  mutate(n = n(), beta_adj = (beta * (n - 1) + 0.5) / n)

m_beta_comp <- betareg(beta_adj ~ comp, data = beta_comp_df)
m_beta_null <- betareg(beta_adj ~ 1,    data = beta_comp_df)
table2_beta <- bbmle::AICtab(m_beta_comp, m_beta_null, weights = TRUE, sort = FALSE)

# ── Individual-species data ───────────────────────────────────────────────────
av_zoop_bm <- av_zoop_body_size_new %>% rownames_to_column("Taxon")

env_abund_sp <- env_abund %>%
  left_join(av_zoop_bm, by = c("Species_Name" = "Taxon")) %>%
  dplyr::rename(Taxon = Species_Name) %>%
  filter(
    lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
    HA >= 0.5, lake_max_depth > 3,
    !Taxon %in% c("Collotheca_spp.", "Eurycercus_lamellatus",
                  "Lecane_spp.", "Monostyla_spp.",
                  "Polyarthra_vulgaris", "Simocephalus_serrulatus",
                  "nauplii"),
    !is.na(Body_mass_ug)
  )

# ── Figure S2a: species log density by fish, ordered by body size ─────────────
fig_s2a <- env_abund_sp %>%
  ggplot(aes(x = reorder(Taxon, Body_mass_ug, FUN = median),
             y = log(zoop_density + 1), fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  ylab("Zooplankton Log Density + 1") + xlab("Zooplankton Taxa") + ggtitle("a)") +
  theme_ms +
  theme(axis.text.x           = element_text(angle = 60, hjust = 1),
        legend.position       = c(0.9, 0.85),
        legend.background     = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# ── Figure S2b: proportion of lakes occupied ──────────────────────────────────
n_by_fish <- env_div %>%
  filter(!is.na(actual_fish_presence)) %>%
  dplyr::group_by(actual_fish_presence) %>%
  dplyr::summarise(total_sites = dplyr::n(), .groups = "drop")

occupancy_df <- env_abund_sp %>%
  group_by(Taxon, Body_mass_ug, actual_fish_presence) %>%
  summarise(n_occ = sum(!is.na(zoop_density) & zoop_density > 0, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(n_by_fish, by = "actual_fish_presence") %>%
  mutate(occupancy = n_occ / total_sites) %>%
  dplyr::rename(Fish = actual_fish_presence)

fig_s2b <- occupancy_df %>%
  ggplot(aes(x = reorder(Taxon, Body_mass_ug), y = occupancy,
             fill = Fish, group = Fish)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence") +
  xlab("Zooplankton Taxa") + ylab("Proportion of Lakes Occupied") + ggtitle("b)") +
  theme_ms + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                   legend.position = "none")

# ── Figure 3: relative change vs body mass — counts + occupancy filter ────────
# Counts (not density) avoids sample_vol/tow_depth/tow_number confounding with
# fish presence. 5% minimum occupancy in BOTH groups ensures stable mean
# estimates and removes species represented by only 1–4 sites in one group.

n_fish_sites <- env_abund_sp %>% ungroup() %>%
  filter(actual_fish_presence == "Yes") %>%
  summarise(n = n_distinct(paste(lake_id, survey_date))) %>% pull(n)
n_nofish_sites <- env_abund_sp %>% ungroup() %>%
  filter(actual_fish_presence == "No") %>%
  summarise(n = n_distinct(paste(lake_id, survey_date))) %>% pull(n)

occ_fish <- env_abund_sp %>% ungroup() %>%
  filter(actual_fish_presence == "Yes", Counts > 0) %>%
  group_by(Taxon) %>%
  summarise(n_fish_occ = n_distinct(paste(lake_id, survey_date)), .groups = "drop")
occ_nofish <- env_abund_sp %>% ungroup() %>%
  filter(actual_fish_presence == "No", Counts > 0) %>%
  group_by(Taxon) %>%
  summarise(n_nofish_occ = n_distinct(paste(lake_id, survey_date)), .groups = "drop")

occ_5pct <- occ_fish %>%
  inner_join(occ_nofish, by = "Taxon") %>%
  mutate(pct_fish   = n_fish_occ   / n_fish_sites   * 100,
         pct_nofish = n_nofish_occ / n_nofish_sites * 100) %>%
  filter(pmin(pct_fish, pct_nofish) >= 5) %>%
  pull(Taxon)

density_change <- env_abund_sp %>%
  filter(Taxon %in% occ_5pct) %>%
  group_by(Taxon, actual_fish_presence, Body_mass_ug) %>%
  summarise(Mean_val = mean(log(Counts + 1), na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = actual_fish_presence, values_from = Mean_val) %>%
  filter(is.finite(No), is.finite(Yes), No > 0, Yes > 0) %>%
  mutate(relative_change = Yes / No)

fig3 <- density_change %>%
  ggplot(aes(x = log(Body_mass_ug + 1), y = relative_change)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ylab("Relative Change in Zooplankton Density") +
  xlab(expression("Zooplankton Log Body Mass (" * mu * "g)")) +
  theme_ms

m_density_bm   <- glm(relative_change ~ log(Body_mass_ug + 1),
                       family = gaussian, data = density_change)
m_density_null <- glm(relative_change ~ 1, family = gaussian, data = density_change)
table3_density <- bbmle::AICtab(m_density_bm, m_density_null,
                                 weights = TRUE, sort = FALSE)

# ── Figure S3 + Table S3: body mass of absent species ────────────────────────
# Use a broader source — only apply elevation/depth/HA filters, NOT the
# Polyarthra_vulgaris exclusion (those are the species we want to identify here)
env_abund_for_absent <- env_abund %>%
  left_join(av_zoop_bm, by = c("Species_Name" = "Taxon")) %>%
  dplyr::rename(Taxon = Species_Name) %>%
  filter(
    lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
    HA >= 0.5, lake_max_depth > 3,
    !Taxon %in% c("Collotheca_spp.", "Eurycercus_lamellatus",
                  "Lecane_spp.", "Monostyla_spp.",
                  "Simocephalus_serrulatus", "nauplii"),
    !is.na(Body_mass_ug)
  )

fish_spp     <- env_abund_for_absent %>%
  filter(actual_fish_presence == "Yes", !is.na(zoop_density), zoop_density > 0) %>%
  pull(Taxon) %>% unique()
fishless_spp <- env_abund_for_absent %>%
  filter(actual_fish_presence == "No",  !is.na(zoop_density), zoop_density > 0) %>%
  pull(Taxon) %>% unique()

absent_df <- av_zoop_bm %>%
  filter(Taxon %in% union(fish_spp, fishless_spp)) %>%
  mutate(absent_from = case_when(
    !Taxon %in% fish_spp     ~ "Fish",
    !Taxon %in% fishless_spp ~ "Fishless",
    TRUE                      ~ NA_character_
  )) %>%
  filter(!is.na(absent_from))

fig_s3 <- absent_df %>%
  ggplot(aes(x = absent_from, y = log(Body_mass_ug + 1), fill = absent_from)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Species Absent From") + ylab("Zooplankton Log Body Mass (µg)") +
  theme_ms + theme(legend.position = "none")

m_absent_fish <- glm(log(Body_mass_ug + 1) ~ absent_from, family = gaussian, data = absent_df)
m_absent_null <- glm(log(Body_mass_ug + 1) ~ 1,           family = gaussian, data = absent_df)
tableS3_absent <- bbmle::AICtab(m_absent_fish, m_absent_null,
                                 weights = TRUE, sort = FALSE)

# ── Figure 4 + S4: CWM (Table 4) ─────────────────────────────────────────────
av_for_cwm <- av_zoop_body_size_new %>% dplyr::select(Body_mass_ug)

cwm_wide <- clean_zoopzz %>%
  pivot_wider(id_cols    = c(lake_id, survey_date),
              names_from  = "Species_Name",
              values_from = "zoop_density") %>%
  unite("id_date", lake_id, survey_date, remove = TRUE) %>%
  column_to_rownames("id_date") %>%
  replace(is.na(.), 0) %>%
  dplyr::select(sort(names(.)))

shared_cwm <- intersect(colnames(cwm_wide), rownames(av_for_cwm))
if (length(shared_cwm) == 0) stop("No shared species between trait and abundance matrices.")
cwm_wide   <- cwm_wide[,    shared_cwm, drop = FALSE]
av_for_cwm <- av_for_cwm[shared_cwm, , drop = FALSE]

zero_spp  <- which(colSums(cwm_wide) == 0)
if (length(zero_spp) > 0) {
  cwm_wide   <- cwm_wide[,  -zero_spp, drop = FALSE]
  av_for_cwm <- av_for_cwm[-zero_spp, , drop = FALSE]
}
zero_sites <- which(rowSums(cwm_wide) == 0)
if (length(zero_sites) > 0) cwm_wide <- cwm_wide[-zero_sites, , drop = FALSE]

tres_bm <- dbFD(av_for_cwm, cwm_wide, corr = "lingoes",
                stand.FRic = TRUE, calc.FDiv = TRUE)

cwm_df <- data.frame(CWM     = tres_bm$CWM$Body_mass_ug,
                     id_date = rownames(cwm_wide)) %>%
  separate(id_date, sep = "_", into = c("lake_id", "survey_date"), extra = "merge") %>%
  mutate(lake_id = as.integer(lake_id))

env_cwm <- env %>%
  left_join(cwm_df, by = c("lake_id", "survey_date")) %>%
  filter(!lake_id %in% c(11505, 42219, 71257, 71282),
         lake_elevation_nbr > 1800, lake_elevation_nbr < 3500,
         HA >= 0.5, lake_max_depth > 3,
         !is.na(CWM))

# Figure 4: CWM boxplot by fish presence
fig4 <- env_cwm %>%
  ggplot(aes(x = actual_fish_presence, y = CWM * 0.01, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence",
                     labels = c("No", "Yes")) +
  xlab("Fish Presence") + ylab("Zooplankton CWM") +
  theme_ms + theme(legend.position = "none")

# Figure S4: CWM vs elevation, faceted by fish presence
fig_s4 <- env_cwm %>%
  ggplot(aes(x = lake_elevation_nbr, y = CWM, colour = actual_fish_presence)) +
  geom_point() + geom_smooth(method = "lm") +
  scale_color_viridis_d(name = "Fish Presence", labels = c("No", "Yes")) +
  facet_grid(~actual_fish_presence) +
  xlab("Elevation (m)") + ylab("Zooplankton CWM") +
  theme_ms +
  theme(legend.position        = c(0.8, 0.9),
        legend.background      = element_blank(),
        legend.box.background  = element_rect(colour = "black"))

m_CWM_interact <- glm(CWM ~ lake_elevation_nbr * actual_fish_presence,
                       family = gaussian, data = env_cwm)
m_CWM_elev     <- glm(CWM ~ lake_elevation_nbr,
                       family = gaussian, data = env_cwm)
m_CWM_fish     <- glm(CWM ~ actual_fish_presence,
                       family = gaussian, data = env_cwm)
m_CWM_null     <- glm(CWM ~ 1, family = gaussian, data = env_cwm)
table4_CWM <- bbmle::AICtab(m_CWM_interact, m_CWM_elev, m_CWM_fish, m_CWM_null,
                              weights = TRUE, sort = FALSE)

# ── Save all objects ──────────────────────────────────────────────────────────
save(
  env_div, env_cwm,
  fig1a, fig1b, fig1c, fig1d,
  fig2,
  fig3,
  fig4,
  fig_s1a, fig_s1b,
  fig_s2a, fig_s2b,
  fig_s3,  fig_s4,
  m_N1_elev, m_N1_fish, m_N1_interact, m_N1_null,
  m_LCBD_elev, m_LCBD_fish, m_LCBD_interact, m_LCBD_null,
  table1_N1, table1_LCBD,
  m_beta_comp, m_beta_null, table2_beta,
  m_density_bm, m_density_null, table3_density,
  m_CWM_interact, m_CWM_elev, m_CWM_fish, m_CWM_null, table4_CWM,
  m_absent_fish, m_absent_null, tableS3_absent,
  file = "data/new_lake_analytics.RData"
)

# ── Save manuscript figures (while data frames are still in scope) ────────────
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)
dir.create("Newfigs", showWarnings = FALSE)

ggsave("Newfigs/Figure1.png",
       plot = plot_grid(fig1a, fig1b, fig1c, fig1d, nrow = 2),
       width = 25, height = 20, dpi = 600, units = "cm")
ggsave("Newfigs/Figure2.png",  plot = fig2,
       width = 15, height = 12, dpi = 600, units = "cm")
ggsave("Newfigs/Figure3.png",  plot = fig3,
       width = 15, height = 12, dpi = 600, units = "cm")
ggsave("Newfigs/Figure4.png",  plot = fig4,
       width = 12, height = 12, dpi = 600, units = "cm")
ggsave("Newfigs/FigureS1.png",
       plot = plot_grid(fig_s1a, fig_s1b, nrow = 1),
       width = 25, height = 12, dpi = 600, units = "cm")
ggsave("Newfigs/FigureS2.png",
       plot = plot_grid(fig_s2a, fig_s2b, nrow = 2),
       width = 30, height = 30, dpi = 600, units = "cm")
ggsave("Newfigs/FigureS3.png", plot = fig_s3,
       width = 12, height = 12, dpi = 600, units = "cm")
ggsave("Newfigs/FigureS4.png", plot = fig_s4,
       width = 25, height = 12, dpi = 600, units = "cm")

message("new_lake_analytics.R complete: all manuscript figures and models saved.")

# ── Export Manuscript Tables ──────────────────────────────────────────────────
dir.create("Tables", showWarnings = FALSE)

# R² helper: pseudo-R² for GLM (McFadden) or betareg summary
r2_mod <- function(m) {
  if (inherits(m, "betareg")) {
    return(round(summary(m)$pseudo.r.squared, 2))
  }
  if (!is.null(m$null.deviance) && m$null.deviance != 0) {
    return(round((m$null.deviance - m$deviance) / m$null.deviance, 2))
  }
  NA_real_
}

# Format a single AICtab into a data.frame with model labels + R²
fmt_table <- function(aic_tab, models_list, model_labels, response_label = NULL) {
  dtab      <- as.data.frame(aic_tab)
  row_count <- nrow(dtab)
  out <- data.frame(
    Response = if (!is.null(response_label)) rep(response_label, row_count) else NA,
    Model    = model_labels,
    dAIC     = round(dtab$dAIC, 1),
    df       = as.integer(dtab$df),
    wi       = round(dtab$weight, 3),
    R2       = sapply(models_list, r2_mod),
    stringsAsFactors = FALSE
  )
  out
}

# ── Table 1: Shannon diversity and LCBD ──────────────────────────────────────
t1_N1 <- fmt_table(
  table1_N1,
  list(m_N1_elev, m_N1_fish, m_N1_interact, m_N1_null),
  c("~Elevation", "~Fish", "~Fish*Elevation", "~1"),
  "Shannon Diversity"
)
t1_LCBD <- fmt_table(
  table1_LCBD,
  list(m_LCBD_elev, m_LCBD_fish, m_LCBD_interact, m_LCBD_null),
  c("~Elevation", "~Fish", "~Fish*Elevation", "~1"),
  "LCBD"
)
write.csv(rbind(t1_N1, t1_LCBD),
          "Tables/Table1_diversity_models.csv", row.names = FALSE)

# ── Table 2: Beta diversity component ────────────────────────────────────────
t2 <- fmt_table(
  table2_beta,
  list(m_beta_comp, m_beta_null),
  c("~Beta diversity component", "~1"),
  "Lake beta diversity"
)
write.csv(t2, "Tables/Table2_beta_diversity.csv", row.names = FALSE)

# ── Table 3: Relative change in species density ──────────────────────────────
t3 <- fmt_table(
  table3_density,
  list(m_density_bm, m_density_null),
  c("~Body Mass", "~1"),
  "Relative Change in Zooplankton Counts (>=5% occupancy filter)"
)
write.csv(t3, "Tables/Table3_density_change.csv", row.names = FALSE)

# ── Table 4: CWM ─────────────────────────────────────────────────────────────
t4 <- fmt_table(
  table4_CWM,
  list(m_CWM_interact, m_CWM_elev, m_CWM_fish, m_CWM_null),
  c("~Fish*Elevation", "~Elevation", "~Fish", "~1"),
  "Zooplankton CWM"
)
write.csv(t4, "Tables/Table4_CWM_models.csv", row.names = FALSE)

# ── Table S1: elevation ~ fish presence (ANOVA) ──────────────────────────────
m_s1 <- lm(lake_elevation_nbr ~ actual_fish_presence,
            data = env_div %>% filter(!is.na(actual_fish_presence)))
s1_coef <- as.data.frame(summary(m_s1)$coefficients)
s1_coef <- cbind(
  Factor = c("(Intercept)", "Fish"),
  round(s1_coef, 3)
)
colnames(s1_coef) <- c("Factor", "Estimate", "Std.Error", "t.value", "P.value")
write.csv(s1_coef, "Tables/TableS1_fish_elevation.csv", row.names = FALSE)

# ── Table S2: species absent from fish or fishless sites ─────────────────────
tableS2_df <- absent_df %>%
  dplyr::select(Taxon, absent_from) %>%
  dplyr::rename(Taxa = Taxon, "Site_Absent_From" = absent_from) %>%
  dplyr::arrange(Site_Absent_From, Taxa)
write.csv(tableS2_df, "Tables/TableS2_absent_species.csv", row.names = FALSE)

# ── Table S3: body mass ~ absent from fish/fishless ──────────────────────────
tS3 <- fmt_table(
  tableS3_absent,
  list(m_absent_fish, m_absent_null),
  c("~Fish", "~1"),
  "Log Body Mass"
)
write.csv(tS3, "Tables/TableS3_absent_bodymass.csv", row.names = FALSE)

message("Tables exported to Tables/")


