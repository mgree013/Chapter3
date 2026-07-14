# new_lake_models.R
# Fit models on processed lake data.

Packages <- c(
  "tidyverse", "bbmle", "betareg", "performance"
)
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) install.packages(missing_packages)
lapply(Packages, library, character.only = TRUE)

load("data/new_lake_processed.RData")

env_div <- env %>%
  left_join(local_diversity, by = c("lake_id", "survey_date")) %>%
  filter(!lake_id %in% c(11505, 42219, 71257, 71282)) %>%
  mutate(
    Com.Size = log(Com.Size + 1),
    Biomass = log(Sum.Biomass + 1)
  ) %>%
  filter(lake_elevation_nbr > 1800, lake_elevation_nbr < 3500, HA >= 0.5, lake_max_depth > 3)

models <- list(
  N1_fish = glm(N1 ~ actual_fish_presence, family = gaussian(link = "identity"), data = env_div),
  N1_elev = glm(N1 ~ lake_elevation_nbr, family = gaussian(link = "identity"), data = env_div),
  N1_interact = glm(N1 ~ lake_elevation_nbr * actual_fish_presence, family = gaussian(link = "identity"), data = env_div),
  betas_fish = betareg(betas.LCBD ~ actual_fish_presence, data = env_div),
  betas_elev = betareg(betas.LCBD ~ lake_elevation_nbr, data = env_div),
  betas_interact = betareg(betas.LCBD ~ lake_elevation_nbr * actual_fish_presence, data = env_div),
  ComSize_fish = glm(Com.Size ~ actual_fish_presence, family = gaussian(link = "identity"), data = env_div),
  ComSize_elev = glm(Com.Size ~ lake_elevation_nbr, family = gaussian(link = "identity"), data = env_div),
  ComSize_interact = glm(Com.Size ~ lake_elevation_nbr * actual_fish_presence, family = gaussian(link = "identity"), data = env_div),
  Biomass_fish = glm(Biomass ~ actual_fish_presence, family = gaussian(link = "identity"), data = env_div),
  Biomass_elev = glm(Biomass ~ lake_elevation_nbr, family = gaussian(link = "identity"), data = env_div),
  Biomass_interact = glm(Biomass ~ lake_elevation_nbr * actual_fish_presence, family = gaussian(link = "identity"), data = env_div)
)

model_summaries <- map(models, summary)
model_r2 <- map(models, r2)

save(models, model_summaries, model_r2, file = "data/new_lake_models.RData")
message("new_lake_models.R complete: models saved to data/new_lake_models.RData")
