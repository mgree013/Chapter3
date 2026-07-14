# new_lake_analysis.R
# Download and process lake zooplankton and environmental data.

Packages <- c(
  "tidyverse", "reshape", "reshape2", "viridis", "vegan", "FD", "adespatial"
)
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) install.packages(missing_packages)
lapply(Packages, library, character.only = TRUE)

if (!file.exists("code/Download_Slip_Data.R")) {
  stop("Run this script from the project root where code/Download_Slip_Data.R exists.")
}

source("code/Download_Slip_Data.R", echo = FALSE)

if (!dir.exists("data")) dir.create("data")

av_zoop_body_size_new <- read.csv("data/length_mass_regress_zoop.csv", row.names = 1) %>%
  dplyr::select(Body_mass_ug)

zoopzz <- dt10 %>% left_join(dt12, by = "SpeciesID")

dt8.1 <- dt8 %>%
  dplyr::select(c(lake_id, survey_date, lake_max_depth, zoo_tow_number, zoo_tow_type, zoo_tow_depth))

clean_zoopzz <- zoopzz %>%
  left_join(dt8.1, by = c("lake_id", "survey_date")) %>%
  left_join(
    dt13 %>%
      group_by(lake_id, survey_date) %>%
      summarise(sample_vol = first(sample_vol), .groups = "drop"),
    by = c("lake_id", "survey_date")
  ) %>%
  group_by(lake_id, survey_date) %>%
  mutate(Max.Subsample = max(Subsample)) %>%
  ungroup() %>%
  group_by(lake_id, survey_date, Species_Name) %>%
  mutate(
    Counts = sum(Number) / Max.Subsample,
    zoop_density = Counts * sample_vol / (zoo_tow_depth * zoo_tow_number),
    Species_Name = str_replace(Species_Name, " ", "_"),
    Species_Name = str_replace(Species_Name, "/", "_")
  ) %>%
  filter(Species_Name != "standard_measurement") %>%
  dplyr::select(lake_id, survey_date, Species_Name, zoop_density) %>%
  distinct(lake_id, survey_date, Species_Name, zoop_density)

zoop_body_mass <- av_zoop_body_size_new %>% rownames_to_column("Taxon")

clean_zoopz_biomass <- zoopzz %>%
  left_join(dt8.1, by = c("lake_id", "survey_date")) %>%
  group_by(lake_id, survey_date) %>%
  mutate(Max.Subsample = max(Subsample)) %>%
  ungroup() %>%
  group_by(lake_id, survey_date, Species_Name) %>%
  mutate(
    Counts = sum(Number),
    Species_Name = str_replace(Species_Name, " ", "_"),
    Species_Name = str_replace(Species_Name, "/", "_")
  ) %>%
  dplyr::rename(Taxon = Species_Name) %>%
  left_join(zoop_body_mass, by = "Taxon") %>%
  group_by(lake_id, survey_date, Taxon) %>%
  mutate(Biomass = Counts * Body_mass_ug) %>%
  dplyr::select(c(lake_id, survey_date, Taxon, Biomass))

site.sp.quad <- cast(clean_zoopzz, lake_id + survey_date ~ Species_Name, value = 'zoop_density')
site.sp.quad <- as.data.frame(site.sp.quad)
site.sp.quad[is.na(site.sp.quad)] <- 0

species <- as.data.frame(site.sp.quad[, 3:ncol(site.sp.quad)])
diversity <- species %>%
  transmute(
    N0 = rowSums(species > 0),
    H = diversity(species),
    N1 = exp(H),
    N2 = diversity(species, "inv"),
    J = H / log(N0),
    E10 = N1 / N0,
    E20 = N2 / N0,
    Com.Size = rowSums(species),
    betas.LCBD = beta.div(species, method = "hellinger", sqrt.D = TRUE)$LCBD
  )

local_diversity <- bind_cols(diversity, site.sp.quad[, 1:2])

clean_zoopz_biomass_site <- clean_zoopz_biomass %>%
  filter(Taxon != "standard_measurement", Biomass > 0) %>%
  group_by(lake_id, survey_date) %>%
  summarise(Sum.Biomass = sum(Biomass), .groups = "drop")

local_diversity <- local_diversity %>%
  left_join(clean_zoopz_biomass_site, by = c("lake_id", "survey_date"))

gps <- read.csv("data/sierra_lakes.csv")

dt5.1 <- dt5 %>%
  pivot_wider(names_from = littoral_type, names_glue = "{littoral_type}_{.value}", values_from = littoral_amount)

dt6.1 <- dt6 %>%
  pivot_wider(names_from = shoreline_type, names_glue = "{shoreline_type}_{.value}", values_from = shoreline_amount)

env <- dt5.1 %>%
  left_join(dt6.1, by = c("lake_id", "survey_date")) %>%
  left_join(dt4, by = "lake_id") %>%
  left_join(dt8, by = c("lake_id", "survey_date")) %>%
  filter(zoo_sample_ind == "Yes") %>%
  dplyr::select(-c(
    zoo_sample_ind, zoo_sample_time, zoo_tow_number, zoo_tow_type, zoo_tow_depth,
    benthic_sample_ind, benthic_sample_percent, nbr_benthic_sweeps,
    lake_fairy_shrimp_ind, lake_shrimp_collection, pool_fairy_shrimp_ind, pool_shrimp_collection,
    amphib_survey_starttime, amphib_survey_endtime, amphib_survey_duration, amphib_survey_desc, amphib_survey_fish_presence,
    fish_survey_type, fish_net_location_type, fish_net_set_datetime, fish_net_pull_datetime
  )) %>%
  left_join(gps, by = "lake_id")

sp_abund_env <- left_join(site.sp.quad, env, by = c("lake_id", "survey_date")) %>%
  filter(lake_id != "70534")

pivot_clean_zoopzz <- clean_zoopzz %>%
  pivot_wider(names_from = "Species_Name", values_from = "zoop_density")

env_abund <- left_join(clean_zoopzz, env, by = c("lake_id", "survey_date")) %>%
  filter(actual_fish_presence %in% c("Yes", "No"))

save(
  clean_zoopzz,
  clean_zoopz_biomass,
  site.sp.quad,
  diversity,
  local_diversity,
  clean_zoopz_biomass_site,
  av_zoop_body_size_new,
  zoop_body_mass,
  env,
  sp_abund_env,
  pivot_clean_zoopzz,
  env_abund,
  file = "data/new_lake_processed.RData"
)

message("new_lake_analysis.R complete: data saved to data/new_lake_processed.RData")
