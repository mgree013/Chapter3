# new_lake_analytics.R
# Load processed data and compute analytics / figure objects.

Packages <- c(
  "tidyverse", "vegan", "viridis", "FD", "betareg", "bbmle", "performance"
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

supp.b <- env_div %>%
  ggplot(aes(x = actual_fish_presence, fill = actual_fish_presence)) +
  geom_bar(stat = "count") +
  stat_count(geom = "text", colour = "black", size = 3.5, aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.25) +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ggtitle("b)") +
  ylab("Number of Lake Sites") +
  theme_minimal() +
  theme(legend.position = "none")

supp.a <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = lake_elevation_nbr, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ggtitle("a)") +
  ylab("Elevation (m)") +
  theme_minimal() +
  theme(legend.position = "none")

fig3a <- env_div %>%
  ggplot(aes(x = lake_elevation_nbr, y = N1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") +
  ylab("Species (Shannon) Diversity") +
  ggtitle("c)") +
  theme_minimal()

fig3b <- env_div %>%
  ggplot(aes(x = lake_elevation_nbr, y = betas.LCBD)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") +
  ylab("Beta-diversity (LCBD)") +
  ggtitle("d)") +
  theme_minimal()

fig3c <- env_div %>%
  ggplot(aes(x = lake_elevation_nbr, y = Com.Size)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") +
  ylab("Community Size") +
  ggtitle("c)") +
  theme_minimal()

fig3d <- env_div %>%
  ggplot(aes(x = lake_elevation_nbr, y = Biomass)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Elevation (m)") +
  ylab("Community Biomass") +
  ggtitle("d)") +
  theme_minimal()

fig2a <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = N1, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ylab("Species (Shannon) Diversity") +
  ggtitle("a)") +
  theme_minimal() +
  theme(legend.position = "none")

fig2b <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = betas.LCBD, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ylab("Beta-Diversity (LCBD)") +
  ggtitle("b)") +
  theme_minimal() +
  theme(legend.position = "none")

fig2c <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = Com.Size, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ylab("Community Size") +
  ggtitle("c)") +
  theme_minimal() +
  theme(legend.position = "none")

fig2d <- env_div %>%
  ggplot(aes(x = actual_fish_presence, y = Biomass, fill = actual_fish_presence)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish Presence", labels = c("No", "Yes")) +
  xlab("Fish Presence") +
  ylab("Community Biomass") +
  ggtitle("d)") +
  theme_minimal() +
  theme(legend.position = "none")

save(
  env_div,
  supp.a,
  supp.b,
  fig2a,
  fig2b,
  fig2c,
  fig2d,
  fig3a,
  fig3b,
  fig3c,
  fig3d,
  file = "data/new_lake_analytics.RData"
)

message("new_lake_analytics.R complete: analytics saved to data/new_lake_analytics.RData")
