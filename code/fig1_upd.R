#Figure 1 FIx


fig3a<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=N1))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("a)") +
  xlab("Elevation (m)")+ylab("Species (Shannon) Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig3b<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=betas.LCBD))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("c)") +
  xlab("Elevation (m)")+ ylab(expression(paste(beta, "-Diversity (LCBD)"))) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +  # Round y-axis labels to 3 decimals
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig2a<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=N1, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Species (Shannon) Diversity")+
  ggtitle("b)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2b<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=betas.LCBD, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence") + ylab(expression(paste(beta, "-Diversity (LCBD)"))) +  # Use the beta symbol on y-axis label
  ggtitle("d)") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +  # Round y-axis labels to 3 decimals
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

library(cowplot)
plot_grid(fig3a,fig2a,fig3b,fig2b, ncol = 2)

############################################################################################################################################
#Intrapsecific Body Size
setwd("/Users/matthewgreen/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3/data")
av_zoop_body_size_new<-read.csv("length_convert.csv")

av_zoop_body_size_new<-av_zoop_body_size_new%>%dplyr::rename(Species_Name = Taxon)

fish_present=env%>%select(c(lake_id,actual_fish_presence, lake_elevation_nbr,HA,lake_max_depth))

#dt11/dt12
zoop_body_size_news<-dt11%>%
  left_join(dt12, by = "SpeciesID")%>%
  mutate(Species_Name = str_replace(Species_Name, " ", "_"),
         Species_Name = str_replace(Species_Name, "/", "_"))%>%
  left_join(av_zoop_body_size_new, by= "Species_Name")%>%
  mutate(ln_length=log(Length),
    new_weight=(lna+b*ln_length),
         new_weight_ln=log(new_weight))%>%
  filter(Species_Name != "standard_measurement")%>%
  left_join(fish_present, by="lake_id")%>%filter(actual_fish_presence=="Yes" |actual_fish_presence=="No")%>%
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3200)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)%>%
  filter(Species_Name != "Ascomorpha_spp.")%>%  filter(Species_Name != "Asplanchna_spp.")%>%  filter(Species_Name != "Daphnia_pulex_pulicaria")%>%
  filter(Species_Name != "Polyarthra_vulgaris")%>%  filter(Species_Name != "Polyphemus_pediculus")%>%  filter(Species_Name != "Synchaeta_spp.")%>% filter(Species_Name !="Trichotria_spp.")%>%
  filter(Species_Name != "Harpacticoida")%>%  filter(Species_Name != "Synchaeta_spp.")%>%  filter(Species_Name != "Polyphemus_pediculus")%>% filter(Species_Name !="Trichotria_spp.")%>% 
  filter(Species_Name !="Alonella_excisa")%>%filter(Species_Name !="Alonella_spp.")%>%filter(Species_Name !="Alona_spp.")
  
species_fish_presence <- zoop_body_size_news %>%
  group_by(Species_Name) %>%
  summarise(
    fish_presence = paste(unique(actual_fish_presence), collapse = ", ")
  )


zoop_body_size_news%>%
  ggplot(aes(x=Species_Name, y=new_weight_ln,  fill=as.factor(actual_fish_presence )))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Body Weight")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #ggtitle("b)") +
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

zoop_body_size_news%>%
  #filter(lake_elevation_nbr>3200)%>%
  ggplot(aes(x=Species_Name, y=new_weight_ln,  fill=as.factor(actual_fish_presence )))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Zooplankton Taxa")+ylab(expression("Zooplankton Body Mass (" * mu * "g/L, Log"[10] * ")"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #ggtitle("b)") +
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

dat<-zoop_body_size_news%>%group_by(Species_Name,actual_fish_presence)%>%summarise(av_weight=mean(new_weight_ln))

# Perform ANOVA
anova_model <- aov(new_weight_ln ~ Species_Name, data = zoop_body_size_news)
zoop_body_size_news%>%select(Species_Name,new_weight_ln,actual_fish_presence)
# Display ANOVA Table
summary(anova_model)
# Tukey's HSD Test
tukey_result <- TukeyHSD(anova_model)
print(tukey_result)

# Plot Tukey's HSD
plot(tukey_result)

wide_df <- pivot_wider(zoop_body_size_news, names_from = actual_fish_presence, values_from = new_weight_ln)

# Paired t-test
# Perform paired t-test for each species

test_results <- zoop_body_size_news %>%
  group_by(Species_Name) %>%
  filter(sum(actual_fish_presence == "Yes") > 1 & sum(actual_fish_presence == "No") > 1) %>%  # Ensure enough data points
  summarise(
    t_test = list(t.test(new_weight_ln ~ actual_fish_presence)),
    p_value = sapply(t_test, function(x) x$p.value),
    t_statistic = sapply(t_test, function(x) x$statistic),
    conf_low = sapply(t_test, function(x) x$conf.int[1]),
    conf_high = sapply(t_test, function(x) x$conf.int[2])
  )%>% filter(p_value < 0.05)

# Print the extracted t-test results
test_results

library(ggplot2)
library(dplyr)
library(multcomp)
# Fit a GLM for each species
# Fit a GLM for each species
glm_results <- zoop_body_size_news %>%
  group_by(Species_Name) %>%
  do(model = glm(new_weight_ln ~ actual_fish_presence, data = ., family = gaussian()))

# Extract GLM summaries (coefficients and p-values)
glm_summaries <- glm_results %>%
  summarise(
    Species_Name = first(Species_Name),
    p_value = summary(model)$coefficients[2, 4]  # p-value for fish presence
  )

# Visualize with boxplots or violin plots
ggplot(zoop_body_size_news, aes(x = actual_fish_presence, y = new_weight_ln, fill = actual_fish_presence)) +
  geom_boxplot() +
  facet_wrap(~ Species_Name) +
  theme_minimal() +
  labs(title = "Boxplot of Species Weight by Fish Presence")

# Check residuals for GLM assumptions
residuals_df <- glm_results %>%
  rowwise() %>%
  mutate(
    fitted_values = list(fitted(model)),
    residuals = list(resid(model))
  ) %>%
  unnest(cols = c(fitted_values, residuals))

# 1. Histogram of Residuals
ggplot(residuals_df, aes(x = residuals)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ Species_Name) +
  theme_minimal() +
  labs(title = "Histogram of Residuals for Each Species")

# 2. Residuals vs. Fitted Plot
ggplot(residuals_df, aes(x = fitted_values, y = residuals)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Species_Name) +
  theme_minimal() +
  labs(title = "Residuals vs. Fitted Values Plot")

# 3. Q-Q Plot for Normality of Residuals
ggplot(residuals_df, aes(sample = residuals)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~ Species_Name) +
  theme_minimal() +
  labs(title = "Q-Q Plot for Residuals")

# Print GLM summaries
print(glm_summaries)




# Extracting and reporting only the key statistics
glm_results_reduced <- zoop_body_size_news %>%
  group_by(Species_Name) %>%
  do({
    model <- glm(new_weight_ln ~ actual_fish_presence, data = ., family = gaussian)
    
    summary_model <- summary(model)
    p_value <- coef(summary_model)[2, 4]  # p-value for actual_fish_presence
    coef_estimate <- coef(summary_model)[2, 1]  # Coefficient estimate
    conf_int <- confint(model)[2, ]  # Confidence interval for the coefficient
    r_squared <- 1 - (summary_model$deviance / summary_model$null.deviance)  # Pseudo R²
    n_obs <- nrow(.)  # Number of observations
    aic <- AIC(model)  # AIC
    
    # Return a simplified results table
    tibble(Species_Name = unique(.$Species_Name),
           p_value = p_value,
           coef_estimate = coef_estimate,
           conf_int_lower = conf_int[1],
           conf_int_upper = conf_int[2],
           r_squared = r_squared,
           n_obs = n_obs,
           aic = aic)
  })

# View reduced results
glm_results_reduced



