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
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)
  
species_fish_presence <- zoop_body_size_news %>%
  group_by(Species_Name) %>%
  summarise(
    fish_presence = paste(unique(actual_fish_presence), collapse = ", ")
  )


zoop_body_size_news%>%
  filter(Species_Name != "Ascomorpha_spp.")%>%  filter(Species_Name != "Asplanchna_spp.")%>%  filter(Species_Name != "Daphnia_pulex_pulicaria")%>%
  filter(Species_Name != "Polyarthra_vulgaris")%>%  filter(Species_Name != "Polyphemus_pediculus")%>%  filter(Species_Name != "Synchaeta_spp.")%>% filter(Species_Name !="Trichotria_spp.")%>%
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
  filter(lake_elevation_nbr>3200)%>%
  ggplot(aes(x=Species_Name, y=new_weight_ln,  fill=as.factor(actual_fish_presence )))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Body Weight")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  #ggtitle("b)") +
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

env_abundz_filter<-env_abundzzz%>%
  left_join(av_zoop_body_size_newer, by="Taxon")%>%
  group_by(Taxon,actual_fish_presence,Body_mass_ug)%>%
  summarise(Mean_density=mean(log(zoop_density+1)))%>%
  pivot_wider(names_from=actual_fish_presence,values_from =Mean_density )%>%
  replace(is.na(.), 0)%>% 
  mutate(change_density=No-Yes, change_1_density=Yes-No,relative_change=Yes/No,abs.change=abs(No-Yes),
         new_change=(Yes-No)/(.5*(Yes+No)))


#Changes in body size and turnover
