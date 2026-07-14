#Title: The influence of non-native fish on stream macroinvertebrate and lake zooplankton communities along elevational gradients
#Authors: Matthew Green, David Herbst, and Kurt Anderson

#Analysis of Lake Communities


#Link to EDI data portal: https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=577&revision=2
#Package ID: edi.577.2 Cataloging System:https://pasta.edirepository.org.
#Data set title: The Sierra Lakes Inventory Project: Non-Native fish and community composition of lakes and ponds in the Sierra Nevada, California.
########################################################################################################################################################################
#Load Packages
Packages <- c("tidyverse","betareg" ,"ggplot2", "vegan", "reshape2","reshape", "adespatial", "sf", "viridis", "FD","multcomp","semPlot","lavaan", "performance")
# Install missing packages
missing_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages)
}
lapply(Packages, library, character.only = TRUE)

distance_df <- function(dist_obj, colname = c("site1", "site2", "distance")) {
  m <- as.matrix(dist_obj)
  combs <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(
    site1 = rownames(m)[combs[, 1]],
    site2 = rownames(m)[combs[, 2]],
    distance = m[combs],
    stringsAsFactors = FALSE
  ) %>% setNames(colname)
}

#Load all data: Run Script "Download_Slip_Data.R"
setwd("~/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3/")
source("code/Download_Slip_Data.R", echo = TRUE)

setwd("/Users/matthewgreen/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3/data")
av_zoop_body_size_new <- read.csv("length_mass_regress_zoop.csv", row.names = 1) %>%
  dplyr::select(Body_mass_ug)
########################################################################################################################################################################
#Part 1) Organize Species and Environmental Data Set: There are a lot of separate files that need to be cleaned and then merged together

#########################################################################
#1A) Correct Sub-sample and uneven sampling data issues
zoopzz<-dt10%>%
  left_join(dt12, by = "SpeciesID") # link species names with unique ID's

dt8.1<-dt8%>%
  dplyr::select(c(lake_id,survey_date,lake_max_depth,zoo_tow_number,zoo_tow_type,zoo_tow_depth)) #remove unnecessary info

#Correct sub sample issue to make equal comparison among sites
clean_zoopzz<-zoopzz%>%
  left_join(dt8.1, by=c("lake_id", "survey_date"))%>%
  left_join(
    dt13 %>%
      group_by(lake_id, survey_date) %>%
      summarise(sample_vol = first(sample_vol), .groups = "drop"),
    by = c("lake_id", "survey_date")
  )%>%
  group_by(lake_id,survey_date)%>%
  mutate(Max.Subsample=max(Subsample))%>%
  ungroup()%>%
  group_by(lake_id,survey_date,Species_Name)%>%
  mutate(
    Counts = sum(Number)/Max.Subsample,
    # Zooplankton density: scale sampled count up by sample_vol, then divide by filtered volume
    zoop_density = Counts * sample_vol / (zoo_tow_depth * zoo_tow_number),
    Species_Name = str_replace(Species_Name, " ", "_"),
    Species_Name = str_replace(Species_Name, "/", "_")
  )%>%
  filter(Species_Name != "standard_measurement")%>%
  dplyr::select(lake_id, survey_date, Species_Name, zoop_density)%>%
  distinct(lake_id, survey_date, Species_Name, zoop_density)

#Biomass
zoop_body_mass <- av_zoop_body_size_new %>%
  rownames_to_column("Taxon")

clean_zoopz_biomass<-zoopzz%>%
  left_join(dt8.1,by=c("lake_id", "survey_date"))%>%
  group_by(lake_id,survey_date)%>%
  mutate(Max.Subsample=max(Subsample))%>%
  ungroup()%>%
  group_by(lake_id,survey_date,Species_Name)%>%
  mutate(Counts=sum(Number),
         Species_Name = str_replace(Species_Name, " ", "_"), #replace spaces in names with "_"
         Species_Name = str_replace(Species_Name, "/", "_"))%>% #replace / in names with "_")%>%
  dplyr::rename(Taxon=Species_Name)%>%
  left_join(zoop_body_mass, by="Taxon")%>%
  group_by(lake_id,survey_date,Taxon)%>%
  mutate(Biomass=Counts*Body_mass_ug)%>%
  dplyr::select(c(lake_id,survey_date,Taxon,Biomass))


#########################################################################
#1B) Pivot Data Set long to wide format for Species abundance Matrix
site.sp.quad <- cast(clean_zoopzz, lake_id+survey_date ~ Species_Name, value='zoop_density') 
site.sp.quad <- as.data.frame(site.sp.quad)  #set data set as data frame
site.sp.quad[is.na(site.sp.quad)] <- 0       #Replace NA's with 0's


#Species Diversity as Response Varible: calculate multiple species diversity metrics
species<-as.data.frame(site.sp.quad[,3:39])   #Select species and abundances only
diversity<-species%>%
  transmute(N0=rowSums(species > 0),
            H= diversity(species),
            N1 =exp(H),
            N2 =diversity(species, "inv"),
            J= H/log(N0),
            E10= (N1/N0),
            E20= (N2/N0),
            Com.Size=rowSums(species), #Total number individuals per site
            betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD) #LCBD (Beta-diversity or variability among sites)

#Add back site info and treatments to diversity data
local_diversity<-cbind(diversity, site.sp.quad[,1:2])

clean_zoopz_biomass_site<-clean_zoopz_biomass%>%
  filter(Taxon!="standard_measurement")%>%
  filter(Biomass>0)%>%
  group_by(lake_id, survey_date)%>%
  summarise(Sum.Biomass=sum(Biomass))

local_diversity<-local_diversity%>%left_join(clean_zoopz_biomass_site,by=c("lake_id", "survey_date"))

str(local_diversity)
str(clean_zoopz_biomass_site)
#########################################################################
#Part 1C) Organize Environmental Data Set and merge with Species Data

#Add in GPS coordinates from LAKEID variable: Using sieera_lakes arc GIS export to excel
gps<-read.csv("sierra_lakes.csv")

dt5.1<-dt5%>%
  pivot_wider(names_from = littoral_type, names_glue = "{littoral_type}_{.value}",values_from=littoral_amount)

dt6.1<-dt6%>% 
  pivot_wider(names_from = shoreline_type,names_glue = "{shoreline_type}_{.value}", values_from=shoreline_amount)

#Join multiple envrionmental data sets with species data: Link by unique site id and date sampled
env<-dt5.1%>%
  left_join(dt6.1, by= c("lake_id", "survey_date"))%>%
  left_join(dt4, by= "lake_id")%>%
  left_join(dt8, by= c("lake_id", "survey_date"))%>%
  filter(zoo_sample_ind=="Yes")%>%
  dplyr::select(-c(zoo_sample_ind,zoo_sample_time,zoo_tow_number,zoo_tow_type,zoo_tow_number,zoo_tow_depth,
                   benthic_sample_ind,benthic_sample_percent,nbr_benthic_sweeps,
                   lake_fairy_shrimp_ind,lake_shrimp_collection,pool_fairy_shrimp_ind,pool_shrimp_collection,
                   amphib_survey_starttime,amphib_survey_endtime,amphib_survey_duration,amphib_survey_desc,amphib_survey_fish_presence,
                   fish_survey_type, fish_net_location_type,fish_net_set_datetime,fish_net_pull_datetime))%>%
  left_join(gps, by="lake_id")

sp_abund_env<-left_join(site.sp.quad,env, by=c("lake_id", "survey_date"))%>%filter(lake_id!="70534")
summary(sp_abund_env)

#########################################################################
#1D) Isolated individual Species Abundance
#Fish presence  and environmental effects on abundance individuals species
pivot_clean_zoopzz<-clean_zoopzz%>%pivot_wider(names_from = "Species_Name",values_from="zoop_density")

env_abund<-left_join(clean_zoopzz,env, by=c("lake_id","survey_date"))%>%filter(actual_fish_presence=="Yes" |actual_fish_presence=="No")

new_pivot<-pivot_clean_zoopzz%>%
  rownames_to_column("Site")%>%
  pivot_longer(Alona_spp.:Trichotria_spp., names_to = "Taxon", values_to= "Density")%>%
  separate("Site",sep="_" ,into=c("lake_id", "survey_date"))

new_pivot$lake_id<-as.integer(new_pivot$lake_id)

new_pivot_env<-new_pivot%>%
  left_join(env, by=c("lake_id","survey_date"))%>%filter(actual_fish_presence=="Yes" |actual_fish_presence=="No")

str(env)
###################################################################################################################################################################################################################################################################################################################################
#Part 2) Analysis and Visualization


##########################################################################
#2A)MAP Spatial Visualization of species diversity metrics

#mapviewOptions(fgb = FALSE)

#new_env_div<-env_div%>%filter(N0>0)

#Diversity Metrics
#snw_sf <- st_as_sf(new_env_div, coords = c("Lon", "Lat"), crs=4326, remove = FALSE)
#mapview(snw_sf, zcol="betas.LCBD", layer.name="LCBD")
#mapview(snw_sf, zcol="N0", layer.name="Species Richness")
#mapview(snw_sf, zcol="N1", layer.name="Species Diversity")

#fish Presence
#mapview(snw_sf, zcol="actual_fish_presence", layer.name="Fish Presence")
#########################################################################

#2B) Explore Relationships among Diversity as a function of environmental variables: Visualization and Stats

env_div<-left_join(env,local_diversity, by=c("lake_id", "survey_date"))%>%filter(lake_id !="11505" & lake_id !="42219" &lake_id !="71257" &lake_id !="71282" )%>%
  mutate(Com.Size=log(Com.Size+1), Biomass=log(Sum.Biomass+1))%>%
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)

supp.b<-env_div%>%
  ggplot(aes(x=actual_fish_presence, fill=actual_fish_presence))+
  geom_bar(stat="count")+
  stat_count(geom = "text", colour = "black", size = 3.5,
             aes(label = ..count..),position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("b)") +
  ylab("Number of Lake Sites")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

supp.a<-env_div%>%
  ggplot(aes(x = actual_fish_presence, y = lake_elevation_nbr, fill=actual_fish_presence))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("a)") +
  ylab("Elevation (m)")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

env_div%>%
  group_by(actual_fish_presence)%>%count()

env_div%>%
  gather(N0,  N1,  Biomass, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=lake_elevation_nbr, y=value, colour=var))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis_d()+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

env_div%>%
  gather(N0, N1,  Biomass, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=lake_elevation_nbr, y=value, colour=actual_fish_presence))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis_d()+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

env_div%>%
  gather(N0, N1,  E10, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=lake_elevation_nbr, y=value, colour=actual_fish_presence))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis_d()+
  xlab("Elevation (m)")+
  facet_grid(var~actual_fish_presence, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())


fig3a<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=N1))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("c)") +
  xlab("Elevation (m)")+ylab("Species (Shannon) Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig3b<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=betas.LCBD))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("d)") +
  xlab("Elevation (m)")+ylab("Beta-diversity (LCBD)")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig3c<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=Com.Size))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("c)") +
  xlab("Elevation (m)")+ylab("Community Size")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig3d<-env_div%>%
  ggplot(aes(x=lake_elevation_nbr, y=Biomass))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  ggtitle("d)") +
  xlab("Elevation (m)")+ylab("Community Biomass")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

dog<-lm(N0~lake_elevation_nbr, data=env_div)
dog<-glm(N0~lake_elevation_nbr,family = poisson(link = "log"), data=env_div)
summary(dog)
r2(dog)

dog<-lm(N1~lake_elevation_nbr, data=env_div)
summary(dog)
dog<-glm(N1~lake_elevation_nbr, family = gaussian(link="identity"),data=env_div)
summary(dog)

dog<-lm(betas.LCBD~lake_elevation_nbr, data=env_div)
summary(dog)

dog<-lm(Com.Size~lake_elevation_nbr, data=env_div)
summary(dog)

dog<-glm(betas.LCBD~actual_fish_presence, family = gaussian(link="identity"),data=env_div)
dog<-glm(N0~actual_fish_presence, family = poisson(link = "log"),data=env_div)
summary(dog)
r2(dog)


#########################################################################
#Diversity as a function of fish presence/absence
env_div%>%
  gather(N0, N1,  Biomass, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=actual_fish_presence, y=value, fill=actual_fish_presence))+
  geom_boxplot()+
  xlab("Fish Presence")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

env_div_av<-env_div%>%
  filter(!is.na(N0))%>%
  group_by(actual_fish_presence)%>%
  summarise(mean_n1=mean(N1), meanN0=mean(N0), mean.beta=mean(betas.LCBD), mean.Com=mean(Com.Size))

mod<-glm(N0~actual_fish_presence, family=poisson(link="log"),env_div)
summary(mod)

mod<-glm(N1~actual_fish_presence, family=gaussian(link="identity"),env_div)
summary(mod)

mod<-glm(Com.Size~actual_fish_presence, family=gaussian(link="identity"),env_div)
summary(mod)

mod<-betareg(betas.LCBD~actual_fish_presence,env_div)
summary(mod)

env_div%>%
  filter(lake_elevation_nbr>2799, lake_elevation_nbr<3602)%>%
  gather( N1, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=actual_fish_presence, y=value, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig2a<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=N1, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Species (Shannon) Diversity")+
  ggtitle("a)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2b<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=betas.LCBD, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Beta-Diversity (LCBD)")+
  ggtitle("b)") +
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2c<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=Com.Size, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+ 
  xlab("Fish Presence")+ylab("Community Size")+
  ggtitle("c)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2d<-env_div%>%
  ggplot(aes(x=actual_fish_presence, y=Biomass, fill=as.factor(actual_fish_presence)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Community Biomass")+
  ggtitle("d)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

#remove unwanted columns for analysis due to missing data
env_divz<-env_div%>%
  dplyr::select(c(N0, H, N1, N2, E10, E20, J, Com.Size, betas.LCBD,actual_fish_presence, lake_drainage_name,Jurisdicti))

env_divz<-as.data.frame(env_divz)


#########################################################################
#2C) Effects on individual species: Analysis and Viz

env_abund%>%
  ggplot(aes(x=actual_fish_presence,y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Log Density")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~Species_Name, scales="free")+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.border = element_blank(),panel.background = element_blank())

env_abund%>%
  ggplot(aes(x=Species_Name,y=(zoop_density),fill=actual_fish_presence))+
  geom_bar(stat = "identity")+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Density")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

env_abund%>%
  ggplot(aes(x=reorder(Species_Name, zoop_density, FUN = median),y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Density")+xlab("Zooplankton Species")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

env_abund%>%
  ggplot(aes(x=reorder(Species_Name, zoop_density, FUN = median),y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Density")+xlab("Zooplankton Species")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

env_abund%>%
  ggplot(aes(x=lake_elevation_nbr,y=log(zoop_density+1),color=Species_Name))+
  geom_point()+
  geom_smooth(method = "lm")+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_color_viridis_d()+
  ylab("Zooplankton Density")+
  facet_wrap(~Species_Name, scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank(), legend.position = "none")


env_abund$Species_Name<-as.factor(env_abund$Species_Name)
clean_zoopzz$lake_id<-as.integer(clean_zoopzz$lake_id)
env$lake_id<-as.integer(env$lake_id)

env_abundz<-clean_zoopzz%>%
  pivot_wider(names_from = "Species_Name",values_from="zoop_density")%>%
  dplyr::select(-c(Lecane_spp.,Monostyla_spp.,Simocephalus_serrulatus))%>%
  left_join(env, by=c("lake_id", "survey_date"))

av_zoop_body_size_news<-av_zoop_body_size_new%>%rownames_to_column("Taxon")
env_abundzzz<-env_abund%>%dplyr::rename(Taxon=Species_Name)%>%left_join(av_zoop_body_size_news, by="Taxon")%>%
  dplyr::rename(Mean_body_size_mm=Body_mass_ug)%>%filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)

env_abundzzz%>%
  filter()%>%
  ggplot(aes(x=reorder(Taxon, Mean_body_size_mm, FUN = median),y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Density")+xlab("Zooplankton Species")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

env_abundzzz%>%
  ggplot(aes(x=lake_elevation_nbr,y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_point()+geom_smooth()+
  facet_wrap(~Taxon, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Density")+xlab("Zooplankton Species")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

##########
#Just Select Species that occur in fish and fishless sites
av_zoop_body_size_newer<-av_zoop_body_size_new%>%rownames_to_column("Taxon")

env_abundz_filter<-env_abundzzz%>%
  left_join(av_zoop_body_size_newer, by="Taxon")%>%
  group_by(Taxon,actual_fish_presence,Body_mass_ug)%>%
  summarise(Mean_density=mean(log(zoop_density+1)))%>%
  pivot_wider(names_from=actual_fish_presence,values_from =Mean_density )%>%
  replace(is.na(.), 0)%>% 
  mutate(change_density=No-Yes, change_1_density=Yes-No,relative_change=Yes/No,abs.change=abs(No-Yes),
         new_change=(Yes-No)/(.5*(Yes+No)))


env_abundz_filter%>%
  ggplot(aes(x=log(Body_mass_ug+1),y=change_density))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_hline(yintercept=0, linetype='dotted', col = 'black')+
  ylab("Change Zooplankton Density")+xlab("Zooplankton Body Size")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

dog<-lm(change_density~log(Body_mass_ug+1),env_abundz_filter)
summary(dog)

env_abundz_filter%>%
  filter(change_1_density>0)%>%
  ggplot(aes(x=log(Body_mass_ug+1),y=new_change))+
  geom_point()+
  geom_hline(yintercept=0, linetype='dotted', col = 'black')+
  geom_smooth(method = "lm")+
  ylab("Cahnge Zooplankton Density")+xlab("Zooplankton Body Size")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())
env_abundz_filter_plot <- env_abundz_filter %>%
  pivot_longer(cols = c(No, Yes), names_to = "actual_fish_presence", values_to = "Mean_density")
env_abundz_filter_plot %>%
  filter(!is.na(relative_change)) %>%
  ggplot(aes(x=reorder(Taxon, Body_mass_ug, FUN = median), y = relative_change)) +
  geom_boxplot()+
  geom_hline(yintercept=1, linetype='dotted', col = 'black')+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Relative Change Zooplankton Density")+xlab("Zooplankton Species")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

new.fig1a<-env_abundz_filter_plot %>%
  filter(!is.na(relative_change))%>%
  ggplot(aes(x= log(Body_mass_ug+1),y=relative_change))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_hline(yintercept=1, linetype='dotted', col = 'black')+
  #ggtitle("a)") +
  ylab("Relative Change in Zooplankton Density")+xlab(expression("Zooplankton Body Mass ( " * mu * "g/L, Log"[10] * ")"))+ 
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

new<-env_abundzzz%>%
  left_join(av_zoop_body_size_newer, by="Taxon")%>%
  group_by(Taxon,actual_fish_presence,Body_mass_ug)%>%
  summarise(Mean_density=mean(log(zoop_density+1)), .groups = "drop")%>%
  pivot_wider(names_from=actual_fish_presence, values_from=Mean_density)

env_abundz_filter_bm<-new%>%
  filter(Taxon %in% c("Alona_spp.", "Alonella_excisa", "Ascomorpha_spp.", "Asplanchna_spp.", "Polyarthra_vulgaris", "Polyphemus_pediculus", "Trichotria_spp."))%>%
  pivot_longer(cols=Yes:No, names_to = "Fish", values_to = "Mean_density")



env_abundz_filtered<-env_abundzzz%>%
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500) %>%
  filter(HA>=0.5) %>%
  filter(lake_max_depth>3) %>%
  filter(!Taxon %in% c("Collotheca_spp.", "Eurycercus_lamellatus", "Lecane_spp.", "Monostyla_spp.", "Polyarthra_vulgaris", "Simocephalus_serrulatus"))

fig1a<-env_abundz_filtered%>%
  ggplot(aes(x=reorder(Taxon, Mean_body_size_mm, FUN = median),y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  ggtitle("a)") +
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Log Density + 1")+xlab("Zooplankton Taxa")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position = c(0.9, 0.85),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))


env_abundz_filtered%>%
  ggplot(aes(x=actual_fish_presence,y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  ggtitle("a)") +
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Log Density + 1")+xlab("Zooplankton Taxa")+
  facet_wrap(~Taxon)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())
env_abundzzz_new<-env_abundzzz%>%
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)%>%
  dplyr::select(c(lake_id,survey_date,Taxon,zoop_density,actual_fish_presence,lake_elevation_nbr,Mean_body_size_mm))%>%
  filter(Taxon !="nauplii")

env_abundzzz_new%>%
  ggplot(aes(x=actual_fish_presence,y=log(zoop_density+1),fill=actual_fish_presence))+
  geom_boxplot()+
  ggtitle("a)") +
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Zooplankton Log Density + 1")+xlab("Zooplankton Taxa")+
  facet_wrap(~Taxon)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

##########################################################################################
#Occupancy 

env_abundz_body<-env_abundz_filter%>%
  dplyr::select(c(Taxon,Body_mass_ug))

#1) Occupnacy of Lakes its found in
env_abundzzz_new_1<-env_abundzzz_new%>%
  pivot_wider(names_from = actual_fish_presence,values_from = zoop_density)%>%
  dplyr::mutate(Fishless.Occupancy=if_else(No>0, 1,0),Fish.Occupancy=if_else(Yes>0, 1,0))%>%
  replace(is.na(.), 0)%>% 
  group_by(Taxon)%>%
  dplyr::summarise(n=n(),Fish.total.occupancy=sum(Fish.Occupancy),Fishless.total.occupancy=sum(Fishless.Occupancy),
                   Yes=Fish.total.occupancy/n, No=Fishless.total.occupancy/n)%>%
  pivot_longer(cols=Yes:No,names_to = "Fish", values_to="occupancy")%>%
  left_join(env_abundz_body, by="Taxon")



new.prop.a<-env_abundzzz_new_1%>%
  ggplot(aes(x=reorder(Taxon, +Body_mass_ug),y=occupancy, fill=Fish, group=Fish))+
  geom_col( size = 0.5, position='dodge')+
  ggtitle("b)") +
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  xlab("Zooplankton Taxa")+ylab("Proportion of Lakes Occupied")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

env_abundzzz_new_1%>%
  ggplot(aes(x=Fish,y=occupancy, fill=Fish))+
  geom_col()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  xlab("Zooplankton Taxa")+ylab("Proportion of Lakes Occupied")+
  facet_wrap(~Taxon)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

env_abundzzz_new_1%>%
  filter(n>8)%>%
  ggplot(aes(x=log(Body_mass_ug+1),y=occupancy, colour=Fish))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_colour_viridis(discrete = TRUE,name = "Fish Presence")+
  facet_grid(~Fish)+
  xlab("Zooplankton Body Size")+ylab("Proportion of Lakes Occupied")+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                           panel.border = element_blank(),panel.background = element_blank())


#2) Occupancy of all lakes
unique(env_abundzzz_new$lake_id)

env_abundzzz_new_2<-env_abundzzz_new%>%
  pivot_wider(names_from = actual_fish_presence,values_from = zoop_density)%>%
  dplyr::mutate(Fishless.Occupancy=if_else(No>0, 1,0),Fish.Occupancy=if_else(Yes>0, 1,0))%>%
  replace(is.na(.), 0)%>% 
  group_by(Taxon)%>%
  dplyr::summarise(n=n(),Fish.total.occupancy=sum(Fish.Occupancy),Fishless.total.occupancy=sum(Fishless.Occupancy),
                   Yes=Fish.total.occupancy/602, No=Fishless.total.occupancy/602, .groups = "drop")%>%
  pivot_longer(cols=Yes:No,names_to = "Fish", values_to="occupancy")%>%
  left_join(env_abundz_body, by="Taxon")



new.prop.a2<-env_abundzzz_new_2%>%
  ggplot(aes(x=reorder(Taxon, +Body_mass_ug),y=occupancy, fill=Fish, group=Fish))+
  geom_col( size = 0.5, position='dodge')+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  xlab("Zooplankton Taxa")+ylab("Percentage of Lakes Occupied")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())
##########################################################################################
#Analysis: 
#remove unwanted columns for analysis due to missing data
env_abundzz<-env_abundz%>%
  unite("id_date",lake_id:survey_date)%>%
  remove_rownames()%>%
  column_to_rownames('id_date')%>%
  dplyr::select(c(Leptodiaptomus_signicauda,nauplii,Daphnia_melanica,Daphnia_dentifera,Conochilus_unicornis,Keratella_spp.,                       
                  Keratella_quadrata,Hesperodiaptomus_shoshone,Polyarthra_dolichoptera,Cyclopoida,Lepadella_spp.,Chydorus_sphaericus,
                  Holopedium_gibberum,Bosmina_longirostris,Trichocerca_capucina,Harpacticoida,Kellicottia_spp.,Polyarthra_vulgaris,Asplanchna_spp.,
                  Trichotria_spp.,Hesperodiaptomus_eiseni,Alona_spp.,Ceriodaphnia_laticaudata,Alonella_excisa,Notholca_spp.,Synchaeta_spp.,Collotheca_spp.,
                  Ascomorpha_spp.,Scapholeberis_mucronata,Diaphanosoma_brachyurum, Chaoborus_trivitattus,Polyphemus_pediculus, Daphnia_pulex_pulicaria,
                  Eurycercus_lamellatus, actual_fish_presence))

env_abundzz<-as.data.frame(env_abundzz)

regional.names<-env_abundzz %>%
  dplyr::select(Leptodiaptomus_signicauda,nauplii,Daphnia_melanica,Daphnia_dentifera,Conochilus_unicornis,Keratella_spp.,                       
                Keratella_quadrata,Hesperodiaptomus_shoshone,Polyarthra_dolichoptera,Cyclopoida,Lepadella_spp.,Chydorus_sphaericus,
                Holopedium_gibberum,Bosmina_longirostris,Trichocerca_capucina,Harpacticoida,Kellicottia_spp.,Polyarthra_vulgaris,Asplanchna_spp.,
                Trichotria_spp.,Hesperodiaptomus_eiseni,Alona_spp.,Ceriodaphnia_laticaudata,Alonella_excisa,Notholca_spp.,Synchaeta_spp.,Collotheca_spp.,
                Ascomorpha_spp.,Scapholeberis_mucronata,Diaphanosoma_brachyurum, Chaoborus_trivitattus,Polyphemus_pediculus, Daphnia_pulex_pulicaria,
                Eurycercus_lamellatus)


#########################################################################
#2D) Body size effects of fish: Community Weighted Mean (CWM)
# Load zooplankton body mass data for CWM calculations
av_zoop_body_size_new <- read.csv("length_mass_regress_zoop.csv", row.names = 1) %>%
  dplyr::select(-Mean_body_size_mm)

pivot_clean_zoopzz<-clean_zoopzz%>%pivot_wider(names_from = "Species_Name",values_from="zoop_density")%>%
  unite("id_date",lake_id:survey_date)%>%
  dplyr::select(-c(Lecane_spp.,Monostyla_spp.,Simocephalus_serrulatus))%>%
  remove_rownames()%>%
  column_to_rownames('id_date')%>% 
  replace(is.na(.), 0)%>%
  dplyr::select(sort(names(.)))

trait_species <- rownames(av_zoop_body_size_new)
comm_species <- intersect(colnames(pivot_clean_zoopzz), trait_species)

if (length(comm_species) == 0) {
  stop("error: no shared species between trait and abundance matrices")
}

extra_abund_species <- setdiff(colnames(pivot_clean_zoopzz), trait_species)
if (length(extra_abund_species) > 0) {
  message("Dropping ", length(extra_abund_species), " abundance species with no trait data.")
}

pivot_clean_zoopzz <- pivot_clean_zoopzz[, comm_species, drop = FALSE]
av_zoop_body_size_new <- av_zoop_body_size_new[comm_species, , drop = FALSE]

zero_total_species <- which(colSums(pivot_clean_zoopzz) == 0)
if (length(zero_total_species) > 0) {
  message("Removing ", length(zero_total_species), " zero-total-abundance species before CWM.")
  pivot_clean_zoopzz <- pivot_clean_zoopzz[, -zero_total_species, drop = FALSE]
  av_zoop_body_size_new <- av_zoop_body_size_new[colnames(pivot_clean_zoopzz), , drop = FALSE]
}

zero_sum_sites <- which(rowSums(pivot_clean_zoopzz) == 0)
if (length(zero_sum_sites) > 0) {
  message("Removing ", length(zero_sum_sites), " zero-sum communities before CWM.")
  pivot_clean_zoopzz <- pivot_clean_zoopzz[-zero_sum_sites, , drop = FALSE]
}

#Calculate CWM
if (dim(pivot_clean_zoopzz)[2] != dim(av_zoop_body_size_new)[1]) {
  stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")
}

tres_bm = dbFD(av_zoop_body_size_new, pivot_clean_zoopzz, corr = ("lingoes"),
               stand.FRic = TRUE, calc.FDiv = TRUE)

#Combine CWM with env and organize data
cwm=tres_bm$CWM
cwm<-cwm%>%
  rownames_to_column("id_date")%>%
  separate("id_date",sep="_" ,into=c("lake_id", "survey_date"))%>%
  dplyr::rename(CWM=Body_mass_ug) #Body_mass_ug or mean_body_size

cwm$lake_id<-as.integer(cwm$lake_id)
env_cwm<-left_join(env,cwm, by=c("lake_id", "survey_date"))%>%dplyr::filter(lake_id !="11505" & lake_id !="42219" &lake_id !="71257" &lake_id !="71282")%>%
  filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%filter(HA>=0.5)%>%filter(lake_max_depth>3)

str(env_cwm)
str(env)
# view(env_cwm)

env$lake_elevation_nbr
#########################################################################
#Visualize influence of fish and elevation

env_cwm%>%
  ggplot(aes(x=lake_drainage_name,y=CWM, fill=actual_fish_presence))+
  geom_boxplot()+
  xlab("Fish Presence")+
  scale_fill_viridis(discrete = TRUE,name = "Lake Network")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig5a<-env_cwm%>%
  #filter(lake_elevation_nbr>2799, lake_elevation_nbr<3602)%>%
  ggplot(aes(x=actual_fish_presence,y=CWM*.01, fill=actual_fish_presence))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Zooplankton CWM")+
  #ggtitle("a)") +
  #geom_signif(comparisons = list(c("No", "Yes")), map_signif_level=TRUE)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")


dog<-aov(CWM~actual_fish_presence,env_cwm)
summary(dog)
r2(dog)

env_cwm%>%
  filter(lake_elevation_nbr>1800, lake_elevation_nbr<3500)%>%
  count(actual_fish_presence)

fig6a<-env_cwm%>%
  ggplot(aes(x=lake_elevation_nbr,y=CWM,colour=actual_fish_presence))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_viridis_d(name = "Fish Presence",labels = c("No", "Yes"))+
  facet_grid(~actual_fish_presence, scales="free")+
  xlab("Elevation (m)")+  
  #ggtitle("a)") +
  ylab("Zooplankton CWM")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.8, 0.9),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

env_cwm%>%
  ggplot(aes(x=log(lake_area_nbr+1),y=CWM,colour=actual_fish_presence))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_viridis_d(name = "Fish Presence")+
  facet_grid(~actual_fish_presence)+
  xlab("Elevation (m)")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fishless<-env_cwm%>%
  filter(actual_fish_presence=="No")

fish<-env_cwm%>%
  filter(actual_fish_presence=="Yes")

mod1<-lm(CWM~lake_elevation_nbr,  data=fishless)
summary(mod1)

mod1<-lm(CWM~lake_elevation_nbr,  data=fish)
summary(mod1)

mod0<-glm(N0~actual_fish_presence*lake_elevation_nbr*lake_max_depth,family=poisson(link ="log" ),env_div)
summary(mod0)
r2(mod0)

#Analaysis: GLM's
mod1<-glm(CWM~lake_elevation_nbr*actual_fish_presence,family=gaussian(link = "identity"),  data=env_cwm)
mod2<-glm(CWM~lake_elevation_nbr,family=gaussian(link = "identity"),  data=env_cwm)
mod3<-glm(CWM~actual_fish_presence,family=gaussian(link = "identity"),  data=env_cwm)
null<-glm(CWM~1,family=gaussian(link = "identity"),  data=env_cwm)

reported.table2 <- bbmle::AICtab(mod1,mod2,mod3, null,weights = TRUE, sort = FALSE)
reported.table2
r2(mod3)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR1

#Conclusion: Lake Elevation and fish presence interactively explain variation in zooplankton community body size

#########################################################################
#2E) Multivariate Ordinations (NMDS)

set.seed(99)

nmds_species <- sp_abund_env %>%
  dplyr::select(c(Alona_spp.:Trichotria_spp.)) %>%
  mutate(across(everything(), as.numeric)) %>%
  replace(is.na(.), 0)

keep <- rowSums(nmds_species) > 0
if (sum(keep) < 3) {
  message("Skipping NMDS: fewer than 3 lake communities with species data.")
} else {
  nmds_species <- nmds_species[keep, , drop = FALSE]
  env_nmds <- sp_abund_env[keep, ]

  dune.rel <- decostand(nmds_species, "total")
  if (any(!is.finite(as.matrix(dune.rel)))) {
    stop("Non-finite values found in NMDS community matrix after standardization.")
  }

  dune.bray <- vegdist(dune.rel)
  dune.nmds <- metaMDS(dune.rel, k = 2, try = 10, trace = FALSE)
  stressplot(dune.nmds)

  plot(dune.nmds, type = "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
  points(dune.nmds$points[, 1], dune.nmds$points[, 2], pch = 1)
  ordihull(dune.nmds, env_nmds$actual_fish_presence, display = "sites", label = FALSE, lwd = 2, col = c("blue", "orange"))
  ordisurf(dune.nmds, env_nmds$lake_elevation_nbr, labcex = 0.9, add = TRUE, col = "forestgreen")

    adonis2(dune.bray ~ actual_fish_presence + lake_elevation_nbr + lake_drainage_name,
      data = env_nmds,
      permutations = 99,
      method = "bray")
    dune.betad <- betadiver(dune.rel, "z")
    adonis2(dune.betad ~ actual_fish_presence + lake_elevation_nbr,
      data = env_nmds,
      permutations = 200)

    adonis2(dune.bray ~ actual_fish_presence,
      data = env_nmds,
      permutations = 99,
      method = "bray")

  dune.envfit <- envfit(dune.nmds, env = env_nmds$lake_elevation_nbr, perm = 999)
  env.scores.dune <- as.data.frame(scores(dune.envfit, display = "vectors"))
  env.data <- cbind(env_nmds$lake_elevation_nbr)
  mds.data.envfit <- envfit(dune.nmds, env.data)
  plot(mds.data.envfit, col = "black", labels = c("Elevation"), lwd = 2)

  mod <- betadisper(dune.bray, env_nmds$actual_fish_presence)
  anova(mod)
  print(mod)
  permutest(mod)
  boxplot(mod)
}


####################################################################################################################################################################################################################################################################################################
#Fish by Elevation
lake.elev.fish<-env_div%>%
  ggplot(aes(x = actual_fish_presence, y = lake_elevation_nbr, fill=actual_fish_presence))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("a)") +
  ylab("Elevation (m)")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")
####################################################################################################################################################################################################################################################################################################
library(betapart)

# Beta diversity for lakes filtered by elevation, habitat, and depth
sp_abund_env_filter <- sp_abund_env %>%
  filter(lake_elevation_nbr > 1800,
         lake_elevation_nbr < 3500,
         HA >= 0.5,
         lake_max_depth > 3) %>%
  unite("site_id", lake_id, survey_date, remove = FALSE)

species <- sp_abund_env_filter %>%
  dplyr::select(c(Alona_spp.:Trichotria_spp.)) %>%
  mutate(across(everything(), as.numeric)) %>%
  replace(is.na(.), 0)

keep <- rowSums(species) > 0
species <- species[keep, , drop = FALSE]
sp_abund_env_filter <- sp_abund_env_filter[keep, ]
rownames(species) <- sp_abund_env_filter$site_id

if (nrow(species) < 3) {
  stop("Not enough lake communities with species data for beta diversity analysis.")
}

lake.dist <- beta.pair.abund(species, index.family = "bray")
lake.dist.bray.df <- distance_df(lake.dist$beta.bray, colname = c("site1", "site2", "dist.bray"))
lake.dist.bal.df <- distance_df(lake.dist$beta.bray.bal, colname = c("site1", "site2", "dist.bal"))
lake.dist.gra.df <- distance_df(lake.dist$beta.bray.gra, colname = c("site1", "site2", "dist.gra"))

spatial.dist <- sp_abund_env_filter %>%
  dplyr::select(site_id, lake_elevation_nbr, Lon, Lat) %>%
  as.data.frame()
rownames(spatial.dist) <- spatial.dist$site_id
spatial.dist$site_id <- NULL
spatial.dist <- replace(spatial.dist, is.na(spatial.dist), 0)
space.df <- distance_df(vegdist(as.matrix(spatial.dist), "euclidean"), colname = c("site1", "site2", "dist.space"))

cwm.dist <- env_cwm %>%
  unite("site_id", lake_id, survey_date, remove = FALSE) %>%
  filter(site_id %in% sp_abund_env_filter$site_id) %>%
  dplyr::select(site_id, CWM) %>%
  replace(is.na(.), 0) %>%
  as.data.frame()
rownames(cwm.dist) <- cwm.dist$site_id
cwm.dist$site_id <- NULL
cwm.df <- distance_df(vegdist(as.matrix(cwm.dist), "euclidean"), colname = c("site1", "site2", "dist.cwm"))

site.fish.1 <- sp_abund_env_filter %>%
  dplyr::select(site_id, actual_fish_presence) %>%
  dplyr::rename(site1 = site_id, fish1 = actual_fish_presence)
site.fish.2 <- sp_abund_env_filter %>%
  dplyr::select(site_id, actual_fish_presence) %>%
  dplyr::rename(site2 = site_id, fish2 = actual_fish_presence)

net.fish.1 <- sp_abund_env_filter %>%
  dplyr::select(lake_drainage_name, site_id) %>%
  dplyr::rename(network1 = lake_drainage_name, site1 = site_id)
net.fish.2 <- sp_abund_env_filter %>%
  dplyr::select(lake_drainage_name, site_id) %>%
  dplyr::rename(network2 = lake_drainage_name, site2 = site_id)

lake.df <- lake.dist.bray.df %>%
  left_join(lake.dist.bal.df, by = c("site1", "site2")) %>%
  left_join(lake.dist.gra.df, by = c("site1", "site2")) %>%
  left_join(space.df, by = c("site1", "site2")) %>%
  left_join(cwm.df, by = c("site1", "site2")) %>%
  left_join(site.fish.1, by = "site1") %>%
  left_join(site.fish.2, by = "site2") %>%
  left_join(net.fish.1, by = "site1") %>%
  left_join(net.fish.2, by = "site2")

lake.df.filter <- lake.df %>%
  mutate(
    Fish.turnover = case_when(
      fish1 == "Yes" & fish2 == "Yes" ~ "Fish2fish",
      fish1 == "Yes" & fish2 == "No" ~ "Fish2nofish",
      fish1 == "No" & fish2 == "No" ~ "Nofish2noFish",
      TRUE ~ "NoFish2fish"
    ),
    Network = if_else(network1 == network2, "Same", "Different")
  )

lake.df.filtered <- lake.df.filter %>%
  pivot_longer(cols = dist.bray:dist.gra, names_to = "comp", values_to = "beta")

beta.a <- lake.df.filtered %>%
  ggplot(aes(x = comp, y = beta, fill = comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Beta component") +
  scale_x_discrete(limits = c("dist.bal", "dist.gra", "dist.bray"), labels = c(expression(beta[bal]), expression(beta[gra]), expression(beta[bray]))) +
  xlab("Turnover and nestedness components") +
  ylab(expression("Zooplankton " * beta * "-diversity")) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), legend.position = "none")

beta.stream.a <- lake.df.filter %>%
  filter(Network == "Same") %>%
  ggplot(aes(x = dist.space, y = dist.cwm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Spatial dissimilarity") +
  ylab("CWM distance") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), legend.position = "none")

beta.stream.b <- lake.df.filter %>%
  filter(Network == "Same") %>%
  ggplot(aes(x = Fish.turnover, y = dist.cwm, fill = Fish.turnover)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, name = "Fish turnover") +
  xlab("Fish turnover") +
  ylab("CWM distance") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), legend.position = "none")

dog <- lm(dist.cwm ~ dist.space * Fish.turnover, lake.df.filter)
summary(dog)

dog <- lm(dist.cwm ~ dist.space, lake.df.filter)
summary(dog)

dog <- glm(dist.cwm ~ Fish.turnover, family = gaussian(link = "identity"), lake.df.filter)
summary(dog)
r2(dog)
