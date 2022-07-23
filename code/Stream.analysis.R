#Title: The influence of non-native fish on stream macroinvertebrate and lake zooplankton communities along elevational gradients
#Authors: Matthew Green, David Herbst, and Kurt Anderson

#Analysis of Stream Communities

###########################################
#Load Packages
Packages <- c("tidyverse","betapart","otuSummary", "ggplot2", "vegan", "reshape2","reshape", "betareg","adespatial", "sf", "mapview", "viridis", "FD","multcomp","semPlot","lavaan", "performance")
lapply(Packages, library, character.only = TRUE)

setwd("~/Dropbox/Manuscipts/Chapter 3/Chapter3/data/")

###########################################
#Load Data
#species<-read.csv(file = "sp.density.update.12.28.19.csv")
species<-read.csv(file = "sp.density.update.12.28.19_traits.csv")
summary(species)

traits<-read.csv("data/Full_full_fn_trait.csv")
summary(traits)

envs <-read.csv("dave.matt.env.full.12.29.19.csv")
summary(envs)

species_biomass<-species%>%rownames_to_column("Site")%>%pivot_longer(-Site,names_to = "Taxon", values_to= "Density")%>%mutate(Number_ind=Density*.412)
species_biomass

all_biomass_density<-left_join(species_biomass,traits, by="Taxon")%>%
  mutate(Biomass=Number_ind*M)%>%dplyr::select(c(Taxon,Site,Density,Number_ind,M,Biomass,Feed_prim_abbrev))
str(all_biomass_density)
###########################################
#organize data
envs<-envs%>%mutate(Euc.dist.lake=log(Euc.dist.lake+1),River.dist.lake=log(River.dist.lake+1),Head.river.dist=log(Head.river.dist+1))
#envs<-envs%>%dplyr::select(c(Site,O.NET))
species_all<-species%>%pivot_longer(-Site, names_to="Taxon", values_to="abundance")%>%filter(abundance>0)
traits_mass<-traits%>%dplyr::rename(Body_mass_mg=M)%>%dplyr::select(c(Taxon, Body_mass_mg))
species_mass_data<-left_join(species_all,traits_mass, by="Taxon")

species_mass_data_env<-left_join(species_mass_data,envs, by="Site")%>%filter(O.NET != "YOUNG")%>%
  filter(Site !=	"Outlet.11007.fishless.2003")%>%
  filter(Site !=	"Outlet.11007.fishless.2004")%>%
  filter(Site !=	"Outlet.10487.trt.2003")%>%
  filter(Site !=	"Outlet.10487.trt.2004")%>%
  filter(Site !=	"Outlet.10477.trt.2003")%>%
  filter(Site !=	"Outlet.10477.trt.2004")%>%
  filter(Site !=	"Outlet.10494.trt.2012")%>%
  filter(Site !=	"	Outlet.Vidette.below.2003")%>%
  filter(Site !=	"	Outlet.Vidette.below.2004")%>%
  filter(Site !=	"	Outlet.Vidette.below.20012")%>%filter(Elevation >3200)

species_mass_data_env

########################################################################################################################
species_mass_data_env%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Fish,y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Log Density")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~Taxon, scales="free")+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Taxon,y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Log Density")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~Taxon, scales="free")+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=reorder(Taxon, Body_mass_mg, FUN = mean),y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Log Density")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Fish,y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  xlab("Fish Presence")+ylab("Log Density")+
  facet_wrap(~Taxon)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

###########
#Just Select Species that occur in fish and fishless sites
species_mass_data_env_sum<-species_mass_data_env%>%
  group_by(Taxon,Fish)%>%
  summarise(Total_density=sum(abundance))%>%
  pivot_wider(names_from = Taxon,values_from = Total_density)%>%
  filter(Fish !="NA")%>%
  pivot_longer(cols=Acari:Zapada,names_to = "Taxon")

no_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="No")%>% filter(is.na(value))

yes_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="Yes")%>% filter(is.na(value))

all_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(is.na(value))

all_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  left_join(traits, by="Taxon")%>%
  filter(is.na(value))%>%
  filter(Taxon !=" Cultus", Taxon != "Sialis", Taxon !="Tabanus", Taxon !="Pleuroceridae", Taxon !="Despaxia", Taxon !="Cenocorixa", Taxon != "Chrysops",
         Taxon !="Capniidae", Taxon !="Rhabdomastix", Taxon != "Hyrda", Taxon !="Monohelea", Taxon != "Forcipomyia", Taxon != "Hemerodromia", Taxon !="Euhirudinea",
         Taxon != "Cultus", Taxon !="Wormaldia")#, Taxon !="Arctopsyche", Taxon !="Chyranda", Taxon !="Neothremma")

#write.csv(all_species_mass_data_env_sum,"all_species_mass_data_env_sum.csv")

fig1d<-all_species_mass_data_env_sum%>%
  ggplot(aes(x=Fish,y=log(M+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  ggtitle("b)") +
  xlab("Species Absent From")+ylab("Macroinvertebrate Body Mass (mg)")+
  #facet_grid(~Feed_prim_abbrev)+
  scale_fill_viridis(discrete = TRUE,name = "Species Absent From", labels = c("Fishless", "Fish"))+
  scale_x_discrete(labels=c("Fishless", "Fish"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.5, 0.9),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))



yes_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="Yes")%>% 
  na.omit()

no_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="No")%>% 
  na.omit()

species_mass_data_env_sum_2<-species_mass_data_env_sum%>%
  na.omit()%>%
  left_join(traits_mass, by="Taxon")

species_mass_data_env_sum_2%>%
  ggplot(aes(x=Fish,y=log(Body_mass_mg+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  ggtitle("d)") +
  xlab("Species Absent From")+ylab("Macroinvertebrate Body Mass (mg)")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("Fishless", "Fish"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.75, 0.85),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
#write.csv(all_species_mass_data_env_sum,"all_species_mass_data_env_sum.csv")

#No:Aedes,Alloperla,Allotrichoma,Blephariceridae,Brachycentrus,Callibaetis,Calliperla,Capniidae,Centroptilum,Cheumatopsyche,Chrysops
#yes: Arctocorisa,Arctopsyche,Atrichopogon,Capniidae,Cenocorixa,Ceratopogon

species_mass_data_env_filter<-species_mass_data_env%>%
  filter(Taxon !="Aedes" ,Taxon !="Alloperla" ,Taxon !="Allotrichoma" ,Taxon !="Blephariceridae" ,Taxon !="Brachycentrus"
         ,Taxon !="Callibaetis" ,Taxon !="Calliperla" ,Taxon !="Capniidae" ,Taxon !="Centroptilum" ,Taxon !="Cheumatopsyche" ,Taxon !="Chrysops" ,
         Taxon !="Claassenia" ,Taxon !="Cleptelmis" ,Taxon !="Cultus" ,Taxon !="Despaxia" ,Taxon !="Deuterophlebia" ,Taxon !="Ephemerella" ,
         Taxon !="Euhirudinea" ,Taxon !="Forcipomyia" ,Taxon !="Glutops" ,Taxon !="Haploperla" ,Taxon !="Hemerodromia" ,Taxon !="Hesperoperla" ,
         Taxon !="Hexatoma",Taxon !="Hirudinea" ,Taxon !="Hydra" ,Taxon !="Hydropsyche" ,Taxon !="Lepidostoma" ,Taxon !="Limnophila" ,
         Taxon !="Malenka" ,Taxon !="Megarcys" ,Taxon !="Monophilus" ,Taxon !="Narpus" ,Taxon !="Nemertea" ,Taxon !="Ochrotrichia" ,Taxon !="Oreodytes" ,
         Taxon !="Orohermes" ,Taxon !="Pedicia" ,Taxon !="Perlinodes" ,Taxon !="Planorbidae" ,Taxon !="Pleuroceridae" ,Taxon !="Polycentropus" ,Taxon !="Rhabdomastix" ,
         Taxon !="Rhithrogena" ,Taxon !="Rhizelmis" ,Taxon !="Sialis" ,Taxon !="Siphlonurus" ,Taxon !="Skwala" ,Taxon !="Stictotarsus" ,Taxon !="Tabanus" ,
         Taxon !="Arctocorisa" ,Taxon !="Arctopsyche" ,Taxon !="Atrichopogon" ,Taxon !="Cenocorixa" ,Taxon !="Ceratopogon" ,Taxon !="Chrysops" ,Taxon !="Chyranda" ,
         Taxon !="Culiseta" ,Taxon !="Cultus" ,Taxon !="Despaxia" ,Taxon !="Dolichopus" ,Taxon !="Graptocorixa" ,Taxon !="Hemerodromia" ,Taxon !="Hexatoma" ,
         Taxon !="Hydroporus" ,Taxon !="Lednia" ,Taxon !="Limnophila" ,Taxon !="Limonia" ,Taxon !="Metacnephia" ,Taxon !="Monohelea" ,Taxon !="Neothremma" ,
         Taxon !="Nixe" ,Taxon !="Ormosia" ,Taxon !="Pleuroceridae" ,Taxon !="Rhabdomastix" ,Taxon !="Sciara" ,Taxon !="Soyedina" ,Taxon !="Stictotarsus" ,
         Taxon !="Tabanus", Taxon !="Wormaldia", Taxon !="Wiedeman", Taxon !="Helodon", Taxon !="Tipula", Taxon !="Optioservus")

fig1b<-species_mass_data_env_filter%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=reorder(Taxon, Body_mass_mg, FUN = mean),y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  ggtitle("b)") +
  xlab("Macroinvertebrate Taxa")+ylab("Macroinvertebrate Log Density + 1")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.9, 0.85),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

species_mass_data_env_filter%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Fish,y=log(abundance+1),fill=as.factor(Fish)))+
  geom_boxplot()+
  xlab("Taxon")+ylab("Log Density")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~Taxon)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env_filter%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Elevation,y=log(abundance+1),color=as.factor(Fish)))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Taxon")+ylab("Log Density")+
  scale_color_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  facet_wrap(~Taxon,scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env_filter%>%
  filter(Fish !="NA")%>%
  ggplot(aes(x=Elevation,y=log(abundance+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Taxon")+ylab("Log Density")+
  facet_wrap(~Taxon,scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())



#Change in density by body size
species_mass_data_env$Fish<-as.factor(species_mass_data_env$Fish)
species_mass_data_env_sum<-species_mass_data_env%>%
  group_by(Taxon,Fish, Body_mass_mg)%>%
  summarise(Mean_density=mean(log(abundance+1)))%>%
  filter(Fish !="NA")%>%
  pivot_wider(names_from = Fish,values_from = Mean_density)%>%
  replace(is.na(.), 0)%>% 
  mutate(change_density=No-Yes, change_1_density=Yes-No,relative_change=Yes/No,abs.change=abs(No-Yes),
         new_change=(Yes-No)/(1/2*(Yes+No)))%>%
  left_join(traits_FFG, by="Taxon")

traits_FFG<-traits%>%dplyr::select(c(Taxon, Feed_prim_abbrev))
species_mass_data_env_sum%>%
  filter(change_density>0)%>%
  ggplot(aes(x=log(Body_mass_mg+1),y=change_density))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("Change Macro Density")+xlab("Maco Body Size")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

new.fig1b<-species_mass_data_env_sum%>%
  filter(relative_change != "NA")%>%
  filter(relative_change > 0.000001)%>%
  filter(relative_change < 5)%>%
  ggplot(aes(x=log(Body_mass_mg+1),y=relative_change))+
  geom_point()+
  ggtitle("b)") +
  geom_smooth(method = "lm",linetype="dashed" )+
  geom_hline(yintercept=1, linetype='dotted', col = 'black')+
  ylab("Relative Change Macroinvertebrate Density")+xlab("Macroinvertebrate Log Body Mass (mg)")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env_sum%>%
  filter(relative_change != "NA")%>%
  ggplot(aes(x=Elevation,y=relative_change))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_hline(yintercept=1, linetype='dotted', col = 'black')+
  ylab("Relative Change Macroinvertebrate Density")+xlab("Macroinvertebrate Body Mass")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

dog<-lm(relative_change~log(Body_mass_mg+1),species_mass_data_env_sum)
summary(dog)

species_mass_data_env_sum%>%
  filter(relative_change != "NA")%>%
  ggplot(aes(x=reorder(Taxon, Body_mass_mg, FUN = median),y=relative_change))+
  geom_boxplot()+
  geom_hline(yintercept=1, linetype='dotted', col = 'black')+
  #facet_wrap(~actual_fish_presence, scales="free")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  ylab("Relative Change Macroinvertebrate Density")+xlab("Macroinvertebrate Species")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())


#All sites Occupancy
species_mass_data_env_filter_2<-species_mass_data_env%>%
  pivot_wider(names_from = Fish,values_from = abundance)%>%
  dplyr::mutate(Fishless.Occupancy=if_else(No>0, 1,0),Fish.Occupancy=if_else(Yes>0, 1,0))%>%
  replace(is.na(.), 0)%>% 
  group_by(Taxon)%>%
  dplyr::summarise(Fish.total.occupancy=sum(Fish.Occupancy),Fishless.total.occupancy=sum(Fishless.Occupancy),n=n(),#Fishless.total.occupancy+Fish.total.occupancy,
                   Yes=Fish.total.occupancy/62, No=Fishless.total.occupancy/62)%>%
  pivot_longer(cols=Yes:No,names_to = "Fish", values_to="occupancy")%>%
  left_join(traits_mass, by="Taxon")


new.prop.b2<-species_mass_data_env_filter_2%>%
  ggplot(aes(x=reorder(Taxon, +Body_mass_mg),y=occupancy, fill=Fish, group=Fish))+
  geom_col( size = 0.5, position='dodge')+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  xlab("Macroinvertebrate Taxa")+ylab("Proportion of Stream Sites Occupied")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env_filter_2%>%
  ggplot(aes(x=Fish,y=occupancy, fill=Fish, group=Fish))+
  geom_col( size = 0.5, position='dodge')+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  facet_wrap(~Taxon,scales="free")+
  xlab("Macroinvertebrate Taxa")+ylab("Proportion of Stream Sites Occupied")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())
species_mass_data_env_filter_1<-species_mass_data_env%>%
  pivot_wider(names_from = Fish,values_from = abundance)%>%
  dplyr::mutate(Fishless.Occupancy=if_else(No>0, 1,0),Fish.Occupancy=if_else(Yes>0, 1,0))%>%
  replace(is.na(.), 0)%>% 
  group_by(Taxon)%>%
  dplyr::summarise(Fish.total.occupancy=sum(Fish.Occupancy),Fishless.total.occupancy=sum(Fishless.Occupancy),n=n(),#Fishless.total.occupancy+Fish.total.occupancy,
                   Yes=Fish.total.occupancy/n, No=Fishless.total.occupancy/n)%>%
  pivot_longer(cols=Yes:No,names_to = "Fish", values_to="occupancy")%>%
  left_join(traits_mass, by="Taxon")


new.prop.b<-species_mass_data_env_filter_1%>%
  ggplot(aes(x=reorder(Taxon, +Body_mass_mg),y=occupancy, fill=Fish, group=Fish))+
  geom_col( size = 0.5, position='dodge')+
  ggtitle("d)") +
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence")+
  xlab("Macroinvertebrate Taxa")+ylab("Proportion of Stream Sites Occupied")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")
############################################################################################################################################
#Diversity
species<-read.csv(file = "sp.density.update.12.28.19_traits.csv", row.names = 1)
species<-species%>%dplyr::select(-c(Chironomidae,Acari,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))

diversity<-species%>%
  transmute(N0=rowSums(species > 0),H= diversity(species),N1 =exp(H),N2 =diversity(species, "inv"),J= H/log(N0),E10= (N1/N0),E20= (N2/N0),Com.Size=log(rowSums(species)+1),betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD  )

#Add Community Biomass
all_biomass_density_div<-all_biomass_density%>%filter(Density>0)%>%
  dplyr::select(-c(Density,M,Number_ind))%>%
  group_by(Site)%>%
  dplyr::summarize(Biomass=sum(Biomass))%>%
  mutate(Biomass=log(Biomass+1))

diversity.data<-diversity%>%rownames_to_column("Site")

all<-diversity.data%>%left_join(envs,by="Site")%>%left_join(all_biomass_density_div,by="Site")%>%filter(Elevation >3200)
diversity.env<-all%>%filter(Fish!="NA")%>%filter(O.NET!="YOUNG")%>%
  filter(Site !=	"Outlet.11007.fishless.2003")%>%
  filter(Site !=	"Outlet.11007.fishless.2004")%>%
  filter(Site !=	"Outlet.10487.trt.2003")%>%
  filter(Site !=	"Outlet.10487.trt.2004")%>%
  filter(Site !=	"Outlet.10477.trt.2003")%>%
  filter(Site !=	"Outlet.10477.trt.2004")%>%
  filter(Site !=	"Outlet.10494.trt.2012")%>%
  filter(Site !=	"	Outlet.Vidette.below.2003")%>%
  filter(Site !=	"	Outlet.Vidette.below.2004")%>%
  filter(Site !=	"	Outlet.Vidette.below.20012")%>%  filter(Elevation >3200)

diversity.env%>%
  gather(  N0, N1,  Biomass, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=Elevation, y=value, colour=var))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

fig3e<-diversity.env%>%
  ggplot(aes(x=Elevation, y=N1))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  ggtitle("f)") +
  xlab("Elevation (m)")+ylab("Species (Shannon) Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

fig3f<-diversity.env%>%ggplot(aes(x=Elevation, y=betas.LCBD))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  xlab("Elevation (m)")+ylab("Beta-diversity (LCBD)")+
  ggtitle("h)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

fig3g<-diversity.env%>%
  ggplot(aes(x=Elevation, y=Com.Size))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  ggtitle("g)") +
  xlab("Elevation (m)")+ylab("Community Size")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")

fig3h<-diversity.env%>%  ggplot(aes(x=Elevation, y=(Biomass)))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  xlab("Elevation (m)")+ylab("Community Biomass")+
  ggtitle("h)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")


diversity.env%>%
  filter(Elevation >3100 ,Elevation <3500 )%>%
  gather(N1, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=as.factor(Fish), y=value, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig2e<-diversity.env%>%
  ggplot(aes(x=as.factor(Fish), y=N1, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Species (Shannon) Diversity")+
  ggtitle("e)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")+
  theme(legend.position = c(0.70, 0.89),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))


fig2f<-diversity.env%>%
  ggplot(aes(x=as.factor(Fish), y=betas.LCBD, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Beta-Diversity (LCBD)")+
  ggtitle("g)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2g<-diversity.env%>%
  ggplot(aes(x=as.factor(Fish), y=Com.Size, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Community Size")+
  ggtitle("g)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

fig2h<-diversity.env%>%
  ggplot(aes(x=as.factor(Fish), y=Biomass, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+ylab("Community Biomass")+
  ggtitle("h)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")


diversity.env%>%
  filter(Elevation >3100 ,Elevation <3500 )%>%
  gather(N0, N1,  Biomass, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=Elevation, y=value, colour=as.factor(Fish)))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

diversity.env%>%
  filter(Elevation >3100 ,Elevation <3500 )%>%
  gather(N0, N1,  E10, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=Elevation, y=value, colour=as.factor(Fish)))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Elevation (m)")+
  facet_grid(var~Fish, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())


env_div_av<-diversity.env%>%
  filter(!is.na(N0))%>%
  group_by(Fish)%>%
  summarise(mean_n1=mean(N1), meanN0=mean(N0), mean.beta=mean(betas.LCBD), mean.Com=mean(Com.Size))

mod<-glm(N0~Fish, family=poisson(link="log"),diversity.env)
summary(mod)

mod<-glm(N1~Fish, family=gaussian(link="identity"),diversity.env)
summary(mod)

mod<-glm(Com.Size~Fish, family=gaussian(link="identity"),diversity.env)
summary(mod)

mod<-betareg(betas.LCBD~Fish,diversity.env)
summary(mod)

diversity.env%>%
  filter(Elevation >3100 ,Elevation <3500 )%>%
  gather( N1, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=as.factor(Fish), y=value, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

############################################################################################################################################
#NMDS
species_env<-species%>%rownames_to_column("Site")%>%left_join(envs, by="Site")%>%filter(Fish!="NA")
species_2<-species_env%>%filter(Elevation>3200)%>%column_to_rownames("Site")%>%dplyr::select(c(Acentrella:Zapada))

sset.seed(99)
species<-species_2
dune.rel<-decostand(species,"total") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.rel, k=2, try=1000) #NMDS code
dune.nmds
stressplot(dune.nmds) #this tells us if our plot is going to work, and it looks good

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
#text(dune.nmds$species[,1], dune.nmds$species[,2], rownames(dune.nmds$species), cex=0.7, col ="black")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) 
ordihull(dune.nmds, species_env$Fish, display="sites", label=F,lwd=2, col=c("blue","orange"))
#ordihull(dune.nmds, groups=sp_abund_env$lake_drainage_name, draw="polygon", label=T)
ordisurf(dune.nmds, species_env$Elevation, prioirty=,labcex=0.9, add = T,col="forestgreen")

#PERMANOVA analysis-Whats driving variation we see above?
adonis2(dune.bray ~ species_env$Fish+species_env$Elevation+species_env$O.NET, permutations = 999, method = "bray")
betad <- betadiver(dune.rel, "z")
adonis(betad ~ species_env$Fish+species_env$Elevation, data=species, perm=200)

adonis2(dune.bray ~ species_env$Fish, permutations = 99, method = "bray")
betad <- betadiver(dune.rel, "z")
adonis(betad ~ species_env$Fish, data=species, perm=200)

dune.envfit <- envfit(dune.nmds, env = species_env$Elevation, perm = 999) #standard envfit
dune.envfit
env.scores.dune <- as.data.frame(scores(dune.envfit, display = "vectors")) #extracts relevant scores from envifit
env.data = cbind(species_env$Elevation)
mds.data.envfit = envfit(dune.nmds, env.data)

plot(mds.data.envfit, col = "black", labels = c( "Elevation"), lwd = 2)

mod <- betadisper(dune.bray, species_env$Fish)
anova(mod)
print(mod)
permutest(mod)
boxplot(mod)

############################################################################################################################################
#CWM

species.traits<-read.csv(file = "sp.density.update.12.28.19_traits.csv", row.names = 1)
tr.traits<-read.csv("Full_full_fn_trait.csv")

species.traits<-species.traits%>%dplyr::select(-c(Chironomidae))


row.traits<-tr.traits%>%filter(Taxon !="Benthic.producers" , Taxon !="Pelagic.producers", Taxon !="Chironomidae")
traitsy<-tr.traits%>%dplyr::rename(Body_mass_mg=M)%>%filter(Taxon !="Benthic.producers" , Taxon !="Pelagic.producers", Taxon !="Chironomidae")%>%dplyr::select(c(Body_mass_mg))
str(traitsy)
#rownames(traitsy)<-rownames(row.traits)
rownames(traitsy)<-row.traits$Taxonomic_name

tres_bm = dbFD(traitsy,species.traits, corr = ("lingoes"),
               stand.FRic = TRUE, calc.FDiv = TRUE)

cwm=tres_bm$CWM
head(cwm)
cwm<-cwm%>%rownames_to_column("Site")
datas<-left_join(envs,cwm, by="Site")

datasz<-datas%>%
  #filter(River.dist.lake>0.1)%>%
  filter(O.NET != "YOUNG")%>%
  filter(Elevation >3200)
  #filter(Elevation >3150 ,Elevation <3500 )
  #filter(Head.river.dist>2.3)%>%
  #filter(Fish != "NA")%>%
  #filter(O.NET != "CASCADE",O.NET != "YOUNG")%>%
  #filter(Body_mass_mg<15)

datasz%>%
  filter(Elevation>2790)%>%
  ggplot(aes(x = Elevation, y = Body_mass_mg))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #geom_smooth(method = "lm")+
  stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  theme_bw()+ylab("CWM")+xlab("Elevation (m)")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig6b<-datasz%>%
  ggplot(aes(x = Elevation, y = Body_mass_mg, colour=Fish))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  scale_color_viridis_d(name = "Fish Presence",labels = c("No", "Yes"))+
  ylab("Macroinvertebrate CWM")+xlab("Elevation (m)")+  ggtitle("b)") +
  facet_grid(~Fish,scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.75, 0.82),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

sdatasz%>%
  gather(Head.river.dist,River.dist.lake,key = "var", value = "value") %>% #PC2,PC3,PC4,River.dist.lake,
  ggplot(aes(x = value, y = Body_mass_mg))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #geom_smooth(method = "lm")+
  stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

fig5b<-datasz%>%
  #filter(Fish != "NA")%>%
  ggplot(aes(x = as.factor(Fish), y = Body_mass_mg, fill=as.factor(Fish)))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("no", "yes"))+
  xlab("Fish Presence")+ylab("Macroinvertebrate CWM")+
  labs(fill='Fish Presence') +
  ggtitle("b)") +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.5, 0.9),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

datasz%>%
  filter(Fish != "NA")%>%
  ggplot(aes(x = as.factor(Fish), y = Body_mass_mg, fill=as.factor(Fish)))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("no", "yes"))+
  xlab("Fish Presence")+ylab("CWM")+
  labs(fill='Fish Presence') +
  ggtitle("b)") +
  facet_grid(~O.NET)+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+
  theme(legend.position = c(0.45, 0.9),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

dog<-aov(Body_mass_mg~Fish,datasz)
summary(dog)
r2(dog)
#########################################################################
#Visualize influence of fish and elevation


############################################################################################################################################

#Fish  By Elevation
supp.c<-stream.elev%>%
  filter(Elevation >3200)%>%
  ggplot(aes(x = as.factor(Fish), y = Elevation, fill=as.factor(Fish)))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("c)") +
  ylab("Elevation (m)")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")+ theme(legend.position = "none")

stream.elev<-diversity.env%>%
  filter(Elevation >3200)#,Elevation <3700 )


dog<-aov(Elevation~Fish,stream.elev)
summary(dog)


supp.d<-stream.elev%>%
  ggplot(aes(x = as.factor(Fish), fill=as.factor(Fish)))+ 
  geom_bar(stat="count")+
  stat_count(geom = "text", colour = "black", size = 3.5,
             aes(label = ..count..),position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("d)") +
  ylab("Number of Stream Sites")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")+ theme(legend.position = "none")

species_mass_data_env_no<-species_mass_data_env%>%
  group_by(Taxon)%>%
  count()

diversity.env%>%
  filter(Elevation >3200)%>%
  ggplot(aes(x = as.factor(Fish), fill=as.factor(Fish)))+ 
  geom_bar(stat="count")+
  stat_count(geom = "text", colour = "black", size = 3.5,
             aes(label = ..count..),position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("b)") +
  ylab("Number of Stream Sites")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")+
  theme(legend.position = c(0.25, 0.85),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

############################################################################################################################################
#Betapart
species_beta<-new_species_envs%>%dplyr::select(c(Site, Acari:Zapada))%>%column_to_rownames("Site")
env_beta<-new_species_envs%>%dplyr::select(c(Site:Fish))
dune.rel<-decostand(species_beta,"total")


sum.dist<-beta.multi.abund(species_beta,index.family="bray")

dist<-beta.pair.abund(species_beta,index.family="bray")
dist.bal<-dist$beta.bray.bal
dist.gra<-dist$beta.bray.gra
dist.bray<-dist$beta.bray

dist.bray.df<-matrixConvert(dist.bray, colname = c("site1", "site2", "dist.bray"))
dist.bal.df<-matrixConvert(dist.bal, colname = c("site1", "site2", "dist.bal"))
dist.gra.df<-matrixConvert(dist.gra, colname = c("site1", "site2", "dist.gra"))

site.fish.1<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site1="Site")%>%dplyr::rename(fish1="Fish")
site.fish.2<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site2="Site")%>%dplyr::rename(fish2="Fish")
net.fish.1<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network1="O.NET")%>%dplyr::rename(site1="Site")
net.fish.2<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network2="O.NET")%>%dplyr::rename(site2="Site")

env_beta_dist<-env_beta%>%column_to_rownames("Site")%>%dplyr::select(c(Elevation,Lon,Lat))
spatial.dist<-vegdist(env_beta_dist, "euclidean")
dist.df<-matrixConvert(spatial.dist, colname = c("site1", "site2", "spatial.dist"))

cwm_dist<-datas%>%column_to_rownames("Site")%>%dplyr::select(c(Body_mass_mg))
cwm.dist<-vegdist(cwm_dist, "euclidean")
cwm.df<-matrixConvert(cwm.dist, colname = c("site1", "site2", "cwm.dist"))

stream.df<-dist.bray.df%>%left_join(dist.bal.df,  by=c("site1", "site2"))%>%left_join(dist.gra.df,  by=c("site1", "site2"))%>%
  left_join(site.fish.1, by="site1")%>%left_join(site.fish.2, by="site2")%>%left_join(net.fish.1, by="site1")%>%left_join(net.fish.2, by="site2")%>%left_join(dist.df,by=c("site1", "site2"))%>%left_join(cwm.df,by=c("site1", "site2"))

stream.df.filter<-stream.df%>%mutate(Fish.turnover=if_else(fish1== "Yes" & fish2== "Yes", 
                                                           "Fish2fish",if_else(fish1== "Yes" & fish2== "No",
                                                                               "Fish2nofish",if_else(fish1== "No" & fish2== "No",
                                                                                                     "Nofish2noFish","NoFish2fish"))))%>%
  mutate(Network=if_else(network1== network2,"Same", "Different"))

stream.df.filter_mod<-stream.df.filter%>%pivot_longer(cols=dist.bray:dist.gra, names_to="beta_div_comp", values_to ="beta_value" )
stream.df.filter_mod%>%group_by(beta_div_comp)%>%summarise(mean.beta=mean(beta_value))

stream.df.filter%>%
  gather(dist.bray,dist.bal,dist.gra,key = "var", value = "value")%>% 
  ggplot(aes(x=Fish.turnover, y=value,fill=Fish.turnover))+
  geom_boxplot()+
  #facet_grid(~Fish.turnover)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("Spatial Dissimilarity")+
  ggtitle("c)") +
  ylab("Beta Diversity")+
  facet_grid(~var)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

stream.df.filter%>%
  filter(Network=="Same")%>%
  gather(dist.bray,dist.bal,dist.gra,key = "var", value = "value")%>% 
  ggplot(aes(x=cwm.dist, y=value,colour=var))+
  geom_point()+geom_smooth(method="lm")+
  scale_colour_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("Spatial Dissimilarity")+
  ggtitle("c)") +
  ylab("Beta Diversity")+
  facet_grid(Fish.turnover~var)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

stream.df.filtered<-stream.df.filter%>%
  pivot_longer(dist.bray:dist.gra, names_to ="comp" ,values_to="beta" )

beta.b<-stream.df.filtered%>%
  ggplot(aes(x=comp, y=beta,fill=comp))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Turnover")+
  scale_x_discrete(limits = c("dist.bal","dist.gra","dist.bray"), labels=c(expression(beta[bal]),expression(beta[gra]),expression(beta[bray])))+
  xlab("Turnover and Nestedness Components")+
  ggtitle("b)") +
 # ylab("Macroinvertebrate Beta Diversity")+
  labs(y=(("Macroinvertebrate \u03B2-diversity")))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

#beta presence/abs
dune.rel[dune.rel > 0] <- 1
dist<-beta.multi(dune.rel,index.family="sorensen")


dist<-beta.pair(dune.rel,index.family="sorensen")
dist.sim<-dist$beta.sim
dist.sne<-dist$beta.sne
dist.sor<-dist$beta.sor

dist.sim.df<-matrixConvert(dist.sim, colname = c("site1", "site2", "dist.sim"))
dist.sne.df<-matrixConvert(dist.sne, colname = c("site1", "site2", "dist.sne"))
dist.sor.df<-matrixConvert(dist.sor, colname = c("site1", "site2", "dist.sor"))

site.fish.1<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site1="Site")%>%dplyr::rename(fish1="Fish")
site.fish.2<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site2="Site")%>%dplyr::rename(fish2="Fish")
net.fish.1<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network1="O.NET")%>%dplyr::rename(site1="Site")
net.fish.2<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network2="O.NET")%>%dplyr::rename(site2="Site")


stream.df<-dist.sim.df%>%left_join(dist.sne.df,  by=c("site1", "site2"))%>%left_join(dist.sor.df,  by=c("site1", "site2"))%>%
  left_join(site.fish.1, by="site1")%>%left_join(site.fish.2, by="site2")%>%left_join(net.fish.1, by="site1")%>%left_join(net.fish.2, by="site2")

stream.df.filter<-stream.df%>%mutate(Fish.turnover=if_else(fish1== "Yes" & fish2== "Yes", 
                                                           "Fish2fish",if_else(fish1== "Yes" & fish2== "No",
                                                                               "Fish2nofish",if_else(fish1== "No" & fish2== "No",
                                                                                                     "Nofish2noFish","NoFish2fish"))))%>%
  mutate(Network=if_else(network1== network2,"Same", "Different"))

stream.df.filter_mod<-stream.df.filter%>%pivot_longer(cols=dist.sim:dist.sor, names_to="beta_div_comp", values_to ="beta_value" )
stream.df.filter_mod%>%group_by(beta_div_comp)%>%summarise(mean.beta=mean(beta_value))


stream.df.filter%>%
  gather(dist.sim,dist.sne,dist.sor,key = "var", value = "value")%>% 
  ggplot(aes(x=Fish.turnover, y=value,fill=Fish.turnover))+
  geom_boxplot()+
  #facet_grid(~Fish.turnover)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("Spatial Dissimilarity")+
  ggtitle("c)") +
  ylab("Beta Diversity")+
  facet_grid(~var)+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

beta.b<-stream.df.filter%>%
  gather(dist.sim,dist.sne,dist.sor,key = "var", value = "value")%>% 
  ggplot(aes(x=var, y=value,fill=var))+
  geom_boxplot()+
  #facet_grid(~Fish.turnover)+
  scale_fill_viridis(discrete = TRUE,name = "Fish Turnover")+
  scale_x_discrete(labels=c(expression(beta[sim]),expression(beta[sne]),expression(beta[sor])))+
  xlab("Turnover and Nestedness Components")+
  ggtitle("b)") +
  ylab("Macroinvertebrate Beta Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")
############################################################################################################################################
#Beta
species<-read.csv(file = "data/sp.density.update.12.28.19_traits.csv")
envs <-read.csv("data/dave.matt.env.full.12.29.19.csv")

new_species_envs<-envs%>%left_join(species, by="Site")%>%filter(Elevation >3200)

species_beta<-new_species_envs%>%dplyr::select(c(Site, Acari:Zapada))%>%column_to_rownames("Site")
env_beta<-new_species_envs%>%dplyr::select(c(Site:Fish))


dune.rel<-decostand(species_beta,"total") #standardize community data
dune.bray<-vegdist(dune.rel)
str(species)
str(env_elevation)

#distance<-matrix(c(env_beta$Lat,env_beta$Lon, env_beta$Elevation),nrow=70, ncol=3)
#dist<-dist(distance, method="euclidean")

env_beta_dist<-env_beta%>%column_to_rownames("Site")%>%dplyr::select(c(Elevation,Lon,Lat))
spatial.dist<-vegdist(env_beta_dist, "euclidean")
dist.df<-matrixConvert(spatial.dist, colname = c("site1", "site2", "spatial.dist"))

cwm_dist<-datas%>%column_to_rownames("Site")%>%dplyr::select(c(Body_mass_mg))
cwm.dist<-vegdist(cwm_dist, "euclidean")
cwm.df<-matrixConvert(cwm.dist, colname = c("site1", "site2", "cwm.dist"))

plot(dist,dune.bray, pch=16 ,xlab="Spatial Dissimilairty (euclidean)",ylab="Macro dissimilarity (bray-curtis)", ylim=c(0,1))

dune.bray.df<-matrixConvert(dune.bray, colname = c("site1", "site2", "cmtny.dist"))
dist.df<-matrixConvert(spatial.dist, colname = c("site1", "site2", "spatial.dist"))
cwm.df<-matrixConvert(cwm.dist, colname = c("site1", "site2", "cwm.dist"))


site.fish.1<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site1="Site")%>%dplyr::rename(fish1="Fish")
site.fish.2<-envs%>%dplyr::select(c(Site,Fish))%>%dplyr::rename(site2="Site")%>%dplyr::rename(fish2="Fish")
net.fish.1<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network1="O.NET")%>%dplyr::rename(site1="Site")
net.fish.2<-envs%>%dplyr::select(c(O.NET,Site))%>%dplyr::rename(network2="O.NET")%>%dplyr::rename(site2="Site")


stream.df<-dune.bray.df%>%left_join(dist.df,  by=c("site1", "site2"))%>%left_join(cwm.df,  by=c("site1", "site2"))%>%
  left_join(site.fish.1, by="site1")%>%left_join(site.fish.2, by="site2")%>%left_join(net.fish.1, by="site1")%>%left_join(net.fish.2, by="site2")

stream.df.filter<-stream.df%>%mutate(Fish.turnover=if_else(fish1== "Yes" & fish2== "Yes", 
                                             "Fish2fish",if_else(fish1== "Yes" & fish2== "No",
                                                              "Fish2nofish",if_else(fish1== "No" & fish2== "No",
                                                                                    "Nofish2noFish","NoFish2fish"))))%>%
  mutate(Network=if_else(network1== network2,"Same", "Different"))%>%
  filter(Network=="Same")
  
  


beta.stream.c<-stream.df.filter%>%
  filter(Network=="Same")%>%
  ggplot(aes(x=spatial.dist, y=cmtny.dist))+ #,colour=Fish.turnover))+
  geom_point()+
  geom_smooth(method="lm")+
  #facet_grid(~Fish.turnover)+
  scale_colour_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("Spatial Dissimilarity")+
  ggtitle("c)") +
  ylab("Beta Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

beta.stream.d<-stream.df.filter%>%
  filter(Network=="Same")%>%
  ggplot(aes(x=Fish.turnover, y=cmtny.dist, fill=Fish.turnover))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("Fish Turnover")+
  ggtitle("d)") +
  ylab("Beta Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

stream.df.filter%>%
  filter(Network=="Same")%>%
  ggplot(aes(x=log(cwm.dist+1), y=cmtny.dist))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~Fish.turnover)+
  scale_colour_viridis(discrete = TRUE,name = "Fish Turnover")+
  xlab("CWM Dissimilarity")+
  ggtitle("c)") +
  ylab("Beta Diversity")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")

dog<-lm(cmtny.dist ~ spatial.dist*Fish.turnover,stream.df.filter)
summary(dog)
############################################################################################################################################