#Mass Abundance Curves by Site

#Matthew and Katherine
#September 17,2020

#Papers
#Relationships between body size and abundance in ecology

###########################################
#Load Packages
Packages <- c("tidyverse", "ggplot2", "vegan", "reshape2","reshape", "betareg","adespatial", "sf", "mapview", "viridis", "FD","multcomp","semPlot","lavaan", "performance")
lapply(Packages, library, character.only = TRUE)

setwd("~/Dropbox/Manuscipts/L-S Food web/Food-Web-Structure-Sierra/data")

###########################################
#Load Data
#species<-read.csv(file = "sp.density.update.12.28.19.csv")
species<-read.csv(file = "sp.density.update.12.28.19_traits.csv")
summary(species)

traits<-read.csv("Full_full_fn_trait.csv")
summary(traits)

env <-read.csv("dave.matt.env.full.12.29.19.csv")
summary(env)

###########################################
#organize data
env<-env%>%mutate(Euc.dist.lake=log(Euc.dist.lake+1),River.dist.lake=log(River.dist.lake+1),Head.river.dist=log(Head.river.dist+1))
envs<-env%>%dplyr::select(c(Site,O.NET))
species_all<-species%>%pivot_longer(-Site, names_to="Taxon", values_to="abundance")%>%filter(abundance>0)
traits_mass<-traits%>%dplyr::rename(Body_mass_mg=M)%>%dplyr::select(c(Taxon, Body_mass_mg))
species_mass_data<-left_join(species_all,traits_mass, by="Taxon")

species_mass_data_env<-left_join(species_mass_data,env, by="Site")%>%filter(O.NET != "YOUNG")%>%
  filter(Site !=	"Outlet.11007.fishless.2003")%>%
  filter(Site !=	"Outlet.11007.fishless.2004")%>%
  filter(Site !=	"Outlet.10487.trt.2003")%>%
  filter(Site !=	"Outlet.10487.trt.2004")%>%
  filter(Site !=	"Outlet.10477.trt.2003")%>%
  filter(Site !=	"Outlet.10477.trt.2004")%>%
  filter(Site !=	"Outlet.10494.trt.2012")%>%
  filter(Site !=	"	Outlet.Vidette.below.2003")%>%
  filter(Site !=	"	Outlet.Vidette.below.2004")%>%
  filter(Site !=	"	Outlet.Vidette.below.20012")

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

###########
#Just Select Species that occur in fish and fishless sites
species_mass_data_env_sum<-species_mass_data_env%>%
  group_by(Taxon,Fish)%>%
  summarise(Total_density=sum(abundance))%>%
  pivot_wider(names_from = Taxon,values_from = Total_density)%>%
  filter(Fish !="NA")%>%
  pivot_longer(cols=Acentrella:Zapada,names_to = "Taxon")

no_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="No")%>% filter(is.na(value))

yes_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(Fish=="Yes")%>% filter(is.na(value))

all_species_mass_data_env_sum<-species_mass_data_env_sum%>%
  filter(is.na(value))

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
  xlab("Macroinvetebrate Taxa")+ylab("Macroinvetebrate Log Density + 1")+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

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

model1<-glm(log(abundance+1)~Fish,family=gaussian(link="identity"),species_mass_data_env_filter)
model2<-glm(log(abundance+1)~Elevation,family=gaussian(link="identity"),species_mass_data_env_filter)
model3<-glm(log(abundance+1)~Body_mass_mg,family=gaussian(link="identity"),species_mass_data_env_filter)
model4<-glm(log(abundance+1)~Fish*Elevation,family=gaussian(link="identity"),species_mass_data_env_filter)
model5<-glm(log(abundance+1)~Fish*Body_mass_mg,family=gaussian(link="identity"),species_mass_data_env_filter)
model6<-glm(log(abundance+1)~Elevation*Body_mass_mg,family=gaussian(link="identity"),species_mass_data_env_filter)
model7<-glm(log(abundance+1)~Fish*Elevation*Body_mass_mg,family=gaussian(link="identity"),species_mass_data_env_filter)
model8<-glm(log(abundance+1)~1,family=gaussian(link="identity"),species_mass_data_env_filter)
reported.table2 <- bbmle::AICtab(model1,model2,model3,model4,model5,model6,model7,model8,weights = TRUE, sort = FALSE)
reported.table2
r2(model7)
pseudoR0 <- ((model7$null.deviance-model7$deviance)/model7$null.deviance)
pseudoR0

#Doroneuria
species_mass_data_env_filter_Doroneuria<-species_mass_data_env_filter%>%
  filter(Taxon=="Doroneuria")

model1<-glm(log(abundance+1)~Fish,family=gaussian(link="identity"),species_mass_data_env_filter_Doroneuria)
model2<-glm(log(abundance+1)~Elevation,family=gaussian(link="identity"),species_mass_data_env_filter_Doroneuria)
model4<-glm(log(abundance+1)~Fish*Elevation,family=gaussian(link="identity"),species_mass_data_env_filter_Doroneuria)
model8<-glm(log(abundance+1)~1,family=gaussian(link="identity"),species_mass_data_env_filter_Doroneuria)
reported.table2 <- bbmle::AICtab(model1,model2,model4,model8,weights = TRUE, sort = FALSE)
reported.table2
r2(model1)

pseudoR0 <- ((model1$null.deviance-model1$deviance)/model1$null.deviance)
pseudoR0


#Change in density by body size
species_mass_data_env$Fish<-as.factor(species_mass_data_env$Fish)
species_mass_data_env_sum<-species_mass_data_env%>%
  group_by(Taxon,Fish, Body_mass_mg)%>%
  summarise(Mean_density=mean(abundance))%>%
  filter(Fish !="NA")%>%
  pivot_wider(names_from = Fish,values_from = Mean_density)%>%
  mutate(change_density=No-Yes, change_1_density=Yes-No)


species_mass_data_env_sum%>%
  filter(change_density>0)%>%
  ggplot(aes(x=log(Body_mass_mg+1),y=log(change_density+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("Change Macro Density")+xlab("Maco Body Size")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

species_mass_data_env_sum%>%
  filter(change_1_density>0)%>%
  ggplot(aes(x=log(Body_mass_mg+1),y=log(change_1_density+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("Change Macro Density")+xlab("Maco Body Size")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),panel.background = element_blank())

############################################################################################################################################
#Diversity
species<-read.csv(file = "sp.density.update.12.28.19_traits.csv", row.names = 1)
species<-species%>%dplyr::select(-c(Chironomidae,Acari,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))

diversity<-species%>%
  transmute(N0=rowSums(species > 0),H= diversity(species),N1 =exp(H),N2 =diversity(species, "inv"),J= H/log(N0),E10= (N1/N0),E20= (N2/N0),Com.Size=log(rowSums(species)+1),betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD  )
#Calculate Beta Diversity (LCBD) for each network
combine<-cbind(species,env$O.NET)

kern<-subset(combine, `env$O.NET`=="KERN")
kern.var<-kern[, colSums(kern != 0) > 0]
kern.beta<-beta.div(kern.var[1:61], method="hellinger",sqrt.D=TRUE)
kern.beta.comp<-beta.div.comp(kern.var[1:61], coef = "S", quant=T)

casc<-subset(combine, `env$O.NET`=="CASCADE")
casc.var<-casc[, colSums(casc != 0) > 0]
casc.beta<-beta.div(casc.var[1:44], method="hellinger",sqrt.D=TRUE)
casc.beta.comp<-beta.div.comp(casc.var[1:44], coef = "S", quant=T)

evo<-subset(combine, `env$O.NET`=="EVO")
evo.var<-evo[, colSums(evo != 0) > 0]
evo.beta<-beta.div(evo.var[1:36], method="hellinger",sqrt.D=TRUE)
evo.beta.comp<-beta.div.comp(evo.var[1:36], coef = "S", quant=T)

bubb<-subset(combine, `env$O.NET`=="BUBBS")
bubb.var<-bubb[, colSums(bubb != 0) > 0]
bubb.beta<-beta.div(bubb.var[1:83], method="hellinger",nperm=999,sqrt.D=TRUE)
bubb.beta.comp<-beta.div.comp(bubb.var[1:83], coef = "S", quant=T)

young<-subset(combine, `env$O.NET`=="YOUNG")
young.var<-young[, colSums(young != 0) > 0]
young.beta<-beta.div(young.var[1:29], nperm=999,method="hellinger",sqrt.D=TRUE)
young.beta.comp<-beta.div.comp(young.var[1:29], coef = "S", quant=T)

rock<-subset(combine, `env$O.NET`=="ROCK")
rock.var<-rock[, colSums(rock != 0) > 0]
rock.beta<-beta.div(rock.var[1:60], nperm=999,method="hellinger",sqrt.D=TRUE)
rock.beta.comp<-beta.div.comp(rock.var[1:60], coef = "S", quant=T)

betas.LCBD<-c(kern.beta$LCBD,casc.beta$LCBD,evo.beta$LCBD,bubb.beta$LCBD,young.beta$LCBD,rock.beta$LCBD)
betas.LCBD<-as.data.frame(betas.LCBD)
betas.LCBD.data<-betas.LCBD%>%rownames_to_column("Site")
diversity.data<-diversity%>%rownames_to_column("Site")

#all<-diversity.data%>%left_join(env,by="Site")%>%left_join(betas.LCBD.data,by="Site")#%>%filter(Elevation >2790)
all<-diversity.data%>%left_join(env,by="Site")%>%filter(Elevation >2790)
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
  filter(Site !=	"	Outlet.Vidette.below.20012")%>%filter(Elevation >2800)

diversity.env%>%
  filter(Elevation >2800)%>%
  gather(  N0, N1,  E10, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=Elevation, y=value, colour=var))+
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  scale_color_viridis_d()+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),legend.position = "none")
lm<-lm(betas.LCBD~Elevation,diversity.env)
summary(lm)

mod<-betareg(betas.LCBD~Elevation,diversity.env)
mod1<-betareg(betas.LCBD~Fish,diversity.env)
mod2<-betareg(betas.LCBD~Elevation*Fish,diversity.env)
mod3<-betareg(betas.LCBD~1,diversity.env)
reported.table2 <- bbmle::AICtab(mod,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2

mod<-glm(N1~Elevation,family=gaussian(link = "identity"),diversity.env)
mod1<-glm(N1~Fish,family=gaussian(link = "identity"),diversity.env)
mod2<-glm(N1~Elevation*Fish,family=gaussian(link = "identity"),diversity.env)
mod3<-glm(N1~1,family=gaussian(link = "identity"),diversity.env)
reported.table2 <- bbmle::AICtab(mod,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod2)

mod<-glm(Com.Size~Elevation,family=gaussian(link = "identity"),diversity.env)
mod1<-glm(Com.Size~Fish,family=gaussian(link = "identity"),diversity.env)
mod2<-glm(Com.Size~Elevation*Fish,family=gaussian(link = "identity"),diversity.env)
mod3<-glm(Com.Size~1,family=gaussian(link = "identity"),diversity.env)
reported.table2 <- bbmle::AICtab(mod,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod2)

diversity.env%>%
  filter(Elevation >2790)%>%
  gather(N0, N1,  E10, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=as.factor(Fish), y=value, fill=as.factor(Fish)))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())


diversity.env%>%
  filter(Elevation >2790)%>%
  gather(N0, N1,  E10, Com.Size, betas.LCBD,key = "var", value = "value")%>% 
  ggplot(aes(x=Elevation, y=value, colour=as.factor(Fish)))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_color_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Elevation (m)")+
  facet_wrap(~var, scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

diversity.env%>%
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
species_env<-species%>%rownames_to_column("Site")%>%left_join(env, by="Site")%>%filter(Fish!="NA")
species_2<-species_env%>%column_to_rownames("Site")%>%dplyr::select(c(Acentrella:Zapada))
set.seed(99)
species<-species_2
dune.rel<-decostand(species,"total") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.rel, k=2, try=1000) #NMDS code
dune.nmds
stressplot(dune.nmds) #this tells us if our plot is going to work, and it looks good

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
#text(dune.nmds$species[,1], dune.nmds$species[,2], rownames(dune.nmds$species), cex=0.7, col ="black")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) 
ordihull(dune.nmds, groups=as.factor(species_env$Fish), draw="polygon", label=T)
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


row.traits<-tr.traits%>%filter(Taxon !="Benthic.producers" & Taxon !="Pelagic.producers")
traitsy<-tr.traits%>%dplyr::rename(Body_mass_mg=M)%>%filter(Taxon !="Benthic.producers" & Taxon !="Pelagic.producers")%>%dplyr::select(c(Body_mass_mg))
str(traitsy)
#rownames(traitsy)<-rownames(row.traits)
rownames(traitsy)<-row.traits$Taxonomic_name

tres_bm = dbFD(traitsy,species.traits, corr = ("lingoes"),
               stand.FRic = TRUE, calc.FDiv = TRUE)

cwm=tres_bm$CWM
FDis=as.data.frame(tres_bm$FDis)
FEve=as.data.frame(tres_bm$FEve)
FRic=as.data.frame(tres_bm$FRic)

head(cwm)

datas<-cbind(env,cwm,FRic,FDis,FEve)
datasz<-datas%>%
  filter(River.dist.lake>0.1)%>%
  #filter(Head.river.dist>2.3)%>%
  #filter(Fish != "NA")%>%
  filter(O.NET != "YOUNG")%>%
  filter(Body_mass_mg<15)

datasz%>%
  filter(Elevation>2790)%>%
  ggplot(aes(x = Elevation, y = Body_mass_mg))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #geom_smooth(method = "lm")+
  stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  theme_bw()+ylab("CWM")+xlab("Elevation (m)")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

datasz%>%
  filter(Fish !="NA")%>%
  filter(Elevation>2790)%>%
  ggplot(aes(x = Elevation, y = Body_mass_mg, colour=as.factor(Fish)))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #geom_smooth(method = "lm")+
  stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  scale_color_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  theme_bw()+ylab("CWM")+xlab("Elevation (m)")+
  facet_wrap(~as.factor(Fish),scales = "free")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

datasz%>%
  gather(Head.river.dist,River.dist.lake,key = "var", value = "value") %>% #PC2,PC3,PC4,River.dist.lake,
  ggplot(aes(x = value, y = Body_mass_mg))+ #, colour=as.factor(Fish))) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  #geom_smooth(method = "lm")+
  stat_smooth(method = glm, method.args = list(family = gaussian(link="identity")))+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

datasz%>%
  filter(Fish != "NA")%>%
  ggplot(aes(x = as.factor(Fish), y = Body_mass_mg, fill=as.factor(Fish)))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("no", "yes"))+
  xlab("Fish Presence")+
  labs(fill='Fish Presence') +
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())

#########################################################################
#Visualize influence of fish and elevation


############################################################################################################################################

#Fish  By Elevation
stream.elev.fish<-diversity.env%>%
  ggplot(aes(x = as.factor(Fish), y = Elevation, fill=as.factor(Fish)))+ 
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE,name = "Fish Presence", labels = c("No", "Yes"))+
  xlab("Fish Presence")+
  ggtitle("b)") +
  ylab("Elevation (m)")+
  #scale_fill_discrete(name = "Fish Presence", labels = c("no", "yes"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank())+ theme(legend.position = "none")+
  theme(legend.position = c(0.25, 0.15),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))


############################################################################################################################################

