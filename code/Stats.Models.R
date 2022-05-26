#Stats Code

######################################################################################################################################################
#Diversity

#lakes
mod0<-glm(N1~actual_fish_presence,family = gaussian(link="identity"), data=env_div)
mod1<-glm(N1~lake_elevation_nbr,family = gaussian(link="identity"), data=env_div)
mod2<-glm(N1~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_div)
mod3<-glm(N1~1,family =gaussian(link="identity"), data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
check_collinearity(mod2)
multicollinearity(mod2)
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-betareg(betas.LCBD~actual_fish_presence, data=env_div)
mod1<-betareg(betas.LCBD~lake_elevation_nbr, data=env_div)
mod2<-betareg(betas.LCBD~lake_elevation_nbr*actual_fish_presence, data=env_div)
mod3<-betareg(betas.LCBD~1, data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-glm(Com.Size~actual_fish_presence,family = gaussian(link="identity"), data=env_div)
mod1<-glm(Com.Size~lake_elevation_nbr,family = gaussian(link="identity"), data=env_div)
mod2<-glm(Com.Size~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_div)
mod3<-glm(Com.Size~1,family =gaussian(link="identity"), data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-glm(Biomass~actual_fish_presence,family = gaussian(link="identity"), data=env_div)
mod1<-glm(Biomass~lake_elevation_nbr,family = gaussian(link="identity"), data=env_div)
mod2<-glm(Biomass~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_div)
mod3<-glm(Biomass~1,family =gaussian(link="identity"), data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

###########################################################################
#Streams
mod0<-glm(N1~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(N1~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(N1~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(N1~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-betareg(betas.LCBD~Fish, data=diversity.env)
mod1<-betareg(betas.LCBD~Elevation, data=diversity.env)
mod2<-betareg(betas.LCBD~Elevation*Fish, data=diversity.env)
mod3<-betareg(betas.LCBD~1, data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-glm(Com.Size~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(Com.Size~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(Com.Size~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(Com.Size~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

mod0<-glm(Biomass~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(Biomass~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(Biomass~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(Biomass~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

######################################################################################################################################################
#Permanova/Dispersion: NMDS

#lakes
a<-adonis2(lake.bray ~ sp_abund_env_filter$actual_fish_presence+sp_abund_env_filter$lake_elevation_nbr, permutations = 99, method = "bray")
summary(a)
betad <- betadiver(dune.rel, "z")
adonis(betad ~ sp_abund_env$actual_fish_presence+sp_abund_env$lake_elevation_nbr, data=species, perm=200)
mod <- betadisper(lake.bray, sp_abund_env_filter$actual_fish_presence)
anova(mod)
print(mod)
permutest(mod)
boxplot(mod)

#Stream
adonis2(dune.bray ~ species_env$Fish+species_env$Elevation, permutations = 99, method = "bray")
betad <- betadiver(dune.bray, "z")
adonis(betad ~ species_env$Fish+species_env$Elevation, data=species, perm=200)
mod <- betadisper(dune.bray, species_env$Fish)
anova(mod)
print(mod)
permutest(mod)
boxplot(mod)



######################################################################################################################################################
#CWM

#lakes
mod0<-glm(CWM~actual_fish_presence,family = gaussian(link="identity"), data=env_cwm)
mod1<-glm(CWM~lake_elevation_nbr,family = gaussian(link="identity"), data=env_cwm)
mod2<-glm(CWM~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_cwm)
mod3<-glm(CWM~1,family =gaussian(link="identity"), data=env_cwm)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

#Streams
mod0<-glm(Body_mass_mg~Fish,family = gaussian(link="identity"), data=datasz)
mod1<-glm(Body_mass_mg~Elevation,family = gaussian(link="identity"), data=datasz)
mod2<-glm(Body_mass_mg~Elevation*Fish,family =gaussian(link="identity"), data=datasz)
mod3<-glm(Body_mass_mg~1,family =gaussian(link="identity"), data=datasz)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

######################################################################################################################################################
#Beta-diversity

#Lakes Beta
analysis<-lake.df.filter%>%
  filter(cmtny.dist>0, cmtny.dist<1)
mod0<-betareg(cmtny.dist~Fish.turnover, data=analysis)
mod1<-betareg(cmtny.dist~spatial.dist, data=analysis)
mod2<-betareg(cmtny.dist~spatial.dist*Fish.turnover, data=analysis)
mod3<-betareg(cmtny.dist~1, data=analysis)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)


#Streams Beta
mod0<-betareg(cmtny.dist~Fish.turnover, data=stream.df.filter)
mod1<-betareg(cmtny.dist~spatial.dist, data=stream.df.filter)
mod2<-betareg(cmtny.dist~spatial.dist*Fish.turnover, data=stream.df.filter)
mod3<-betareg(cmtny.dist~1, data=stream.df.filter)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
r2(mod3)

######################################################################################################################################################
#Species Responses

#Relative Change
env_abundz_filter_plot_1<-env_abundz_filter_plot%>%
  filter(relative_change != "NA")%>%
  filter(relative_change > 0.000001)%>%
  filter(relative_change < 5)

mod0<-glm(relative_change~log(Body_mass_ug+1),family=gaussian(link="identity"),env_abundz_filter_plot_1)
modnull<-glm(relative_change~1,family=gaussian(link="identity"),env_abundz_filter_plot_1)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)

species_mass_data_env_sum_1<-species_mass_data_env_sum%>%
  filter(relative_change != "NA")%>%
  filter(relative_change > 0.000001)%>%
  filter(relative_change < 5)

mod0<-glm(relative_change~log(Body_mass_mg+1),family=gaussian(link="identity"),species_mass_data_env_sum_1)
modnull<-glm(relative_change~1,family=gaussian(link="identity"),species_mass_data_env_sum_1)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)


######################################################################################################################################################
#Supplemental

#lakes
mod0<-glm(lake_elevation_nbr~actual_fish_presence,family= gaussian(link = "identity") ,env_div)
summary(mod0)

#streams
mod0<-glm(Elevation~Fish,family= gaussian(link = "identity") ,stream.elev)
summary(mod0)

######################################################################################################################################################
#Ind Species Models

#All sp
mod0<-glm(zoop_density~actual_fish_presence*Taxon,family = gaussian(link="identity"), data=env_abundzzz_new)
mod1<-glm(zoop_density~lake_elevation_nbr*Taxon,family = gaussian(link="identity"), data=env_abundzzz_new)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence*Taxon,family =gaussian(link="identity"), data=env_abundzzz_new)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=env_abundzzz_new)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

#"Leptodiaptomus_signicauda" "nauplii"                   "Daphnia_melanica"          "Daphnia_dentifera"        
# "Conochilus_unicornis"      "Keratella_spp."            "Keratella_quadrata"        "Hesperodiaptomus_shoshone"
# "Polyarthra_dolichoptera"   "Chydorus_sphaericus"       "Holopedium_gibberum"       "Lepadella_spp."           
# "Bosmina_longirostris"      "Cyclopoida"                "Trichocerca_capucina"      "Harpacticoida"            
# "Kellicottia_spp."          "Polyarthra_vulgaris"       "Asplanchna_spp."           "Trichotria_spp."          
# "Hesperodiaptomus_eiseni"   "Alona_spp."                "Ceriodaphnia_laticaudata"  "Alonella_excisa"          
# "Ascomorpha_spp."           "Synchaeta_spp."            "Chaoborus_trivitattus"     "Daphnia_pulex_pulicaria"  
# "Diaphanosoma_brachyurum"   "Polyphemus_pediculus"  
#A) Lakes
Bosmina<-env_abundzzz_new%>%filter(Taxon=="Bosmina_longirostris")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Bosmina)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Bosmina)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Bosmina)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Bosmina)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Ceriodaphnia_laticaudata<-env_abundzzz_new%>%filter(Taxon=="Ceriodaphnia_laticaudata")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Ceriodaphnia_laticaudata)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Ceriodaphnia_laticaudata)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Ceriodaphnia_laticaudata)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Ceriodaphnia_laticaudata)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Chaoborus_trivitattus<-env_abundzzz_new%>%filter(Taxon=="Chaoborus_trivitattus")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Chaoborus_trivitattus)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Chaoborus_trivitattus)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Chaoborus_trivitattus)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Chaoborus_trivitattus)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Chydorus_sphaericus<-env_abundzzz_new%>%filter(Taxon=="Chydorus_sphaericus")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Chydorus_sphaericus)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Chydorus_sphaericus)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Chydorus_sphaericus)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Chydorus_sphaericus)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Conochilus_unicornis<-env_abundzzz_new%>%filter(Taxon=="Conochilus_unicornis")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Conochilus_unicornis)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Conochilus_unicornis)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Conochilus_unicornis)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Conochilus_unicornis)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Cyclopoida<-env_abundzzz_new%>%filter(Taxon=="Cyclopoida")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Cyclopoida)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Cyclopoida)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Cyclopoida)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Cyclopoida)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Daphnia_dentifera<-env_abundzzz_new%>%filter(Taxon=="Daphnia_dentifera")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Daphnia_dentifera)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Daphnia_dentifera)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Daphnia_dentifera)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Daphnia_dentifera)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Daphnia_melanica<-env_abundzzz_new%>%filter(Taxon=="Daphnia_melanica")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Daphnia_melanica)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Daphnia_melanica)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Daphnia_melanica)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Daphnia_melanica)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Daphnia_pulex_pulicaria<-env_abundzzz_new%>%filter(Taxon=="Daphnia_pulex_pulicaria")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Daphnia_pulex_pulicaria)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Daphnia_pulex_pulicaria)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Daphnia_pulex_pulicaria)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Daphnia_pulex_pulicaria)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Diaphanosoma_brachyurum<-env_abundzzz_new%>%filter(Taxon=="Diaphanosoma_brachyurum")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Diaphanosoma_brachyurum)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Diaphanosoma_brachyurum)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Diaphanosoma_brachyurum)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Diaphanosoma_brachyurum)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Harpacticoida<-env_abundzzz_new%>%filter(Taxon=="Harpacticoida")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Harpacticoida)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Harpacticoida)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Harpacticoida)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Harpacticoida)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Hesperodiaptomus_eiseni<-env_abundzzz_new%>%filter(Taxon=="Hesperodiaptomus_eiseni")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Hesperodiaptomus_eiseni)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Hesperodiaptomus_eiseni)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Hesperodiaptomus_eiseni)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Hesperodiaptomus_eiseni)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Hesperodiaptomus_shoshone<-env_abundzzz_new%>%filter(Taxon=="Hesperodiaptomus_shoshone")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Hesperodiaptomus_shoshone)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Hesperodiaptomus_shoshone)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Hesperodiaptomus_shoshone)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Hesperodiaptomus_shoshone)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Hesperodiaptomus_shoshone_1<-Hesperodiaptomus_shoshone%>%
  pivot_wider(names_from = actual_fish_presence,values_from = zoop_density)%>%
  dplyr::mutate(Fishless.Occupancy=if_else(No>0, 1,0),Fish.Occupancy=if_else(Yes>0, 1,0))%>%
  replace(is.na(.), 0)%>% 
  group_by(Taxon,)%>%
  dplyr::summarise(n=n(),Fish.total.occupancy=sum(Fish.Occupancy),Fishless.total.occupancy=sum(Fishless.Occupancy),
                   fish.prop.occupancy=Fish.total.occupancy/n, fishless.prop.occupancy=Fishless.total.occupancy/n)

Holopedium_gibberum<-env_abundzzz_new%>%filter(Taxon=="Holopedium_gibberum")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Holopedium_gibberum)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Holopedium_gibberum)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Holopedium_gibberum)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Holopedium_gibberum)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Kellicottia_spp.<-env_abundzzz_new%>%filter(Taxon=="Kellicottia_spp.")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Kellicottia_spp.)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Kellicottia_spp.)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Kellicottia_spp.)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Kellicottia_spp.)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Keratella_quadrata<-env_abundzzz_new%>%filter(Taxon=="Keratella_quadrata")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Keratella_quadrata)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Keratella_quadrata)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Keratella_quadrata)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Keratella_quadrata)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Keratella_spp.<-env_abundzzz_new%>%filter(Taxon=="Keratella_spp.")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Keratella_spp.)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Keratella_spp.)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Keratella_spp.)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Keratella_spp.)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Lepadella_spp.<-env_abundzzz_new%>%filter(Taxon=="Lepadella_spp.")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Lepadella_spp.)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Lepadella_spp.)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Lepadella_spp.)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Lepadella_spp.)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Leptodiaptomus_signicauda<-env_abundzzz_new%>%filter(Taxon=="Leptodiaptomus_signicauda")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Leptodiaptomus_signicauda)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Leptodiaptomus_signicauda)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Leptodiaptomus_signicauda)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Leptodiaptomus_signicauda)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Polyarthra_dolichoptera<-env_abundzzz_new%>%filter(Taxon=="Polyarthra_dolichoptera")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Polyarthra_dolichoptera)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Polyarthra_dolichoptera)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Polyarthra_dolichoptera)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Polyarthra_dolichoptera)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Trichocerca_capucina<-env_abundzzz_new%>%filter(Taxon=="Trichocerca_capucina")
mod0<-glm(zoop_density~actual_fish_presence,family = gaussian(link="identity"), data=Trichocerca_capucina)
mod1<-glm(zoop_density~lake_elevation_nbr,family = gaussian(link="identity"), data=Trichocerca_capucina)
mod2<-glm(zoop_density~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=Trichocerca_capucina)
mod3<-glm(zoop_density~1,family =gaussian(link="identity"), data=Trichocerca_capucina)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)
######################################################################################################################################################
#Streams
species_mass_data_env_filter

Ostracoda<-species_mass_data_env_filter%>%filter(Taxon=="Ostracoda")
mod0<-glm(abundance~Fish,family = gaussian(link="identity"), data=Ostracoda)
mod1<-glm(abundance~Elevation,family = gaussian(link="identity"), data=Ostracoda)
mod2<-glm(abundance~Elevation*Fish,family =gaussian(link="identity"), data=Ostracoda)
mod3<-glm(abundance~1,family =gaussian(link="identity"), data=Ostracoda)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

Acari<-species_mass_data_env_filter%>%filter(Taxon=="Acari")
mod0<-glm(abundance~Fish,family = gaussian(link="identity"), data=Acari)
mod1<-glm(abundance~Elevation,family = gaussian(link="identity"), data=Acari)
mod2<-glm(abundance~Elevation*Fish,family =gaussian(link="identity"), data=Acari)
mod3<-glm(abundance~1,family =gaussian(link="identity"), data=Acari)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)
r2(mod1)
r2(mod2)

######################################################################################################################################################