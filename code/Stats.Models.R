#Title: The influence of non-native fish on stream macroinvertebrate and lake zooplankton communities along elevational gradients
#Authors: Matthew Green, David Herbst, and Kurt Anderson

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
#Body Size species missing

dog<-glm(log(Body_mass_mg+1)~Fish,family=gaussian(link="identity"),all_species_mass_data_env_sum)
summary(dog)
r2(dog)

mod0<-glm(log(Body_mass_mg+1)~Fish,family=gaussian(link="identity"),all_species_mass_data_env_sum)
modnull<-glm(log(Body_mass_mg+1)~1,family=gaussian(link="identity"),all_species_mass_data_env_sum)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2


mod0<-glm(log(Body_mass_ug+1)~Fish,family=gaussian(link="identity"),env_abundz_filter_bm)
modnull<-glm(log(Body_mass_ug+1)~1,family=gaussian(link="identity"),env_abundz_filter_bm)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2

######################################################################################################################################################
#Supplemental

#lakes
mod0<-glm(lake_elevation_nbr~actual_fish_presence,family= gaussian(link = "identity") ,env_div)
summary(mod0)

#streams
mod0<-glm(Elevation~Fish,family= gaussian(link = "identity") ,stream.elev)
summary(mod0)

######################################################################################################################################################


