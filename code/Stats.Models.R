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

######################################################################################################################################################
#Species Responses

#Relative Change

mod0<-glm(relative_change~log(Body_mass_ug+1),family=gaussian(link="identity"),env_abundz_filter_plot)
modnull<-glm(relative_change~1,family=gaussian(link="identity"),env_abundz_filter_plot)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)


mod0<-glm(relative_change~log(Body_mass_mg+1),family=gaussian(link="identity"),species_mass_data_env_sum)
modnull<-glm(relative_change~1,family=gaussian(link="identity"),species_mass_data_env_sum)
reported.table2 <- bbmle::AICtab(mod0,modnull,weights = TRUE, sort = FALSE)
reported.table2
r2(mod0)

#Ind Species Models

######################################################################################################################################################
#Supplemental

#lakes
mod0<-glm(lake_elevation_nbr~actual_fish_presence,family= gaussian(link = "identity") ,env_div)
summary(mod0)

#streams
mod0<-glm(Elevation~Fish,family= gaussian(link = "identity") ,stream.elev)
summary(mod0)

######################################################################################################################################################