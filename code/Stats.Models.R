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
r2(mod1)

mod0<-betareg(betas.LCBD~actual_fish_presence, data=env_div)
mod1<-betareg(betas.LCBD~lake_elevation_nbr, data=env_div)
mod2<-betareg(betas.LCBD~lake_elevation_nbr*actual_fish_presence, data=env_div)
mod3<-betareg(betas.LCBD~1, data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

mod0<-glm(Com.Size~actual_fish_presence,family = gaussian(link="identity"), data=env_div)
mod1<-glm(Com.Size~lake_elevation_nbr,family = gaussian(link="identity"), data=env_div)
mod2<-glm(Com.Size~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_div)
mod3<-glm(Com.Size~1,family =gaussian(link="identity"), data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

mod0<-glm(Biomass~actual_fish_presence,family = gaussian(link="identity"), data=env_div)
mod1<-glm(Biomass~lake_elevation_nbr,family = gaussian(link="identity"), data=env_div)
mod2<-glm(Biomass~lake_elevation_nbr*actual_fish_presence,family =gaussian(link="identity"), data=env_div)
mod3<-glm(Biomass~1,family =gaussian(link="identity"), data=env_div)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

###########################################################################
#Streams
mod0<-glm(N1~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(N1~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(N1~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(N1~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

mod0<-betareg(betas.LCBD~Fish, data=diversity.env)
mod1<-betareg(betas.LCBD~Elevation, data=diversity.env)
mod2<-betareg(betas.LCBD~Elevation*Fish, data=diversity.env)
mod3<-betareg(betas.LCBD~1, data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

mod0<-glm(Com.Size~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(Com.Size~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(Com.Size~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(Com.Size~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

mod0<-glm(Biomass~Fish,family = gaussian(link="identity"), data=diversity.env)
mod1<-glm(Biomass~Elevation,family = gaussian(link="identity"), data=diversity.env)
mod2<-glm(Biomass~Elevation*Fish,family =gaussian(link="identity"), data=diversity.env)
mod3<-glm(Biomass~1,family =gaussian(link="identity"), data=diversity.env)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2,mod3,weights = TRUE, sort = FALSE)
reported.table2
r2(mod1)

######################################################################################################################################################
#Permanova/Dispersion: NMDS

######################################################################################################################################################
#CWM


######################################################################################################################################################
#Beta-diversity

######################################################################################################################################################
#Species Responses

######################################################################################################################################################