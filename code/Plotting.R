#Title: The influence of non-native fish on stream macroinvertebrate and lake zooplankton communities along elevational gradients
#Authors: Matthew Green, David Herbst, and Kurt Anderson

#plotting Code

library(cowplot)
library(tidyverse)

###########################################################################
# Figure 1

#plot_grid(new.fig1a,new.fig1b,fig1c,fig1d, nrow=2)
plot_grid(new.fig1a,new.fig1b, nrow=1)
plot_grid(fig1c,fig1d, nrow=1)


plot_grid(fig1a,fig1b, nrow=2)
plot_grid(new.prop.a,new.prop.b,nrow=2)
plot_grid(fig1a,fig1b,new.prop.a,new.prop.b,nrow=4)

#plot_grid(new.prop.a2,new.prop.b2,nrow=2)

###########################################################################
# Figure 2
#plot_grid(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f, nrow=2)
plot_grid(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f,fig2g,fig2h, nrow=2)

###########################################################################
# Figure 3
#plot_grid(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f, nrow=2)
plot_grid(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f,fig3g,fig3h, nrow=2)

plot_grid(fig2a,fig2b,fig3a,fig3b,fig2e,fig2f,fig3e,fig3f,nrow=2)
plot_grid(fig2a,fig3a,fig2b,fig3b,fig2e,fig3e,fig2f,fig3f,nrow=2)

###########################################################################
# Figure 4
par(mfrow=c(1,2))    # set the plotting area into a 1*2 array


sp_abund_env_filter<-sp_abund_env%>%filter(lake_elevation_nbr >1800, lake_elevation_nbr <3500)%>%
  filter(HA>=0.5)%>%filter(lake_max_depth>3)
#lake
set.seed(99)
species<-sp_abund_env_filter[,3:39]
lake.rel<-decostand(species,"total") #standardize community data
lake.bray<-vegdist(lake.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
lake.nmds=metaMDS(lake.bray, k=2, try=1000) #NMDS code
lake.nmds
#stressplot(dune.nmds)

plot(lake.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
#text(dune.nmds$species[,1], dune.nmds$species[,2], rownames(dune.nmds$species), cex=0.7, col ="black")
points(lake.nmds$points[,1], lake.nmds$points[,2],  pch = 16, col=c("#440154FF","#FDE725FF")[sp_abund_env_filter$actual_fish_presence]) 
#ordihull(lake.nmds, sp_abund_env_filter$actual_fish_presence, display="sites", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordiellipse(lake.nmds, groups=sp_abund_env_filter$actual_fish_presence, display="sites",kind="sd",conf=.9,draw="lines", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordisurf(lake.nmds, sp_abund_env_filter$lake_elevation_nbr, prioirty=,labcex=0.9, add = T,col="forestgreen")

line = 1
cex = 1
side = 3
adj=-0.05
mtext("a)", side=side, line=line, cex=cex, adj=adj)

#stream
env_elevation<-envs%>%filter(Elevation >3200)%>%dplyr::select(c(Elevation,Site))
species_env<-species_env%>%filter(Elevation >3200)%>%dplyr::select(c(Elevation,Site,Fish))
species_env$Fish<-as.factor(species_env$Fish)

set.seed(99)
species<-species_2%>%rownames_to_column("Site")%>%left_join(env_elevation, by="Site")%>%
  filter(Elevation >3200)%>%dplyr::select(-c(Elevation))%>%column_to_rownames("Site")#%>%dplyr::select(-c(Chironomidae,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))
dune.rel<-decostand(species,"hellinger") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.bray, k=2, try=1000) #NMDS code
dune.nmds
#stressplot(dune.nmds) 

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 16, col=c("#440154FF","#FDE725FF")[species_env$Fish]) 
#ordilabel(dune.nmds, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
#ordihull(dune.nmds, species_env$Fish, display="sites", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordiellipse(dune.nmds, groups=species_env$Fish, display="sites",kind="sd",conf=.9,draw="lines", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordisurf(dune.nmds, species_env$Elevation, prioirty=,labcex=0.9, add = T,col="forestgreen")
line = 1
cex = 1
side = 3
adj=-0.05
mtext("b)", side=side, line=line, cex=cex, adj=adj)

# calculating means and SE for each province
dune.nmds.df<-as.data.frame(dune.nmds$points)
dune.nmds.df.full<-dune.nmds.df%>%rownames_to_column("Site")%>%
  left_join(species_env, by="Site")

`calcSE` <-
  function(data, index) sd(data[index])

meanNMDS1 = as.data.frame(aggregate(dune.nmds.df.full$MDS1~dune.nmds.df.full$Fish, FUN=mean))
names(meanNMDS1) = c("fish", "NMDS1")
meanNMDS2 = as.data.frame(aggregate(dune.nmds.df.full$MDS2~dune.nmds.df.full$Fish, FUN=mean))
names(meanNMDS2) = c("fish", "NMDS2")
seNMDS1 = as.data.frame(aggregate(dune.nmds.df.full$MDS1~dune.nmds.df.full$Fish, FUN=calcSE))
names(seNMDS1) = c("fish", "NMDS1se")
seNMDS2 = as.data.frame(aggregate(dune.nmds.df.full$MDS2~dune.nmds.df.full$Fish, FUN=calcSE))
names(seNMDS2) = c("fish", "NMDS2se")


##merging all the pieces
df = merge(meanNMDS1, meanNMDS2 )
df = merge(df, seNMDS1)
df = merge(df, seNMDS2)

## plot of mean and SE of scores  ##
ggplot(data = df,aes(x = NMDS1,y = NMDS2,shape = fish)) + 
  geom_point(aes(size=5)) + 
  ylim(-0.3, 0.3)+
  theme_bw() +
  theme(axis.line.x.top = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.x.bottom = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y.left = element_line(colour = 'black', size=0.5, linetype='solid'))
        


dev.off()
###########################################################################
# Figure 5:CWM  by Fish
plot_grid(fig5a,fig5b,nrow=1)


###########################################################################
# Figure 6:CWM  by Elevation
plot_grid(fig6a,fig6b,nrow=1)


###########################################################################
# Figure 7
plot_grid(beta.stream.a,beta.stream.b,beta.stream.c,beta.stream.d,nrow=2)

a<-plot_grid(beta.a,beta.b,nrow=1)
file4 <- tempfile("a", fileext = ".pdf")
save_plot(file4,base_width = 12, base_height = 9,
          dpi = 300)
ggsave("plotname.png", plot = a, width = 25, height = 10, dpi=600, units = "cm")


###########################################################################
###########################################################################
#Supplemental Figure S1
#plot_grid(lake.elev.fish,stream.elev.fish, nrow=1)

plot_grid(supp.a,supp.b,supp.c,supp.d, nrow=2)

