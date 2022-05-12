#plotting Code

library(cowplot)

###########################################################################
# Figure 1
plot_grid(fig1a,fig1b, nrow=2)


###########################################################################
# Figure 2
plot_grid(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f, nrow=2)

###########################################################################
# Figure 3
plot_grid(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f, nrow=2)

###########################################################################
# Figure 4
par(mfrow=c(1,2))    # set the plotting area into a 1*2 array

#lake
set.seed(99)
species<-sp_abund_env[,3:39]
dune.rel<-decostand(species,"total") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.rel, k=2, try=10) #NMDS code
dune.nmds
#stressplot(dune.nmds)

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
#text(dune.nmds$species[,1], dune.nmds$species[,2], rownames(dune.nmds$species), cex=0.7, col ="black")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) 
ordihull(dune.nmds, sp_abund_env$actual_fish_presence, display="sites", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordisurf(dune.nmds, sp_abund_env$lake_elevation_nbr, prioirty=,labcex=0.9, add = T,col="forestgreen")

line = 1
cex = 1
side = 3
adj=-0.05
mtext("a)", side=side, line=line, cex=cex, adj=adj)

#stream
set.seed(99)
species<-species_2
dune.rel<-decostand(species,"total") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.rel, k=2, try=1000) #NMDS code
dune.nmds
#stressplot(dune.nmds) 

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) 
ordihull(dune.nmds, species_env$Fish, display="sites", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordisurf(dune.nmds, species_env$Elevation, prioirty=,labcex=0.9, add = T,col="forestgreen")
line = 1
cex = 1
side = 3
adj=-0.05
mtext("b)", side=side, line=line, cex=cex, adj=adj)


dev.off()
###########################################################################
# Figure 5:CWM  by Fish
plot_grid(fig5a,fig5b,nrow=1)


###########################################################################
# Figure 6:CWM  by Elevation
plot_grid(fig6a,fig6b,nrow=1)


###########################################################################
# Figure 7


###########################################################################
###########################################################################
#Supplemental Figure S1
plot_grid(lake.elev.fish,stream.elev.fish, nrow=1)
