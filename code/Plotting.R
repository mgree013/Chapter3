#plotting Code

library(cowplot)

###########################################################################
# Figure 1
plot_grid(fig1a,fig1b, nrow=2)

plot_grid(new.fig1a,new.fig1b, nrow=1)

###########################################################################
# Figure 2
#plot_grid(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f, nrow=2)
plot_grid(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f,fig2g,fig2h, nrow=2)

###########################################################################
# Figure 3
#plot_grid(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f, nrow=2)
plot_grid(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f,fig3g,fig3h, nrow=2)

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
points(lake.nmds$points[,1], lake.nmds$points[,2],  pch = 1) 
ordihull(lake.nmds, sp_abund_env_filter$actual_fish_presence, display="sites", label=F,lwd=2, col=c("#440154FF","#FDE725FF"))
ordisurf(lake.nmds, sp_abund_env_filter$lake_elevation_nbr, prioirty=,labcex=0.9, add = T,col="forestgreen")

line = 1
cex = 1
side = 3
adj=-0.05
mtext("a)", side=side, line=line, cex=cex, adj=adj)

#stream
env_elevation<-envs%>%filter(Elevation >3200)%>%dplyr::select(c(Elevation,Site))
species_env<-species_env%>%filter(Elevation >3200)%>%dplyr::select(c(Elevation,Site,Fish))

set.seed(99)
species<-species_2%>%rownames_to_column("Site")%>%left_join(env_elevation, by="Site")%>%
  filter(Elevation >3200)%>%dplyr::select(-c(Elevation))%>%column_to_rownames("Site")#%>%dplyr::select(-c(Chironomidae,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))
dune.rel<-decostand(species,"hellinger") #standardize community data
dune.bray<-vegdist(dune.rel) #calculate dissimilarity among sites (i.e. dissimilarity matrix)
dune.nmds=metaMDS(dune.bray, k=2, try=1000) #NMDS code
dune.nmds
#stressplot(dune.nmds) 

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) 
#ordilabel(dune.nmds, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
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
plot_grid(beta.stream.a,beta.stream.b,beta.stream.c,beta.stream.d,nrow=2)


###########################################################################
###########################################################################
#Supplemental Figure S1
#plot_grid(lake.elev.fish,stream.elev.fish, nrow=1)

plot_grid(supp.a,supp.b,supp.c,supp.d, nrow=2)

