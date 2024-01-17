#5_plotting: Plotting analysis resutls from bee individual - plant species 
  #networks from FFAR 2019 field data. Loads in a few objects created in 3_objects
  #then generates plots. Requires packages within 1_initiation


#load in objects 
load('data/quant_distSBNHB.RData')
load('data/indivNet_rbcl1.RData')
source('2_functions.R')


#Figure 1: visualizing networks with high (top) to low (bottom) hub scores
par(mfrow=c(3,1),mai=c(.1,.1,.1,.1))
lapply(indivNet_rbcl1[c(2,5,9)],function(x){plotweb(x)})
#lapply(indivNet_rbcl1[c(5)],function(x){plotweb(x, text.rot=90)})
par(mfrow=c(1,1))


#Figure 2 A: rbcL (plant) variation by hub score
ggplot(quant_distSB5NHB[quant_distSB5NHB$label=='rbcl',],
       aes(x=degree,y=avgDist,col=GenusSpecies))+
  scale_color_viridis(discrete=TRUE)+
  theme_ms() +
  geom_point() +
  geom_smooth(method=lm,se=F)+
  #  geom_smooth(method=lm,aes=(x=quant_dist5[quant_dist5$label=='rbcl',]$degree, quant_dist5[quant_dist5$label=='rbcl',]$y=avgDist))+
  ylab("Average distance to rbcL centroid (within species)")+
  xlab("Z-score of maximum plant degree")


#Figure 2 B: microbial variation by hub score
ggplot(quant_distSB5NHB[quant_distSB5NHB$label=='micro',],
       aes(x=degree,y=avgDist,col=GenusSpecies))+
  scale_color_viridis(discrete=TRUE)+
  theme_ms() +
  geom_point() +
  geom_smooth(method=lm,se=F)+
  #  geom_smooth(method=lm,aes=(x=quant_dist5[quant_dist5$label=='rbcl',]$degree, quant_dist5[quant_dist5$label=='rbcl',]$y=avgDist))+
  ylab("Average distance to 16s centroid (within species)")+
  xlab("Z-score of maximum plant degree")
