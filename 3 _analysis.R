#Analysis of bee individual - plant species networks from FFAR 2019
  #field data

  #Requires packages found in 1_packages, and loads functions from 2_functions


#load functions
source('2_functions.R')

#load network data
load('data/rbclNets.RData')
load('data/paraNets.RData')
load('data/microNets.RData')
load('data/spec_RBCL_16s.RData')
load('data/trees.RData')
load('data/covarmatrix_community.RData')


####*****************************************

#Build objects, explore data 

####*****************************************

##Plants - rbcl

#remove NAs
indivNet_rbcl1 <- lapply(indivNet_rbcl,function(x){
  keep.ls <- unlist(lapply(c(1:length(rownames(x))),function(y){
    all(!is.na(x[y,]))
  }))
  return(x[keep.ls,])
})

#make a binary version of the matrices
bin_rbcl <- lapply(indivNet_rbcl, function(x){
  x[x>0] <- 1
  return(x)
})

#estimate centroids and calculate distances via vegan::betadisper
rbcl_dis <- netToDisper(indivNet_rbcl1,spec,'rbcl')

bin_rbcl_dis <- netToDisper(bin_rbcl,spec,'rbcl')


##microbes
bin_micro <- lapply(indivNet_micro, function(x){
  x[x>0] <- 1
  return(x)
})

#estimate centroids and calculate distances via vegan::betadisper
micro_dis <- netToDisper(indivNet_micro,spec,'micro')
bin_micro_dis <- netToDisper(bin_micro,spec,'micro')

##visualize relative abun of bacterial sp
#plot(rowSums(indivNet_micro[[1]],na.rm=T))


## --- number of microbes/plants present in diff percentages of the sample (Table3)

#put all sites together into 1 matrix
bin_micro_all <- do.call(cbind,bin_micro) 

#number of total bee individuals: 1686
length(colnames(bin_micro_all))

# number of microbes found in at least 5% (84/1686) of samples
length(rownames(bin_micro_all[rowSums(bin_micro_all,na.rm=T)>84,]))
#38/758, or 5.01% of ASVs

# number of microbes found in at least 10% of samples 168
length(rownames(bin_micro_all[rowSums(bin_micro_all,na.rm=T)>168,]))
#23, or 3.03% of microbial ASVs

# number of microbes found in at least 20% of samples 337
length(rownames(bin_micro_all[rowSums(bin_micro_all,na.rm=T)>337,]))
#11, or 1.45% of microbial ASVs

# number of microbes found in at least 50% of samples 843
length(rownames(bin_micro_all[rowSums(bin_micro_all,na.rm=T)>843,]))
#3, or 0.395% of microbial ASVs

#plants
bin_rbcl_all <- do.call(cbind,bin_rbcl) #put all sites together into 1 matrix

length(colnames(bin_rbcl_all))
#1178

# number of plants found in at least 5% (59/1178) of samples
length(rownames(bin_rbcl_all[rowSums(bin_rbcl_all,na.rm=T)>59,]))
#21/123, or 17.1% of plant ASVs

# number of microbes found in at least 10% of samples 
length(rownames(bin_rbcl_all[rowSums(bin_rbcl_all,na.rm=T)>118,]))
#13, or 10.5% of plant ASVs

# number of microbes found in at least 20% of samples 
length(rownames(bin_rbcl_all[rowSums(bin_rbcl_all,na.rm=T)>236,]))
#8, or 6.50% of plant ASVs

# number of microbes found in at least 50% of samples 
length(rownames(bin_rbcl_all[rowSums(bin_rbcl_all,na.rm=T)>589,]))
#0, or 0% of plant ASVs


## -- presence of floral specialist microbes in the bee microbiome samples

#floral associated microbes (from Vannette annual review 2020 and Alvarez 2012)

flr_bact_list <- c("Acinetobacter", "Rosenbergiella", "Metschnikowia", "Pseudomonas", "Pantoea")

flr_bacts <- grep(paste(flr_bact_list, collapse="|"),rownames(bin_micro_all))
  #just three rows: 146, 148, 160



####*****************************************

#Is variability correlated between interaction types?

####*****************************************

#reorganize the micro and rbcl dataframes
micro.df <- data.frame('obs' = paste(micro_dis$GenusSpecies, micro_dis$site, sep="_"),
                       'microDist' = micro_dis$avgDist,
                       'microN' = micro_dis$n)
rbcl.df <- data.frame('obs' = paste(rbcl_dis$GenusSpecies, rbcl_dis$site, sep="_"),
                      'rbclDist' = rbcl_dis$avgDist,
                      'rbclN' = rbcl_dis$n)

#combine them into a single df, add a couple organizer columns
comp.df <- full_join(micro.df,rbcl.df)
comp.df$site <- word(comp.df$obs,2,sep="_")
comp.df$GenusSpecies <- word(comp.df$obs,1,sep="_")

#cor.lme <- lme(microDist~rbclDist,data=comp.df, random=~1|GenusSpecies,na.action=na.omit)
#summary(cor.lme)
##yes, variability in rbcl predicts variability in microbes

#plot(comp.df$microDist~comp.df$rbclDist)
##ah, maybe b/c singletons: several data points w/ 0s for both

#Filter out species with >4 observations, honeybees, and line 67 (can't remember why???)
comp.filt <- droplevels(comp.df[comp.df$microN>4 & comp.df$rbclN>4,])
comp.filt <- comp.filt[-67,]
comp.filtNHB <- subset(comp.filt,comp.filt$GenusSpecies!="Apis mellifera")


#RESULT 1: rbcl variation corr w/ microbe var, quant nets no zeros 
cor.lme.filt <- lme(microDist~rbclDist,data=comp.filtNHB, random=~1|GenusSpecies,na.action=na.omit)
summary(cor.lme.filt)
r.squaredGLMM(cor.lme.filt)
  #RESULT PRESENTED

#repeating with honeybees included: no qualitative change
cor.lme.filtHB <- lme(microDist~rbclDist,data=comp.filt, random=~1|GenusSpecies,na.action=na.omit)
summary(cor.lme.filtHB)


####--------
#Repeating the above with Core baceria only 
####--------

#currently based on Graystock, Rehan, and McFrederick 2017, conservation genetics
core.sp <- c("Rosenbergiella", "Pseudomonas", "Gilliamella", "Lactobacillus", "Caulobacter", "Snodgrassella", "Acinetobacter", "Corynebacterium", "Sphingomonas", "Commensalibacter", "Methylobacterium", "Massilia","Stenotrophomonas")

core.ls <- do.call(c,lapply(core.sp,
                            function(x){
                              grep(x,names(spec))
                            }))

core.micro <- illumSplit(spec,"Site",core.ls)

#had to get rid of 0s here, since vegdist got angry if i didn't
core.micro <- lapply(core.micro,function(x){
  filt <- x[,colSums(x,na.rm=T)>0]
  return(filt)
})

micro_core <- netToDisper(core.micro,spec,'micro')

core.df <- data.frame('obs' = paste(micro_core$GenusSpecies, micro_core$site, sep="_"),
                      'microDist' = micro_core$avgDist,
                      'microN' = micro_core$n)

comp.core <- full_join(core.df,rbcl.df)
comp.core$site <- word(comp.core$obs,2,sep="_")
comp.core$GenusSpecies <- word(comp.core$obs,1,sep="_")

comp.core.filt <- droplevels(comp.core[comp.core$microN>4 & comp.core$rbclN>4,])


#Result: no qual change
core.lme.filt <- lme(microDist~rbclDist,data=comp.core.filt, random=~1|GenusSpecies,na.action=na.omit,subset=comp.core.filt$rbclDist>0)
summary(core.lme.filt)
#rbcl dist correlated with core micro dist


####*****************************************

#Hubbyness metrics

####*****************************************

#takes list of networks, calculates the different potential hubbyness metrics for each
huby <- netTohub(indivNet_rbcl1,names(indivNet_rbcl1))
  #note: bmotif::mcount seems to be deprecated and are no longer used below 
      #(but are still in the functions). Degree provided a better hubbyness metric
      #than motif composition.
save(huby,file='data/huby.RData')
load('data/huby.RData')

#quant_dist$M44 <- huby$M44[match(quant_dist$site,huby$site)]
#quant_dist$M38 <- huby$M38[match(quant_dist$site,huby$site)]
#quant_dist$M44f <- huby$M44f[match(quant_dist$site,huby$site)]
#quant_dist$clustering <- huby$clustering[match(quant_dist$site,huby$site)]
quant_dist$degree <- huby$degree[match(quant_dist$site,huby$site)]

quant_dist5$degree <- huby$degree[match(quant_dist5$site,huby$site)]
quant_dist5NHB$degree <- huby$degree[match(quant_dist5NHB$site,huby$site)]


#Figure 1: visualizing networks with high (top) to low (bottom) hub scores
par(mfrow=c(3,1),mai=c(.1,.1,.1,.1))
lapply(indivNet_rbcl1[c(2,5,9)],function(x){plotweb(x)})
#lapply(indivNet_rbcl1[c(5)],function(x){plotweb(x, text.rot=90)})
par(mfrow=c(1,1))

##to plot them individually (to better see species labels):
#plotweb(indivNet_rbcl1[[2]],text.rot=90)

##to visualize all networks:
#par(mfrow=c(4,4))
#lapply(indivNet_rbcl1,function(x){plotweb(x)})
#par(mfrow=c(1,1))

#if par is funny and doesnt have axis labels after above:
dev.off()

##To visualize degree hub score
#barplot(huby$degree,ylab='Max Degree z-scrore'


####
#Effect of Hubbyness on variability:
####

#Model without phylogenetic correction (not presented)
huby5.lme4 <- lmer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|site),
                   data=quant_dist5NHB,
                   control = lmerControl(optimizer = 'bobyqa')
)
#no main effect of netDeg BUT
#sig effect of rbcl, such that its higher than micro (rbcl more variable)
#sig interaction of rbcl with netDeg, such that low netDeg <-> high avgDist 
#getting rid of HB (by removing "NHB" from the df) also doesn't change it qual

summary(huby5.lme4)
#lmerTest does bootstrapping for lme4

#Adding in bee phylogeny to the model
#Almer one needs the bee tree rather than covar matrix
A <- Matrix::Matrix(ape::vcv(phylo),sparse=T)
colnames(A) <- rownames (A) <- colnames(co.var.mat)
huby5.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|site),
                     data=quant_dist5NHB,
                     control = lmerControl(optimizer = 'bobyqa'),
                     A=list(GenusSpecies=A)
)
summary(huby5.almer)
  #Qualitatively the same as the models w/o phylogeny. 
  #RESULT PRESENTED
  

##another representation of the pattern: flat for micro, downsloping for rbcl
#bargraph.CI(response=avgDist,x.factor=label,group=degree,data=quant_dist5)


#Figure 2 A: rbcL (plant) variation by hub score
ggplot(quant_dist5NHB[quant_dist5NHB$label=='rbcl',],aes(x=degree,y=avgDist,col=GenusSpecies))+
  scale_color_viridis(discrete=TRUE)+
  theme_ms() +
  geom_point() +
  geom_smooth(method=lm,se=F)+
  #  geom_smooth(method=lm,aes=(x=quant_dist5[quant_dist5$label=='rbcl',]$degree, quant_dist5[quant_dist5$label=='rbcl',]$y=avgDist))+
  ylab("Average distance to rbcL centroid (within species)")+
  xlab("Z-score of maximum plant degree")


#Figure 2 B: microbial variation by hub score
ggplot(quant_dist5NHB[quant_dist5NHB$label=='micro',],aes(x=degree,y=avgDist,col=GenusSpecies))+
  scale_color_viridis(discrete=TRUE)+
  theme_ms() +
  geom_point() +
  geom_smooth(method=lm,se=F)+
  #  geom_smooth(method=lm,aes=(x=quant_dist5[quant_dist5$label=='rbcl',]$degree, quant_dist5[quant_dist5$label=='rbcl',]$y=avgDist))+
  ylab("Average distance to 16s centroid (within species)")+
  xlab("Z-score of maximum plant degree")

#  xlab("Z-score of maximum plant degree")







