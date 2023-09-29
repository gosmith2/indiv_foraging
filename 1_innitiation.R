## 1. Innitiation
#First script of analysis for the individual foraging project


rm(list=ls())
library(bipartite)
library(tidyverse) 
library(parallel) 
library(vegan)
library(stringr)
library(nlme)
library(fossil)
library(reshape2)

source('functions.R')


#load the rbcl reads
load('c:/Users/14132/Dropbox/sunflower_saved/data/RBCL2018.RData')

#load in the compiled data with info on each sample
complete <- read.csv('c:/Users/14132/Dropbox/sunflower_saved/data/relational/relational/traditional/specimens-complete.csv')
complete <- dat.clean(complete)
complete$SiteYr <- paste(complete$Site, substring(complete$Date,7,10), sep = "_")
complete$SpSiteYr <- paste(complete$GenusSpecies, complete$SiteYr, sep = "_")

#subset the larger dataset to just samples that were processed in rbcl
sampled <- complete[complete$UniqueID %in% RBCL2018$UniqueID,]


#get species visitation abundances from complete to weight species within networks
sites <- complete[complete$SiteYr %in% sampled$SiteYr,]

abun.ls <- do.call(rbind, lapply(unique(sites$SpSiteYr), function(x){
  
  sp <- data.frame('SpSiteYr' = x)
  sp$abun <- length(sites$SpSiteYr[sites$SpSiteYr==x])
  return(sp)
}))



#quick data cleaning, adding some organizing columns, then merging the dataframes
samp.df <- dat.clean(sampled)
samp.df$SiteYr <- paste(samp.df$Site, substring(samp.df$Date,7,10), sep = "_")
samp.df$SpSiteYr <- paste(samp.df$GenusSpecies, samp.df$SiteYr, sep = "_")

sampAll.df <- left_join(samp.df,RBCL2018, by = 'UniqueID')

save(sampAll.df,file='data/sampAll.RData')

load('data/sampAll.RData')

#length object to isolate just the reads
len <- c((length(samp.df)+1):length(sampAll.df))

siteSplit.ls <- lapply(unique(sampAll.df$SiteYr),function(x){
  site <- sampAll.df[sampAll.df$SiteYr==x,]
})
names(siteSplit.ls) <- unique(sampAll.df$SiteYr)


#make an object of just the read numbers for community stuff

nets.ls <- lapply(siteSplit.ls,function(x){
  nets <- illumSplit(x,"SpSiteYr",len)
})

weightNets.ls <- lapply(nets.ls,function(x){
  sps <- lapply(names(x),function(y){
    sp <- x[[y]]
    sp1 <- sp * (abun.ls$abun[(abun.ls$SpSiteYr %in% y)]/length(colnames(sp)))
    colnames(sp1) <- paste(colnames(sp1),sub("\\_.*", "", y),sep="_")
    return(round(sp1))
  })
  names(sps) <- names(x)
  return(sps)
})

save(nets.ls,file='data/nets.RData')

save(weightNets.ls,file='weightNets.RData')
#nets <- illumSplit(siteSplit.ls,"SpSiteYr",len)


#to randomize indivs w/in species (submatrices), lapply vaznull on nets, then cbind sp together

cores <- 1

vazNets100 <- vaznullFromSp(nets.ls,100,5, FALSE, cores)
save(vazNets100,file='data/vazNets100.RData')

vazNets100abun <- vaznullFromSp(weightNets.ls,100,5, FALSE, cores)

vazNetsBin100 <- vaznullFromSp(nets.ls,100,5, cores, binary = TRUE)


splvlLS100<- lapply(vazNets100, function(x){
  mclapply(x,function(y){
    mets <- specieslevel(y)
  },mc.cores=cores)
})
  #takes whole network (siteYr), with colunns as indivs. calculates node level stats

splvlLS100abun <- lapply(vazNets100abun,function(x){
  mclapply(x,function(y){
    mets <- specieslevel(y)
  },mc.cores=cores)
})


splvlLSBin<- lapply(vazNetsBin100, function(x){
  mclapply(x,function(y){
    mets <- specieslevel(y)
  },mc.cores=cores)
})


save(splvlLS100, file='data/splvlLS100.RData')
save(splvlLSBin100, file='data/splvlLSBin100.RData')
save(splvlLS100abun,file='data/splvlLS100abun.RData')

metrics <- names(splvlLS100abun[[1]]$BM_2018$`higher level`)

#netlvl100 <- lapply(vazNets100, function(x){
#  mclapply(x,function(y){
#    mets <- networklevel(y, metrics)
#  },mc.cores=cores)
#})


splvlComp100 <- splvl.dfer(splvlLS100)

splvlComp100abun <- splvl.dfer(splvlLS100abun)

splvlCompBin100 <- splvl.dfer(splvlLSBin100)


splvlZ100 <- splvl.Zer(splvlComp100,"SpSiteYr")

splvlZ100abun <- splvl.Zer(splvlComp100abun,'SpSiteYr')


splvlZBin100 <- splvl.Zer(splvlCompBin100,"SpSiteYr")


spNetsabun <- build.spNet(weightNets.ls)
  
spNetsSum <- sum.spNet(spNetsabun)

splvlZ100abun$spdeg <- spNetsSum$degree[match(splvlZ100abun$interaction,spNetsSum$SpSiteYr)]
splvlZ100abun$meanIndivDeg <- lapply(splvlZ100abun$interaction,function(x){
  sp <- subset(splvlComp100abun,splvlComp100abun$SpSiteYr==x)
  meandeg <- mean(sp$degree)
})
  #yea, species have much higher degree than indivs (duh)

splvlZBin100$spdeg <- spNetsSum$degree[match(splvlZBin100$interaction,spNetsSum$SpSiteYr)]
splvlZBin100$meanIndivDeg <- lapply(splvlZBin100$interaction,function(x){
  sp <- subset(splvlCompBin100,splvlCompBin100$SpSiteYr==x)
  meandeg <- mean(sp$degree)
})



save(splvlZ100abun,file='data/splvlZ100abun.RData')


#Result 1: 
  #individuals have much lower degree than species
  #AGREE for both quant and binomial


  #extras: 
  #weighted betweenness lower than expected in many (only a few sig)
  #proportional gen and similarity much higher than expected
    #all sig, but low variability in Z scores makes me think its an artifact
  #d' lower than expected, most sig

#######################

## netlvl

#######################


load('data/vazNets.RData')

metrics.net <- c('NODF',
                 'weighted NODF',
                 'H2',
                 'robustness',
                 'vulnerability',
                 'niche overlap')

netlvl100.df <- build.netlvl(vazNets100,metrics.net) 
save(netlvl100.df,file='data/netlvl100.RData')

netlvlBin100.df <- build.netlvl(vazNetsBin100,metrics.net, weighted=FALSE)
save(netlvlBin100.df,file='data/netlvlBin100.RData')

netlvl100abun.df <- build.netlvl(vazNets100abun, metrics.net)
save(netlvl100abun.df,file='data/netlvl100abun.RData')


  #takes whole siteYr, but columns are indivs. calculates network-level metrics

netlvlZ100 <- netlvl.Zer(netlvl100.df)
  #so NODF is higher at 3 sites, h2 MUCH higher at all sites, overlapHL lower at all sites, vuln much lower at all sites 
  #tur found NODF lower, h2 lower
  #so... even for specialists in presence of burst of their resource, indivs are specialized relative to sp
    #BUT diff areas / types of networks vary in overall specialization, nestedness


netlvlZBin100 <- netlvl.Zer(netlvlBin100.df)
  #Nodf higher at all, overlap higher at all, (h2 all 0s)
  #vulnerability higher at all

netlvlZ100abun <- netlvl.Zer(netlvl100abun.df)


#modularity
plotModuleWeb(computeModules(vazNets100abun[[1]]$BM_2018)) #5 mods
  #interestingly, only 4 mods? w/ the abun weighted network. likely
plotModuleWeb(computeModules(vazNets100abun[[1]]$Turk_2018), labsize=0.75) # 5 mods, 6 w/ abun
plotModuleWeb(computeModules(vazNets100abun[[1]]$Barrios_2018),labsize=0.8) #4 mods, 4 w/ abun
plotModuleWeb(computeModules(vazNets100abun[[1]]$Turk28_2018)) # 5 mods for both
  #and yea, doing this on the binomial nets makes a ton of weird modules


#actual webs:

bmCol <- colLister(colnames(vazNets100abun[[1]]$BM_2018))
plotweb(vazNets100abun[[1]]$BM_2018, text.rot=90, col.high=bmCol)

turkCol <- colLister(colnames(vazNets100abun[[1]]$Turk_2018))
plotweb(vazNets100abun[[1]]$Turk_2018, text.rot=90, col.high=turkCol)

barCol <- colLister(colnames(vazNets100abun[[1]]$Barrios_2018))
plotweb(vazNets100abun[[1]]$Barrios_2018, text.rot=90, col.high=barCol)

turk28Col <- colLister(colnames(vazNets100abun[[1]]$Turk28_2018))
plotweb(vazNets100abun[[1]]$Turk28_2018, text.rot=90, col.high=turk28Col)

#RESULT 2:
  #Tur found i-sp less nested than expected, h2 lower than expected
  #DISAGREE
    #whether binomial or quant, NODF higher than nulls by a lot
    #h2 higher in quant, no variation in binomial

#######################

## netlvl w/in species

#######################

metrics.net <- c('NODF',
                 'weighted NODF',
                 'H2',
                 'robustness',
                 'vulnerability',
                 'niche overlap')

intraNets100 <- build.intraNetlvl(nets.ls,100,5,metrics.net,binary=FALSE)

save(intraNets100,file='data/intraNets100.RData')

intraZ <- calc.Z(intraNets,'SpSiteYr')
  #so trend seems to be v low HL niche overlap, 
  #higher NODF (though only a couple of these pass the threshold, all Zs are positive)
  #high h2 (most above 2z)

#RESULT 3:
  #tur found that w/in species, networks are not that nested
  #Mostly agree: 6/14 sig nested, compared to their 6/21

  #extras: niche overlap HL very neg, h2 very positive, vuln LL v neg
