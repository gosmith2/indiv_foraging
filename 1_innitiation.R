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
load('c:/Users/thebo/Dropbox/sunflower_saved/data/RBCL2018.RData')

#load in the compiled data with info on each sample
complete <- read.csv('c:/Users/thebo/Dropbox/sunflower_saved/data/relational/relational/traditional/specimens-complete.csv')

#subset the larger dataset to just samples that were processed in rbcl
sampled <- complete[complete$UniqueID %in% RBCL2018$UniqueID,]


#quick data cleaning, adding some organizing columns, then merging the dataframes
samp.df <- dat.clean(sampled)
samp.df$SiteYr <- paste(samp.df$Site, substring(samp.df$Date,7,10), sep = "_")
samp.df$SpSiteYr <- paste(samp.df$GenusSpecies, samp.df$SiteYr, sep = "_")

sampReads.df <- left_join(samp.df,RBCL2018, by = 'UniqueID')
#samp.df <- dat.dates(samp.df)


#length object to isolate just the reads
len <- c((length(samp.df)+1):length(sampReads.df))

siteSplit.ls <- lapply(unique(sampReads.df$SiteYr),function(x){
  site <- sampReads.df[sampReads.df$SiteYr==x,]
})
names(siteSplit.ls) <- unique(sampReads.df$SiteYr)


#make an object of just the read numbers for community stuff

nets.ls <- lapply(siteSplit.ls,function(x){
  nets <- illumSplit(x,"SpSiteYr",len)
})

save(nets.ls,file='data/nets.RData')
#nets <- illumSplit(siteSplit.ls,"SpSiteYr",len)


#to randomize indivs w/in species (submatrices), lapply vaznull on nets, then cbind sp together

cores <- 1

vazNets <- vaznullFromSp(nets.ls,6,5, cores)
save(vazNets,file='data/vazNets.RData')

vazNetsBin <- vaznullFromSp(nets.ls,6,5, cores, binary = TRUE)


splvlLS<- lapply(vazNets, function(x){
  mclapply(x,function(y){
    mets <- specieslevel(y)
  },mc.cores=cores)
})

splvlLSBin<- lapply(vazNetsBin, function(x){
  mclapply(x,function(y){
    mets <- specieslevel(y)
  },mc.cores=cores)
})


save(splvlLS, file='data/splvlLS.RData')

metrics <- names(splvlLS[[1]]$BM_2018$`higher level`)

netlvl <- lapply(vazNets, function(x){
  mclapply(x,function(y){
    mets <- networklevel(y, metrics)
  },cores=cores)
},mc.cores=cores)


splvlComp <- splvl.dfer(splvlLS)

splvlCompBin <- splvl.dfer(splvlLSBin)


splvlZ <- splvl.Zer(splvlComp,"SpSiteYr")

splvlZBin <- splvl.Zer(splvlCompBin,"SpSiteYr")


spNets <- build.spNet(nets.ls)
  
spNetsSum <- sum.spNet(spNets)

splvlZ$spdeg <- spNetsSum$degree[match(splvlZ$interaction,spNetsSum$SpSiteYr)]
  #yea, species have much higher degree than indivs (duh)



#######################

## netlvl

#######################


load('data/vazNets.RData')

metrics.net <- c('NODF',
                'H2',
                'robustness',
                'vulnerability',
                'niche overlap')

netlvl.df <- build.netlvl(vazNets,metrics.net) 

netlvlZ <- netlvl.Zer(netlvl.df)
  #so NODF seems higher, h2 seems higher, overlapHL seems lower, and 
  #tur found NODF lower, h2 lower
  #so... even for specialists in presence of burst of their resource, indivs are specialized relative to sp
    #BUT diff areas / types of networks vary in overall specialization, nestedness


#######################

## netlvl w/in species

#######################

metrics.net <- c('NODF',
                 'H2',
                 'robustness',
                 'vulnerability',
                 'niche overlap')

intraNets <- build.intraNetlvl(nets.ls,6,5,metrics.net,binary=FALSE)

intraZ <- calc.Z(intraNets,'SpSiteYr')
  #so trend seems to be v low HL niche overlap, 
  #higher NODF (though only a couple of these pass the threshold, all Zs are positive)
  #high h2 (most above 2z)

