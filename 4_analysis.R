#4_analysis: Analysis of bee individual - plant species networks from FFAR 2019
#field data. Loads in a objects created in 3_objects, then performs the main 
#analyses
rm(list=ls())

#load in objects 
load('data/spec_RBCL_16s_bloom.RData')
load('data/quant_distsSB.RData')
load('data/covarmatrix_community.Rdata')
load('data/hubySB.RData')
load('data/dists_micro_rbcl.Rdata')

source('2_functions.R')

## ================================================
## Are micro and rbcl community compositions correlated?
## ================================================

# Remove 2 observations of this one species that was not
#processed
spec <- spec[spec$GenusSpecies != "Lasioglossum sp. d",]

sitelistSB <- unique(spec$SiteBloom)

## run an MRM on each site, then add descriptive columns to make the
## MRM output readable
MRMresultsSB <- MRMrunnerSBbee(sitelistSB,spec, dist.rbcl,
                            dist.phylo.microbes,
                            bee.dist)

MRMresultsSB$hub <- hubySB$HubDegree[match(MRMresultsSB$site, hubySB$site)]
MRMresultsSB$loc <- word(MRMresultsSB$site, 1, sep="_")
MRMresultsSB$bloom <- word(MRMresultsSB$site, 2, sep="_")
write.csv(MRMresultsSB, 'MRMtable.csv')

## Repeating above, excluding sunflower
MRMresultsSB.NOSFbee <- MRMrunnerSBbee(sitelistSB,
                                       spec,dist.rbcl.NOSF,dist.phylo.microbes,
                                       bee.dist)
MRMresultsSB.NOSFbee$hub <-
  hubySB$HubDegree[match(MRMresultsSB.NOSFbee$site,
                         hubySB$site)]
MRMresultsSB.NOSFbee$loc <- word(MRMresultsSB.NOSFbee$site, 1, sep="_")
MRMresultsSB.NOSFbee$bloom <- word(MRMresultsSB.NOSFbee$site, 2, sep="_")
write.csv(MRMresultsSB.NOSFbee, 'MRMtable_noSF.csv')

##================================================
## Hub score analysis
##================================================
## phylogenetically corrected lmm of hub score on centroid distance
huby5SB.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSB5NHB,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SB.almer)

#Removing the random slope effect of bee species (i.e., changing)
#the first random effect term to "(1|GenusSpecies)" also does not effect the 
#output.
huby5SBNoSlope.almer <- Almer(avgDist~degree*label+(1|GenusSpecies)+(1|round/loc),
                       data=quant_distSB5NHB,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SBNoSlope.almer)


#Repeating this analysis with quant_distSB5 (i.e., including honey bees), 
  #quant_distSBCore5, or quant_distSBCore5NHB as the dataframe does not 
  #qualitatively change the output. These alternatives (listed as supplemental 
  #material) are below. 
huby5SB.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSB5,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SB.almer)


huby5SBCore.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSBCore5,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SBCore.almer)

huby5SBCoreNHB.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSBCore5NHB,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SBCoreNHB.almer)
