#4_analysis: Analysis of bee individual - plant species networks from FFAR 2019
#field data. Loads in a objects created in 3_objects, then performs the main 
#analyses


#load in objects 
load('data/spec_RBCL_16s.RData')
load('data/quant_distSB5.RData')
load('data/compSB.filtNHB.RData')
load('data/covarmatrix_community.Rdata')
load('data/hubySB.RData')
load('data/indiv_16s.Rdata')
load("data/rbcl_dists.Rdata")
load('data/bee_dist.RData')

source('2_functions.R')

#================================================
#Correlation between micro and rbcl composition
#================================================

###Are micro and rbcl community compositions correlated?

#add organizer columns, generate a list of site+bloom combos
spec$BloomRound <- lapply(1:length(spec$SFBloomStatus), function (x){
  if(spec$SFBloomStatus[x] %in% (c("peak","starting to bloom", "ending bloom"))){
    return("bloom")
  } else {
    if(spec$SFBloomStatus[x] == "before bloom"){
      return("before bloom")
    }else{
      if(spec$SFBloomStatus[x] == "after bloom"){
        return("after bloom")
      }else{
        if(spec$SampleRound[x] %in% c(1,2)){
          return("before bloom")
        } else{
          if(spec$SampleRound[x] %in% c(3,4)){
            return("bloom")
          } else{
            return("after bloom")
          }
        }
      }
    }
  }})
spec$SiteBloom <- paste(spec$Site, spec$BloomRound, sep="_")

#NEW line removing 2 observations of this one species that was not processed, 
#which broke MRM for one of the networks
spec <- spec[spec$GenusSpecies != "Lasioglossum sp. d",]

sitelistSB <- unique(spec$SiteBloom)

#run an MRM on each site, then add descriptive columns to make the MRM output
#readable
MRMresultsSB <- MRMrunnerSB(sitelistSB,spec,dist.rbcl,dist.phylo.microbes,dist.beeSB)
MRMresultsSB$hub <- hubySB$HubDegree[match(MRMresultsSB$site,hubySB$site)]
MRMresultsSB$loc <- word(MRMresultsSB$site,1,sep="_")
MRMresultsSB$bloom <- word(MRMresultsSB$site,2,sep="_")
write.csv(MRMresultsSB,'MRMtable.csv')
  #Above chunk is unchanged. NOTE: dist.SB (bee community composition distance) is in there, 
  #but isn't actually called. "bee.dist" is the new species-level phylogenetic distance matrix for bees. 

#No SF
MRMresultsSB.NOSF <- MRMrunnerSB(sitelistSB,spec,dist.rbcl.NOSF,dist.phylo.microbes,dist.beeSB)
MRMresultsSB.NOSF$hub <- hubySB$HubDegree[match(MRMresultsSB.NOSF$site,hubySB$site)]
MRMresultsSB.NOSF$loc <- word(MRMresultsSB.NOSF$site,1,sep="_")
MRMresultsSB.NOSF$bloom <- word(MRMresultsSB.NOSF$site,2,sep="_")

#No SF, bee phylo included
MRMresultsSB.NOSFbee <- MRMrunnerSBbee(sitelistSB,spec,dist.rbcl.NOSF,dist.phylo.microbes,bee.dist)
MRMresultsSB.NOSFbee$hub <- hubySB$HubDegree[match(MRMresultsSB.NOSFbee$site,hubySB$site)]
MRMresultsSB.NOSFbee$loc <- word(MRMresultsSB.NOSFbee$site,1,sep="_")
MRMresultsSB.NOSFbee$bloom <- word(MRMresultsSB.NOSFbee$site,2,sep="_")


#================================================
#Hub score analysis
#================================================
#phylogenetically corrected lmm of hub score on centroid distance
huby5SB.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSB5NHB,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=co.var.mat)
)
summary(huby5SB.almer)
  #repeating this analysis with quant_distSB5NHB, quant_distSBCore5, or 
  #quant_distSBCore5NHB as the dataframe does not qualitatively change the output


#repeating with a non-phylogenetic lmm gives qualitatively identical results
hubySB5.lme4 <- lmer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                     data=quant_distSB5NHB,
                     control = lmerControl(optimizer = 'bobyqa')
)
summary(hubySB5.lme4)
