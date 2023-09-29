#4_analysis: Analysis of bee individual - plant species networks from FFAR 2019
#field data. Loads in a objects created in 3_objects, then performs the main 
#analyses


#load in objects 
load('data/spec_RBCL_16s.RData')
load('data/quant_distSB5.RData')
load('data/spec_RBCL_16s.RData')
load('data/compSB.filtNHB.RData')
load('data/hubySB.RData')

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
sitelistSB <- unique(spec$SiteBloom)

#run an MRM on each site, then add descriptive columns to make the MRM output
#readable
MRMresultsSB <- MRMrunnerSB(sitelistSB,spec,dist.rbcl,dist.phylo.microbes,dist.beeSB)
MRMresultsSB$hub <- hubySB$HubDegree[match(MRMresultsSB$site,hubySB$site)]
MRMresultsSB$loc <- word(MRMresultsSB$site,1,sep="_")
MRMresultsSB$bloom <- word(MRMresultsSB$site,2,sep="_")
write.csv(MRMresultsSB,'MRMtable.csv')

#================================================
#Hub score analysis
#================================================
#phylogenetically corrected lmm of hub score on centroid distance
huby5SB.almer <- Almer(avgDist~degree*label+(1+degree|GenusSpecies)+(1|round/loc),
                       data=quant_distSB5NHB,
                       control = lmerControl(optimizer = 'bobyqa'),
                       A=list(GenusSpecies=A)
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
