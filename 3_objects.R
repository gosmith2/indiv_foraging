## 3_objects: Construction of objects used in analysis of bee
## individual - plant species networks from FFAR 2019 field data

## Requires packages found in 1_initiation, sources in functions from
## 2_functions, and loads some pre-constructed objects from "data"
## directory (networks, phylogenetic trees)

rm(list=ls())
## load functions
source('2_functions.R')

## load network data
load('data/rbclNets.RData')
load('data/microNets.RData')
load('data/spec_RBCL_16s.RData')
load('data/trees.RData')
load('data/covarmatrix_community.RData')

## ================================================
## create bloom status groupings
## ================================================

## Add a bloom status (pre, during, and post-bloom) column to spec for
## easier grouping
spec$BloomRound <- lapply(1:length(spec$SFBloomStatus), function (x){
  if(spec$SFBloomStatus[x] %in%
     c("peak","starting to bloom", "ending bloom")){
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
save(spec, file='data/spec_RBCL_16s_bloom.RData')

## remove rows of all NAs
indivNet_rbcl1 <- lapply(indivNet_rbcl, function(x){
  keep.ls <- unlist(lapply(c(1:length(rownames(x))),function(y){
    all(!is.na(x[y,]))
  }))
  return(x[keep.ls,])
})
save(indivNet_rbcl1,file="data/indivNet_rbcl1.RData")

## make a binary version of the matrices
bin_rbcl <- lapply(indivNet_rbcl, function(x){
  x[x > 0] <- 1
  return(x)
})
bin_micro <- lapply(indivNet_micro, function(x){
  x[x > 0] <- 1
  return(x)
})

## save these objects for later use
save(bin_rbcl,file="data/bin_rbcl.RData")
save(bin_rbcl,file="data/bin_micro.RData")

## split site-level matrixes by bloom status
indivNet_rbclSB <- BloomSplit(indivNet_rbcl1, spec)
indivNet_microSB <- BloomSplit(indivNet_micro, spec)

indivNet_rbclSB_noHB<- BloomSplit(indivNet_rbcl1,
                       spec[spec$GenusSpecies != "Apis mellifera",])
indivNet_microSB_noHB <- BloomSplit(indivNet_micro,
                        spec[spec$GenusSpecies != "Apis mellifera",])

## ================================================
## prep matrices for MRMs
## ================================================

## rbcl distance matrices, all sites
rbcl <- colnames(spec)[grepl("RBCL", colnames(spec))]
comm.rbcl.indiv <- makeIndivComm(spec, rbcl)
dist.rbcl <- as.matrix(vegan::vegdist(comm.rbcl.indiv,
                                      "altGower"))

## make a version of the distance matrix without sunflower
comm.rbcl.indiv.NOSF <- makeIndivComm(spec,
                        rbcl[rbcl !=  "RBCL:Asteraceae_Helianthus_annuus"])
dist.rbcl.NOSF <- as.matrix(vegan::vegdist(comm.rbcl.indiv.NOSF,
                                           "altGower"))

save(dist.rbcl, dist.rbcl.NOSF, file="data/rbcl_dists.Rdata")

## microbe distance matrices, all sites
microbes <- colnames(spec)[grepl("16s", colnames(spec))]
comm.microbes.indiv <- makeIndivComm(spec, microbes)

## TAKES A LONGGGGG TIME
dist.phylo.microbes <- mod.unifrac(comm.microbes.indiv*100,
                                   tree.16s)
dist.phylo.microbes <- as.matrix(dist.phylo.microbes)

## list of apis specimens
no.apis <- spec$UniqueID[spec$GenusSpecies == "Apis mellifera"]

## drop the honeybee specimens
dist.phylo.microbes.NHB <- dist.phylo.microbes[
  !rownames(dist.phylo.microbes) %in% no.apis,
  !colnames(dist.phylo.microbes) %in% no.apis]

save(dist.phylo.microbes, dist.phylo.microbes.NHB,
     file="data/indiv_16s.Rdata")

##================================================
## Calculate centroid distances
##================================================
## uses matrix split by bloom status
## with honey bees
rbclSB_dis <- netToDisper(indivNet_rbclSB, spec,'rbcl')
microSB_dis <- netToDisper(indivNet_microSB, spec,'micro')

## ================================================
## Calculate hubScore, merge micro + rbcl
## ================================================

## calculate hub score
hubySB <- netToDeg(indivNet_rbclSB,names(indivNet_rbclSB))
## limit to observations with more than 4 individuals
hubySB <- hubySB[hubySB$polN > 4,] 
save(hubySB, file="data/hubySB.RData")

## combine distance objects, add some organizer columns
quant_distSB <- rbind(microSB_dis, rbclSB_dis)
quant_distSB$degree <-
  hubySB$HubDegree[match(quant_distSB$site,hubySB$site,
                         nomatch=NA_integer_)]

quant_distSB$loc <- word(quant_distSB$site,1, sep="_")
quant_distSB$round <- word(quant_distSB$site,2, sep="_")

## exclude observations with fewer than 5 individuals
quant_distSB5 <- quant_distSB[quant_distSB$n > 4,]
save(quant_distSB5,file="data/quant_distSB5.RData")  

## exclude honeybees as a dataframe for testing
quant_distSB5NHB <- quant_distSB5[quant_distSB5$GenusSpecies != "Apis mellifera",]
save(quant_distSB5NHB,file="data/quant_distSB5NHB.RData")  

## ================================================
## Repeat object construction with "core" microbes
## ================================================
## define "core", get full column names of microbial ASVs in these genera
core.sp <- c("Rosenbergiella", "Pseudomonas", "Gilliamella",
             "Lactobacillus", "Caulobacter", "Snodgrassella",
             "Acinetobacter", "Corynebacterium", "Sphingomonas",
             "Commensalibacter", "Methylobacterium",
             "Massilia","Stenotrophomonas")

core.ls <- do.call(c,lapply(core.sp,function(x){
                              grep(x,names(spec))
                            }))

## split by new list
core.micro.sitebloom <- illumSplit(spec,"SiteBloom",core.ls)
core.micro.SB <- lapply(core.micro.sitebloom,function(x){
  filt <- x[,colSums(x,na.rm=T)>0]
  return(filt)
})

## calculate centroid distance
microSB_dis_core <- netToDisper(core.micro.SB,spec,'micro')

## combine core microbial distances with rbcl distances
quant_distSBCore <- rbind(rbclSB_dis,microSB_dis_core)
quant_distSBCore$degree <-
  hubySB$HubDegree[match(quant_distSBCore$site,hubySB$site,
                         nomatch=NA_integer_)]
quant_distSBCore$loc <- word(quant_distSBCore$site,1,sep="_")
quant_distSBCore$round <- word(quant_distSBCore$site,2,sep="_")

## exclude observations with fewer than 5 individuals
quant_distSBCore5 <- quant_distSBCore[quant_distSBCore$n>4,]
save(quant_distSBCore5,file="data/quant_distSBCore5.RData")  

## exclude honeybees as an optional dataframe for testing
quant_distSBCore5NHB <- quant_distSBCore5[
  quant_distSBCore5$GenusSpecies!="Apis mellifera",]
save(quant_distSBCore5NHB,file="data/quant_distSBCore5NHB.RData")  
