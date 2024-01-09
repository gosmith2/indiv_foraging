## 3_objects: Construction of objects used in analysis of bee
## individual - plant species networks from FFAR 2019 field data

## Requires packages found in 1_initiation, sources in functions from
## 2_functions, and loads some pre-constructed objects from "data"
## directory (networks, phylogenetic trees)

## load functions
source('2_functions.R')

## load network data
load('data/rbclNets.RData')
load('data/microNets.RData')
load('data/spec_RBCL_16s.RData')
load('data/trees.RData')
load('data/covarmatrix_community.RData')


#================================================
## Tweak objects for use in analysis
#================================================

## Add a bloom status (pre, during, and post-bloom) column to spec for easier grouping
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

## remove NAs, save for plotting
indivNet_rbcl1 <- lapply(indivNet_rbcl,function(x){
  keep.ls <- unlist(lapply(c(1:length(rownames(x))),function(y){
    all(!is.na(x[y,]))
  }))
  return(x[keep.ls,])
})
save(indivNet_rbcl1,file="data/indivNet_rbcl1.RData")

## make a binary version of the matrices
bin_rbcl <- lapply(indivNet_rbcl, function(x){
  x[x>0] <- 1
  return(x)
})
bin_micro <- lapply(indivNet_micro, function(x){
  x[x>0] <- 1
  return(x)
})

## save these objects for later use
save(bin_rbcl,file="data/bin_rbcl.RData")
save(bin_rbcl,file="data/bin_micro.RData")

## split site-level matrixes by bloom status
indivNet_rbclSB <- BloomSplit(indivNet_rbcl1,spec)
indivNet_microSB <- BloomSplit(indivNet_micro,spec)

#================================================
## prep matrices for MRMs
#================================================
#rbcl distance----------------------
rbcl <- colnames(spec)[grepl("RBCL", colnames(spec))]
comm.rbcl.indiv <- makeIndivComm(spec, rbcl)

dist.rbcl <- as.matrix(vegan::vegdist(comm.rbcl.indiv,
                                      "altGower"))

##Optional for Lauren: remove honeybees from rbcl dist
#no.apidae <- spec$UniqueID[spec$GenusSpecies=="Apis mellifera"]

#dist.rbcl.NHB <- dist.rbcl[!rownames(dist.rbcl) %in% no.apidae, !colnames(dist.rbcl) %in% no.apidae]


#make a version of the above without sunflower
comm.rbcl.indiv.NOSF <- makeIndivComm(spec, rbcl[-37])
dist.rbcl.NOSF <- as.matrix(vegan::vegdist(comm.rbcl.indiv.NOSF,
                                      "altGower"))

save(dist.rbcl, dist.rbcl.NOSF, file="data/rbcl_dists.Rdata")


#microbe distance------------------------
microbes <- colnames(spec)[grepl("16s", colnames(spec))]

comm.microbes.indiv <- makeIndivComm(spec, microbes)

## TAKES A LONGGGGG TIME
dist.phylo.microbes <- mod.unifrac(comm.microbes.indiv*100,
                                   tree.16s)
dist.phylo.microbes <- as.matrix(dist.phylo.microbes)

#list of apidae specimens
no.apidae <- spec$UniqueID[spec$GenusSpecies=="Apis mellifera"]

## drop the honeybee specimens
dist.phylo.microbes <- dist.phylo.microbes[
  !rownames(dist.phylo.microbes) %in% no.apidae,
  !colnames(dist.phylo.microbes) %in% no.apidae]

save(dist.phylo.microbes, dist.microbes, file="saved/distmats/indiv_16s.Rdata")

#bee community distance between sites------------------------
bee.commSB <- makeCommStructSB(spec.wild, "GenusSpecies")
dist.beeSB <- as.matrix(vegdist(bee.commSB$comm,"gower"))
save(dist.beeSB,file='data/dist.beeSB.RData')


#================================================
#Calculate centroid distances
#================================================
#rbcl_dis <- netToDisper(indivNet_rbcl1,spec,'rbcl')
rbclSB_dis <- netToDisper(indivNet_rbclSB,spec,'rbcl')

#micro_dis <- netToDisper(indivNet_micro,spec,'micro')
microSB_dis <- netToDisper(indivNet_microSB,spec,'micro')

#bin_rbcl_dis <- netToDisper(bin_rbcl,spec,'rbcl')
#bin_micro_dis <- netToDisper(bin_micro,spec,'micro')


#================================================
#Calculate hubScore, merge micro + rbcl
#================================================

#calculate hub score
hubySB <- netToDeg(indivNet_rbclSB,names(indivNet_rbclSB))
hubySB <- hubySB[hubySB$polN>4,] #limit to observations with more than 4 individuals
save(hubySB,file="data/hubySB.RData")

#combine distance objects, add some organizer columns
quant_distSB <- rbind(microSB_dis,rbclSB_dis)
quant_distSB$degree <- hubySB$HubDegree[match(quant_distSB$site,hubySB$site,nomatch=NA_integer_)]
quant_distSB$loc <- word(quant_distSB$site,1,sep="_")
quant_distSB$round <- word(quant_distSB$site,2,sep="_")

#exclude observations with fewer than 5 individuals
quant_distSB5 <- quant_distSB[quant_distSB$n>4,]
#save(quant_distSB5,file="data/quant_distSB5.RData")  

#exclude honeybees as a dataframe for testing
quant_distSB5NHB <- quant_distSB5[quant_distSB5$GenusSpecies!="Apis mellifera",]
save(quant_distSB5NHB,file="data/quant_distSB5NHB.RData")  


#================================================
#Repeat object construction with "core" microbes
#================================================
#define "core", get full column names of microbial ASVs in these genera
core.sp <- c("Rosenbergiella", "Pseudomonas", "Gilliamella", "Lactobacillus", "Caulobacter", "Snodgrassella", "Acinetobacter", "Corynebacterium", "Sphingomonas", "Commensalibacter", "Methylobacterium", "Massilia","Stenotrophomonas")
core.ls <- do.call(c,lapply(core.sp,function(x){
                              grep(x,names(spec))
                            }))

#split by new list
core.micro.sitebloom <- illumSplit(spec,"SiteBloom",core.ls)
core.micro.SB <- lapply(core.micro.sitebloom,function(x){
  filt <- x[,colSums(x,na.rm=T)>0]
  return(filt)
})

#calculate centroid distance
microSB_dis_core <- netToDisper(core.micro.SB,spec,'micro')

#combine core microbial distances with rbcl distances
quant_distSBCore <- rbind(rbclSB_dis,microSB_dis_core)
quant_distSBCore$degree <- hubySB$HubDegree[match(quant_distSBCore$site,hubySB$site,nomatch=NA_integer_)]
quant_distSBCore$loc <- word(quant_distSBCore$site,1,sep="_")
quant_distSBCore$round <- word(quant_distSBCore$site,2,sep="_")

#exclude observations with fewer than 5 individuals
quant_distSBCore5 <- quant_distSB[quant_distSBCore$n>4,]
#save(quant_distSBCore5,file="data/quant_distSBCore5.RData")  

#exclude honeybees as an optional dataframe for testing
quant_distSBCore5NHB <- quant_distSBCore5[quant_distSBCore5$GenusSpecies!="Apis mellifera",]
#save(quant_distSBCore5NHB,file="data/quant_distSBCore5NHB.RData")  

## **********************************************************************
## Generate bee phylogenetic distance matrix from dated tree
## **********************************************************************
## load the compiled tree
dated.tree <- read.tree("data/calibees-207t-95p-spruceup-chronos-strict-clock-125Ma.tre")

## cleaning so the tips match the species in the dataset

dated.tree$tip.label <- gsub("_BLX\\d+", "", dated.tree$tip.label)

dated.tree$tip.label <- gsub("_", " ", dated.tree$tip.label)

## dated.tree <- drop.tip(dated.tree, dated.tree$tip.label[duplicated(dated.tree$tip.label)])

dated.tree$tip.label[dated.tree$tip.label == "Lasioglossum spE"] <-
  "Lasioglossum sp. e"
dated.tree$tip.label[dated.tree$tip.label == "Sphecodes sp"] <-
  "Sphecodes spp."
dated.tree$tip.label[dated.tree$tip.label == "Xylocopa sp"] <-
  "Xylocopa varipuncta"
dated.tree$tip.label[dated.tree$tip.label == "Ashmeadiella aridula astragali"] <-
  "Ashmeadiella sp."
dated.tree$tip.label[dated.tree$tip.label == "Hylaeus mesillae"] <-
  "Hylaeus morphoD"
dated.tree$tip.label[dated.tree$tip.label == "Hylaeus rudbeckiae"] <-
  "Hylaeus morphoC"
dated.tree$tip.label[dated.tree$tip.label == "Hylaeus bisinuatus"] <-
  "Hylaeus morphoA"
dated.tree$tip.label[dated.tree$tip.label == "Hylaeus calvus"] <-
  "Hylaeus morphoB"
dated.tree$tip.label[dated.tree$tip.label == "Triepeolus heterurus"] <-
  "Triepeolus sp. a"
dated.tree$tip.label[dated.tree$tip.label == "Lasioglossum tegulare"] <-
  "Lasioglossum tegulariforme"

dated.tree$tip.label[dated.tree$tip.label == "Apis mellifera HC22"] <-
  "Apis mellifera"

no.tips <- unique(spec$GenusSpecies[!spec$GenusSpecies %in%
                                      dated.tree$tip.label])

phylo <- keep.tip(phy=dated.tree,
                  tip= dated.tree$tip.label[dated.tree$tip.label
                                            %in%
                                              unique(spec$GenusSpecies)])

#phylogenetic relatedness matrix
co.var.mat <- ape::vcv.phylo(phylo)

#phylogenetic distance matrix
bee.dist <- cophenetic.phylo(phylo)

save(co.var.mat, phylo, file="data/covarmatrix_community.Rdata")
save(bee.dist, phylo, file = "data/bee_dist.RData")
