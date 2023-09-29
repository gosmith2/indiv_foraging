#6_supplement: Descriptive statisitcs of bee individual - plant species networks 
  #from FFAR 2019
  #field data. Loads in a few objects created in 3_objects, then calculates some
  #descriptors and lists that help describe the dataset


#load in objects 
load('data/spec_RBCL_16s.RData')
load('data/bin_micro.RData')
load('data/bin_rbcl.RData')


#================================================
#Supplemental descritpors of the dataset
#================================================

#---------------------------
#Number of sunflower-captured individuals carrying non-sunflower pollen 
#---------------------------
sfonly <- spec %>% filter(TransectType == "SF", Pollen == 1); 872:994
sfonly.mat <- as.matrix(sfonly[,872:944])
sfonly.mat1 <- sfonly.mat1[rowSums(sfonly.mat1,na.rm=T)>0,]
sfonly.mat1n <- sfonly.mat1[,colnames(sfonly.mat1)!="RBCL:Asteraceae_Helianthus_annuus"]
length(rownames(sfonly.mat1n[rowSums(sfonly.mat1n,na.rm=T)>0,]))
#220/361 individuals from SF transects had non-sf pollen

#---------------------------
#Number of margin-captured individuals carrying sunflower pollen 
#---------------------------
sfno <- spec %>% filter(TransectType != "SF", Pollen == 1,SFBloomStatus %in% c("peak","starting to bloom","ending bloom"))
sfno.mat <- as.matrix(sfno[,872:994])
sfno.mat[sfno.mat>0] <- 1
sum(sfno.mat[,colnames(sfno.mat)=="RBCL:Asteraceae_Helianthus_annuus"],na.rm=T)
#26/48 individuals from margin transects during sf bloom had SF pollen
###split site-level networks into site+sample round networks

#---------------------------
#Number of microbes/plants present in diff percentages of the sample (Table 3)
#---------------------------

#microbes-------------------
#put all sites together into 1 matrix
bin_micro_all <- do.call(cbind,bin_micro) 
#Save this object for Supplementary Table S2.1
write.csv(rowSums(bin_micro_all,na.rm=T),'microASVlist.csv')

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

#plants--------------------------
bin_rbcl_all <- do.call(cbind,bin_rbcl) #put all sites together into 1 matrix
#Save this object for Supplementary Table S2.2
write.csv(rowSums(bin_rbcl_all,na.rm=T),'rbclASVlist.csv')

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

#---------------------------
#Presence of floral specialist microbes in the bee microbiome samples
#---------------------------
#floral associated microbes (from Vannette annual review 2020 and Alvarez 2012)

flr_bact_list <- c("Acinetobacter", "Rosenbergiella", "Metschnikowia", "Pseudomonas", "Pantoea")

flr_bacts <- grep(paste(flr_bact_list, collapse="|"),rownames(bin_micro_all))
#just three rows: 146, 148, 160


#---------------------------
#Bee species by flower capture location (Supplementary Table S1)
#---------------------------
tableCaps <- do.call(rbind,lapply(unique(spec$PlantGenusSpecies), function(x){
  proc <- spec[spec$Pollen==T|spec$Gut==T,]
  print(x)
  plant <- proc[proc$PlantGenusSpecies==x,]
  pols <- do.call(cbind,lapply(unique(proc$GenusSpecies),function(y){
    #browser()
    plant.df <- data.frame("n" = length(plant$UniqueID[plant$GenusSpecies==y]))
    rownames(plant.df) <- x
    colnames(plant.df) <- y
    #print(y)
    return(plant.df)
  }))
  return(pols)
}))
tableCaps <- tableCaps[order(rowSums(tableCaps),decreasing=T),order(colSums(tableCaps),decreasing=T)]
write.csv(tableCaps,"CaptureTable.csv")

#---------------------------
#(Supplementary Table S1)
#---------------------------

