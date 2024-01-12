# This file contains functions used in the script 3_analysis

'%!in%' <- function(x,y)!('%in%'(x,y))


illumSplit <- function(data,split,len){
  ## splits illumina data (which has already been attached to
  ## character data) by a split column, then converts the data into a
  ## list of network objects
  
splitList <- lapply(unique(data[,split]),function(x){
  site <- subset(data,data[,split]==x)
  site.mx <- site[,len]
  rownames(site.mx) <- site$UniqueID
  site.mx <- t(site.mx)
  
})
  names(splitList) <- unique(data[,split])
  return(splitList)
}

zscorer <- function(data, score){ #
  ## sub-function; calculates z scores
  max <- max(data[,score],na.rm=T)
  mean <- mean(data[,score],na.rm=T)
  sd <- sd (data[,score],na.rm=T)
  return((max-mean)/sd)
}


makeIndivComm <- function(spec, col.sp.names){ #
  ## converts the whole dataset (spec) into a community matrix based
  ## on a list of column names
  
  comm.indiv <- spec[, col.sp.names]
  rownames(comm.indiv) <- spec$UniqueID
  ## find those that were not screened to drop them
  not.screened <- apply(comm.indiv, 1, function(x) all(is.na(x)))
  comm.indiv <- comm.indiv[!not.screened, ]
  
  comm.indiv[is.na(comm.indiv)] <- 0
  comm.indiv  <- bipartite::empty(comm.indiv)
  return(comm.indiv)
}

netToDisper <- function(data, meta, label){ 
  ## calculates dispersion (in this case centroid distance) from
  ## network objects meta refers to a larger dataset containing
  ## metadata for the specimens, not currently used label refers to a
  ## tag for which the dispersions are calculated (e.g., rbcl reads)
  
  sites <- lapply(names(data), function(x){
    site.net <- data[[x]] 
    ids <- colnames(site.net)    
    sp.ls <- meta$GenusSpecies[meta$UniqueID %in% ids]

    if(length(sp.ls) > 0){
      dis <- betadisper(vegdist(t(site.net),na.rm=T),sp.ls,type='centroid')

      ## extract averages and number of individuals for each species
      avgs <- lapply(unique(dis$group), function(y){
        avgDist <- mean(dis$distances[dis$group == y])
        len <- length(dis$distances[dis$group == y])
        site.data <- data.frame('GenusSpecies' = y,
                                'avgDist' = avgDist,
                                'n' = len)
        return(site.data)
      })

      ## combine averages into a dataframe
      site.dat <- do.call(rbind, avgs)
      site.dat$site <- x
      site.dat$label <- label
      print(x)
      
    } else {
      site.dat <- data.frame('GenusSpecies'=NA,
                             'avgDist' = NA,
                             'n' = NA,
                             "site" = x,
                             "label" = label)
    }
    print(x)
    return(site.dat)
  })
  
  return(do.call(rbind, sites))
  
}


netToDeg <- function(data,sitelist){ #
 #calculates hub score metric (relative degree) for each network within a "sitelist" object
  sites <- lapply(sitelist,function(x){
    mets <- specieslevel(data[[x]],index=c('degree'))
    plants <- mets$'lower level'
    deg <- zscorer(plants,'degree')
    site <- data.frame('site' = x,"HubDegree" = deg,"plantN" = dim(data[[x]])[[1]],"polN" = dim(data[[x]])[[2]])
    
    
    return(site)
  })
  do.call(rbind,sites)
}


#Call a specific set of ggplot theme elements
theme_ms <- function(base_size=14, base_family="sans") { #
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5,
                                      margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = 'white'),
            plot.background = element_rect(colour = NA, fill = 'white'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'black'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'black'),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            ##legend.position = "none",
            legend.background = element_rect(fill ='white'),
            legend.text = element_text(color = 'black'),
            legend.key = element_rect(colour = NA, fill = 'white'),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
                                        #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'black'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="white"),
            strip.text = element_text(face="bold", colour =
                                                     'black')
            ))
}



BloomSplit <- function(network.list,meta=spec){ #!
  #splits site-level networks into site+bloom status level networks (i.e., 
  #splits "Turk" into "Turk_before bloom", "Turk_bloom" etc.)
  
  newlist <- do.call(c,lapply(names(network.list), function(x){
    ## narrow down to just the right site-level network
    site <- network.list[[x]] 
    
    rounds <- lapply(unique(meta$SiteBloom[meta$Site==x]), function(y){
      
      ## subset that site network to just specimens collected during a
      ## given sample round
      bloomround <- site[,colnames(site) %in%
                          meta$UniqueID[meta$SiteBloom==y]]
      print(y)
      return(bloomround)
    })
    
    ##name the objects in the rounds list
    names(rounds) <- unique(meta$SiteBloom[meta$Site==x])
    
    ## remove any networks with dim 0 (i.e., no pollinators)
    clean.ls <- unlist(lapply(rounds, function(z){
      if(length(colnames(z))==0){
        return(FALSE)
      }else{
        if(dim(z)[[2]]==0){
          return(FALSE)
        }else{
          return(TRUE)
        }
      }
      
    }))
    rounds <- rounds[clean.ls]                     
    
    return(rounds)
  }))
  return(newlist)
}

MRMrunner <- function(sitelist, spec.data, rbcl, micro, bees){
  #Runs MRM tests for each network in the site list
  #Depricated version
  siteSR.list <- do.call(rbind,lapply(sitelist,function(x){
    print(x)
    #browser()
    
    ID.list <- spec.data$UniqueID[spec.data$SiteSample==x]
    
    if(length(rownames(rbcl)[rownames(rbcl)%in%ID.list])<3){
      
      res <- data.frame("site"=x,'coef' = NA, "p" = NA, "r2"=NA)
      
    }else{
      
      rbcl.filt <- rbcl[rownames(rbcl)%in%ID.list,colnames(rbcl)%in%ID.list]
      micro.filt <- micro[rownames(micro)%in%ID.list,colnames(micro)%in%ID.list]
      #bees.filt <- bees[x,x]
      
      summary <- MRM(as.dist(micro.filt) ~ as.dist(rbcl.filt), nperm=10^4)
      
      res <- data.frame("site"=x,'coef' = summary$coef[2], "p" = summary$coef[4], "r2" = summary$r.squared[[1]])
      #browser()
    }
    
    return(res)
  }))
  
  return(siteSR.list)
}

MRMrunnerSB <- function(sitelist, spec.data, rbcl, micro, bees){ #
  #Runs MRM tests for each network in the site list
  #This version is the old version used for the last submission, which now gives unequal matrix errors

  siteSB.list <- do.call(rbind,lapply(sitelist,function(x){
    print(x)
    #browser()
    
    ID.list <- spec.data$UniqueID[spec.data$SiteBloom==x]

    #skip networks that are too small
    if(length(rownames(rbcl)[rownames(rbcl)%in%ID.list])<3){
      
      res <- data.frame("site"=x,'coef' = NA, "p" = NA, "r2"=NA)
      
    }else{
      #filter to just specimens at the given network
      rbcl.filt <- rbcl[rownames(rbcl)%in%ID.list,colnames(rbcl)%in%ID.list]
      micro.filt <- micro[rownames(micro)%in%ID.list,colnames(micro)%in%ID.list]
      #bees.filt <- bees[x,x]

      #Run an MRM on the filtered matrices
      summary <- MRM(as.dist(micro.filt) ~ as.dist(rbcl.filt), nperm=10^4)

      #extract relevant summary data into a 1 row dataframe to combine with the rest
      res <- data.frame("site"=x,'coef' = summary$coef[2], "p" = summary$coef[4], "r2" = summary$r.squared[[1]])
      #browser()
    }
    
    return(res)
  }))
  
  return(siteSB.list)
}

MRMrunnerSB <- function(sitelist, spec.data, rbcl, micro, bees){
  #New version resolving the unequal matrix size issue. Doesn't add bees to the MRM (the one below this does)
  
  siteSB.list <- do.call(rbind,lapply(sitelist,function(x){
    print(x)
    
    ID.list <- spec.data$UniqueID[spec.data$SiteBloom==x]
    
    #browser()
    #site <- spec.data[spec.data$SiteBloom==x,]
    
    if(length(rownames(rbcl)[rownames(rbcl)%in%ID.list])<3){
      
      res <- data.frame("site"=x,'coef' = NA, "p" = NA, "r2"=NA)
      print("skipped1")
      
    }else{
      print("started")

      rbcl.filt <- rbcl[rownames(rbcl)%in%ID.list,colnames(rbcl)%in%ID.list]
      micro.filt <- micro[rownames(micro)%in%ID.list,colnames(micro)%in%ID.list]
      
      micro.filt.sm <- micro.filt[rownames(micro.filt)%in%rownames(rbcl.filt),colnames(micro.filt)%in%colnames(rbcl.filt)]
      rbcl.filt.sm <- rbcl.filt[rownames(rbcl.filt)%in%rownames(micro.filt),colnames(rbcl.filt)%in%colnames(micro.filt)]
      
      if(length(rownames(rbcl.filt.sm))<3){
        
        res <- data.frame("site"=x,'coef' = NA, "p" = NA, "r2"=NA)
        print("skipped2")
        
      } else {
        
        #print("filtered")
        #browser()
        
        summary <- MRM(as.dist(micro.filt.sm) ~ as.dist(rbcl.filt.sm), nperm=10^4)
        
        res <- data.frame("site"=x,'coef' = summary$coef[2], "p" = summary$coef[4], "r2" = summary$r.squared[[1]])
        print("tested")
      }
    }
    
    return(res)
  }))
  
  return(siteSB.list)
}

MRMrunnerSBbee <- function(sitelist, spec.data, rbcl, micro, bees){
  ##Runs MRM tests for each network in the site list
  
  ## Version using the new filtering and including bee phylo in the model (when there are multiple bee species present)
  
  siteSB.list <- do.call(rbind,lapply(sitelist,function(x){
    print(x)
    
    ID.list <- spec.data$UniqueID[spec.data$SiteBloom==x]
        
    if(length(rownames(rbcl)[rownames(rbcl)%in%ID.list])<3){
      
      res <- data.frame("site"=x,'rbcl_coef' = NA, "rbcl_p" = NA, 'bee_coef' = NA, "bee_p" = NA,"r2"=NA)
      print("skipped1")
      
    }else{
      print("started")
      
      rbcl.filt <- rbcl[rownames(rbcl)%in%ID.list,colnames(rbcl)%in%ID.list]
      micro.filt <- micro[rownames(micro)%in%ID.list,colnames(micro)%in%ID.list]
      
      micro.filt.sm <- micro.filt[rownames(micro.filt)%in%rownames(rbcl.filt),colnames(micro.filt)%in%colnames(rbcl.filt)]
      rbcl.filt.sm <- rbcl.filt[rownames(rbcl.filt)%in%rownames(micro.filt),colnames(rbcl.filt)%in%colnames(micro.filt)]
      
      if(length(rownames(rbcl.filt.sm))<3){
        
        res <- data.frame("site"=x,'rbcl_coef' = NA, "rbcl_p" = NA, 'bee_coef' = NA, "bee_p" = NA,"r2"=NA)
        print("skipped2")
        
      } else {

        bees.filt <- as.matrix(do.call(rbind, lapply(rownames(rbcl.filt.sm), function(x){
          comp.list <- unlist(lapply(colnames(rbcl.filt.sm), function(y){

            x.sp <- spec.data$GenusSpecies[spec.data$UniqueID==x]
            y.sp <- spec.data$GenusSpecies[spec.data$UniqueID==y]
            return(bees[x.sp,y.sp])
          }))
        })))
        rownames(bees.filt) <- rownames(rbcl.filt.sm)
        colnames(bees.filt) <- rownames(rbcl.filt.sm)

        if(all(bees.filt[1,1] == bees.filt)){
          summary <- MRM(as.dist(micro.filt.sm) ~ as.dist(rbcl.filt.sm), nperm=10^4)
          print("w/o bees")
          res <- data.frame("site"=x,'rbcl_coef' = summary$coef[2], "rbcl_p" = summary$coef[4], 'bee_coef' = NA, "bee_p" = NA, "r2" = summary$r.squared[[1]])
        } else {
          summary <- MRM(as.dist(micro.filt.sm) ~ as.dist(rbcl.filt.sm) + as.dist(bees.filt), nperm=10^4)
          print("w/ bees")
          res <- data.frame("site"=x,'rbcl_coef' = summary$coef[2], "rbcl_p" = summary$coef[5], 'bee_coef' = summary$coef[3], "bee_p" = summary$coef[6], "r2" = summary$r.squared[[1]])
          print("tested")
        }
       
      }
    }
    
    return(res)
  }))
  
  return(siteSB.list)
}

samp2site.spp <- function(site,spp,abund) { #
  ## sub function, helps conversion from sample to site level
  
    x <- tapply(abund, list(site=site,spp=spp), sum)
    x[is.na(x)] <- 0
    return(x)
}

mod.unifrac <- function (comm, tree){
  ## Calculates Unifrac distances from a community matrix and a
  ## phylogenetic tree
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute UniFrac")
  }
  comm <- as.matrix(comm)
  s <- nrow(comm)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(comm)
  colnames(phylodist) <- rownames(comm)
  comm_comb <- matrix(NA, s * (s - 1)/2, ncol(comm))
  colnames(comm_comb) <- colnames(comm)
  i <- 1
  for (l in 1:(s - 1)) {
    for (k in (l + 1):s) {
      comm_comb[i, ] <- comm[l, ] + comm[k, ]
      i <- i + 1
    }
  }
  pdcomm <- pd.mod(comm, tree)
  pdcomm_comb <- pd.mod(comm_comb, tree)
  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdcomm[l, "PD"]
    for (k in (l + 1):s) {
      pdk <- pdcomm[k, "PD"]
      pdcomb <- pdcomm_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] = (pdcomb - pdsharedlk)/pdcomb
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}

pd.mod <- function (samp, tree, include.root = TRUE){
  ## Sub function used in mod.unifrac, from picante
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (include.root) {
    tree <- node.age(tree)
  }
  species <- colnames(samp)
  SR <- rowSums(ifelse(samp > 0, 1, 0))
  nlocations <- dim(samp)[1]
  nspecies <- dim(samp)[2]
  PDs <- rep(NA, nlocations)
  for (i in 1:nlocations) {
    present <- species[samp[i, ] > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in%
                                         present))]
    if (length(present) == 0) {
      PDs[i] <- 0
    }
    else if (length(present) == 1) {
      if (!include.root) {
        warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.")
        PDs[i] <- NA
      }
      else {
        PDs[i] <- tree$ages[which(tree$edge[, 2] == which(tree$tip.label ==
                                                          present)[1])]
      }
    }
    else if (length(treeabsent) == 0) {
      PDs[i] <- sum(tree$edge.length)
    }
    else {
      sub.tree <- drop.tip(tree, treeabsent)
      if (include.root) {
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[,2] %in%
                                               which(tree$tip.label %in% present))])
        PDs[i] <- sum(sub.tree$edge.length) + (orig.tree.depth -
                                               sub.tree.depth)
      }
      else {
        PDs[i] <- sum(sub.tree$edge.length)
      }
    }
  }
  PDout <- data.frame(PD = PDs, SR = SR)
  rownames(PDout) <- rownames(samp)
  return(PDout)
}
