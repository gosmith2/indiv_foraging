# This file contains functions used in the script 3_analysis

'%!in%' <- function(x,y)!('%in%'(x,y))

fix.white.space <- function(d) {
  #cleans character data by removing white space
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)
  
  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))
  
  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}

dat.clean <- function(spec.dat) {
  #cleans data of common issues from and formatting and blanks
  spec.dat$GenusSpecies <- fix.white.space(paste(spec.dat$Genus,
                                                 spec.dat$Species,
                                                 spec.dat$SubSpecies))
  
  spec.dat$PlantGenusSpecies <-  fix.white.space(paste(spec.dat$PlantGenus,
                                                       spec.dat$PlantSpecies,
                                                       spec.dat$PlantVar,
                                                       spec.dat$PlantSubSpecies))
  
  spec.dat$Int <-  fix.white.space(paste(spec.dat$GenusSpecies,
                                         spec.dat$PlantGenusSpecies))
  spec.dat$IntGen <-  fix.white.space(paste(spec.dat$Genus,
                                            spec.dat$PlantGenus))
  return(spec.dat)
}

dat.dates <- function(spec.dat) {
  #cleans date data by reformatting
  spec.dat$Date <- as.Date(spec.dat$Date, format='%m/%d/%y')
  spec.dat$Doy <- as.numeric(strftime(spec.dat$Date, format='%j'))
  spec.dat$Year <- as.numeric(format(spec.dat$Date,'%Y'))
  return(spec.dat)
}

dat.rm.blanks <- function(spec.dat) {
  #cleans character data by removing blanks
  spec.dat <- spec.dat[spec.dat$PlantGenusSpecies != "",]
  spec.dat <- spec.dat[spec.dat$GenusSpecies != "",]
  return(spec.dat)
}


illumSplit <- function(data,split,len){
  ##splits illumina data (which has already been attached
  #to character data) by a split column, then converts the 
  #data into a list of network objects
  
  splitList <- lapply(unique(data[,split]),function(x){
    
    site <- subset(data,data[,split]==x)
    
    site.mx <- site[,len]
    
    rownames(site.mx) <- site$UniqueID
    
    site.mx <- t(site.mx)
    
    
  })
  
  names(splitList) <- unique(data[,split])
  return(splitList)
}

netSums <- function(data, binary=FALSE){
  #takes a list of network objects (e.g., from 
  #illumSplit) and calculates summary stats using 
  #specieslevel. can be calculated either with reads
  #acting as observations, or by converting the reads
  #to presence/absence (binary = TRUE)
  
  ls <- lapply(names(data),function(x){
    if(binary==FALSE){
      mets <- specieslevel(data[[x]])
      
      df <- data.frame("UniqueID"=rownames(mets$`higher level`))
      df <- cbind(df,mets$`higher level`)
      df$Site <- x
      return(df)
    } else {
      
      bin <- ifelse(data[[x]]>0,1,0)
      
      mets <- specieslevel(bin)
      
      df <- data.frame("UniqueID"=rownames(mets$`higher level`))
      df <- cbind(df,mets$`higher level`)
      df$Site <- x
      return(df)
      
    }
  })
  
  do.call(rbind,ls)
}

vaznullFromSp <- function(data, iterations, threshold, binary=FALSE, cores){
  #data should be a nested list: each element in main list should be a list of networks (i.e., species networks within site within a single "data" list object)
  #iterations is the number of complete null networks
  #threshold is the minimum number of individuals to have of a given species w/in a site to randomize them.
  
  iter <- c(1:iterations)
  
  #subset down to site
  nets <- mclapply(data,function(x){
    vazSub1(x,threshold=threshold,binary= binary, obs=TRUE)
  },mc.cores=cores)
  
  vazes.ls <- lapply(iter,function(x){
    vazes <- mclapply(data, function(y){      
      vazSub1(y, threshold=threshold,binary= binary, obs=FALSE)
    },mc.cores=cores)
  })
  
  nets.ls <- list(nets)
  nets.ls <- c(nets.ls,vazes.ls)
  return(nets.ls)
}


vazSub <- function(data, threshold, obs = FALSE, binary){
  #data should be a single species network
  #threshold is the minimum number of individuals to randomize
  #obs is whether or not nulls should be created, or whether the data is just being reorganized
  if(obs == FALSE & length(colnames(data))>=threshold){
    if(binary == FALSE){
      sp <- data[(rowSums(data)>0),]
      vaz <- vaznull(1,sp)
      rownames(vaz[[1]]) <- rownames(sp)
      colnames(vaz[[1]]) <- colnames(sp)
      return(vaz[[1]])
      
    } else {
      sp <- data[(rowSums(data)>0),]
      bin <- ifelse(sp>0,1,0)
      vaz <- vaznull(1,bin)
      rownames(vaz[[1]]) <- rownames(sp)
      colnames(vaz[[1]]) <- colnames(sp)
      return(vaz[[1]])
    }
    
  } else if(all(dim(data)>1)){
    if(binary==FALSE){
      sp <- data[(rowSums(data)>0),]
      return(sp)
    } else{
      sp <- data[(rowSums(data)>0),]
      bin <- ifelse(sp>0,1,0)
      return(bin)
    }
    
  } else {
    if(binary==FALSE){
      sp <- data
      return(sp)
    } else{
      sp <- data
      bin <- ifelse(data>0,1,0)
      return(bin)
    }
  }
}


vazSub1 <- function(data, threshold, obs, binary, visitweight = FALSE) { 
  ##sub-function
  
  #down to species
  net <- lapply(data,function(y){
    
    #if(null=="vaz"){
    sp <- vazSub(y, threshold=threshold,obs=obs, binary=binary)
    #} else {
    #  sp <- basSub(y, threshold=threshold,obs=obs, binary=binary)
    #}
    
    if(length(rownames(sp))<length(rownames(y))){
      extras <- matrix(nrow = (length(rownames(y))-length(rownames(sp))), 
                       ncol = length(colnames(sp)))
      rownames(extras) <- rownames(y)[!(rownames(y) %in% rownames(sp))]
      colnames(extras) <- colnames(sp)
      extras[anyNA(extras)] <- 0
      
      sp1 <- rbind(sp,extras)
      sp1 <- sp1[row.names(y),,drop=FALSE]
      return(sp1)
      
    } else {
      return(sp)
    }
    
  })
  
  do.call(cbind,net)
}


splvl.dfer <- function(data){ 
  #builds a df from species level metrics
  
  do.call(rbind,lapply(c(1:length(data)),function(x){
    
    iteration <- data[[x]]
    
    #narrow down to site
    site <- lapply(iteration,function(y){
      pol <- y$'higher level'
      pol$UniqueID <- rownames(pol)
      
      left_join(pol,samp.df, by='UniqueID')
      
    })
    
    sites <- do.call(rbind, site)
    sites$iteration <- x
    
    return(sites)
    
  }))
}

splvl.Zer <- function(data, split){
  #generates a summary z score df from species level network metrics
  
  do.call(rbind,lapply(unique(data[,split]),function(x){
    inter <- subset(data,data[,split]==x)
    
    means <- do.call(cbind,lapply(metrics,function(y){
      met <- unlist(lapply(unique(inter$iteration),function(z){
        iter <- inter[inter$iteration==z,]
        mean(iter[,y])
      }))
      met1 <- data.frame(y = met)
      names(met1)<- c(y)
      return(met1)
    }))
    
    zscores <- do.call(cbind,lapply(names(means), function(x){
      zscore <- (means[1,x] - mean(means[,x]))/
        (sd(means[,x])+10^-10)
      
      zscore1 <- data.frame(x = zscore)
      names(zscore1) <- c(x)
      return(zscore1)
    })
    )
    
    zscores$interaction <- x
    zscores$indivs <-length(inter$UniqueID)/length(unique(inter$iteration))
    return(zscores)
  })
  )
}


zscorer <- function(data, score){
  #sub-function; calculates z scores
  max <- max(data[,score],na.rm=T)
  mean <- mean(data[,score],na.rm=T)
  sd <- sd (data[,score],na.rm=T)
  return((max-mean)/sd)
}

build.spNet <- function(data){
  #sumarizes species level metrics and builds a df
  lapply(data,function(x){
    sps <- lapply(x,function(y){
      sp <- as.data.frame(rowSums(y))
    })
    sites <- do.call(cbind, sps)
    names(sites) <- names(x)
    return(sites)
  })
}

sum.spNet <- function(data){
  #summarizes / reorganizes species level network metrics
  do.call(rbind,lapply(spNets,function(x){
    #node stats for each sp*site*yr combo
    site <- specieslevel(x)
    hl <- site$'higher level'
    hl$SpSiteYr <- rownames(hl)
    rownames(hl) <- NULL
    return(hl)
  }))
}


build.netlvl<- function(data,metrics, weighted=TRUE){
  #Builds a df of network level metrics from raw observation data 
  #data = list of sites, each of which has observations intended to become a network
  #metrics = list of metrics you  want to calculate
  #weighted = whether you want appropriate metrics to be weighted
  
  ls <- lapply(data, function(x){
    nets <- mclapply(x,function(y){
      mets <- as.data.frame(networklevel(y, metrics, weighted=weighted))
    },mc.cores=cores)
    nets1 <- do.call(cbind,nets)
    names(nets1) <- names(x)
    net.ls <- t(nets1)
  })
  
  df <- do.call(rbind,lapply (c(1:length(ls)),function(x){
    iter <- as.data.frame(ls[[x]])
    iter$SiteYr <- rownames(iter)
    rownames(iter) <- NULL
    iter$iteration <- x
    return(iter)
  }))
  
  return(df)
} 

netlvl.Zer <- function(data){
  #generates a summary z score df from network level metrics
  do.call(rbind,lapply(unique(data$SiteYr),function(x){
    site <- subset(data, data$SiteYr==x)
    metrics <- lapply(names(site[c(1:(length(site)-2))]),function(y){
      zscore <- (site[1,y] - mean(site[,y]))/
        (sd(site[,y])+10^-10)
      zscore1 <- data.frame(y = zscore)
      names(zscore1) <- c(y)
      return(zscore1)
    })
    
    site.df <- do.call(cbind,metrics)
    site.df$Site <- x
    return(site.df)
  }))
}


build.intraNetlvl <- function(data, iterations, threshold, metrics, binary, weighted=TRUE){
  #returns single dataframe of summary stats for each spsiteyr with more indivs than the threshold  
  
  iter <- c(1:iterations)
  
  observed <- lapply(data,function(x){
    sp <- lapply(x,function(y){
      if(length(colnames(y))>=threshold){
        vazSub(y, obs=TRUE,threshold=threshold,binary=binary)
      }
    })
    sp <- sp %>% discard(is.null)
  })
  
  #create list w/ randomized iterations
  net.ls <- lapply(iter,function(x){
    thresh <- lapply(data,function(y){
      sp <- lapply(y,function(z){
        if(length(colnames(z))>=threshold){
          vazSub(z, obs=FALSE,threshold=threshold,binary=binary)
        }
        
      })
      sp <- sp %>% discard(is.null)
    })
  })
  
  nets.ls <- list(observed)
  nets.ls <- c(nets.ls,net.ls)
  
  #get summary stats for nets.ls
  
  sums <- lapply(nets.ls,function(x){
    site <- lapply(x,function(y){
      sps <- lapply(y,function(z){
        sp <- as.data.frame(networklevel(z,metrics, weighted=weighted))
      })
      sps1 <- do.call(cbind,sps) 
      names(sps1) <- names(y)
      sps.df <- as.data.frame(t(sps1))
      sps.df$SpSiteYr <- rownames(sps.df)
      rownames(sps.df) <- NULL
      return(sps.df)
    })
    
    sites.df <- do.call(rbind,site)
    
  })
  
  len <- c(1:length(sums))
  
  sims <- lapply(len,function(x){
    sums[[x]]$iteration <- x
    return(sums[[x]])
  })
  
  all <- do.call(rbind,sims)
  
  return(all)
}



calc.Z <- function(data, split){
  #calculates z scores
  do.call(rbind,lapply(unique(data[,split]),function(x){
    site <- data[(data[,split]==x),]
    metrics <- lapply(names(site[c(1:(length(site)-2))]),function(y){
      zscore <- (site[1,y] - mean(site[,y]))/
        (sd(site[,y])+10^-10)
      zscore1 <- data.frame(y = zscore)
      names(zscore1) <- c(y)
      return(zscore1)
    })
    
    site.df <- do.call(cbind,metrics)
    site.df$Site <- x
    return(site.df)
  }))
}


nullweb = function(ref)
  #taken from ESM
{
  ref <- empty(ref)
  ref[ref>0] <- 1
  margin <- ifelse(ncol(ref)<nrow(ref),2,1)
  currentweb <- matrix(0,ncol=ncol(ref),nrow=nrow(ref))
  pc <- colMeans(ref)
  pr <- rowMeans(ref)
  if(margin==2)
  {
    for(i in 1:ncol(ref))
    {
      currentweb[,i] <- (pc[i]+pr)/2
    }
  } else {
    for(i in 1:nrow(ref))
    {
      currentweb[i,] <- (pr[i]+pc)/2
    }
  }
  return(apply(currentweb,margin,function(x)rbinom(length(x),1,x)))
}



colLister <- function(data){
  #categorizes bee species for this dataset into colors
  cols <- lapply(data,function(x){
    if(str_detect(x,"Melissodes agilis")==TRUE){
      col <- "blue"
    } else if(str_detect(x,"Svastra obliqua")==TRUE){
      col <- 'green'
    } else if(str_detect(x,"Halictus ligatus")==TRUE){
      col <- 'red'
    } else if(str_detect(x,"Melissodes lupina")==TRUE){
      col <- 'gray'
    } else if(str_detect(x,"Lasioglossum incompletum")==TRUE){
      col <- 'orange'
    } else if(str_detect(x,"Melissodes robustior")==TRUE){
      col <- 'yellow'
    } else if(str_detect(x,"Melissodes stearnsi")==TRUE){
      col <- 'pink'
    } else if(str_detect(x,"Apis mellifera")==TRUE){
      col <- 'purple'
    } else if(str_detect(x,"Halictus tripartitus")==TRUE){
      col <- 'black'
    } else {
      col <- 'light gray'
    }
    return(col)
  })
  
  return(unlist(cols))
  
}




basSub <- function(data, threshold, obs = FALSE, binary){
  #data should be a single species network
  #threshold is the minimum number of individuals to randomize
  #obs is whether or not nulls should be created, or whether the data is just being reorganized
  if(obs == FALSE & length(colnames(data))>=threshold){
    if(binary == FALSE){
      sp <- data[(rowSums(data)>0),]
      bas <- t(nullweb(sp))
      rownames(bas) <- rownames(sp)
      colnames(bas) <- colnames(sp)
      print('1')
      return(bas)
      
    } else {
      sp <- data[(rowSums(data)>0),]
      bin <- ifelse(sp>0,1,0)
      bas <- t(nullweb(bin))
      rownames(bas) <- rownames(sp)
      colnames(bas) <- colnames(sp)
      print('2')
      return(bas)
    }
    
  } else if(all(dim(data)>1)){
    if(binary==FALSE){
      sp <- data[(rowSums(data)>0),]
      print('3)')
      return(sp)
    } else{
      sp <- data[(rowSums(data)>0),]
      bin <- ifelse(sp>0,1,0)
      print('4')
      return(bin)
    }
    
  } else {
    if(binary==FALSE){
      sp <- data
      print('5')
      return(sp)
    } else{
      sp <- data
      print('6')
      bin <- ifelse(data>0,1,0)
      return(bin)
    }
  }
}


vaznull.fast <- function(web) {
  #a faster calculation of vazquez nulls for networks
  rs.p <- rowSums(web)/sum(web)
  cs.p <- colSums(web)/sum(web)
  P <- P1 <- tcrossprod(rs.p, cs.p)
  finalmat <- matrix(0, nrow(web), ncol(web))
  n.int.finalmat <- 0
  while (n.int.finalmat < sum(dim(web))) {
    sel <- sample(1:length(web), 1, prob = P, replace = TRUE)
    selc <- floor((sel - 1)/(dim(web)[1])) + 1
    selr <- ((sel - 1)%%dim(web)[1]) + 1
    if (sum(finalmat[, selc]) == 0 | sum(finalmat[selr,
    ]) == 0) {
      finalmat[sel] <- 1
      P[sel] <- 0
    }
    n.int.finalmat <- sum(rowSums(finalmat) > 0) + sum(colSums(finalmat) >
                                                         0)
  }
  conn.remain <- sum(web > 0) - sum(finalmat > 0)
  if (conn.remain > 0) {
    if (length(which(finalmat == 0)) == 1) {
      add <- which(finalmat == 0)
    }
    else {
      add <- sample(which(finalmat == 0), conn.remain,
                    prob = P1[finalmat == 0])
    }
    finalmat[add] <- 1
  }
  int.remain <- sum(web) - sum(finalmat)
  if (int.remain > 0) {
    add <- sample(which(finalmat > 0),
                  int.remain, prob = P1[which(finalmat >
                                                0)], replace = TRUE)
    finalmat[as.numeric(names(table(add)))] <-
      finalmat[as.numeric(names(table(add)))] +
      table(add)
  }
  finalmat
}

netToDisper <- function(data, meta, label){
  #calculates dispersion (in this case centroid distance) from network objects
  #meta refers to a larger dataset containing metadata for the specimens, not currently used
  #label refers to a tag for which the dispersions are calculated (e.g., rbcl reads)
  
  sites <- lapply(names(data), function(x){
    #print(x)
    #browser()
    site.net <- data[[x]] 
    ids <- colnames(site.net)
    
    sp.ls <- meta$GenusSpecies[meta$UniqueID %in% ids]
    
    if(length(sp.ls)>0){
      dis <- betadisper(vegdist(t(site.net),na.rm=T),sp.ls,type='centroid')
      
      avgs <- lapply(unique(dis$group), function(y){
        avgDist <- mean(dis$distances[dis$group==y])
        len <- length(dis$distances[dis$group==y])
        #browser()
        site.data <- data.frame('GenusSpecies'=y,
                                'avgDist' = avgDist,
                                'n' = len)
        return(site.data)
      })
      
      site.dat <- do.call(rbind,avgs)
      
      site.dat$site <- x
      site.dat$label <- label
      print(x)
      return(site.dat)
      
    } else {
      site.dat <- data.frame('GenusSpecies'=NA,
                             'avgDist' = NA,
                             'n' = NA,
                             "site" = x,
                             "label" = label)
      return(site.dat)
    }
    
    print(x)
    #if(length(dis)>0){
    #  print("yup")
    #}
    
  })
  
  return(do.call(rbind,sites))
  
}

netToDisperAll <- function(data, label){
  #calculates dispersion (in this case centroid distance) from network objects
  #label refers to a tag for which the dispersions are calculated (e.g., rbcl reads)
  
  sites <- lapply(names(data), function(x){
    site.net <- data[[x]] 
    col.length <- length(colnames(site.net))
    
    dis <- betadisper(vegdist(t(site.net),na.rm=T),
                      group=as.factor(rep('sp',times=col.length)),
                      type='centroid')
    
    site.data <- data.frame('site' = x,
                            'avgDist' = mean(dis$distances),
                            'n' = length(dis$distances),
                            'label' = label)
    return(site.data)
  })
  
  return(do.call(rbind,sites))
  
}

netTohub <- function(data,sitelist){
  #calculates hub score metrics for each network within a "sitelist" object
  #metrics include degree, closeness, and a number of network motifs describing
  #highly asymmetrical node arrangements
  #NOTE: only degree is used in analyses; the motifs (bmotif::mcount) seems to 
  #depricated
  sites <- lapply(sitelist,function(x){
    mets <- specieslevel(data[[x]],index=c('degree','closeness'))
    plants <- mets$'lower level'
    site <- data.frame('site' = x,'clustering' = max(clustering_local_tm(data[[x]])$lc,na.rm=T))
    
    
    df.met <- do.call(cbind,lapply(c('degree','weighted.closeness'),function(y){
      met <- data.frame(y = zscorer(plants,y))
      names(met) <- y
      return(met)
    }))
    
    site <- cbind(site,df.met)
    
    #motifs <- mcount(data[[x]],
    #                 six_node = T,
    #                 normalisation = T,
    #                 mean_weight = T,
    #                 standard_dev = F)
    
    #site$M44 <- (motifs$mean_weight[44] - mean(motifs$mean_weight[18:44]))/sd(motifs$mean_weight[18:44])
    #site$M38 <- (motifs$mean_weight[38] - mean(motifs$mean_weight[18:44]))/sd(motifs$mean_weight[18:44])
    #site$M17 <- (motifs$mean_weight[17] - mean(motifs$mean_weight[8:17]))/sd(motifs$mean_weight[18:44])
    #site$M13 <- (motifs$mean_weight[13] - mean(motifs$mean_weight[8:17]))/sd(motifs$mean_weight[18:44])
    
    #site$M44f <- (motifs$frequency[44] - mean(motifs$frequency[18:44]))/sd(motifs$frequency[18:44])
    #site$M38f <- (motifs$frequency[38] - mean(motifs$frequency[18:44]))/sd(motifs$frequency[18:44])
    #site$M17f <- (motifs$frequency[17] - mean(motifs$frequency[8:17]))/sd(motifs$frequency[18:44])
    #site$M13f <- (motifs$frequency[13] - mean(motifs$frequency[8:17]))/sd(motifs$frequency[18:44])
    
    return(site)
  })
  do.call(rbind,sites)
}

netToDeg <- function(data,sitelist){
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

## GGPLOT THEMES

#Call a specific set of ggplot theme elements
theme_dark_black <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = 'black'),
            plot.background = element_rect(colour = NA, fill = 'black'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'white'),
            axis.line.x = element_line(colour="white"),
            axis.line.y = element_line(colour="white"),
            axis.ticks = element_line(colour="white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            legend.background = element_rect(fill ='black'),
            legend.text = element_text(color = 'white'),
            legend.key = element_rect(colour = NA, fill = 'black'),
            ## legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'white'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="black"),
            strip.text = element_text(face="bold", colour =
                                        'white')
    ))
}

#Call a specific set of ggplot theme elements
theme_ms <- function(base_size=14, base_family="sans") {
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

#split site-level networks into site+sample round networks
SRsplit <- function(network.list,meta=spec){
  newlist <- do.call(c,lapply(names(network.list),function(x){
    #browser()
    #narrow down to just the right site-level network
    site <- network.list[[x]] 
    
    rounds <- lapply(unique(meta$SiteSample[meta$Site==x]), function(y){
      
      #subset that site network to just specimens collected during a given sample round
      siteround <- site[,colnames(site) %in% meta$UniqueID[meta$SiteSample==y]]
      
      #drop empty plant levels (but don't break for 1 weird case). EDIT: commented out, not necessary. doesn't change calculations
      #if(length(colnames(siteround))>0){
      #  siteround <- siteround[rowSums(siteround)>0,]
      #}
      
      #browser()
      print(y)
      #print(dim(siteround))
      return(siteround)
    })
    #name the objects in the rounds list
    #browser()
    names(rounds) <- unique(meta$SiteSample[meta$Site==x])
    
    #remove any networks with dim 0 (i.e., no pollinators)
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
BloomSplit <- function(network.list,meta=spec){
  newlist <- do.call(c,lapply(names(network.list),function(x){
    #browser()
    #narrow down to just the right site-level network
    site <- network.list[[x]] 
    
    rounds <- lapply(unique(meta$SiteBloom[meta$Site==x]), function(y){
      
      #subset that site network to just specimens collected during a given sample round
      bloomround <- site[,colnames(site) %in% meta$UniqueID[meta$SiteBloom==y]]
      
      #drop empty plant levels (but don't break for 1 weird case). EDIT: commented out, not necessary. doesn't change calculations
      #if(length(colnames(siteround))>0){
      #  siteround <- siteround[rowSums(siteround)>0,]
      #}
      
      #browser()
      print(y)
      #print(dim(siteround))
      return(bloomround)
    })
    #name the objects in the rounds list
    #browser()
    names(rounds) <- unique(meta$SiteBloom[meta$Site==x])
    
    #remove any networks with dim 0 (i.e., no pollinators)
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

dis_combiner <- function(dfmicro,dfplant){
  #reorganize the dataframes
  micro.df <- data.frame('obs' = paste(dfmicro$GenusSpecies, dfmicro$site, sep=":"),
                         'microDist' = dfmicro$avgDist,
                         'microN' = dfmicro$n)
  rbcl.df <- data.frame('obs' = paste(dfplant$GenusSpecies, dfplant$site, sep=":"),
                        'rbclDist' = dfplant$avgDist,
                        'rbclN' = dfplant$n)
  
  
  #combine them into a single df, add a couple organizer columns
  comp <- full_join(micro.df,rbcl.df)
  #browser()
  comp$site <- word(comp$obs,2,sep=":")
  comp$GenusSpecies <- word(comp$obs,1,sep=":")
  return(comp)
}


MRMrunner <- function(sitelist, spec.data, rbcl, micro, bees){
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

MRMrunnerSB <- function(sitelist, spec.data, rbcl, micro, bees){
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
  #Version using the new filtering and including bee phylo in the model (when there are multiple bee species present)
  
  siteSB.list <- do.call(rbind,lapply(sitelist,function(x){
    print(x)
    browser()
    #site <- spec.data[spec.data$SiteBloom==x,]
    
    ID.list <- spec.data$UniqueID[spec.data$SiteBloom==x]
    
    #thing <- spec.data[spec.data$UniqueID%in%ID.list,c("Site","SFBloomStatus","Pollen","Gut","GenusSpecies")]
    
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
       # browser()
        bees.filt <- as.matrix(do.call(rbind, lapply(rownames(rbcl.filt.sm), function(x){
          comp.list <- unlist(lapply(colnames(rbcl.filt.sm), function(y){
            #browser()
            x.sp <- spec.data$GenusSpecies[spec.data$UniqueID==x]
            y.sp <- spec.data$GenusSpecies[spec.data$UniqueID==y]
            return(bees[x.sp,y.sp])
          }))
        })))
        rownames(bees.filt) <- rownames(rbcl.filt.sm)
        colnames(bees.filt) <- rownames(rbcl.filt.sm)
        
        #print("filtered")
        #browser()
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

samp2site.spp <- function(site,spp,abund) {
    x <- tapply(abund, list(site=site,spp=spp), sum)
    x[is.na(x)] <- 0
    return(x)
}

makeCommStruct <- function(spec.dat, type){
    ## prep site by species matrix
    prep.comm <- aggregate(spec.dat[, type],
                           list(site= spec.dat$Site,
                                status= spec.dat$SiteType,
                                sp= spec.dat[, type]), length)

    comm <-  samp2site.spp(site= prep.comm$site,
                           spp= prep.comm$sp, abund=
                                                  prep.comm$x)
    sites <- rownames(comm)
    site.type <- spec.dat$SiteType[match(rownames(comm),
                                         spec.dat$Site)]
    adjsf <- spec.dat$AdjSF[match(rownames(comm),
                                  spec.dat$Site)]

    comm <- bipartite::empty(comm)

    return(list(comm=comm,
                ## sites=sites,
                site.type = site.type,
                adjsf = adjsf))
}

makeCommStructSR <- function(spec.dat, type){
  ## prep site by species matrix
  prep.comm <- aggregate(spec.dat[, type],
                         list(site= spec.dat$SiteSample,
                              status= spec.dat$SiteType,
                              sp= spec.dat[, type]), length)
  
  comm <-  samp2site.spp(site= prep.comm$site,
                         spp= prep.comm$sp, abund=
                           prep.comm$x)
  sites <- rownames(comm)
  site.type <- spec.dat$SiteType[match(rownames(comm),
                                       spec.dat$SiteSample)]
  adjsf <- spec.dat$AdjSF[match(rownames(comm),
                                spec.dat$SiteSample)]
  comm <- bipartite::empty(comm)
  
  return(list(comm=comm,
              ## sites=sites,
              site.type = site.type,
              adjsf = adjsf))
}
makeCommStructSB <- function(spec.dat, type){
  ## prep site by species matrix
  prep.comm <- aggregate(spec.dat[, type],
                         list(site= spec.dat$SiteBloom,
                              status= spec.dat$SiteType,
                              sp= spec.dat[, type]), length)
  #print(1)
  comm <-  samp2site.spp(site= prep.comm$site,
                         spp= prep.comm$sp, abund=
                           prep.comm$x)
  sites <- rownames(comm)
  site.type <- spec.dat$SiteType[match(rownames(comm),
                                       spec.dat$SiteBloom)]
  adjsf <- spec.dat$AdjSF[match(rownames(comm),
                                spec.dat$SiteBloom)]
  comm <- bipartite::empty(comm)
  
  return(list(comm=comm,
              ## sites=sites,
              site.type = site.type,
              adjsf = adjsf))
}


makeStructFromComm <- function(spec, spp){
    spp.pre.comm <- spec[, c("Site", spp)]
    spp.pre.comm <- spp.pre.comm[!apply(spec[, spp], 1,
                                                  function(x) all(is.na(x))),]

    spp.pre.comm <- spp.pre.comm  %>%
        group_by(Site) %>%
        summarise_each(list(mean))

    site.type <- spec$SiteType[match(spp.pre.comm$Site,
                                     spec$Site)]

    adjsf <- spec$AdjSF[match(spp.pre.comm$Site,
                             spec$Site)]

    sites <- spp.pre.comm$Site

    comm <- spp.pre.comm
    comm$Site <- NULL
    comm <- as.matrix(comm)
    rownames(comm) <- spp.pre.comm$Site
    comm[is.na(comm)] <- 0
    comm <- bipartite::empty(comm)

    list(comm=comm,
         ## sites=sites,
         site.type=site.type,
         adjsf = adjsf)
}


getParComm <- function(parasite){
    parasite <- aggregate(list(Parasite=spec[, parasite]),
                          list(GenusSpecies=spec$GenusSpecies,
                               Site=spec$AltSiteName,
                               SiteType=spec$SiteType),
                          function(x) sum(x) / length(x))

    parasite.comm <- samp2site.spp(parasite$Site,
                                   parasite$GenusSpecies,
                                   parasite$Parasite)
    return(parasite.comm)
}

mod.unifrac <- function (comm, tree)
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute UniFrac")
  }
  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for UniFrac calculation")
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

pd.mod <- function (samp, tree, include.root = TRUE)
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (include.root) {
    if (!is.rooted(tree)) {
      stop("Rooted tree required to calculate PD with include.root=TRUE argument")
    }
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
      if (!is.rooted(tree) || !include.root) {
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
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
        }
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[,
                                                         2] %in% which(tree$tip.label %in% present))])
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
