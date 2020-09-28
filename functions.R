#Functions


fix.white.space <- function(d) {
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
  spec.dat$Date <- as.Date(spec.dat$Date, format='%m/%d/%y')
  spec.dat$Doy <- as.numeric(strftime(spec.dat$Date, format='%j'))
  spec.dat$Year <- as.numeric(format(spec.dat$Date,'%Y'))
  return(spec.dat)
}

dat.rm.blanks <- function(spec.dat) {
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

#not currently used
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
  if(obs == FALSE & length(colnames(data))>threshold){
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

vazSub1 <- function(data, threshold, obs, binary) { 

  #down to species
  net <- lapply(data,function(y){
    
    sp <- vazSub(y, threshold=threshold,obs=obs, binary=binary)
    
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


build.spNet <- function(data){
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
  do.call(rbind,lapply(spNets,function(x){
  #node stats for each sp*site*yr combo
  site <- specieslevel(x)
  hl <- site$'higher level'
  hl$SpSiteYr <- rownames(hl)
  rownames(hl) <- NULL
  return(hl)
}))
}

build.netlvl<- function(data,metrics){
  ls <- lapply(data, function(x){
    nets <- mclapply(x,function(y){
      mets <- as.data.frame(networklevel(y, metrics))
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


build.intraNetlvl <- function(data, iterations, threshold, metrics, binary){
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
        sp <- as.data.frame(networklevel(z,metrics))
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
  
  netlvl.Zer <- function(data){
    
    
  
  len <- c(1:length(data))
  
  sims <- lapply(len,function(x){
    data[[x]]$iteration <- x
  })
  
  all <- do.call(rbind,sims)
  
  
  
} 












vaznull.fast <- function(web) {
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