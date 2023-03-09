


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


zscorer <- function(data, score){
  max <- max(data[,score],na.rm=T)
  mean <- mean(data[,score],na.rm=T)
  sd <- sd (data[,score],na.rm=T)
  return((max-mean)/sd)
}

netTohub <- function(data,sitelist){
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

    motifs <- mcount(data[[x]],
                     six_node = T,
                     normalisation = T,
                     mean_weight = T,
                     standard_dev = F)
    
    site$M44 <- (motifs$mean_weight[44] - mean(motifs$mean_weight[18:44]))/sd(motifs$mean_weight[18:44])
    site$M38 <- (motifs$mean_weight[38] - mean(motifs$mean_weight[18:44]))/sd(motifs$mean_weight[18:44])
    site$M17 <- (motifs$mean_weight[17] - mean(motifs$mean_weight[8:17]))/sd(motifs$mean_weight[18:44])
    site$M13 <- (motifs$mean_weight[13] - mean(motifs$mean_weight[8:17]))/sd(motifs$mean_weight[18:44])
    
    site$M44f <- (motifs$frequency[44] - mean(motifs$frequency[18:44]))/sd(motifs$frequency[18:44])
    site$M38f <- (motifs$frequency[38] - mean(motifs$frequency[18:44]))/sd(motifs$frequency[18:44])
    site$M17f <- (motifs$frequency[17] - mean(motifs$frequency[8:17]))/sd(motifs$frequency[18:44])
    site$M13f <- (motifs$frequency[13] - mean(motifs$frequency[8:17]))/sd(motifs$frequency[18:44])
    
    return(site)
  })
  do.call(rbind,sites)
}

netToDisper <- function(data, meta, label){
  
  sites <- lapply(names(data), function(x){
    site.net <- data[[x]] 
    ids <- colnames(site.net)
    
    sp.ls <- meta$GenusSpecies[meta$UniqueID %in% ids]
    
    dis <- betadisper(vegdist(t(site.net),na.rm=T),sp.ls,type='centroid')

    avgs <- lapply(unique(dis$group), function(y){
      avgDist <- mean(dis$distances[dis$group==y])
      n <- length(dis$distances[dis$group==y])
      site.data <- data.frame('GenusSpecies'=y,
                              'avgDist' = avgDist,
                              'n' = n)
      return(site.data)
    })
    
    site.dat <- do.call(rbind,avgs)

    site.dat$site <- x
    site.dat$label <- label
    print(x)
    return(site.dat)
  })
  
  return(do.call(rbind,sites))
  
}

netToDisperAll <- function(data, label){
  
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




plotCommDistbyGroup  <- function(dist.mat,
                                 groups,
                                 species.type,
                                 f.path="figures/pcoas",
                                 group.name){
  
  f.pcoa <- function(){
    #dist.mat.na <- apply(dist.mat, 1, function(x) any(is.na(x)))
    #dist.mat <- dist.mat[!dist.mat.na, !dist.mat.na]
    pcoa.comm <- cmdscale(dist.mat)
    dist.m <- as.dist(dist.mat)
    print("betadisper result")
    print(anova(betadisper(dist.m,groups)))
    
    pcoa.mod <- adonis(dist.mat~groups)
    plot(NA, asp=1,  cex=1.5,
         ylim=range(pcoa.comm[,2]),
         xlim=range(pcoa.comm[,1]),
         xlab='',
         ylab='',
         xaxt='n',
         yaxt='n',
         cex.lab=1.5)
    
    cols <- brewer.pal(length(unique(groups)), "Set3")
    names(cols) <- unique(groups)
    for(group in unique(groups)){
      ## all points sitting on top of eachother so triple jitter
      points(pcoa.comm[groups == group,],
             col=cols[group], pch=16, cex=1.5)
      points(pcoa.comm[groups == group,],
             col="black", pch=1, cex=1.5)
    }
    ordihull(pcoa.comm, groups)
    
    legend("topright", legend= unique(groups),
           bty="n", col=cols[unique(groups)],
           pch=16, cex=1)
    legend("topright", legend= unique(groups),
           bty="n", col="black", pch=1,
           cex=1)
    
    mtext('PCoA1', 1, line=2, cex=1.5)
    mtext('PCoA2', 2, line=2, cex=1.5)
    return(pcoa.mod)
  }
  ## function for plotting PcoA axes
  pdf.f(f.pcoa,
        file= file.path(f.path,
                        sprintf("%s_%s_pcoa.pdf",
                                species.type, group.name)),
        width=7, height=7)
}


mcount1 <- function (M, six_node = FALSE, normalisation, mean_weight, standard_dev) 
{
  if (inherits(M, "matrix") != TRUE) {
    stop("'M' must be an object of class 'matrix'")
  }
  print(1)
  if (!all(apply(M, 1:2, is.numeric))) {
    stop("Elements of 'M' must be numeric")
  }
  print(2)
  if (!all(apply(M, 1:2, function(x) length(x) > 0))) {
    stop("Elements of 'M' cannot have length zero")
  }
  print(3)
  if (!all(apply(M, 1:2, function(x) x >= 0))) {
    stop("Elements of 'M' must be greater than or equal to zero")
  }
  print(4)
  if (inherits(normalisation, "logical") != TRUE) {
    stop("'normalisation' must be of class 'logical' i.e. TRUE or FALSE")
  }
  print(5)
  if (inherits(six_node, "logical") != TRUE) {
    stop("'six_node' must be of class 'logical' i.e. TRUE or FALSE")
  }
  print(6)
  if (mean_weight == FALSE & standard_dev == TRUE) {
    stop("Cannot have standard_dev = TRUE and mean_weight = FALSE. If you want the standard deviations, set standard_dev = TRUE and mean_weight = TRUE")
  }
  print(7)
  if (six_node == TRUE & standard_dev == TRUE) {
    warning("Standard deviation values are not available for six node motifs. Standard deviation will only be returned for 2-5 node motifs")
  }
  print(8)
  W <- M
  M[M > 0] <- 1
  dimnames(M) <- NULL
  p <- dim(M)[2]
  z <- dim(M)[1]
  J <- matrix(rep(1, z * p), nrow = z, ncol = p)
  JP <- matrix(rep(1, p * p), nrow = p, ncol = p)
  JZ <- matrix(rep(1, z * z), nrow = z, ncol = z)
  MT <- t(M)
  N <- J - M
  NT <- t(N)
  P <- MT %*% M
  Q <- MT %*% N
  R <- NT %*% M
  Z <- M %*% MT
  Y <- M %*% NT
  X <- N %*% MT
  dP <- apply(M, MARGIN = 2, sum)
  jP <- rep(1, p)
  dZ <- apply(M, MARGIN = 1, sum)
  jZ <- rep(1, z)
  if (six_node == TRUE) {
    if (p < z) {
      J3 <- array(rep(1, p * p * p), c(p, p, p))
      AP <- maketensor(M, M)
      BP <- maketensor(M, N)
      CP <- maketensor(N, M)
      DP <- maketensor(N, N)
      MA <- tensor::tensor(MT, AP, 2, 1)
      MB <- tensor::tensor(MT, BP, 2, 1)
      MC <- tensor::tensor(MT, CP, 2, 1)
      MD <- tensor::tensor(MT, DP, 2, 1)
      Na <- tensor::tensor(NT, AP, 2, 1)
      NB <- tensor::tensor(NT, BP, 2, 1)
      NC <- tensor::tensor(NT, CP, 2, 1)
      K3 <- J3
      for (i in 1:p) {
        for (j in 1:p) {
          K3[i, j, j] <- 0
          K3[j, i, j] <- 0
          K3[j, j, i] <- 0
        }
      }
    }
    if (p >= z) {
      J3 <- array(rep(1, z * z * z), c(z, z, z))
      AP <- maketensor(MT, MT)
      BP <- maketensor(MT, NT)
      CP <- maketensor(NT, MT)
      DP <- maketensor(NT, NT)
      MA <- tensor::tensor(M, AP, 2, 1)
      MB <- tensor::tensor(M, BP, 2, 1)
      MC <- tensor::tensor(M, CP, 2, 1)
      MD <- tensor::tensor(M, DP, 2, 1)
      Na <- tensor::tensor(N, AP, 2, 1)
      NB <- tensor::tensor(N, BP, 2, 1)
      NC <- tensor::tensor(N, CP, 2, 1)
      K3 <- J3
      for (i in 1:z) {
        for (j in 1:z) {
          K3[i, j, j] <- 0
          K3[j, i, j] <- 0
          K3[j, j, i] <- 0
        }
      }
    }
    MA <- MA * K3
    MB <- MB * K3
    MC <- MC * K3
    MD <- MD * K3
    Na <- Na * K3
    NB <- NB * K3
    NC <- NC * K3
  }
  if (six_node == FALSE) {
    if (normalisation == TRUE) {
      out <- data.frame(motif = 1:17, nodes = c(2, rep(3, 
                                                       2), rep(4, 4), rep(5, 10)), frequency = NA, normalise_sum = NA, 
                        normalise_sizeclass = NA, normalise_levelsize = NA, 
                        normalise_nodesets = NA)
    }
    else {
      out <- data.frame(motif = 1:17, nodes = c(2, rep(3, 
                                                       2), rep(4, 4), rep(5, 10)), frequency = NA)
    }
  }
  else {
    if (normalisation == TRUE) {
      out <- data.frame(motif = 1:44, nodes = c(2, rep(3, 
                                                       2), rep(4, 4), rep(5, 10), rep(6, 27)), frequency = NA, 
                        normalise_sum = NA, normalise_sizeclass = NA, 
                        normalise_levelsize = NA, normalise_nodesets = NA)
    }
    else {
      out <- data.frame(motif = 1:44, nodes = c(2, rep(3, 
                                                       2), rep(4, 4), rep(5, 10), rep(6, 27)), frequency = NA)
    }
  }
  if (six_node == FALSE) {
    for (i in 1:17) {
      out[i, "frequency"] <- countmotif(x = M, motif = i, 
                                        z = z, p = p, JP = JP, JZ = JZ, P = P, Q = Q, 
                                        R = R, Z = Z, Y = Y, X = X, dP = dP, jP = jP, 
                                        dZ = dZ, jZ = jZ)
    }
  }
  else {
    for (i in 1:44) {
      out[i, "frequency"] <- countmotif(x = M, motif = i, 
                                        z = z, p = p, JP = JP, JZ = JZ, P = P, Q = Q, 
                                        R = R, Z = Z, Y = Y, X = X, dP = dP, jP = jP, 
                                        dZ = dZ, jZ = jZ, J3 = J3, MA = MA, MB = MB, 
                                        MC = MC, MD = MD, Na = Na, NB = NB, NC = NC)
    }
  }
  if (normalisation == TRUE) {
    out$normalise_sum <- out$frequency/sum(out$frequency)
    out$normalise_sizeclass <- do.call("c", lapply(split(out, 
                                                         out$nodes), function(df) df$frequency/sum(df$frequency)))
    out$normalise_levelsize <- NA
    out$normalise_levelsize[1] <- out$frequency[1]/sum(out$frequency[1])
    out$normalise_levelsize[2] <- out$frequency[2]/sum(out$frequency[2])
    out$normalise_levelsize[3] <- out$frequency[3]/sum(out$frequency[3])
    out$normalise_levelsize[4] <- out$frequency[4]/sum(out$frequency[4])
    out$normalise_levelsize[5:6] <- out$frequency[5:6]/sum(out$frequency[5:6])
    out$normalise_levelsize[7] <- out$frequency[7]/sum(out$frequency[7])
    out$normalise_levelsize[8] <- out$frequency[8]/sum(out$frequency[8])
    out$normalise_levelsize[9:12] <- out$frequency[9:12]/sum(out$frequency[9:12])
    out$normalise_levelsize[13:16] <- out$frequency[13:16]/sum(out$frequency[13:16])
    out$normalise_levelsize[17] <- out$frequency[17]/sum(out$frequency[17])
    if (six_node == TRUE) {
      out$normalise_levelsize[18] <- out$frequency[18]/sum(out$frequency[18])
      out$normalise_levelsize[19:24] <- out$frequency[19:24]/sum(out$frequency[19:24])
      out$normalise_levelsize[25:37] <- out$frequency[25:37]/sum(out$frequency[25:37])
      out$normalise_levelsize[38:43] <- out$frequency[38:43]/sum(out$frequency[38:43])
      out$normalise_levelsize[44] <- out$frequency[44]/sum(out$frequency[44])
    }
    sets <- node_sets(M, six_node = six_node)
    out$normalise_nodesets <- out$frequency/sets
    out[do.call(cbind, lapply(out, is.nan))] <- NA
  }
  if (mean_weight) {
    out$mean_weight <- mean_weight(W, mc = out$frequency, 
                                   six_node = six_node)
  }
  if (standard_dev) {
    out[1:17, "standard_dev"] <- motif_sd(W, mc = out$frequency, 
                                          meanw = out$mean_weight)
    if (six_node) {
      out[18:44, "standard_dev"] <- NA
    }
  }
  return(out)
}
