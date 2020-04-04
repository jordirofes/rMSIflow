### This funcion separates all the images

sepMatrix <- function(peakMatrix){
  matrixList <- list()
  if(length(peakMatrix$numPixels) == 1){
    matrixList[[1]] <- peakMatrix
  } else {
  for(i in 1:length(peakMatrix$numPixels)){
    limitDown <- sum(peakMatrix$numPixels[1:(i-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:i])
    if(i == 1){limitDown <- 1}
    int <- limitDown:limitUp
    matrixList[[i]] <- rMSIproc::`[.rMSIprocPeakMatrix`(peakMatrix, int)}
  }
  return(matrixList)
}

### This function plots the intensities of a given mass on the images

initialPlotter <- function(peakMatrix, matrixL, mz, img, norma, sv, file){
  a <- vector("list", length(img))
  # if(is.na(groups)){
  #   groups <-  peakMatrix$names
  #   }
  if (img == 0){
      a <- rMSIproc::plotPeakImageG(peakMatrix, mz, normalization = norma)
      if(sv == T){svPlot(file, a)}
      return(a)
  } else{
    for (j in (1:length(img))){
      
      if(!is.na(norma)){
      limitDown <- sum(peakMatrix$numPixels[1:(j-1)])+1
      limitUp <- sum(peakMatrix$numPixels[1:j])
      if(j == 1){limitDown <- 1}
      int <- limitDown:limitUp} else{int = 1}
      
      a[[j]] <- rMSIproc::plotPeakImageG(matrixL[[img[j]]], mz, normalization = norma[int])
    }
  }
  pl <- lapply(1:length(img2plot), function(i){
    plotly::ggplotly(p = a[[i]])
  })
  if(sv == T){svPlot(file, a)}
  return(subplot(pl))
}

### This function saves given plot on a pdf file  

svPlot <- function(flName, plots2save){
  pdf(flName)
  print(plots2save)
  dev.off()

}

### This function makes a PCA from the peak matrix

pcaJr <- function(peakMatrix, cnt = T, scl = T, norma = NA){
  dataPca <- peakMatrix$intensity
  if(!is.na(norma)){
    dataPca <- dataPca/norma
    }
     colnames(dataPca) <- peakMatrix$mass
  pcaJr <- prcomp(dataPca, center = cnt, scale. = scl) ### Fem la PCA de la matriu de pics
  return(pcaJr)
}

### This function creates a two dimensional plot with the especified principal components

pca2dimPlot <- function(peakMatrix, pcaJr, grpImg, pc1, pc2, sv2, fileN2) {
  if(!is.na(grpImg)){
  listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, grpImg, 1:length(grpImg))
  } else {listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, peakMatrix$names, 1:length(peakMatrix$names)) }
  
  vectorGroups <- unlist(listGroups)
  pcaData <- list()
  pcaData[[1]] <- pcaJr$x[,pc1]
  pcaData[[2]] <- pcaJr$x[,pc2]
  pcaData[[3]] <- vectorGroups
  dataPca <- as.data.frame(pcaData)
  colnames(dataPca) <- c("PC1", "PC2", "Groups")
  rownames(dataPca) <- c(1:length(dataPca$PC1))
  variance <- pcaJr$sdev^2
  xvariance <- round(variance/sum(variance)*100, 1)
  
  plotPca <- ggplot(dataPca,aes(x = PC1, y = PC2, colour = Groups))+ geom_point(alpha = 0.5)+ stat_ellipse() + 
    xlab(paste0("PC", pc1 ,"(", xvariance[pc1],"%)")) + ylab(paste0("PC",pc2, "(", xvariance[pc2],"%)")) + 
    theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_point(aes(x= 0, y =0), color= "black") 
  if(sv2 == T){
    svPlot(fileN2, plotPca)
  }
  return(plotPca)
}

### This function plots a chosen principal component into selected images

pcaPlotImg <- function(peakMatrix, matrixLst, pcaJr, img ,grpImg, pc, sv, fileN){
  a <- vector("list", length(img))
  if(is.na(grpImg)){ grpImg <- peakMatrix$names}
  
  if (img == 0){
    a <- rMSIproc::plotValuesImageG(peakMatrix, pixel_values = pcaJr$x[,pc], plot_labels =  grpImg)
      if(sv == T){svPlot(fileN, a)}
    return(a)
  } else{
    for (j in (1:length(img))){
      
      limitDown <- sum(peakMatrix$numPixels[1:(img[j]-1)])+1
      limitUp <- sum(peakMatrix$numPixels[1:img[j]])
      if(img[j] == 1){limitDown <- 1}
      int <- limitDown:limitUp
      a[[j]] <- rMSIproc::plotValuesImageG(matrixLst[[img[j]]], pixel_values = pcaJr$x[int,pc] , plot_labels = grpImg[img[j]])
    }
  }
  pl <- lapply(1:length(img), function(i){
    plotly::ggplotly(p = a[[i]])
  })
    if(sv == T){svPlot(fileN, a)}
  return(plotly::subplot(pl))
  
}

pcChooseNum <- function(pcaJr, perc = 0.9){
  
  variance <- pcaJr$sdev^2
  pcVar <- sapply(variance, function(dev){
    dev/sum(variance)*100
  })
  
  cumulative <- sapply(1:length(variance), function(x){sum(variance[1:x])})/sum(variance)*100
  
  pcplotData <- data.frame(Cumulative = cumulative, PC = 1:length(pcaJr$sdev))
  bound <- which.min(abs(cumulative-perc*100))
  
  
  plot <- ggplotly(ggplot(pcplotData) + geom_line(aes(x = PC, y = Cumulative)) + geom_point(aes(x = PC, y = Cumulative)) + theme_minimal() +
                     ylab("Cummulative %") + xlab("Principal Component") + geom_segment(aes(x = bound, y = 0, xend = bound, yend = cumulative[bound]), linetype = "dashed") +
                     geom_segment(aes(x = 0, y = cumulative[bound], xend = bound, yend = cumulative[bound]), linetype = "dashed"))
  
  # plot <- ggplotly(ggplot(pcplotData) + geom_line(aes(x = PC, y = Std)) + geom_point(aes(x = PC, y = Std)) + theme_minimal() +
  #                      ylab("Explained %") + xlab("Principal Component") + geom_segment(aes(x = bound, y = 0, xend = bound, yend = Std[bound] )))
  return(plot)
}

### These functions plot the Average Raw Spectrum of selected images in an interactive way

medSpecRaw <- function(directory){
  bad <- list.files(wDir, pattern = "^ramdisk")
  imageName <- list.files(wDir, pattern = "proc.tar")
  imageName <- imageName[!(imageName %in% bad)]

  # Then we run the code so it extracts the image information and calculates the average spectrum
  avSpec <- lapply(imageName, function(name){
    image <- rMSI::LoadMsiData(file.path(wDir, name))
    mass <- image$mass
    aS <- rMSIproc::AverageSpectrum(image) 
    specData <- cbind(mass, aS)
    return(specData)
  })
  return(avSpec)
}

medSpecRawplot <- function(specData, img2plot2){
  
specPlot <- ggplot()
parsed <- lapply(img2plot2, function(x){
  data.frame(mz = specData[[x]][,1], int = specData[[x]][,2], image = rep(paste("Image",x),
             times = nrow(specData[[x]])))
})
parsed <- do.call(rbind, parsed)
  plots2 <- ggplotly(ggplot(parsed) + geom_line(aes(x = mz, y = int, 
                          colour = image)) + theme_minimal() + xlab("Mass") + ylab("Intensity") + 
                          theme(legend.title = element_blank()))
   return(plots2)
}

### These functions plot the Average Spectrum of the peak matrix

medSpecP <- function(peakMatrix, pixels, norma = NA){
    datapm <- peakMatrix$intensity
    if(!is.na(norma)){
      datapm <- datapm/norma
    }
    pixelMatrix <- datapm[pixels,]
    avgSpecpm <- apply(pixelMatrix, 2, mean)
    return(avgSpecpm)
  }

medSpecComp <- function(peakMat, pixels1, pixels2, normalization, name1, name2, sav2, filename2){
  medSpec1 <- medSpecP(peakMat, pixels1, normalization)
  medSpec2 <- medSpecP(peakMat, pixels2, normalization)
  medSpecdf <- as.data.frame(cbind(medSpec1, medSpec2, peakM$mass))
  colnames(medSpecdf) <- c("Spec1", "Spec2", "Mass") 
  
  aveSpec <- ggplot(medSpecdf) + geom_segment(aes(x = Mass,y = 0, xend = Mass, yend = Spec1, colour = name1)) +
                                geom_segment(aes(x = Mass, y = 0, xend = Mass, yend = Spec2, colour = name2))  +
                                xlab("Mass") + ylab("Intensity") + theme_minimal() + theme(legend.title = element_blank())
  
  if(sav2 == T){ svPlot(filename2, aveSpec)}
  aveSpecInteractive <- ggplotly(aveSpec)
  return(aveSpecInteractive)
  
}
  
  
### K-means
elbowMethod <- function(peakMatrix, norma, testClus = 10){
  
  datapm <- peakMatrix$intensity
  if(!is.na(norma)){
    datapm <- datapm/norma
  }
  whitss <- sapply(1:testClus, function(x){
    km <- kmeans(datapm, x)
    km$tot.withinss
  })
  df <- data.frame("NC" = 1:testClus, "TW" =  whitss)
  plotKmeans <- ggplot(df) + geom_point(aes(x = NC, y =TW)) + geom_line(aes(x = NC, y = TW)) +
    ylab("Total Withinss") + xlab("Number of Clusters") + theme_minimal()
  
  return(plotKmeans)
}  

kmeansCluster <- function(peakMatrix, intensities, norma, numCl, together){
  clusData <- intensities
  
  if(!is.na(norma)){
    clusData <- clusData/norma
  }
  
  clData <- list()
  if(together == F){
  for (i in 1:length(peakMatrix$numPixels)){
    limitDown <- sum(peakMatrix$numPixels[1:(i-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:i])
    int <-limitDown:limitUp
    if (i == 1){ limitDown <- 1}

      clData[[i]] <- kmeans(clusData[limitDown:limitUp, ], numCl)

  }
  } else{
  clData[[1]] <- kmeans(clusData, numCl)
  }
  return(clData)
  
}

clusterDataPlotting <- function(peakMatrix, matrixL, img, groups, clusterData, sv, file){
  
  clusterData <- lapply(clusterData, function(x){
    return(x[1]) 
  })
  allClusImages <- unlist(clusterData)
  
  int2plot <- lapply(img, function(y){
    limitDown <- sum(peakMatrix$numPixels[1:(y-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:y])
    if(y == 1){limitDown <- 1}
    int <- limitDown:limitUp
  })
  
   if(img == 0){
    kplot <- rMSIproc::plotClusterImageG(peakMatrix, allClusImages)
    if(sv == T){svPlot(file, kplot)}
    return(kplot)
  } else{
     
     int2plot <- unlist(int2plot)
     matrix2plot <- rMSIproc::`[.rMSIprocPeakMatrix`(peakMatrix, int2plot)
     plot <- rMSIproc::plotClusterImageG(matrix2plot, allClusImages[int2plot])
      return(plot)
  }
  if(sv == T){svPlot(file, kplot)}
  # kplot <- lapply(kplots, function(x){
  #   plots <- ggplotly(kplots)
  #   return(plots)
  # })
  # return(subplot(kplot))
  }

clusterKmeansComparison <- function(peakMatrix, norma, img, groups, cluster, clusterData){

  if(is.na(groups)){               # Assigning groups to image names if vector class equals NA
    groups <- peakMatrix$names 
  }
  clusterData <- lapply(clusterData, function(x){   # Getting the clusterization from the k-means data
    return(x[1])  
  })
clusterData <- unlist(clusterData)

clusSpec <- lapply(img, function(x){          # We calculate the medium spectra from chosen images and chosen clusters
    limitDown <- sum(peakMatrix$numPixels[1:(x-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:x])
    if(x == 1){limitDown <- 1}
    int <- limitDown:limitUp
  lapply(cluster,function(y){
    data.frame(Spectra = medSpecP(peakM, which(clusterData[int] == y), norma), mz = peakM$mass, Image = groups[x], Cluster = y)
    })
  })  

clusSpec <- unlist(clusSpec, recursive = F)
clusSpec <- do.call(rbind, clusSpec)
clusSpec <- clusSpec[which(!is.na(clusSpec$Spectra)),]

plot1 <- ggplotly(ggplot(clusSpec) + geom_segment(aes(x = mz, xend = mz, y = 0, yend = Spectra, colour = paste(Image, "Cluster", Cluster))) + 
  theme_minimal() + theme(legend.title = element_blank()) + ylab("Intensity") + xlab("M/Z"))

return(plot1)
}


### Fold Change and P-values from medium spectrum

compareClusMedSpec <- function(refSpec, compSpec, peakMatrix, clusterData, norma){
  clusterData <- lapply(clusterData, function(x){   # Getting the clusterization from the k-means data
    return(x[1])  
  })
  clusterData <- unlist(clusterData)
  dt <- peakMatrix$intensity/norma
  
  ###################################################################### Mean Spectra
  
  specs2comp <- mapply(function(x, y){
    limitDown <- sum(peakMatrix$numPixels[1:(x-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:x])
    if(x == 1){limitDown <- 1}
    int <- limitDown:limitUp
    Specs <- dt[which(clusterData[int] == y),]
  
  }, c(refSpec[1], compSpec[1]), c(refSpec[2], compSpec[2]), SIMPLIFY = F)

  ###################################################################### P-Values
  
  pvalues <- sapply(1:length(peakMatrix$mass), function(com){
    test <- t.test(specs2comp[[1]][,com], specs2comp[[2]][,com])
    test$p.value
  })
  pvalues <- p.adjust(pvalues, method = "fdr")
  
  
  
  ###################################################################### Fold Change
  means2comp <- lapply(specs2comp, function(x){
    apply(x, 2, mean)
  })
  
  foldChange <- means2comp[[2]]/means2comp[[1]]
  change <- which(foldChange < 1)

# 
#   foldChange[change] <- means2comp[[1]][change]/means2comp[[2]][change]
  foldChange <- log2(foldChange)
  
  ###################################################################### Error estimation
  
  auTList <- c(196.966570, 393.933140, 590.899710, 787.866280, 984.832850)
  
  rmIndex <- sapply(auTList, function(tMass){
    diffMatrix <- (abs(peakMatrix$mass - tMass))
    index <- which.min(diffMatrix)
    return(index) ## This will return the index of the masses where the gold peaks are based on the real mass                                                 ## with the minimum difference with the teorical mass
  })
  auMass <- peakMatrix$mass[rmIndex]
  
  error <- abs((auTList - auMass)/auTList*10^6)
  
  massError <- sapply(peakMatrix$mass, function(mass){
    
    error[which.min(abs(auMass - mass))]
    
  })
  
  
  ###################################################################### Mounting the data-frame
  
  tag = paste("Image", compSpec[1] ,"Cluster",  compSpec[2],"Vs", "Image", refSpec[1], "Cluster", refSpec[2])
  compData <- list()
  compData[[1]] <- data.frame("mz" = peakMatrix$mass, "Log2FC" = foldChange, "pvalues" = pvalues, "ErrorEstimation" = massError)

  names(compData) <- tag
  return(compData)
}

volcanoPlotJ <- function(data2plot){
  
  volcPlot <-ggplot(data2plot[[1]]) + geom_point(aes(x = Log2FC, y = -log10(pvalues), color = mz)) + geom_vline(aes(xintercept = log2(2)), linetype = "dashed") +
    geom_hline(aes(yintercept = log10(0.05), linetype = "1")) + geom_vline(aes(xintercept = -log2(2)),linetype = "dashed") + 
    theme_minimal() + theme(legend.position = "none")
  return(volcPlot)
}