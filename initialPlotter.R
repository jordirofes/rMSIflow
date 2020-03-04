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

initialPlotter <- function(peakMatrix, matrixL, mz, img, groups, norm, sv, file){
  a <- vector("list", length(img))
  if(is.na(groups)){
    groups <-  peakMatrix$names
    }
  if (img == 0){
      a <- rMSIproc::plotPeakImageG(peakMatrix, mz, plot_labels = groups, normalization = norm)
      if(sv == T){svPlot(file, a)}
      return(a)
  } else{
    for (j in (1:length(img))){
      a[[j]] <- rMSIproc::plotPeakImageG(matrixL[[img[j]]], mz, plot_labels = groups[img[j]], normalization = norm)
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

pca <- function(peakMatrix, cnt, scl, norm = NA){
  dataPca <- peakMatrix$intensity
  if(!is.na(norm)){
    dataPca <- dataPca/norm
    }
     colnames(dataPca) <- peakMatrix$mass
  pca <- prcomp(dataPca, center = cnt, scale. = scl) ### Fem la PCA de la matriu de pics
  return(pca)
}

### This function creates a two dimensional plot with the especified principal components

pca2dimPlot <- function(peakMatrix, pca, grpImg, pc1, pc2, sv2, fileN2) {
  if(!is.na(grpImg)){
  listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, grpImg, 1:length(grpImg))
  } else {listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, peakMatrix$names, 1:length(peakMatrix$names)) }
  
  vectorGroups <- unlist(listGroups)
  pcaData <- list()
  pcaData[[1]] <- pca$x[,pc1]
  pcaData[[2]] <- pca$x[,pc2]
  pcaData[[3]] <- vectorGroups
  dataPca <- as.data.frame(pcaData)
  colnames(dataPca) <- c("PC1", "PC2", "Groups")
  rownames(dataPca) <- c(1:length(dataPca$PC1))
  sdev <- round(pca$sdev/sum(pca$sdev)*100)
  
  plotPca <- ggplot(dataPca,aes(x = PC1, y = PC2, colour = Groups))+ geom_point(alpha = 0.5)+ stat_ellipse() + 
    xlab(paste0("PC", pc2plot1 ,"(", sdev[pc2plot1],"%)")) + ylab(paste0("PC",pc2plot2, "(", sdev[pc2plot2],"%)")) + 
    theme_minimal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_point(aes(x= 0, y =0), color= "black") 
  if(sv2 == T){
    svPlot(fileN2, plotPca)
  }
  return(plotPca)
}

### This function plots a chosen principal component into selected images

pcaPlotImg <- function(peakMatrix, matrixLst, pca, img ,grpImg, pc, sv, fileN){
  a <- vector("list", length(img))
  if(is.na(grpImg)){ grpImg <- peakMatrix$names}
  
  if (img == 0){
    a <- rMSIproc::plotValuesImageG(peakMatrix, pixel_values = pca$x[,pc], plot_labels =  grpImg)
      if(sv == T){svPlot(fileN, a)}
    return(a)
  } else{
    for (j in (1:length(img))){
      
      limitDown <- sum(peakMatrix$numPixels[1:(img[j]-1)])+1
      limitUp <- sum(peakMatrix$numPixels[1:img[j]])
      if(img[j] == 1){limitDown <- 1}
      int <- limitDown:limitUp
      a[[j]] <- rMSIproc::plotValuesImageG(matrixLst[[img[j]]], pixel_values = pca$x[int,pc] , plot_labels = grpImg[img[j]])
    }
  }
  pl <- lapply(1:length(img), function(i){
    plotly::ggplotly(p = a[[i]])
  })
    if(sv == T){svPlot(fileN, a)}
  return(plotly::subplot(pl))
  
}

pcChooseNum <- function(pca, perc = 0.9){
  pca$sdev
  
  pcVar <- sapply(pca$sdev, function(dev){
    dev/sum(pca$sdev)*100
  })
  
  cumulative <- sapply(1:length(pca$sdev), function(x){sum(pca$sdev[1:x])})/sum(pca$sdev)*100
  
  pcplotData <- data.frame(Cumulative = cumulative, PC = 1:length(pca$sdev))
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

medSpecP <- function(peakMatrix, pixels, norm = NA){
  data <- peakMatrix$intensity
  if(!is.na(norm)){
    data <- data/norm
  }
  pixelMatrix <- peakMatrix$intensity[pixels,]
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
  
kmeansCluster <- function(peakMatrix, norm, numCl){
  clusData <- peakMatrix$intensity
  
  if(!is.na(norm)){
    clusData <- clusData/norm
  }
  
  clData <- list()
  
  for (i in 1:length(peakMatrix$numPixels)){ 
    limitDown <- sum(peakMatrix$numPixels[1:(i-1)])+1
    limitUp <- sum(peakMatrix$numPixels[1:i])
    int <-limitDown:limitUp
    if (i == 1){ limitDown <- 1}                                                 
    
      clData[[i]] <- kmeans(clusData[limitDown:limitUp, ], numCl)
    
  }
  return(clData)
  
}

clusterDataPlotting <- function(peakMatrix, matrixL, img, groups, clusterData, sv, file){
  
  clusterData <- lapply(clusterData, function(x){
    return(x[1]) 
  })
   if(img == 0){
    allClusImages <- unlist(clusterData)
    kplot <- rMSIproc::plotClusterImageG(peakMatrix, allClusImages)
    if(sv == T){svPlot(file, kplot)}
    return(kplot)
  } else{
  
    allClusImages <- unlist(clusterData, recursive = FALSE)
  kplot <- lapply(img, function(x){
      plots <- rMSIproc::plotClusterImageG(matrixL[[x]], allClusImages[[x]])
      return(plots)
  })
  if(sv == T){svPlot(file, kplot)}
  return(subplot(kplot))
  kplot <- lapply(kplots, function(x){
    plots <- ggplotly(kplots)
    return(plots)
  })
  return(subplot(kplot))
  }
  
}

clusterKmeansComparison <- function(peakMatrix, norm, img, groups, cluster, clusterData, sv, file){

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
    data.frame(Spectra = medSpecP(peakM, which(clusterData[int] == y), norm), mz = peakM$mass, Image = groups[x], Cluster = y)
    })
  })  
clusSpec <- unlist(clusSpec, recursive = F)
clusSpec <- do.call(rbind, clusSpec)
plot1 <- ggplotly(ggplot(clusSpec) + geom_segment(aes(x = mz, xend = mz, y = 0, yend = Spectra, colour = paste(Image, "Cluster", Cluster))) + 
  theme_minimal() + theme(legend.title = element_blank()) + ylab("Intensity") + xlab("M/Z"))

return(plot1)
}
