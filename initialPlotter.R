### This funcion separates all the images

sepMatrix <- function(peakmatrix){
  matrixList <- list()
  
  for(i in 1:length(peakM$numPixels)){
    limitDown <- sum(peakM$numPixels[1:(i-1)])+1
    limitUp <- sum(peakM$numPixels[1:i])
    if(i == 1){limitDown <- 1}
    int <- limitDown:limitUp
    matrixList[[i]] <- rMSIproc::`[.rMSIprocPeakMatrix`(peakM, int)
  }
  return(matrixList)
}

### This function plots the intensities of a given mass on the images

initialPlotter <- function(pM, matrixL, mz, img, groups, norm = NA, sv, file){
  a <- vector("list", length(img))
  if(is.na(groups)){
    groups <-  pM$names}
  
  if (img == 0){
      a <- rMSIproc::plotPeakImageG(pM, mz, plot_labels = groups, normalization = norm)
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

pca <- function(peakMatrix, dataPca , cnt, scl, norm = NA){
  norm <- toupper(norm)
  if(!is.na(norm)){
    if(norm == "TIC"){ dataPca <- dataPca/peakMatrix$normalization$TIC}
    else if(norm == "RMS"){dataPca <- dataPca/peakMatrix$normalization$RMS}
    else if(norm  == "ACQTIC"){dataPca <- dataPca/peakMatrix$normalization$AcqTic}
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
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
  sdev <- round(pca$sdev/sum(pca$sdev)*100, 4)
  
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

medSpecP <- function(peakMat, dataMat ,pixels, normalization = NA){
  data <- peakMat$intensity
  normalization <- toupper(normalization)
  if(!is.na(normalization)){
    if(normalization == "TIC"){data <- data/peakMat$normalization$TIC}
    else if(normalization == "RMS"){data <- data/peakMat$normalization$RMS}
    else if(normalization == "ACQTIC"){data <- data/peakMat$normalization$AcqTic}
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
  }
  pixelMatrix <- peakMat$intensity[pixels,]
  avgSpecpm <- apply(pixelMatrix, 2, mean)
                             
  return(avgSpecpm)
}

# medSpecComp <- function(peakMat, pixels1, pixels2, normalization, name1, name2, sav2, filename2)
  
  
### K-means
  
kmeansCluster <- function(peakmatrix, clusData ,norm4, numCl){
  
  norm <- toupper(norm4)
  if(!is.na(norm)){
    if(norm == "TIC"){ clusData <- clusData/peakmatrix$normalization$TIC}
    else if(norm == "RMS"){clusData <- clusData/peakmatrix$normalization$RMS}
    else if(norm  == "ACQTIC"){clusData <-clusData/peakmatrix$normalization$AcqTic}
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
  }
  clData <- list()
  
  for (i in 1:length(peakmatrix$numPixels)){ 
    if (i == 1){
      clData[[i]] <- kmeans(clusData[1:peakmatrix$numPixels[i], ], numCl)
    } else {
      clData[[i]] <- kmeans(clusData[(sum(peakmatrix$numPixels[1:(i-1)])+1):sum(peakmatrix$numPixels[1:i]), ], numCl)
    }
  }
  clusImages <- lapply(clData, function(x){
    return(x[1]) 
  })
  allClusImages <- unlist(clusImages)
  
  
  return(plotly::ggplotly(rMSIproc::plotValuesImageG(peakmatrix, allClusImages)))
  
}

