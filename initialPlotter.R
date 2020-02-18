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
  return(pl)
}

svPlot <- function(flName, plots2save){
  pdf(flName)
  print(plots2save)
  dev.off()

}

### This function makes a PCA from the peak matrix

pca <- function(peakMatrix, cnt, scl, norm = NA){
  dataRaw <- peakMatrix$intensity
  norm <- toupper(norm)
  if(!is.na(norm)){
    if(norm == "TIC"){dataRaw/peakMatrix$normalization$TIC}
    else if(norm == "RMS"){dataRaw/peakMatrix$normalization$RMS}
    else if(norm  == "ACQTIC"){dataRaw/peakMatrix$normalization$AcqTic}
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
  }
     colnames(dataRaw) <- peakMatrix$mass
  pcaRaw <- prcomp(dataRaw, center = cnt, scale. = scl) ### Fem la PCA de la matriu de pics
  return(pcaRaw)
}

## This function creates a two dimensional plot with the especified principal components

pca2dimPlot <- function(peakMatrix, pca, grpImg, pc1, pc2, sv2, fileN2) {
  listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, grpImg, 1:length(grpImg))
  
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

pcaPlotImg <- function(peakMatrix, matrixLst, pca, img ,grpImg, pc, sv, fileN){
  a <- vector("list", length(img))
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

medSpec <- function(peakMat, pixels, normalization = NA){
  data <- peakMat$intensity
  normalization <- toupper(normalization)
  if(!is.na(normalization)){
    if(normalization == "TIC"){data/peakMat$normalization$TIC}
    else if(normalization == "RMS"){data/peakMat$normalization$RMS}
    else if(normalization == "ACQTIC"){data/peakMat$normalization$AcqTic}
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
  }
  pixelMatrix <- peakMat$intensity[pixels,]
  avgSpecpm <- apply(pixelMatrix, 2, mean)
                             
  return(avgSpecpm)
}

medSpecComp <- function(peakMat, pixels1, pixels2, normalization, name1, name2){
  medSpec1 <- medSpec(peakMat, pixels1, normalization)
  medSpec2 <- medSpec(peakMat, pixels2, normalization)
  medSpecdf <- as.data.frame(cbind(medSpec1, medSpec2, peakM$mass))
  colnames(medSpecdf) <- c("Spec1", "Spec2", "Mass") 

  aveSpec <- plotly::ggplotly(ggplot(medSpecdf) + geom_segment(aes(x = Mass,y = 0, xend = Mass, yend = Spec1, colour = name1)) +
             geom_segment(aes(x = Mass, y = 0, xend = Mass, yend = Spec2, colour = name2))  +
             xlab("Mass") + ylab("Intensity") + theme_minimal() + theme(legend.title = element_blank())
             
             )
  return(aveSpec)
}

