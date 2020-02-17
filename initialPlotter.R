
initialPlotter <- function(pM, matrixL, mz, img, groups, norm = NA, sv, file){
  a <- vector("list", length(img))
  if (img == 0){
      a <- rMSIproc::plotPeakImageG(pM, mz, plot_labels = groups, normalization = norm)
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


pcaPlotter <- function(peakMatrix, pc1, pc2, cnt, scl, grpImg , normalization = NA, sv2, fileN2){
  dataRaw <- peakMatrix$intensity
  toupper(normalization)
  if(!is.na(normalization)){
    if(normalization == "TIC"){dataRaw/peakMatrix$normalization$TIC}
    else if(normalization == "RMS"){dataRaw/peakMatrix$normalization$RMS}
    else if(normalization == "ACQTIC"){dataRaw/peakMatrix$normalization$AcqTic}
    
    else {stop("The normalization name is not one of the list. Check if you've written it correctly or if you don't want the data to be normalized write a NA in normalization")}
  }
     colnames(dataRaw) <- peakMatrix$mass
  pcaRaw <- prcomp(dataRaw, center = cnt, scale. = scl) ### Fem la PCA de la matriu de pics
  
  listGroups <- mapply(function(groups,i){
    rep(groups,  peakMatrix$numPixels[i])
  }, grpImg, 1:length(grpImg))
  
  vectorGroups <- unlist(listGroups)
  pcaData <- list()
  pcaData[[1]] <- pcaRaw$x[,pc1]
  pcaData[[2]] <- pcaRaw$x[,pc2]
  pcaData[[3]] <- vectorGroups
  dataPca <- as.data.frame(pcaData)
  colnames(dataPca) <- c("PC1", "PC2", "Groups")
  rownames(dataPca) <- c(1:length(dataPca$PC1))
  sdev <- round(pcaRaw$sdev/sum(pcaRaw$sdev)*100)
  
  plotPca <- ggplot(dataPca)+ geom_point(aes(x = PC1, y = PC2, colour = Groups), alpha = 0.5) + 
             xlab(paste0("PC1 (", sdev[pc2plot1],"%)")) + ylab(paste0("PC1 (", sdev[pc2plot2],"%)")) + 
             theme_minimal() + geom_abline(slope = 0,intercept = 0) + geom_abline(slope = 1000, intercept = 0) + 
             geom_point(aes(x= 0, y =0), color= "black")
  
  if(sv2 == T){
    svPlot(fileN2, plotPca)
  }
  return(plotPca)
}
 


