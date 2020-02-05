matrixSeparation <- function(pM){
  
  pkl <- vector("list", length = length(pM$numPixels))
  normList <- vector("list", length = length(pM$normalizations))
  
 for(l in (1:length(pM$numPixels))){
    pkl[[l]] <- pM
    
    limitDown <- sum(pM$numPixels[1:(l-1)])+1
    limitUp <- sum(pM$numPixels[1:l])
    if(l == 1){limitDown <- 1}
    int <- limitDown:limitUp
    
    pkl[[l]]$intensity <- pM$intensity[int,]
    pkl[[l]]$SNR <- pM$SNR[int,]
    pkl[[l]]$area <- pM$area[int,]
    pkl[[l]]$pos <- pM$pos[int,]    
    pkl[[l]]$posMotors <- pM$posMotors[int,]

    pkl[[l]]$numPixels <- pM$numPixels[l]
    pkl[[l]]$names <- pM$names[l]
    pkl[[l]]$uuid <- pM$uuid[l]
    
     for(z in 1:length(pM$normalizations)){
       normList[[z]] <- pM$normalizations[[z]][int]
     }
    normList2 <- as.data.frame(normList)
    colnames(normList2) <- colnames(pM$normalizations)
    pkl[[l]]$normalizations <- normList2
 }
  
  return(pkl)
    
  }


initialPlotter <- function(pM, matrixL, mz, img, groups, norm = NA){
  a <- vector("list", length(img))
  if (img == 0){
      a <- rMSIproc::plotPeakImageG(pM, mz, plot_labels = groups, normalization = norm)
      return(a)
  } else{
    for (j in (1:length(img))){
      a[[j]] <- rMSIproc::plotPeakImageG(matrixL[[img[j]]], mz, plot_labels = groups[img[j]], normalization = norm)
    }
  }
  return(a)
  }
 


