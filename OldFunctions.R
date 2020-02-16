### Finding gold peaks

```{r}
## This function searches the peak matrix for the intensities of the gold peaks with teorical mass auTList
auTList <- c(196.966570, 393.933140, 590.899710, 787.866280, 984.832850)

rmIndex <- sapply(auTList, function(tMass){
  diffMatrix <- (abs(peakM$mass - tMass))
  index <- which.min(diffMatrix)
  return(index) ## This will return the masses where the gold peaks are based on the real mass                                                 ## with the minimum difference with the teorical mass
})
auMass <- peakM$mass[rmIndex]

#########################

## 
# rmIndex2 <- sapply(auTList, function(tMass){
#   diffMatrix <- (abs(tMass - peakM$mass)/tMass*10^6)
#   index <- which(diffMatrix < 50)
#   return(index) 
#   })
# auMass50 <- peakM$mass[rmIndex2]

#########################

## All intensities from every Au peak
aupeakMatrix <- peakM$intensity[,rmIndex]
aunames <- c("Au1", "Au2","Au3", "Au4","Au5")
colnames(aupeakMatrix) <- aunames


anotation <- rMSIproc::peakAnnotation(PeakMtx = peakM)

```

### Calculating Au intensity quotients

```{r}
# Mean intensities by image
meanImagePeakAu <- matrix(nrow = length(peakM$numPixels), ncol = length(auMass))

for(i in 1:length(peakM$numPixels)){
  limitDown <- sum(peakM$numPixels[1:(i-1)])+1
  limitUp <- sum(peakM$numPixels[1:i])
  if(i == 1){limitDown <- 1}
  int <- limitDown:limitUp
  for(z in 1:length(auTList)){
    meanImagePeakAu[i,z] <-  mean(aupeakMatrix[int, z])
  }
  colnames(meanImagePeakAu) <- c("Au1", "Au2","Au3", "Au4","Au5")
}

rMSIproc::plotPeakImage(peakMatrix = peakM, mz = 30)

# Calculating all Au quotients
auQ <- c()
quotient <- c()
qnames <- c()
for(a in 1:ncol(aupeakMatrix)){
  for(b in 1:ncol(aupeakMatrix)){
    if(a == b){next}
    quotient <- aupeakMatrix[,a]/aupeakMatrix[,b]
    auQ <- cbind(auQ, quotient)
    qnames <- as.list(cbind(qnames, paste0(colnames(aupeakMatrix)[a],colnames(aupeakMatrix)[b])))
  }
}

colnames(auQ) <- qnames


# Au and Au quotients plots
pdf(file = "Au intensities.pdf")
for(c in 1:ncol(aupeakMatrix)){
  print(plotValuesImageG(peakMatrix = peakM, pixel_values = aupeakMatrix[,c], title_label = colnames(aupeakMatrix)[c],scale_label = auTList[c]))
}
dev.off()

pdf(file = "Au quotients.pdf")
for(i in 1:ncol(auQ)){
  print(plotValuesImageG(peakMatrix = peakM, pixel_values = auQ[,i], title_label = colnames(auQ)[i]))
}
dev.off()


### Percentatge of pixels with Au peaks for image 

perImagePeakAu <- matrix(nrow = length(peakM$numPixels), ncol = length(auMass))
for(i in 1:length(peakM$numPixels)){
  limitDown <- sum(peakM$numPixels[1:(i-1)])+1
  limitUp <- sum(peakM$numPixels[1:i])
  if(i == 1){limitDown <- 1}
  int <- limitDown:limitUp
  for(z in 1:length(auTList)){
    perImagePeakAu[i,z] <-  length(which(aupeakMatrix[int,z] > 1))/length(aupeakMatrix[int,z])*100
  }
  colnames(perImagePeakAu) <- c("Au1", "Au2","Au3", "Au4","Au5")
}




apply(peakM$intensity, 1, max)

```


## Not needed as the last version of rMSIproc contais a peakmatrix subsetting function

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
## Not needed as the last version of rMSIproc contais a peakmatrix subsetting function

