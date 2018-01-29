###############################################################################
#
# NMSIM_V3.R
#
# 2017-02-07
# Damon Joyce (damon_joyce@nps.gov)
#
# These functions will analyze the data from NMSIM tig files.
#
# 
# nmsimTA(data, threshold_in_decibels):
#   Return time above your threshold per point (nSecs)
# nmsimTAUD(data, threshold_in_decibels):
#   Return time audible/detectable per point(nSecs)
# nmsimSEL(data):
#   Return SEL per point
# nmsimRuntime(data):
#   Return total model vehicle runtime
#
###############################################################################
# Pick the R Script for Compression (RData files to csv files) ----------------------------------------
setwd("F:\\THRO_TAR1818\\NMSIM\\Output_Data\\LNATinclu_60m_2sec\\")

LNAT <- 28

LEQdur <- 24 * 3600

# Load in Traffic Scenario (number of cars) - use the LOADING script --------
loadingFactors <- read.table("TrafficScenario_THRO_2014.tsv", header=TRUE)

# Define the datafiles to be used in analysis -----------------------------
rDataFiles <- dir(pattern="*.Rdata")
# Define TimeCompression Equation (makes 3D matrix with time preserved into 2D matrix) -----------------------
tCompress <-function(nOps, singleOpVal, periodMin)
{  
  N <- sum(nOps)
  if(N != 0)
  {
    A <- sum(nOps * singleOpVal)
    output <- 100 * N * (1-exp(-1 * ((N + 1) * A)/(periodMin * N)))/(N + 1)
    if(is.nan(output)) output <- NA
    output
  }else NA
}



# Define fillInGaps function (Fill in the data from analysis into the blank "NA" matrix) --------------------------
fillInGaps <- function(vec)
{
  utmSpacing <- mean(diff(vec), na.rm=T)
  firstDataPt <- which(!is.na(vec))[1]
  utmFixed <- rep(NA, length(vec))
  
  #center
  utmFixed[firstDataPt] <- vec[firstDataPt]
  utmFixed[(firstDataPt+1):length(vec)] <- vec[firstDataPt] + round(cumsum(rep(utmSpacing, length(vec) - firstDataPt)))
  if(firstDataPt != 1)
  {
    utmFixed[(firstDataPt-1):1] <- vec[firstDataPt] - round(cumsum(rep(utmSpacing, firstDataPt - 1)))
  }
  utmFixed
}

# Define Time Audible Function --------------------------------------------
nmsimTAUD <- function(nData, threshold)
{
  temp <- rep(0, length(nData))
  for(i in 1:length(nData))
  {
    nData.s <- nData[[i]]
    if(length(nData.s) > 1)
    {
      timeDiff <- diff(nData.s[,"time_ms"]) / 2
      timeAtPt <- c(timeDiff, 0) + c(0, timeDiff)
      temp[i] <- sum(timeAtPt[nData.s[,"dp"] > (threshold)], na.rm=T) / 1000
    }
  }
  temp
}

# Define Time Above Function ----------------------------------------------
nmsimTA <- function(nData, threshold, splAmb)
{
  temp <- rep(0, length(nData))
  for(i in 1:length(nData))
  {
    nData.s <- nData[[i]]
    if(length(nData.s) > 1)
    {
      if(!is.na(splAmb)) nData.s[,"awt_cb"] <- 100*log10(10^(nData.s[,"awt_cb"]/100) + 10^(splAmb/10))
      timeDiff <- diff(nData.s[,"time_ms"]) / 2
      timeAtPt <- c(timeDiff, 0) + c(0, timeDiff)
      temp[i] <- sum(timeAtPt[nData.s[,"awt_cb"] > (threshold * 10)], na.rm=T) / 1000
    }
  }
  temp
}

# Define SEL Function -----------------------------------------------------
nmsimSEL <- function(nData)
{
  temp <- rep(-Inf, length(nData))
  for(i in 1:length(nData))
  {
    nData.s <- nData[[i]]
    if(length(nData.s) > 1)
    {
      nData.s[which(is.na(nData.s[,"awt_cb"])),"awt_cb"] <- -Inf
      timeDiff <- diff(nData.s[,"time_ms"]) / 2
      timeAtPt <- c(timeDiff, 0) + c(0, timeDiff)
      temp[i] <- round(100*log10(sum(10^(nData.s[,"awt_cb"] / 100) * (timeAtPt / 1000))), 0)
      #centibels
      if(temp[i] == 0) temp[i] <- -Inf
    }
  }
  temp
}

# Load R DataFiles In & Combine Data -----------------------------------------------------
load(rDataFiles[1])

#combData <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))

taudData <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
ta35Data <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
ta40Data <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
ta45Data <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
ta52Data <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
ta60Data <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))
selData <- array(NA, dim=c(dim(loadingFactors)[1], length(nmsimData)))

outData <- NULL
outData$taud <- rep(NA, length(nmsimData))

outData$ta35 <- rep(NA, length(nmsimData))
outData$ta40 <- rep(NA, length(nmsimData))
outData$ta45 <- rep(NA, length(nmsimData))
outData$ta52 <- rep(NA, length(nmsimData))
outData$ta60 <- rep(NA, length(nmsimData))

outData$laeq <- rep(NA, length(nmsimData))
outData$sel <- rep(NA, length(nmsimData))


for(i in 1:dim(loadingFactors)[1])
{
  if(i > 1) load(rDataFiles[grep(loadingFactors[i,1],rDataFiles)])
  taudData[i,] <- nmsimTAUD(nmsimData, 7)
  
  ta35Data[i,] <- nmsimTA(nmsimData, 35, LNAT)
  ta40Data[i,] <- nmsimTA(nmsimData, 40, LNAT)
  ta45Data[i,] <- nmsimTA(nmsimData, 45, LNAT)
  ta52Data[i,] <- nmsimTA(nmsimData, 52, LNAT)
  ta60Data[i,] <- nmsimTA(nmsimData, 60, LNAT)
  
  selData[i,] <- nmsimSEL(nmsimData)
}

for(i in 1:length(nmsimData))
{
  outData$taud[i] <- tCompress(loadingFactors$NUM, taudData[,i], LEQdur)
  
  outData$ta35[i] <- tCompress(loadingFactors$NUM, ta35Data[,i], LEQdur)
  outData$ta40[i] <- tCompress(loadingFactors$NUM, ta40Data[,i], LEQdur)
  outData$ta45[i] <- tCompress(loadingFactors$NUM, ta45Data[,i], LEQdur)
  outData$ta52[i] <- tCompress(loadingFactors$NUM, ta52Data[,i], LEQdur)
  outData$ta60[i] <- tCompress(loadingFactors$NUM, ta60Data[,i], LEQdur)
  
  outData$laeq[i] <- 100*log10(sum(loadingFactors$NUM * 10^(selData[,i]/100)) / LEQdur + 10^(LNAT/10))
  outData$sel[i]  <- 100*log10(sum(loadingFactors$NUM * 10^(selData[,i]/100)) + (LEQdur * 10^(LNAT/10)))
}

utm.X <- fillInGaps(utm.X)
utm.Y <- fillInGaps(utm.Y)
# Calculate Time Audible (d' > 7) --------
##export TAUD layer to datafile 
plotData <- matrix(outData$taud, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TAUD")

write.table(file="2014_TAUD.csv", plotData, sep=',', row.names=F, quote=F)

# Time Above Metrics ------------------------------------------------------

# Calculate Time Above (d' > 35) - To calculate Time above, change d' (to 30, etc.) --------

##export TABV layer to datafile
plotData <- matrix(outData$ta35, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TABV_35")

write.table(file="2014_TABV_35.csv", plotData, sep=',', row.names=F, quote=F)

# Calculate Time Above (d' > 40)  --------

##export TABV layer to datafile
plotData <- matrix(outData$ta40, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TABV_40")

write.table(file="2014_TABV_40.csv", plotData, sep=',', row.names=F, quote=F)

# Calculate Time Above (d' > 45) --------

##export TABV layer to datafile
plotData <- matrix(outData$ta45, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TABV_45")

write.table(file="2014_TABV_45.csv", plotData, sep=',', row.names=F, quote=F)

# Calculate Time Above (d' > 52)  --------

##export TABV layer to datafile
plotData <- matrix(outData$ta52, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TABV_52")

write.table(file="2014_TABV_52.csv", plotData, sep=',', row.names=F, quote=F)

# Calculate Time Above (d' > 60)  --------

##export TABV layer to datafile
plotData <- matrix(outData$ta60, nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","TABV_60")

write.table(file="2014_TABV_60.csv", plotData, sep=',', row.names=F, quote=F)

# Calculate SEL and LEQ ---------------------------------------------------

##export LEQ layer to datafile
plotData <- matrix(round(outData$laeq,0), nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","LEQ24_CB")

write.table(file="2014_LEQ.csv", plotData, sep=',', row.names=F, quote=F, na="-999")

##export SEL to datafile
plotData <- matrix(round(outData$sel,0), nrow=length(utm.X), ncol=length(utm.Y), byrow=T)

tempLims <- data.frame(x=c(615000, 635000), y=c(5270000,5277000))

subIdx.X <- which(utm.X > tempLims$x[1] & utm.X < tempLims$x[2])
subIdx.Y <- which(utm.Y > tempLims$y[1] & utm.Y < tempLims$y[2])

rownames(plotData) <- utm.X
colnames(plotData) <- utm.Y

plotData <- as.data.frame(as.table(plotData))
colnames(plotData) <- c("X","Y","SEL_CB")

write.table(file="2014_SEL.csv", plotData, sep=',', row.names=F, quote=F, na="-999")
