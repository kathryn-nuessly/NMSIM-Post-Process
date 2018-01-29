###############################################################################
#
# NMSIM_V3.R
#
# 2017-02-07
# Damon Joyce (damon_joyce@nps.gov)
#
# These functions will convert and analyze NMSIM tig files.
#
# ConvertTIG2RDATA(tigFile): Convert TIG file into a binary RData list object,
#   with the number of elements equal to the number of points
#
###############################################################################

# Defines Conversion Function & Defines how to read through header --------
ConvertTIG2RDATA <- function(tigFile)
{
  fData <- readLines(tigFile, n=50)
	numI <- 0
	numJ <- 0
	
	fileHeaderLen <- 0
	utmPos <- 0
	ptHeaderLen <- 0
	
	linesPerPoint <- 0
	
	for(i in 1:30)
	{
		if(grepl("^ii",fData[i])){
			numI <- as.numeric(gsub("^.+[[:space:]]+","",fData[i]))
		}else if(grepl("^jj",fData[i])){
			numJ <- as.numeric(gsub("^.+[[:space:]]+","",fData[i]))
		}else if(grepl("---End File Header---",fData[i])){
			fileHeaderLen <- i
		}else if(grepl("^Integ",fData[i])){
			linesPerPoint <- as.numeric(gsub("^.+[[:space:]]+Spec:","",fData[i]))
		}else if(grepl("^UTM",fData[i])){
			utmPos <- i
		}else if(grepl("^SP#",fData[i])){
			ptHeaderLen <- i
		}
	}



# Defines UTM Coordinates from i,j Points ---------------------------------

	utmPos <- utmPos - fileHeaderLen
	ptHeaderLen <- ptHeaderLen - fileHeaderLen + 2
	
	utmData <- NULL
	
	myConn <- file(tigFile, open = "r")
	hdrCount <- sum(nchar(readLines(myConn, n=fileHeaderLen)) + 1 + 1)
	mainCount <- sum(nchar(readLines(myConn, n = (ptHeaderLen + linesPerPoint + 1))) + 1 + 1)
	close(myConn)
	
	nLines <-  (file.size(tigFile) - hdrCount) / mainCount
	
	nmsimData <- as.list(rep(NA,numI * numJ))
	
	utm.X <- rep(NA, numI)
	utm.Y <- rep(NA, numJ)
	
	###
	myConn <- file(tigFile, open = "r", raw=TRUE)
	invisible(readLines(myConn, n=fileHeaderLen))
	ptm <- proc.time()
	for(i in 1:100)
	{
	  hdrData <- readLines(myConn, n = ptHeaderLen)
	  
	  utmData <- c(as.integer(substr(hdrData[utmPos],28,33)), as.integer(substr(hdrData[utmPos],49,55)))
	  
	  iIdx <- as.integer(substring(hdrData[1],7,10))
	  jIdx <- as.integer(substring(hdrData[1],11,14))
	  
	  utm.X[iIdx] <- utmData[1]
	  utm.Y[jIdx] <- utmData[2]
	  
	  ptIdx <- as.integer(((iIdx - 1) * numJ) + jIdx)
	  
	  tempData <- array(scan(file=myConn,nlines=linesPerPoint,quiet=T),dim=c(37,linesPerPoint))[c(2,4,37),]
	  invisible(readLines(myConn, n=1))
	  
	  nmsimData[[ptIdx]] <- cbind(as.integer(1000*tempData[1,]),as.integer(tempData[2,]),as.integer(tempData[3,]))
	}
	
	ptm <- proc.time() - ptm
	ptm <- round(ptm[3] / 100 * nLines,0)
	
	print(paste("Estimated conversion time: ", sprintf("%02d:%02d:%02d", floor(ptm / 3600), floor(ptm / 60), ptm %% 60), sep=""))
	###
	
	myConn <- file(tigFile, open = "r")
	
	invisible(readLines(myConn, n=fileHeaderLen))
	
	myPB <- txtProgressBar(style=3)
	
	ptm <- proc.time()
	
	for(i in 1:nLines)
	{
		hdrData <- readLines(myConn, n = ptHeaderLen)

		utmData <- c(as.integer(substr(hdrData[utmPos],28,33)), as.integer(substr(hdrData[utmPos],49,55)))
		
		iIdx <- as.integer(substring(hdrData[1],7,10))
		jIdx <- as.integer(substring(hdrData[1],11,14))
	
		utm.X[iIdx] <- utmData[1]
		utm.Y[jIdx] <- utmData[2]
		
		ptIdx <- as.integer(((iIdx - 1) * numJ) + jIdx)

		tempData <- array(scan(file=myConn,nlines=linesPerPoint,quiet=T),dim=c(37,linesPerPoint))[c(2,4,37),]
		invisible(readLines(myConn, n=1))

		temp1 <- tempData[1,]
		
		temp2 <- tempData[2,]
		temp2[temp2 == -999] <- NA
		
		temp3 <- tempData[3,]
		temp3[temp3 == -93] <- NA
		
		tempData <- cbind(as.integer(1000*temp1),as.integer(temp2),as.integer(temp3))
		colnames(tempData) <- c("time_ms","awt_cb","dp")
		
		nmsimData[[ptIdx]] <- tempData

		setTxtProgressBar(myPB, i/nLines)
	}
	
	ptm <- proc.time() - ptm
	ptm <- round(ptm[3], 0)
	print(paste("Actual conversion time: ", sprintf("%02d:%02d:%02d", floor(ptm / 3600), floor(ptm / 60), ptm %% 60), sep=""))
	
	close(myConn)
  close(myPB)
  
  save(nmsimData, utm.X, utm.Y, file = gsub(".tig$",".Rdata",tigFile))
}



# Defines Run Time to Preserve Time for future compression calcula --------
nmsimRuntime <- function(nData)
{
  out <- rep(NA, length(nData))
  for(i in 1:length(nData))
  {
    out[i] <- max(nData[[i]]$timestep) - min(nData[[i]]$timestep)
  }
  mean(out)  
}