###############################################################################
#
# NMSIM_V3.R
#
# 2016-09-21
# Damon Joyce (damon_joyce@nps.gov)
#
# These functions will call the script and data files to convert TIG files to
# R Datafiles
#
###############################################################################

# Pick the R Script for Conversion ----------------------------------------
if(!exists("ConvertTIG2RDATA", mode="function"))
{
  myRscript <- choose.files(caption="Choose NMSIM.R script",multi=F, filters=Filters["R",])
  source(myRscript)
}

# Choose the folder location for datafiles for conversion -----------------
myDir <- choose.dir()
myFiles <- list.files(myDir,"*.ti.",recursive=T, full.names=T)


# Run the Conversion Script -----------------------------------------------
for(t in myFiles)
{
  ConvertTIG2RDATA(t)
}
