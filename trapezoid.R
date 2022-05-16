#!/usr/bin/R
# Ryan Preble, Spring 2022
# Boyce Thompson Institute, 533 Tower Road, Ithaca, NY 14853


# CHECK INSTALLED PACKAGES. INSTALL IF NECESSARY
packages <- c("mzR", "Spectra", "MSnbase", "optparse")
list.of.packages <- c("mzR", "Spectra", "MSnbase", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# REQUIRED LIBS
library(mzR)
library(Spectra)
library(MSnbase)
library(optparse)


### KEY FUNCTIONS
trapezoid <- function(d){
  rtime <- d$rtime
  intense <- d$intensity
  
  intense[is.na(intense)] <- 0
  
  area <- 0
  
  for (i in 1:(length(rtime)-1)){
    area <- area + ((0.5)*(intense[i] + intense[i+1])*(rtime[i+1] - rtime[i]))
  }
  
  return( area )
}


### TAKE OPTIONS
optionlist <- list(
  
  make_option(c("-m","--mz"), action="store", type="double", dest="mz", default=NULL, help="Mass number to search for. By default, all mass numbers in sample are considered."),
  make_option(c("-t","--delta"), action="store", type="double", dest="delta", default=0.005, help="The tolerance for the given mass number. Default is 0.005 mass units."),
  make_option(c("-r","--minrt"), action="store", type="double", dest="minrt", default=0, help="The minimum retention time, in seconds, to be used. Default is zero seconds."),
  make_option(c("-R","--maxrt"), action="store", type="double", dest="maxrt", help="The maximum retention time, in seconds, to be used. If this is not supplied, retention times from minrt to the end of the experiment will be considered."),
  make_option(c("-f","--file"), action="store", type="character", dest="file", default=NULL, help="The mzXML/mzML file to be analyzed. Overrides \"-dir\", which uses all files in the given directory."),
  make_option(c("-d","--dir"), action="store", type="character", dest="dir", default=getwd(), help="The directory with mzXML/mzML files to be quantified. If \"-file\" is supplied, that single file will be used instead. Default is the working directoy. Files should have the correct file extension."),
  make_option(c("-o","--outfile"), action="store", type="character", dest="outfile", default=NULL, help="The file that receives output data. Will auto-generate a name by default"),
  make_option(c("-x", "--format"), action="store", type="character", dest="format", default="mzML", help="The file format to be used. Can be mzXML or mzML, mzML by default.")
  
)

parser <- OptionParser(option_list=optionlist)
arguments <- parse_args(parser)

### PROCESS OPTIONS
file <- arguments$file
dir <- arguments$dir
mz <- arguments$mz
delta <- arguments$delta
minrt <- arguments$minrt
maxrt <- arguments$maxrt
outfile <- arguments$outfile
format<- arguments$format

infiles <- 0
outfile_mz <- NULL
outfile_maxrt <- NULL


if(format != "mzML" && format != "mzXML"){
  warning("Bad file type supplied! mzML will be used by default.")
  format <- "mzML"
}
if(!is.null(file)){
  message(paste("File name supplied. Data will be read from ", file, "\n"), sep="")
  infiles <- c(file)
} else {
  message(paste("Reading all ", format ," files in directory: ", dir, "\n"))
  setwd(dir = dir)
  infiles <- list.files(dir, pattern=paste(".", format, sep=""))
}

if(is.null(mz)){
  warning("No mass number supplied. All mass numbers in experiment will be considered.\n")
  outfile_mz <- "all"
} else {
  outfile_mz <- mz
}
if(is.null(maxrt)){
  warning("Retention time window not specified. Retention times to end of experiment will be considered.\n")
  outfile_maxrt <- "max"
} else {
  outfile_maxrt <- maxrt
}

if(is.null(outfile)){
  warning("Outfile not specified. Data will be written to an auto-generated name.")
  outfile <- paste("mz-", outfile_mz, "_rt-", minrt, "-", outfile_maxrt, ".csv", sep="")
}

message("Beginning analysis. This may take some time.")

feature_area <- c()
files <- c()

### MAIN PROCESSING STEP
for (i in 1:length(infiles)) {
  tryCatch({
    ### READ IN FILE
    mzXML <- MSnbase::readMSData(infiles[i], mode = "onDisk")
    files <- c(files, infiles[i])
    
    ### ISOLATE DATA
    if(!is.null(mz)){ #subset by mz if that field was supplied
      mzXML <- filterMz(mzXML, mz=c(mz-delta, mz+delta), msLevel=1)
    }
    if(is.null(maxrt)){ #subset according to min/max rt supplied
      mzXML <- filterRt(mzXML, c(minrt, max(mzXML@featureData@data$retentionTime)))
    } else {
      mzXML <- filterRt(mzXML, c(minrt, maxrt))
    }
    
    ### EXTRACT CHROMATOGRAM
    chromatogram <- chromatogram(mzXML, aggregationFun="sum")
    
    ### ANALYZE CHROMATOGRAM CONTENT
    data <- chromatogram@.Data[1]
    intensities <- as.data.frame(data[[1]])
    feature_area <- c(feature_area, trapezoid(intensities))
  },
  error=function(error_message){
    message("An error occurred. Continuing because this script probably won't break anything.\n")
    message(error_message)
    message("\n")
  })
}

dataframe <- data.frame(files, feature_area)

write.csv(dataframe, file=outfile)

message("Done.\n")






