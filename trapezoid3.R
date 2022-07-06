#!/usr/bin/R
# Ryan Preble, Spring 2022
# Boyce Thompson Institute, 533 Tower Road, Ithaca, NY 14853


# CHECK INSTALLED PACKAGES. INSTALL IF NECESSARY
packages <- c("mzR", "Spectra", "MSnbase", "argparse", "stringr")
list.of.packages <- c("mzR", "Spectra", "MSnbase", "argparse", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# REQUIRED LIBS
library(argparse, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(mzR, quietly = TRUE)
library(Spectra, quietly = TRUE)
library(MSnbase, quietly = TRUE)
library(stringr, quietly= TRUE)

### KEY FUNCTIONS #############################################################
#Precondition: data has already been filtered for retention time
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

#Precondition: data has not been filtered for retention time. d is an 
# OnDiskMSnExp.Retention times are supplied with left/right shifts taken into 
# account. 
#Filters background noise five seconds on either side of the chromatogram
background <- function(d, mymin, mymax){

  left <- mymin
  right <- mymax
  
  tempdata <- chromatogram(filterRt(d, c(left, right)))
  tempdata <- tempdata@.Data[1]
  tempdata <- as.data.frame(tempdata[[1]])

  intense <- tempdata$intensity
  
  intense[is.na(intense)] <- 0
  
  mybackground <- (mymax - mymin)*(intense[1]+intense[length(intense)])/2
  return(mybackground)
  
}

mainfunc <- function(myfile, mymz, mydelt, mymin, mymax, mytrap, mybackground){
  
  tryCatch({
    ### READ IN FILE
    mzXML <- MSnbase::readMSData(myfile, mode = "onDisk")
    
    ### ISOLATE DATA
    if(!is.null(mymz)){ #subset by mz if that field was supplied
      mzXML <- filterMz(mzXML, mz=c(mymz-mydelt, mymz+mydelt), msLevel=1)
    }
    if(is.null(mymax)){ #subset according to min/max rt supplied
      mzXML <- filterRt(mzXML, c(mymin, max(mzXML@featureData@data$retentionTime)))
    } else {
      mzXML <- filterRt(mzXML, c(mymin, mymax))
    }
    
    ### EXTRACT CHROMATOGRAM
    chromatogram <- chromatogram(mzXML, aggregationFun="sum")
    
    ### ANALYZE CHROMATOGRAM CONTENT
    data <- chromatogram@.Data[1]
    intensities <- as.data.frame(data[[1]])
    
    if (is.null(mymax)){
      return(c(myfile, mytrap(intensities)))
    }
    
    return(c(myfile, mytrap(intensities)-mybackground(mzXML, mymin-2, mymax+2)))
  },
  error=function(error_message){
    message("An error occurred. Continuing because this script probably won't break anything.\n")
    message(error_message)
    message("\n")
    return(c(myfile, "ERROR"))
  })
}
###############################################################################



### TAKE OPTIONS ##############################################################
parser <- ArgumentParser()

parser$add_argument("-m","--mz", action="store", type="double", dest="mz", default=NULL, help="Mass number to search for. By default, all mass numbers in sample are considered.")
parser$add_argument("-t","--delta", action="store", type="double", dest="delta", default=0.005, help="The tolerance for the given mass number. Default is 0.005 mass units.")
parser$add_argument("-r","--minrt", action="store", type="double", dest="minrt", default=0, help="The minimum retention time, in seconds, to be used. Default is zero seconds.")
parser$add_argument("-R","--maxrt", action="store", type="double", dest="maxrt", help="The maximum retention time, in seconds, to be used. If this is not supplied, retention times from minrt to the end of the experiment will be considered.")
parser$add_argument("-f","--file", action="store", type="character", dest="file", default=NULL, help="The mzXML/mzML file to be analyzed. Overrides \"-dir\", which uses all files in the given directory.")
parser$add_argument("-d","--dir", action="store", type="character", dest="dir", default=getwd(), help="The directory with mzXML/mzML files to be quantified. If \"-file\" is supplied, that single file will be used instead. Default is the working directoy. Files should have the correct file extension.")
parser$add_argument("-o","--outfile", action="store", type="character", dest="outfile", default=NULL, help="The file that receives output data. Will auto-generate a name by default")
parser$add_argument("-x", "--format", action="store", type="character", dest="format", default="mzML", help="The file format to be used. Can be mzXML or mzML, mzML by default.")
parser$add_argument("-c", "--threads", action="store", type="integer", dest="threads", default=1, help="The number of threads to run in parallel. Default is 1 thread.")

arguments <- parser$parse_args()
###############################################################################



### PROCESS OPTIONS ###########################################################
file <- arguments$file
dir <- arguments$dir
mz <- arguments$mz
delta <- arguments$delta
minrt <- arguments$minrt
maxrt <- arguments$maxrt
outfile <- arguments$outfile
format<- arguments$format
threads <- arguments$threads

if(threads < 1){
  warning("Invalid number of threads! Defaulting to one.")
  threads <- 1
}

registerDoParallel(threads)
parparam <- SerialParam()
register(parparam, default=TRUE)


infiles <- 0
outfile_mz <- NULL
outfile_maxrt <- NULL


if(format != "mzML" && format != "mzXML"){
  warning("Bad file type supplied! mzML will be used by default.")
  format <- "mzML"
} else {
  message(paste("Format: ", format, sep=""))
}
if(!is.null(file)){
  message(paste("File name supplied. Data will be read from ", file, "\n"), sep="")
  infiles <- c(file)
} else {
  message(paste("Reading all ", format ," files in directory: ", dir, "\n"))
  if(format == "mzML"){infiles <- list.files(dir, pattern=".mzML")}
  if(format == "mzXML"){infiles <- list.files(dir, pattern=".mzXML")}
  
  print(infiles)
  message("\n\n")
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
  outfile <- paste("",str_replace_all(str_replace_all(Sys.time()," ","_"), ":", "-"),"_mz-", outfile_mz, "_rt-", minrt, "-", outfile_maxrt, ".csv", sep="")
}


message("Beginning analysis. This may take some time.")
###############################################################################




### MAIN PROCESSING STEP ######################################################
results <- bplapply(X=infiles, FUN=mainfunc, BPPARAM=parparam, mymz=mz, mydelt=delta, mymin=minrt, mymax=maxrt, mytrap=trapezoid, mybackground=background)

data <- c("file,feature_area")
for(i in 1:length(results)){
  data <- c(data, paste(results[[i]][1], results[[i]][2], sep=","))
}
###############################################################################

message(paste("file: ", outfile))
print(data)
write.csv(data, file=outfile, row.names=F)




message("Done.\n")
