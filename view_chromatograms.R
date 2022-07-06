#!/usr/bin/R
# Ryan Preble, Spring 2022
# Boyce Thompson Institute, 533 Tower Road, Ithaca, NY 14853
# Make adjustments as needed for it to work on your system and samples. 

setwd("D:/Mass_Spec/Pawpaw/20220617")

# CHECK INSTALLED PACKAGES. INSTALL IF NECESSARY
packages <- c("mzR", "Spectra", "MSnbase","ggplot2")
list.of.packages <- c("mzR", "Spectra", "MSnbase", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# REQUIRED LIBS
library(mzR)
library(Spectra)
library(MSnbase)
library(grDevices)
library(ggplot2)

# REPLACE THE FOLLOWING LINE WITH MZML/MZXML FILES FOR YOUR SAMPLES
files<- c("1-0.mzML", "1-24.mzML", "2-0.mzML","2-24.mzML","3-0.mzML","3-24.mzML")
#c("1-0.mzML", "1-24.mzML", "2-0.mzML","2-24.mzML","3-0.mzML","3-24.mzML")

### READ IN FILE
mzXML <- MSnbase::readMSData(files, mode = "onDisk")

# CHANGE THE MASS TO YOUR DESIRED MZ
massnum <- 597.4725

### ISOLATE DATA (REPLACE THE MASS NUMBERS, DELTA, AND RETENTION TIME WITH YOUR DESIRED VALUES)
mzXML <- filterMz(mzXML, mz=c(massnum-0.005, massnum+0.005), msLevel=1)

mzXML <- filterRt(mzXML, c(600, 690)) #Change this line for retention time changes


### EXTRACT CHROMATOGRAM
chromatogram <- chromatogram(mzXML, aggregationFun="sum")

data <- chromatogram@.Data
rtime <- c()
intensity <- c()
filename <- c()

for(i in 1:length(data)){
  rtime <- c(rtime, as.data.frame(data[i][[1]])$rtime)
  intensity <- c(intensity, as.data.frame(data[i][[1]])$intensity)
  filename <- c(filename, rep(files[i], times=length(as.data.frame(data[i][[1]])$rtime)))
}

data <- data.frame(rtime, intensity, filename)

### ANALYZE CHROMATOGRAM CONTENT
#plot(chromatogram, col=rainbow(length(files)), main=paste("Intensity of ",massnum," mz Signal", sep=""), xlab="Retention Time (s)", ylab="Intensity")
plot <- ggplot(data = data, aes(x=rtime, y=intensity, color=filename, group=filename)) + 
  geom_line() + 
  ylim(c(0,max(data$intensity + 1000000))) + 
  labs(title=paste("Intensity of ",massnum," mz Signal", sep=""), x="Retention time (s)", y="Intensity", color="")

plot(plot)

message("Done.\n")