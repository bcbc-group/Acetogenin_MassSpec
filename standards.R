#!/usr/bin/R
# Ryan Preble, Spring 2022
# Boyce Thompson Institute, 533 Tower Road, Ithaca, NY 14853
# Make adjustments as needed for it to work on your system and samples. 

setwd("D:/Mass_Spec/Pawpaw/20211029")

# CHECK INSTALLED PACKAGES. INSTALL IF NECESSARY
packages <- c("mzR", "Spectra", "MSnbase")
list.of.packages <- c("mzR", "Spectra", "MSnbase")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# REQUIRED LIBS
library(mzR)
library(Spectra)
library(MSnbase)
library(grDevices)

# REPLACE THE FOLLOWING LINE WITH MZML/MZXML FILES FOR YOUR STANDARDS
files<- c("20211029_standard_1.mzML", "20211029_standard_2.mzML", "20211029_standard_3.mzML", "20211029_standard_4.mzML", "20211029_standard_5.mzML", "20211029_standard_6.mzML", "20211029_standard_7.mzML")

### READ IN FILE
mzXML <- MSnbase::readMSData(files, mode = "onDisk")

### ISOLATE DATA (REPLACE THE MASS NUMBERS, DELTA, AND RETENTION TIME WITH YOUR DESIRED VALUES)
mzXML <- filterMz(mzXML, mz=c(597.4725-0.005, 597.4725+0.005), msLevel=1)

mzXML <- filterRt(mzXML, c(600, 725))


### EXTRACT CHROMATOGRAM
chromatogram <- chromatogram(mzXML, aggregationFun="sum")

### ANALYZE CHROMATOGRAM CONTENT
plot(chromatogram, col=rainbow(length(files)), main="Intensity of 597.4725 mz Signal in Standard Curves", xlab="Retention Time (s)", ylab="Intensity")

message("Done.\n")