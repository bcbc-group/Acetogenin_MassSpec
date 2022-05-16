# Acetogenin_MassSpec
Scripts for analyzing mass spectrometry data from annonaceous plants.

standards.R - This file was used to create figures for the standard curve used to determine the acetogenin content of a plant. This may be adapted for similar purposes with some simple edits.
trapezoid.R - This file iteratively processes mzML and mzXML files and determines the feature area for a given mass number and a given retention time. 
trapezoid2.R - This file accomplishes the same thing as trapezoid.R, but it utilizes multicore parallel processing. It is more likely to run into synchrony issues, so if you encounter errors, try running your analysis using trapezoid.R. It will be slower, but is less likely to fail. 
