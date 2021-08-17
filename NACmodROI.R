# Script to trigger activity detection from calcium imaging recordings
# This script contains functions to detect activity adopted from Prada, J. et al. (2018). An open source tool for automatic spatiotemporal assessment of calcium transients and local ‘signal-close-to-noise’ activity in calcium imaging data. PLoS Comput Biol, 14(3), e1006054. (Neural-Activity-Cubic)
# The collection of scripts does not provide a complete API for function execution, the implementation of the used functions is planned for a new release of Neural-Activity-Cubic.

# Below directories and specific settings for detection can be specified. When settings are completed run (source) the file to start detection

# directories
save_dir <- 'results'
avi_dir <- 'csv' # input directory, for avi AND csv files
roi_dir <- 'roi'
setwd() # set working directory

# settings
input_data <- "csv" # format of input data to analyse, "avi" for avi videos, "csv" for pre-calculated roi data in RData format in one directory
measure_these_days <- list("20160218")
frames_to_analyse <- c(0,1200) # implemented for mode "csv" only! vector specifying start frame and end frame of analysis, if "c(NA,NA)" the all frames will be analysed // in perfusion: c(1200,7200) [perfusion phase], c(0,1200) [baseline phase], perfusion with safety limits c(1500,6900) // in single cell measurements c(0,3000)
SAT <- 3.1 # Signal average threshold (1.3 for single cell, 3.1 for 10x magnification)
SNR_vec <- c(1.75) # Signal to noise ratio (2.25 for single cell, 2.25 for 10x magnification)
WS <- 8 # Window size
MAC <- 1 # Minimum activity count (not implemented)
VarianceEnable <- TRUE # Measure Variance?
VarWindow <- 30 # Size of variance window
session <- 12 # Evaluation Session
for (SNR in SNR_vec) {
  source("NACmodROI_functions.R") # source function .R file
  if(input_data=="avi") {
    source("NACmodROI_loop.R") # starts detection loop
  } else {
    source("NACmodROI_csv.R") # starts detection from csv files
  }
}
