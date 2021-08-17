# libraries
library(stringr)
library(imager)
library(RImageJROI)
library(utils)
require(TTR)
require(IM)
require(tiff)
require(wmtsa)
require(changepoint)
require(cpm)
require(pracma)
require(sp)
require(graphics)
require(MASS)

for(meansfile in list.files(avi_dir)) {
    day_run <- sub(".csv","",meansfile)
    day <- strsplit(day_run,"_")[[1]][1]
    run <- strsplit(day_run,"_")[[1]][2]
    cat("-- Analysing ",day_run,"\n")
    x_dim <- 348
    y_dim <- 260
    roiMidpoints <- NA
    roiArea <- NA
    cat("---- Reading roi means from file \n")
    roiMeans <- read.csv(file=paste0(avi_dir,"/",meansfile))[-1]
    roiMeans <- roiMeans[frames_to_analyse[1]:frames_to_analyse[2],]
    cat("---- ROI means loaded \n")
    detectActivity(roiMeans,roiMidpoints,roiArea,WS,SNR,SAT,MAC,VarianceEnable,VarWindow,day,run,session,x_dim,y_dim,save_dir)
    cat("---- Activity detection completed \n")
}