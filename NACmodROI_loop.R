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

for(day in measure_these_days) {
  day_path <- file.path(roi_dir,day)
  if(file.exists(day_path)) {
    files_roi_avi <- get_roi_avi_files(day,avi_dir,roi_dir)
    cat(length(files_roi_avi[[1]]),"pair of matching ROI and AVI files found for measuring day",day,"\n")
    for(i in 1:length(files_roi_avi[[1]])) {
      run <- files_roi_avi[[1]][i]
      cat("-- Analysing run",run,"\n")
      buf <- get_array_from_avi(files_roi_avi[[3]][i])
      cat("---- AVI file loaded \n")
      roi_coords <- get_roi_coords(files_roi_avi[[2]][i])
      cat("---- ROI coordinates loaded \n")
      x_dim <- dim(buf)[1]
      y_dim <- dim(buf)[2]
      #first_frame <- save_first_frame_annotated(buf[,,1,1],roi_coords,save_dir,run,day)
      pointsInRoi <- get_points_in_roi(roi_coords,x_dim,y_dim)
      roiMidpointsArea <- get_roi_midpoints(pointsInRoi,WS)
      roiMidpoints <- roiMidpointsArea[[1]]
      roiArea <- roiMidpointsArea[[2]]
      cat("---- Reading roi means from frames: \n")
      roiMeans <- get_means_per_roi(buf,pointsInRoi,roiMidpoints,WS)
      cat("---- Means per roi calculated \n")
      rm(buf)
      detectActivity(roiMeans,roiMidpoints,roiArea,WS,SNR,SAT,MAC,VarianceEnable,VarWindow,day,run,session,x_dim,y_dim,save_dir)
      cat("---- Activity detection completed \n")
    }
  }else {
    cat("No ROI directory found for measuring day",day,"\n")
  }
}