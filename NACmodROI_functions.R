# this file includes all functions used in NACmodROI.R

# function to find corresponding roi and avi files for a measurement day
get_roi_avi_files <- function(day, # measurement day identifier
                              avi_dir, # directory of avi files
                              roi_dir # directory of roi files
) {
  # create a vector of the numbers of the roi files found in roi directory
  roi_no <- c()
  measurement_roi_path <- file.path(roi_dir,day)
  for(roi in list.files(measurement_roi_path)) {
    roi <- str_remove(str_remove(roi,pattern=".zip"),pattern = paste0(day,"_"))
    roi_no <- c(roi_no,roi)
  }
  numbers <- c()
  files_roi <- c()
  files_avi <- c()
  # find all avi files in avi directory of measurement day
  measurement_avi_path <- file.path(avi_dir,day)
  for(avi_file in list.files(measurement_avi_path,pattern = ".avi")) {
    avi_no <- strsplit(avi_file,"_")[[1]][2]
    # check if avi number has a corresponding roi number and then add to list of paired files
    if(avi_no %in% roi_no) {
      avi_path <- file.path(measurement_avi_path,avi_file)
      roi_path <- file.path(measurement_roi_path,paste0(day,"_",avi_no,".zip"))
      numbers <- c(numbers, avi_no)
      files_roi <- c(files_roi, roi_path)
      files_avi <- c(files_avi,avi_path)
    }
  }
  return(list(numbers, # run number
              files_roi, # path to roi file
              files_avi # path to avi file
  ))
}

# function to find corresponding roi and csv files for a measurement day, copied from get_roi_avi_files therefore "avi" in variable names
get_roi_csv_files <- function(day, # measurement day identifier
                              avi_dir, # directory of csv files
                              roi_dir # directory of roi files
) {
  # create a vector of the numbers of the roi files found in roi directory
  roi_no <- c()
  measurement_roi_path <- file.path(roi_dir,day)
  for(roi in list.files(measurement_roi_path)) {
    roi <- str_remove(str_remove(roi,pattern=".zip"),pattern = paste0(day,"_"))
    roi_no <- c(roi_no,roi)
  }
  numbers <- c()
  files_roi <- c()
  files_avi <- c()
  # find all avi files in avi directory of measurement day
  measurement_avi_path <- file.path(avi_dir,day)
  for(avi_file in list.files(measurement_avi_path,pattern = ".csv")) {
    avi_no <- str_remove(strsplit(avi_file,"_")[[1]][2],".csv")
    # check if avi number has a corresponding roi number and then add to list of paired files
    if(avi_no %in% roi_no) {
      avi_path <- file.path(measurement_avi_path,avi_file)
      roi_path <- file.path(measurement_roi_path,paste0(day,"_",avi_no,".zip"))
      numbers <- c(numbers, avi_no)
      files_roi <- c(files_roi, roi_path)
      files_avi <- c(files_avi,avi_path)
    }
  }
  return(list(numbers, # run number
              files_roi, # path to roi file
              files_avi # path to avi file
  ))
}

# function to read avi file into R
get_array_from_avi <- function(avi_file # path to single avi file
) {
  # avi file is loaded with imager package and multiplied by 255 to get 8 bit  values
  vid <- load.video(fname=avi_file, maxSize = 5)
  vid <- vid*255
  return(vid)
}

# function to read roi zip files into R and return x,y coordinates of roi shape
get_roi_coords <- function(roi_file # path to single roi zip file
) {
  # roi zip file read using ImageJROI package
  roi_zip <- read.ijzip(roi_file)
  roi_x <- list()
  roi_y <- list()
  for(j in 1:length(roi_zip)) {
    roi_x[[j]] <- roi_zip[[j]]$coords[,1]
    roi_y[[j]] <- roi_zip[[j]]$coords[,2]
  }
  return(list(roi_x,roi_y))
}

# function to detect all pixels enclosed in the rois
get_points_in_roi <- function(roi_coords, # x,y coordinates of all rois in roi zip file
                              x_dim, # width of analysed video
                              y_dim # height of analysed video
) {
  pointsInRoi <- vector(mode="list",length = length((roi_coords[[1]])))
  for(x in 1:x_dim)
    for(y in 1:y_dim)
      for(r in 1:length(roi_coords[[1]]))
        if(point.in.polygon(x,y,roi_coords[[1]][[r]],roi_coords[[2]][[r]])) {
          pointsInRoi[[r]] <- rbind(pointsInRoi[[r]],c(x,y))
        }
  return(pointsInRoi)
}

# function to get midpoints and area of rois (at the moment midpoints not needed anymore)
get_roi_midpoints <- function(pointsInRoi,WS) {
  roiMidpoints = matrix(0,length(pointsInRoi),2)
  roiArea <- c()
  for(r in 1:length(pointsInRoi)) {
    roiMidpoints[r,1] <- mean(pointsInRoi[[r]][,1]) #ceiling(mean(pointsInRoi[[r]][,1])/WS)
    roiMidpoints[r,2] <- mean(pointsInRoi[[r]][,2]) #ceiling(mean(pointsInRoi[[r]][,2])/WS)
    roiArea[r] <- nrow(pointsInRoi[[r]])
  }
  roiMidpoints <- ceiling(roiMidpoints/WS)
  return(list(roiMidpoints,roiArea))
}

# function to get means per frame for each roi
get_means_per_roi <- function(buf, # array of avi video information
                              pointsInRoi, # list containing lists with all points in particular rois
                              roiMidpoints, # midpoint of each roi (no needed)
                              WS # window size of singnal detection
                              ) {
  z_dim <- dim(buf)[3]
  roiMeans <- data.frame(matrix(ncol=length(pointsInRoi),nrow=z_dim))
  for(frame_count in 1:z_dim) {
    this_frame <- buf[,,frame_count,1]
    for(l in 1:length(pointsInRoi)) {
      roiMeans[frame_count,l] <- mean(this_frame[cbind(pointsInRoi[[l]][,1],pointsInRoi[[l]][,2])])
    }
    rm(this_frame)
  }
  return(roiMeans)
}

# funtion to detect activity and variance area in traces of calcium imaging
# the function is compied from Neural Activity cubic (Prada et al, 2018)
detectActivity <- function(roiMeans, # dataframe of means per frame per roi
                           roiMidpoints, # midpoints of rois (not needed)
                           roiArea, # roi area
                           WS, # window size
                           SNR, # signal to noise ratio
                           SAT, # signal average thresholf
                           MAC, # minimum activity count (not implemented)
                           VarianceEnable, # boolean whether variance is calculated
                           VarWindow, # window of variance area calculation
                           rDay, # day of measurements in format YYYYMMDD
                           rRun, # run of measurement of specified day in format XX
                           session, # session of detection
                           x_dim, # width of video
                           y_dim, # height of video
                           save_dir # directory to save csv files
                           ) {
  MNT <- SAT
  tamano <- WS
  WinSize <- VarWindow
  midpoints <- roiMidpoints
  minPeaks <- MAC
  
  numberRois = ncol(roiMeans)
  
  FN = paste(toString(rDay),"_",toString(str_pad(rRun,width=2,pad="0")), sep="")
  SN = toString(str_pad(session,width= 2,pad = "0"))
  
  peak_list=list(list())
  wavtrans_list=list(list())
  wavtree_list=list(list())
  
  for (j in 1:ncol(roiMeans)) {
    series=roiMeans[,j]
    peak_list[[j]]=0
    
    if(mean(series)>MNT){
      
      wav_trans = wavCWT(series)
      wavtrans_list[[j]]=wav_trans
      er1 <- tryCatch(wavCWTTree(wav_trans, n.octave.min=1, tolerance=0.0, type="maxima"),error=function(e) e)
      
      if(inherits(er1, "error")){
        wavtree_list[[j]]=0
        next}
      
      wav_tree = wavCWTTree(wav_trans, n.octave.min=1, tolerance=0.0, type="maxima") 
      wavtree_list[[j]]=wav_tree
      er2 <- tryCatch(wavCWTPeaks(wav_tree, snr.min=SNR, scale.range=NULL, length.min=5, noise.span=NULL, noise.fun="quantile", noise.min=NULL),error=function(e2) e2)
      
      if(!inherits(er2, "error")){
        
        peaks = wavCWTPeaks(wav_tree, snr.min=SNR, scale.range=NULL, length.min=5, noise.span=NULL, noise.fun="quantile", noise.min=NULL)
        peak_list[[j]]=peaks
        
      }
    }
  }
  
  ## Output: CSV file
  
  
  #### Data frame creation for individual recording
  
  dataRoi <- data.frame("ROI" = rep(0,numberRois),"Activity Count"= rep(0,numberRois),"Variance Area"= rep(0,numberRois))
  
  dataAll <- data.frame("Measurement" = rep(0,numberRois),"ROI" = rep(0,numberRois),"ImageJ identifier" = rep(0,numberRois),"Area"=rep(0,numberRois),"Activity Count"= rep(0,numberRois),"Variance Area"= rep(0,numberRois),"WS"=rep(0,numberRois),"SNR"=rep(0,numberRois),"SAT"=rep(0,numberRois),"MAC"=rep(0,numberRois))
  
  
  for (p in 1:ncol(roiMeans)) {
    # activity count
    if (is.list(peak_list[[p]][1])) roiCount <- length(peak_list[[p]][[1]])
    else roiCount <- 0
    
    # variance
    trace = roiMeans[,p]
    if(VarianceEnable){   					
      meanSig = rollapply(trace, WinSize, mean, fill = "extend", align = "right")
      varSig = rollapply(trace, WinSize, var, fill = "extend", align = "right")
      varSig1  = varSig+meanSig
      varSig2 = -varSig+meanSig
      areaSig = round(abs(polyarea(c(1:length(varSig1),length(varSig1):1),c(varSig1,varSig2[length(varSig2):1]))), digits=3)
      
    } else {
      areaSig = "NA"
    }
    
    #dataRoi[p, ] = c(p, roiCount, areaSig)
    dataAll[p, ] = c(FN, p, colnames(roiMeans)[p], roiArea[p], roiCount, areaSig, tamano,SNR, MNT, minPeaks)
  }
  
  # save csv file with counts for individual recording
  #csvROI = paste0(save_dir,"/csv/results_",toString(FN),"_WS",toString(tamano),"_SNR",toString(SNR),"_SAT",toString(MNT),"_MAC",toString(minPeaks),".csv")
  #write.csv(dataRoi, file=csvROI, row.names=F, na="")
  
  # save csv file with roi means for individual recordings
  csvMeans <- paste0(save_dir,"/roidata/means_",toString(FN),"_WS",toString(tamano),"_SNR",toString(SNR),"_SAT",toString(MNT),"_MAC",toString(minPeaks),".csv")
  write.csv(roiMeans,file=csvMeans,row.names = F,na="")
  
  # save csv file with peaks for individual recordings
  datPeaks <- paste0(save_dir,"/roidata/peaks_",toString(FN),"_WS",toString(tamano),"_SNR",toString(SNR),"_SAT",toString(MNT),"_MAC",toString(minPeaks),".rda")
  save(peak_list,file=datPeaks)
  
  # add counts to session summary
  csvAll = paste(save_dir, "/session_",toString(SN),"_results.csv",sep="")
  
  write.table(dataAll, csvAll, col.names = !file.exists(csvAll), append = T, row.names=F)
  
}


# save_first_frame_annotated <- function(frame1, # two dimensional array of the first frame in 8bit grayscale
#                           roi_coords, # coordinates of rois
#                           save_dir, # directory to save image
#                           run, # run of measurement
#                           day # day of measurement
# ) {
#   melt_frame<- melt(frame1) # melt to dataframe
#   roi_coords_melt <- merge(data.frame(id=row.names(melt(roi_coords[[1]])),x=melt(roi_coords[[1]])[,1],roi=melt(roi_coords[[1]])[,2]),data.frame(id=row.names(melt(roi_coords[[2]])),y=melt(roi_coords[[2]])[,1],roi=melt(roi_coords[[2]])[,2]))
#   ggplot() + 
#     #geom_raster(data=melt_frame, mapping=aes(x = Var1, y = Var2, fill=value)) +
#     scale_fill_gradient(low="black", high="white") +
#     theme_bw() + 
#     geom_polygon(data=roi_coords_melt,mapping = aes(x=c(1,2,3),y=c(2,3,4)),group=1)
#   }
#   imager::save.image(frame1,paste0(save_dir,"/images/",day,"_",run,".jpg"),quality = 1)
# }
