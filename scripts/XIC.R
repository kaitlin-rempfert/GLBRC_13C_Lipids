#function for plotting XIC  
raw_XIC_plot <- function(targets, raw) {
  group_colors <- paste0(RColorBrewer::brewer.pal(10, "Spectral"))
  names(group_colors) <- raw@phenoData$catalog_number
  
  for (row in 1:nrow(targets)){
    XIC <- xcms::chromatogram(raw, mz = c(targets$mz_db[row]-.008, targets$mz_db[row]+.008), rt = c(targets$rt_db[row]*60 - .5*targets$rt_window[row], targets$rt_db[row]*60 + .5*targets$rt_window[row])) 
    plot(XIC, col = group_colors[raw$catalog_number], main = targets$compound[row])
  }
}

#function for extracting apex RT, apex intensity, and peak width of XICs from target list
extract_peak_param <- function(targets, raw) {
  output_total = data.frame()
  for (row in 1:nrow(targets)){
    XIC <- xcms::chromatogram(raw, mz = c(targets$mz_db[row]-.008, targets$mz_db[row]+.008), rt = c(targets$rt_db[row]*60 - .5*targets$rt_window[row], targets$rt_db[row]*60 + .5*targets$rt_window[row] )) 
    output <- matrix(ncol=6, nrow=length(XIC))
    output[,1] <- targets$compound[row]
    output[,2] <- raw$catalog_number
    output[,3] <- raw$amendment
    for (sample in 1:length(XIC)){
      samp <- XIC[1, sample]
      samp <- xcms::chromPeaks(samp)
      samp <- as.data.frame(samp)
      samp <- slice(samp, which.max(into))
      output[sample, 4] <- samp$rt[1]
      output[sample, 5] <- samp$maxo[1]
      output[sample, 6] <- samp$rtmax[1] - samp$rtmin[1]
    }
    output <- as.data.frame(output)
    output_total <- rbind(output_total, output)
  }
  colnames(output_total) <- c('compound','catalog_number','amendment', 'RT', 'Intensity', 'RT_width')
  output_total$RT <- as.numeric(as.character(output_total$RT))/60
  return(output_total)
}

#function for extracting apex RT, apex intensity, and peak width of XICs from target list
extract_peak_area <- function(targets, raw) {
  output_total = data.frame()
  for (row in 1:nrow(targets)){
    XIC <- xcms::chromatogram(raw, mz = c(targets$mz_db[row]-.008, targets$mz_db[row]+.008), rt = c(targets$rt_db[row]*60 - .5*targets$rt_window[row], targets$rt_db[row]*60 + .5*targets$rt_window[row] )) 
    output <- matrix(ncol=5, nrow=length(XIC))
    output[,1] <- targets$compound[row]
    output[,2] <- raw$sample_name
    for (sample in 1:length(XIC)){
      samp <- XIC[1, sample]
      samp <- xcms::chromPeaks(samp)
      samp <- as.data.frame(samp)
      samp <- slice(samp, which.max(into))
      output[sample, 3] <- samp$mz[1]
      output[sample, 4] <- samp$rt[1]
      output[sample, 5] <- samp$into[1]
    }
    output <- as.data.frame(output)
    output_total <- rbind(output_total, output)
  }
  colnames(output_total) <- c('Compound','sample_name', 'mz', 'rt_min', 'area')
  output_total$rt_min <- as.numeric(as.character(output_total$rt_min))/60
  return(output_total)
}


#function for extracting apex RT, apex intensity, and peak width of XICs from target list
extract_peak_apex_h <- function(targets, raw) {
  output_total = data.frame()
  for (row in 1:nrow(targets)){
    XIC <- xcms::chromatogram(raw, mz = c(targets$mz_db[row]-.008, targets$mz_db[row]+.008), rt = c(targets$rt_db[row]*60 - .5*targets$rt_window[row], targets$rt_db[row]*60 + .5*targets$rt_window[row] )) 
    output <- matrix(ncol=5, nrow=length(XIC))
    output[,1] <- targets$compound[row]
    output[,2] <- raw$sample_name
    for (sample in 1:length(XIC)){
      samp <- XIC[1, sample]
      samp <- xcms::chromPeaks(samp)
      samp <- as.data.frame(samp)
      samp <- slice(samp, which.max(maxo))
      output[sample, 3] <- samp$mz[1]
      output[sample, 4] <- samp$rt[1]
      output[sample, 5] <- samp$maxo[1]
    }
    output <- as.data.frame(output)
    output_total <- rbind(output_total, output)
  }
  colnames(output_total) <- c('Compound','sample_name', 'mz', 'rt_min', 'area')
  output_total$rt_min <- as.numeric(as.character(output_total$rt_min))/60
  return(output_total)
}
