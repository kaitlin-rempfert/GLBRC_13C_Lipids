# load required packages
library(tidyverse)
library(ggplot2)
library(xcms)

##batch 4 is all negative mode files - create raw data folder with just filenames that contain "Neg"
## Load in raw mzXML files
#path to .mzXML files
dda_folder_b4 = "/people/remp858/EC/raw_data/batch4"

#list all files inside specified folder
mzXMLfiles_b4 <- list.files(dda_folder_b4, recursive = TRUE, full.names = TRUE)

## Create a phenodata data.frame of filenames in order from folder
pd <- data.frame(extension = mzXMLfiles_b4, sample_name = basename(mzXMLfiles_b4))

#import metadata file 
meta <- readxl::read_excel("incubation_metadata_w_qc.xlsx")

#filter just positive
meta_b4 <- meta %>% filter(batch == "4") %>% filter(category != "exclude")
#sort by injection order
meta_b4_ordered <- meta_b4 %>% arrange(injection_order)
#merge with phenodata frame 
pd_b4 <- left_join(meta_b4_ordered, pd, by = "sample_name")
#export ordered file list
samps <- as.vector(pd_b4$extension)

#load data with phenodata
dda_data_b4 <- readMSData(samps, pdata = new("NAnnotatedDataFrame", pd_b4), centroided = TRUE, mode = "onDisk")

#filter out empty spectra 
dda_data_b4_filt <- dda_data_b4 %>% filterEmptySpectra()

## Total Ion Chromatograms
#Plot up total ion chromatograms
#get total ion chromatograms from the raw, centroided data
tis_b4 <- chromatogram(dda_data_b4_filt, aggregationFun = "sum")
#assign colors based on label
group_colors <- paste0(RColorBrewer::brewer.pal(12, "Paired"))
names(group_colors)<- c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "qc", "blank", "BTLE")
#plot total ion chromatograms
#pos
png(file = "tic_b4.png")
p1 <- plot(tis_b4,  col = group_colors[dda_data_b4_filt$group], main = "Batch 1 Total Ion Chromatogram")
legend("topright", legend = c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "QC", "Blank", "BTLE"), fill = group_colors) 
print(p1)
dev.off()

## Base Peak Chromatograms
#get base peak chromatograms from the raw, centroided data
bpis_b4 <- chromatogram(dda_data_b4_filt, aggregationFun = "max")
#plot base peak chromatograms
#pos
png(file = "bpi_b4.png")
p3 <- plot(bpis_b4, col = group_colors[dda_data_b4_filt$group], main = "Batch 1 Base Peak Chromatogram")
legend("topright", legend = c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "QC", "Blank", "BTLE"), fill = group_colors) 
print(p3)
dev.off()

## Total Ion Current
#Plot up total ion current by file
#extract total ion current for each file
tc_b4 <- split(tic(dda_data_b4_filt), f = fromFile(dda_data_b4_filt))
#plot
#pos
png(file = "tc_b4.png")
p4 <-  boxplot(tc_b4, ylab = "intensity", col = group_colors[dda_data_b4_filt$group], main = "Batch 1 Total Ion Current")
legend("topright", legend = c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "QC", "Blank", "BTLE"), fill = group_colors) 
print(p4)
dev.off()

#peak-picking parameters
cwp <- CentWaveParam(A = 4.289723e-07, ppm=1, Instrument=2, peakwidth=c(2.5, 45), snthresh = 10, noise=1000, prefilter=c(3, 10000), firstBaselineCheck = TRUE, integrate=2)

#pick peaks
dda_data_peaks_b4 <- findChromPeaks(dda_data_b4_filt, param = cwp)
saveRDS(dda_data_peaks_b4, "dda_data_stds_peaks_b4.rds")

#read in parameters
resultRetcorGroup <- readRDS("/people/remp858/EC/data/IPO_opt_2021_09_17/grouping_optimized.rds")

density.bw = resultRetcorGroup$best_settings$bw
density.max = resultRetcorGroup$best_settings$max
density.minfrac = resultRetcorGroup$best_settings$minfrac
density.minsamp = resultRetcorGroup$best_settings$minsamp
density.mzwid = resultRetcorGroup$best_settings$mzwid
density.sleep = 0 

loess.missing = resultRetcorGroup$best_settings$missing
loess.extra = resultRetcorGroup$best_settings$extra
loess.family = resultRetcorGroup$best_settings$family

#peak grouping parameters
pdp <- PeakDensityParam(
  sampleGroups = dda_data_peaks_b4@phenoData$group,
  bw = density.bw,
  #will need to experiment with this
  minFraction = 0.1,
  minSamples = density.minsamp,
  binSize = density.mzwid,
  maxFeatures = density.max)

pgp <- PeakGroupsParam(
  minFraction = .95,
  extraPeaks = 2,
  smooth = "loess",
  span = 0.4,
  subset = which(dda_data_peaks_b4@phenoData$group== "qc"),
  subsetAdjust = "average",
  family = "gaussian")

#initial grouping
dda_data_b4_grouped <- groupChromPeaks(dda_data_peaks_b4, param = pdp)
saveRDS(dda_data_b4_grouped, "dda_data_b4_grouped.rds")

#loess retention time correction 
dda_data_b4_rt <-adjustRtime(dda_data_b4_grouped, param = pgp)
saveRDS(dda_data_b4_rt, "dda_data_b4_rt.rds")

#regroup retention time corrected data 
dda_data_b4_rt_grouped  <- groupChromPeaks(dda_data_b4_rt, param = pdp)
saveRDS(dda_data_b4_rt_grouped, "dda_data_b4_rt_grouped.rds")

#Visualize Grouping and Retention Time Adjustments
#visualize alignment 
#pos
png(file = "rt_b4.png")
p5 <- plotAdjustedRtime(dda_data_b4_rt_grouped, col = group_colors[dda_data_b4_rt_grouped$group],
                        peakGroupsCol = "grey", peakGroupsPch = 1)
legend("topright", legend = c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "QC", "Blank", "BTLE"), fill = group_colors) 
print(p5)
dev.off()

#QC on picked peaks- is intensity normalized?
#pos
png(file = "qc_b4.png")
p6 <- boxplot(featureValues(dda_data_b4_rt_grouped, value="into") +1, 
              log="y", las=2, col = group_colors[dda_data_b4_rt_grouped$group], ylab = "intensity", main = "Peak Intensities")
legend("topright", legend = c("MC_2m", "MC_1y", "MS_2m", "MS_1y", "WC_2m", "WC_1y", "WS_2m", "WS_1y", "QC", "Blank", "BTLE"), fill = group_colors) 
print(p6)
dev.off()

## Fill Peaks
#filter to just MS1
dda_data_b4_rt_grouped_ms1 <- dda_data_b4_rt_grouped %>% filterMsLevel(msLevel. = 1)

#fill peaks
fpp <- FillChromPeaksParam(ppm = 5)
dda_data_b4_fill= fillChromPeaks(dda_data_b4_rt_grouped_ms1, param = fpp)

#save filled peaks
saveRDS(dda_data_b4_fill, "dda_data_b4_fill.rds")
