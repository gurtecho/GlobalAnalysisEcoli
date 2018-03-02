
#label peaks as they

setwd("Lab Work/Kosuri Lab/GlobalAnalysisEcoli/scripts/")

library("ggplot2")
library("dplyr")
library("sqldf")

options(stringsAsFactors = T)


plus_peaks <- read.table("../../endo/processed_data/custom_peaks/plus_called_peaks_threshold1.3_merge3_min20.bed")
minus_peaks <- read.table("../../endo/processed_data/custom_peaks/plus_called_peaks_threshold1.3_merge3_min20.bed")