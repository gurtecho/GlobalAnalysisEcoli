setwd("~/Documents/GlobalAnalysisEcoli/scripts/")
setwd("~/Lab Work/Kosuri Lab/GlobalAnalysisEcoli/scripts/")

library("ggplot2")
library("dplyr")
library("wesanderson")
names(wes_palettes)
library("tidyr")
library("ggjoy")r.
library("reshape2")
require(cowplot)
library('RInside')
library("stringr")
library("plotly")
library(GenomicRanges)
#source("https://bioconductor.org/biocLite.R")
#biocLite("metagene")
library("metagene")

#biocLite("metagene")

options(stringsAsFactors = T)
        
#The path (relative or absolute) to the BAM files must be in a vector:
          
bam_files <- list.files("../processed_data/", pattern = "*.bam", full.names = T)#list of file names
bam_files <- c("../processed_data/Plus_RNA.bam")
seqlevels(bam_files)


regions <- list.files("../ref/", pattern = "U00096.2_gene_only.bed", full.names = T)
regions <- c("../ref/U")


names(bed) <- NULL
colnames(bed)<-c("chr", "start", "end", "x")
gr1 <- toGRanges(bed, format="BED", header=F)


#Convert bed into GRange

bed <- read.table("../ref/U00096.2_gene_only.bed",sep = "\t", header = F)
bed$V1 <- "U00096.2"

temp <- 
  lapply(split(bed, bed$V4), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$V4))
  })


#bam_files <- c(system.file("../processed_data/Plus_DNA.bam", package = "metagene"),
#                system.file("../processed_data/Minus_DNA.bam", package = "metagene"),
#                system.file("../processed_data/Plus_RNA.bam", package = "metagene"),
#                system.file("../processed_data/Minus_RNA.bam", package = "metagene"))

#regions <- c(system.file("../ref/U00096.2_gene_only.bed", package="metagene"))      


#Creation of a metagene object for RNA-seq analysis

mg <- metagene$new(regions = temp, 
                   bam_files = bam_files, 
                   assay = 'rnaseq', 
                   paired_end = FALSE)


seqlevelsStyle("../processed_data/Minus_DNA.bam")
seqlevels("../processed_data/Minus_DNA.bam")
