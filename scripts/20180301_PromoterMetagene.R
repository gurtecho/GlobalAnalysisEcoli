setwd("~/Documents/GlobalAnalysisEcoli/scripts/")

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
library("metagene")

#bioclite("metagene")

options(stringsAsFactors = F)
        
#The path (relative or absolute) to the BAM files must be in a vector:
          
bam_files <- list.files("../processed_data/", pattern = "*.bam", full.names = T)#list of file names

regions <- list.files("../ref/", pattern = "U00096.2_gene_only.bed", full.names = T)

bed <- read.table("../ref/U00096.2_gene_only.bed",sep = "\t", header = F)
names(bed) <- NULL
colnames(bed)<-c("chr", "start", "end", "x")
gr1 <- toGRanges(bed, format="BED", header=F)





bam_files <- c(system.file("~/Documents/GlobalAnalysisEcoli/processed_data/Plus_DNA.bam", package = "metagene"),
#                system.file("~/Documents/GlobalAnalysisEcoli/processed_data/Minus_DNA.bam", package = "metagene"),
#                system.file("~/Documents/GlobalAnalysisEcoli/processed_data/Plus_RNA.bam", package = "metagene"),
#                system.file("~/Documents/GlobalAnalysisEcoli/processed_data/Minus_RNA.bam", package = "metagene"))
# regions <- c(system.file("../ref/U00096.2_gene_only.bed", package="metagene"))      


#Creation of a metagene object for RNA-seq analysis

mg <- metagene$new(regions = regions, 
                   bam_files = bam_files, 
                   assay = 'rnaseq', 
                   paired_end = FALSE)
seqlevelsStyle("../processed_data/Minus_DNA.bam")
seqlevels("../processed_data/Minus_DNA.bam")
