
#label peaks as they

setwd("Lab Work/Kosuri Lab/GlobalAnalysisEcoli/scripts/")

library("ggplot2")
library("dplyr")
library("sqldf")
install.packages('compare')
library("compare")
options(stringsAsFactors = F)


plus_peaks <- read.table("../../endo/processed_data/custom_peaks/plus_called_peaks_threshold1.3_merge3_min20.bed",
                         col.names = c("genome", "start", "end", "score", "strand"), header = F)
minus_peaks <- read.table("../../endo/processed_data/custom_peaks/minus_called_peaks_threshold1.3_merge3_min20.bed", 
                          col.names = c('genome', 'start', 'end', 'score', 'strand'), header = F)

plus_gene_annots <- read.table("../ref/U00096.2_gene_only.bed") %>% filter(V6 == "+") %>% mutate(upstream = V2-150, downstream = V3+ 150)
plus_gene_annots$V1 <- "U00096.2"

gene_annots <- read.table("../ref/U00096.2_gene_only.bed")
gene_annots$V1 <- "U00096.2"
names(gene_annots) <- NULL
write.table(gene_annots, "../ref/U00096.2_gene_only.bed", quote = F, row.names = F, sep = "\t")

#Get all plus regions

plus_upstream <- plus_gene_annots %>% select('V1', 'upstream', 'V2', 'V4', 'V5', 'V6')
names(plus_upstream) <- NULL
write.table(plus_upstream, "../processed_data/PromoterGeneLandscape/plus_upstream.bed", quote = F, row.names = F, sep = "\t")

plus_internal_downstream <- plus_gene_annots %>% select('V1', 'V2', 'downstream', 'V4', 'V5', 'V6')
names(plus_internal_downstream) <- NULL
write.table(plus_internal_downstream, "../processed_data/PromoterGeneLandscape/plus_internal_downstream.bed", quote = F, row.names = F, sep = "\t")


#USING THESE FILES, IDENTIFY GENES WITH UPSTREAM PROMOTERS AND INTERGENIC/DOWNSTREAM ANTISENSE PROMOTERS


#Identify genes with dual regulation (sense and antisense)
sense_regulated_plus <-read.table("../processed_data/PromoterGeneLandscape/upstream_sense_regulated__genes.bed")
antisense_regulated_plus <- read.table("../processed_data/PromoterGeneLandscape/antisense_regulated__genes.bed")


dual_regulatedGenes <- semi_join(sense_regulated_plus, antisense_regulated_plus, by = 'V4')
names(dual_regulatedGenes) <- NULL

write.table(dual_regulatedGenes, "../processed_data/PromoterGeneLandscape/dual_regulated_plus.bed", quote = F, row.names = F, sep = "\t")



minus_gene_annots <- read.table("../ref/U00096.2_gene_only.bed") %>% filter(V6 == "-") %>% mutate(upstream = V2-150, downstream = V3+ 150)
minus_gene_annots$V1 <- "U00096.2"


#Find plus promoters that have 5' promoter

sqldf("SELECT a.start, a.end, a.score, b.upstream, b.V2, b.V3, b.downstream
      FROM plus_peaks AS a 
      LEFT JOIN plus_gene_annots AS b 
      ON b.start BETWEEN a.peak_start AND a.peak_end
      OR b.end BETWEEN a.peak_start AND a.peak_end




