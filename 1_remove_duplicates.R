# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# Load files for each genotyping chip
library(dplyr)
mega <- read.delim("MEGA_Chip.bim", header=F, quote="")
neur <- read.delim("Neuro_Chip.bim", header=F, quote="")

# Load call rate for all loci on each chip
cr.mega <- read.table("mega.lmiss",header = TRUE, sep = "")
cr.neur <- read.table("neur.lmiss",header = TRUE, sep = "")

colnames(cr.mega)[2] <- "V2"
colnames(cr.neur)[2] <- "V2"

# Restructure genomic position
mega$position <- paste0(mega$V1,":",mega$V4)
neur$position <- paste0(neur$V1,":",neur$V4)

# Merge genotype data with genomic position
mega_cr <- full_join(mega,cr.mega)
neur_cr <- full_join(neur,cr.neur)

# Merge both data frames into one
merged <- full_join(mega_cr,neur_cr)

# KEEP DUPLICATED POSITION WITH HIGHEST CALL RATE
merged <- merged[order(merged$position, 
                       merged$F_MISS),]

merged_filter <- merged %>% distinct(position, V5, V6, .keep_all= TRUE)

merged_filter <- merged_filter[,c(1,2,3,4,5,6)]

# Write list of loci that will be kept in merged file

write.table(merged_filter$V2, "keep_snps.txt", quote = F, row.names = F, col.names = F, sep = "\t")
