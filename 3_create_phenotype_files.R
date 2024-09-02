##########          STEP 1.0
##########          SETUP NECESSARY FILES FOR GWAS
# Clear working space
rm(list = ls(all.names = T))

# Install packages
if(!require(plyr)){
  install.packages("plyr")
  library(plyr)}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)}

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)}

if(!require(limma)){
  install.packages("limma")
  library(limma)}

if(!require(BiocManager)){
  install.packages("BiocManager")
  library(BiocManager)}

if(!require(qvalue)){
  BiocManager::install("qvalue")
  library(qvalue)}

if(!require(ggplot2)){
  BiocManager::install("ggplot2")
  library(ggplot2)}

if(!require(data.table)){
  install.packages("data.table")
  library(data.table)}

# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# Create merged phenotype data of unrestricted + restricted
# Load files
re <- read.csv("RESTRICTED_mmoerel_11_21_2018_3_14_17.csv", fill = TRUE)
unre <- read.csv("unrestricted_mmoerel_11_21_2018_1_55_15.csv", fill = TRUE)

# Adapt necessary columns
re <- re[,c(1:3,6,15,21)]
unre <- unre[,c(1,4,88,116,156,158,449,450)]

# Merge dataframes
pheno <- merge(re, unre, by = "Subject")

# FEMALE = 0   MALE = 1
# Change column names
# Select genotyped observations and remove HasGT column
names(pheno)[1] <- paste("IID")
names(pheno)[2] <- paste("age")
names(pheno)[7] <- paste("sex")
names(pheno)[4] <- paste("FID")
pheno$sex <- as.character(pheno$sex)
pheno$sex[pheno$sex=="M"] <- 0
pheno$sex[pheno$sex=="F"] <- 1
pheno$sex <- as.factor(pheno$sex)

library(dplyr)
pheno <- pheno[pheno$HasGT=="true",]
pheno <- pheno[,-3]
pheno <- pheno[,c(3,1,2,6,5,4,7,8,9,10,11,12)]

pca <- read.table("merged_plink_qc_pca.eigenvec",header = F, sep = "")
colnames(pca) <- c("FID","IID","pc1","pc2","pc3","pc4","pc5","pc6","pc7",
                   "pc8","pc9","pc10","pc11","pc12","pc13","pc14","pc15",
                   "pc16","pc17","pc18","pc19","pc20")
pheno_pca <- full_join(pheno,pca)

pheno_pca <- pheno_pca[!is.na(pheno_pca$pc1),]

pheno_pca <- pheno_pca[,-c(11,12)]
pheno_pca <- pheno_pca[!is.na(pheno_pca$age),]
pheno_pca <- pheno_pca[!is.na(pheno_pca$sex),]

# RUN THIS SECTION ONLY TO GET CORRELATIONS BETWEEN MEMORY TEST RESIDUALS
# pheno_pca <- pheno_pca[!is.na(pheno_pca$IWRD_TOT),]
# pheno_pca <- pheno_pca[!is.na(pheno_pca$PicSeq_Unadj),]
# pheno_pca <- pheno_pca[!is.na(pheno_pca$ListSort_Unadj),]
# 
# rownames(pheno_pca) <- pheno_pca$IID
# pheno_pca$IWRD_TOT_res <- glm(formula = IWRD_TOT~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = pheno_pca)$residuals
# pheno_pca$PicSeq_Unadj_res <- glm(formula = PicSeq_Unadj~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = pheno_pca)$residuals
# pheno_pca$ListSort_Unadj_res <- glm(formula = ListSort_Unadj~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = pheno_pca)$residuals
# 
# pheno <- pheno_pca[,-c(1:30)]
# 
# cormat <- pheno %>%
#   cor(use = 'everything', method = 'pearson')
# 
# corrplot(cormat, method = 'color', addCoef.col = 'black',
#          type = "lower", tl.col = 'black', cl.pos = 'n',
#          mar=c(0,0,2,0),number.cex = 0.75)
# 
# library(rstatix)
# cor_mat(
#   pheno,
#   vars = NULL,
#   method = "pearson",
#   alternative = "two.sided",
#   conf.level = 0.95)
# 
# cor_pmat(
#   pheno,
#   method = "pearson",
#   alternative = "two.sided",
#   conf.level = 0.95)

iwrd <- pheno_pca[!is.na(pheno_pca$IWRD_TOT),]
pics <- pheno_pca[!is.na(pheno_pca$PicSeq_Unadj),]
list <- pheno_pca[!is.na(pheno_pca$ListSort_Unadj),]

# Calculate residuals
iwrd$IWRD_TOT_res <- glm(formula = IWRD_TOT~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = iwrd)$residuals
pics$PicSeq_Unadj_res <- glm(formula = PicSeq_Unadj~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = pics)$residuals
list$ListSort_Unadj_res <- glm(formula = ListSort_Unadj~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20, data = list)$residuals

# Create qq-plots and lambda values
qqnorm(iwrd$IWRD_TOT_res)
qqline(iwrd$IWRD_TOT_res)
shapiro.test(iwrd$IWRD_TOT_res)

# Create qq-plots and lambda values
qqnorm(pics$PicSeq_Unadj_res)
qqline(pics$PicSeq_Unadj_res)
shapiro.test(pics$PicSeq_Unadj_res)

# Create qq-plots and lambda values
qqnorm(list$ListSort_Unadj_res)
qqline(list$ListSort_Unadj_res)
shapiro.test(list$ListSort_Unadj_res)

# Save the phenotype files
write.table(iwrd[,c(1,2,31)], "pheno_iwrd_residuals_3.txt", sep = "\t", col.names = F,
            row.names = F, quote = F)

write.table(pics[,c(1,2,31)], "pheno_pics_residuals_3.txt", sep = "\t", col.names = F,
            row.names = F, quote = F)

write.table(list[,c(1,2,31)], "pheno_list_residuals_3.txt", sep = "\t", col.names = F,
            row.names = F, quote = F)
