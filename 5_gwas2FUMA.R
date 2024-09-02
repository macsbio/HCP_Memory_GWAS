iwrd <- read.table("geno_assoc_iwrd_residuals.fastGWA", header = T, sep = "\t")
iwrd <- na.omit(iwrd)
iwrd$POS <- as.integer(iwrd$POS)
iwrd$CHR <- as.integer(iwrd$CHR)
iwrd$N <- as.integer(iwrd$N)
iwrd <- iwrd[iwrd$AF1>=0.05 & iwrd$AF1<=0.95,]
iwrd$CHR[iwrd$CHR==24] <- 23
iwrd$CHR[iwrd$CHR==25] <- 23
write.table(iwrd, "FUMA_input_iwrd.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)

pics <- read.table("geno_assoc_pics_residuals.fastGWA", header = T, sep = "\t")
pics <- na.omit(pics)
pics$POS <- as.integer(pics$POS)
pics$CHR <- as.integer(pics$CHR)
pics$N <- as.integer(pics$N)
pics <- pics[pics$AF1>=0.05 & pics$AF1<=0.95,]
pics$CHR[pics$CHR==24] <- 23
pics$CHR[pics$CHR==25] <- 23
write.table(pics, "FUMA_input_pics.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)

list <- read.table("geno_assoc_list_residuals.fastGWA", header = T, sep = "\t")
list <- na.omit(list)
list$POS <- as.integer(list$POS)
list$CHR <- as.integer(list$CHR)
list$N <- as.integer(list$N)
list <- list[list$AF1>=0.05 & list$AF1<=0.95,]
list$CHR[list$CHR==24] <- 23
list$CHR[list$CHR==25] <- 23
write.table(list, "FUMA_input_list.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)

