iwrd <- read.table("iwrd.loco.mlma", header = T, sep = "\t")
iwrd <- na.omit(iwrd)
iwrd <- iwrd[,c(1:5,7:9)]
colnames(iwrd) <- c("chr","rsid","bp","a1","a2","b","se","p")
iwrd$pos <- as.integer(iwrd$bp)
iwrd$chr <- as.integer(iwrd$chr)
write.table(iwrd, "iwrd_imputed_results.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)

pics <- read.table("pics.loco.mlma", header = T, sep = "\t")
pics <- na.omit(pics)
pics <- pics[,c(1:5,7:9)]
colnames(pics) <- c("chr","rsid","bp","a1","a2","b","se","p")
pics$pos <- as.integer(pics$bp)
pics$chr <- as.integer(pics$chr)
write.table(pics, "pics_imputed_results.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)

list <- read.table("list.loco.mlma", header = T, sep = "\t")
list <- na.omit(list)
list <- list[,c(1:5,7:9)]
colnames(list) <- c("chr","rsid","bp","a1","a2","b","se","p")
list$pos <- as.integer(list$bp)
list$chr <- as.integer(list$chr)
write.table(list, "list_imputed_results.txt", col.names = T,
            sep = "\t", quote= F, row.names = F)
