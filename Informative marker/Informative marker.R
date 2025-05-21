setwd("~/Desktop/Cranberry/20250214/InforMarker")
d <- read.csv("Genotype_Hybrid(diploid markers).csv",row.names = 1)
d <- as.matrix(d)

Count <- c()
for (i in 1:nrow(d)) {
 Count[i] <- length(unique(d[i,]))
}
Count_d <- data.frame(SNP=rownames(d),length=Count)
Count_multi <- Count_d[Count_d$length > 1,]
d_dive <- d[which(rownames(d) %in% Count_multi$SNP),]
write.csv(d_dive,"InforMarker_Hybrid(diploid markers).csv")



