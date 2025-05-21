library(tidyverse)
setwd("~/Desktop/Cranberry/20250214/PIC")

d_diplo <- read.csv("InforMarker_diverse(diploid).csv",row.names = 1)
#names(d_diplo)[1] <- "SNP"
d_auto <- read.csv("InforMarker_diverse(autotetraploid).csv",row.names = 1)
#names(d_auto)[1] <- "SNP"
d_diplo_overlap <- d_diplo[which(rownames(d_diplo) %in% intersect(rownames(d_diplo),rownames(d_auto))),]
d_auto_overlap <- d_auto[which(rownames(d_auto)  %in% intersect(rownames(d_diplo),rownames(d_auto))),]
identical(rownames(d_diplo_overlap),rownames(d_auto_overlap))


### Function for allele frequency calculation
calc_allele_frequencies <- function(d_diplo_overlap,ploidy1, d_auto_overlap,ploidy2) {
  # Diploid: 2 copies per individual
  result_2x <- apply(d_diplo_overlap, 1, function(x) {
    allele_sum <- sum(x, na.rm = TRUE)
    allele_count <- sum(!is.na(x))*ploidy1
    c(allele_sum = allele_sum, allele_count = allele_count)
  })
  
  # Autotetraploid: 4 copies per individual
  result_4x <- apply(d_auto_overlap, 1, function(y) {
    allele_sum <- sum(y, na.rm = TRUE)
    allele_count <- sum(!is.na(y)) *ploidy2
    c(allele_sum = allele_sum, allele_count = allele_count)
  })
  
  # Combine both
  total_allele_sum <- result_2x["allele_sum", ] + result_4x["allele_sum", ]
  total_allele_count <- result_2x["allele_count", ] + result_4x["allele_count", ]
  allele_freq <- total_allele_sum / total_allele_count
  
  all_allele_frequencies <- data.frame(
    SNP = rownames(d_diplo_overlap),
    p1 = allele_freq,
    p2 = 1 - allele_freq,
    maf = pmin(allele_freq, 1 - allele_freq)
  )
  
  return(all_allele_frequencies)
}

Fre <-calc_allele_frequencies(d_diplo_overlap,2, d_auto_overlap,4)


### Function for PIC calculation
calc_pic <- function(x) {
  freq_squared <- x^2
  outer_matrix <- outer(freq_squared, freq_squared)
  upper_tri_sum <- sum(outer_matrix[upper.tri(outer_matrix)])
  pic <- 1 - sum(freq_squared) - 2*upper_tri_sum
  return(pic)
}

PIC_results <- apply(Fre[, c("p1", "p2")], 1, calc_pic)
PIC_df <- data.frame(SNP = Fre$SNP, PIC = PIC_results)
summary(PIC_df$PIC)
length(which(PIC_df$PIC>=0.3))
#write.csv(PIC_df,"PIC_acrossDiploTetraplo.csv",row.names = F)

## plot
ggplot(d_final,aes(PIC)) +
  geom_histogram(alpha=0.3, 
                 # position="identity",
                 bins = 30,
                 color = "black",
                 #linetype="dashed",
                 size=0.1)+
  ylab("SNP count")+xlab("PIC acorss the entire sample set ")+
  geom_vline(xintercept = 0.3,linetype = "dashed", col = "red",lwd=0.5)+
  theme_bw()+
  theme(text = element_text(size = 15))

plot(PIC_df$PIC)


