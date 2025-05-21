library(tidyverse)
setwd("~/Desktop/Cranberry/20250214/PIC")
d_diplo_t <- read.csv("InforMarker_diverse(diploid).csv",row.names = 1,check.names = F)

### Function for allele frequency calculation
calc_allele_frequencies <- function(d_diplo_t, ploidy) {
  allele_frequencies <- apply(d_diplo_t, 1, function(x) {
    count_sum <- sum(!is.na(x))  
    allele_sum <- sum(x, na.rm = TRUE) 
    if (count_sum != 0) {allele_sum / (ploidy * count_sum)} else {NA}
  })
  
  all_allele_frequencies <- data.frame(SNP = rownames(d_diplo_t), p1= allele_frequencies, p2= 1-allele_frequencies,maf=pmin(allele_frequencies,(1-allele_frequencies)))
  return(all_allele_frequencies)
}

Fre <-calc_allele_frequencies(d_diplo_t,2)


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
write.csv(PIC_df,"PIC_Hybrid (diploid).csv",row.names = F)

## plot
ggplot(PIC_df,aes(PIC)) +
  geom_histogram(alpha=0.3, position="identity",bins = 20)+
  theme_bw()+
  theme(text = element_text(size = 15))
