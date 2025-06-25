library(tidyverse)
setwd("~/Desktop/Upoad github for cranberry/PIC")
d_diplo_t <- read.csv("InforMarker_oxy4x.csv",row.names = 1,check.names = F)

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

Fre <-calc_allele_frequencies(d_diplo_t,4)


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
write.csv(PIC_df,"PIC_oxy4x.csv",row.names = F)

## plot
macro <- read.csv("PIC_macro2x(diverse).csv")
summary(macro$PIC)
micro <- read.csv("PIC_micro2x.csv")
summary(micro$PIC)
oxy <- read.csv("PIC_oxy4x.csv")
summary(oxy$PIC)

all <- merge(macro,micro,by="SNP",all = T)
all2 <- merge(all,oxy,by="SNP",all = T)
write.csv(all2, "PIC_threeSpecies.csv",row.names = F)

library(RColorBrewer)
d_all <- rbind.data.frame(marco,mirco,oxy)
d_all$Species <- rep(c("V. macrocarpon (2x)","V. microcarpum (2x)","V. oxycoccos (4x)"),c(1815,556,1046))
##histogram
png("PIC_Species.png",height = 2000,width = 3000,res = 300)
ggplot(d_all) +
  geom_histogram(aes(PIC, fill = Species), 
                 alpha=0.8, position="identity")+
  xlim(0,0.4)+ylab("SNP count")+
  scale_fill_manual(values = c("#a6d854", "#fc8d62", "#8da0cb"),
                    labels = c(
                      expression(" "*italic("V. macrocarpon")*" (2x)"),
                      expression(" "*italic("V. microcarpum")*" (2x)"),
                      expression(" "*italic("V. oxycoccos")*" (4x)")
                    ))+theme_bw()+
  theme(text = element_text(size = 20))
dev.off()

ggplot(PIC_df,aes(PIC)) +
  geom_histogram(alpha=0.3, position="identity",bins = 20)+
  theme_bw()+
  theme(text = element_text(size = 15))

#Shared markers
macro_overlap <- macro[which(macro$SNP %in% Reduce(intersect, list(macro$SNP, micro$SNP, oxy$SNP))),]
summary(macro_overlap$PIC)
length(which(macro_overlap$PIC>=0.3))
micro_overlap <- micro[which(micro$SNP %in% Reduce(intersect, list(macro$SNP, micro$SNP, oxy$SNP))),]
summary(micro_overlap$PIC)
length(which(micro_overlap$PIC>=0.3))
oxy_overlap <- oxy[which(oxy$SNP %in% Reduce(intersect, list(macro$SNP, micro$SNP, oxy$SNP))),]
summary(oxy_overlap$PIC)
length(which(oxy$PIC>=0.3))

all <- cbind.data.frame(macro_overlap,micro_overlap[,2],oxy_overlap[,2])
length(which(all$PIC>=0.3 & all$`micro_overlap[, 2]`>=0.3&all$`oxy_overlap[, 2]`>=0.3))
View(all[which(all$PIC >= 0.25 & all$`micro_overlap[, 2]` >= 0.25 & all$`oxy_overlap[, 2]` >= 0.25),])

#Unique markers
macro_unique <- macro[-which(macro$SNP %in% (c(intersect(macro$SNP, micro$SNP),intersect(macro$SNP, oxy$SNP)))),]
summary(macro_unique$PIC)
length(which(macro_unique$PIC>=0.3))
micro_unique <- micro[-which(micro$SNP %in% (c(intersect(macro$SNP, micro$SNP),intersect(micro$SNP, oxy$SNP)))),]
summary(micro_unique$PIC)
length(which(micro_unique$PIC>=0.3))
oxy_unique <- oxy[-which(oxy$SNP %in% (c(intersect(macro$SNP, oxy$SNP),intersect(micro$SNP, oxy$SNP)))),]
summary(oxy_unique$PIC)
length(which(oxy_unique$PIC>=0.3))

