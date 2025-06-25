library(tidyverse)
library(data.table)
library(ggpmisc)
setwd("~/Desktop/Upoad github for cranberry/Marker distribution plot")

#marker distribution 

#Diverse 

macro <- read.csv("PIC_macro2x(diverse).csv")
micro <- read.csv("PIC_micro2x.csv")
oxy <- read.csv("PIC_oxy4x.csv")

length(intersect(macro$SNP,micro$SNP))
length(intersect(macro$SNP,oxy$SNP))
length(intersect(micro$SNP,oxy$SNP))
length(Reduce(intersect, list(macro$SNP, micro$SNP, oxy$SNP)))

macro_micro_oxy_overlap <- macro[which(macro$SNP %in% Reduce(intersect, list(macro$SNP, micro$SNP, oxy$SNP))),]
macro_micro_overlap <- macro[which(macro$SNP %in% intersect(macro$SNP, micro$SNP)),]
macro_oxy_overlap <- macro[which(macro$SNP %in% intersect(macro$SNP, oxy$SNP)),]
micro_oxy_overlap <- micro[which(micro$SNP %in% intersect(micro$SNP, oxy$SNP)),]
macro_unique <- macro[-which(macro$SNP %in% (c(intersect(macro$SNP, micro$SNP),intersect(macro$SNP, oxy$SNP)))),]
micro_unique <- micro[-which(micro$SNP %in% (c(intersect(macro$SNP, micro$SNP),intersect(micro$SNP, oxy$SNP)))),]
oxy_unique <- oxy[-which(oxy$SNP %in% (c(intersect(macro$SNP, oxy$SNP),intersect(micro$SNP, oxy$SNP)))),]

all <- data.frame(SNP=c(macro_micro_oxy_overlap$SNP,macro_micro_overlap$SNP,macro_oxy_overlap$SNP,micro_oxy_overlap$SNP,
                        macro_unique$SNP,micro_unique$SNP, oxy_unique$SNP),
                  Marker=rep(c("Shared marker among all","Shared marker between the macro and micro","Shared marker between the macro and oxy","Shared marker between the micro and oxy","Unique marker in the macro","Unique marker in the micro","Unique marker in the oxy"),c(343,524,955,365,679,10,69)))

all$CHROM <- stringr::str_split_fixed(all$SNP,pattern = "\\_",n=2)[,1]
all$POS <- as.numeric(stringr::str_split_fixed(all$SNP,pattern = "\\_",n=2)[,2])/10^6
all$CHROM <- factor(all$CHROM,levels = sort(unique(all$CHROM)))
all$Marker <- factor(all$Marker,levels = c("Shared marker among all","Shared marker between the macro and micro","Shared marker between the macro and oxy",
                                            "Shared marker between the micro and oxy","Unique marker in the macro","Unique marker in the micro",
                                            "Unique marker in the oxy"))

Size <- fread("sizes.genome.txt",header = F)
Size <- Size[c(1:12),]
names(Size) <- c("chromosome","size")
Size$chromosome <- sort(unique(all$CHROM))
Size$chromosome<- factor(Size$chromosome,levels = Size$chromosome)
library(RColorBrewer)
png("SNP_species.png",height = 3000,width = 4500,res = 300)
ggplot() +
  geom_segment(data = Size,
               aes(x = chromosome, xend = chromosome, y = 0, yend = size/10^6),
               lineend = "round", color = "grey", size = 2) +
 # geom_point(data = Size, aes(x = chromosome, y = 0), color = "grey70", fill = "white", size = 3, shape = 21) +
 # geom_point(data = Size, aes(x = chromosome, y = size), color = "grey70", fill = "white", size = 3, shape = 21) +
  ylab("Genome position (Mbp)") + xlab("") + theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  geom_segment(data = all,
               aes(x = as.numeric(as.factor(CHROM)) - 0.2, 
                   xend = as.numeric(as.factor(CHROM)) + 0.2, 
                   y = POS, yend = POS,colour=Marker),
               size = 0.5) +
 # scale_color_manual(values = c("#999999","red","blue"))+
  scale_color_manual(values =c("black", "blue", "red", "green", "#e6ab02", "pink","purple"))+
 # scale_color_manual(values = brewer.pal(7, "Set1"))+
  coord_flip() +
  scale_y_continuous(labels = scales::comma)
dev.off()

### Hybrid population

d_diplo_marker <- read.csv("PIC_Hybrid (diploid).csv")
d_auto_marker <- read.csv("PIC_Hybrid (autotetraploid).csv") 
all <- rbind.data.frame(d_diplo_marker,d_auto_marker)
all$Markerploidy <- rep(c("Diploid behavior","Tetraploid behavior"),c(936,212))
all$CHROM <- stringr::str_split_fixed(all$SNP,pattern = "\\_",n=2)[,1]
all$POS <- as.numeric(stringr::str_split_fixed(all$SNP,pattern = "\\_",n=2)[,2])
#all$Chromosome <- gsub("^chr","",all$Chromosome)
all$CHROM <- factor(all$CHROM,levels = sort(unique(all$CHROM)))

Size <- fread("sizes.genome.txt",header = F)
Size <- Size[c(1:12),]
names(Size) <- c("chromosome","size")
Size$chromosome <- sort(unique(all$CHROM))
Size$chromosome<- factor(Size$chromosome,levels = Size$chromosome)

png("SNP_plot_hybrid.png",height = 2300,width = 3500,res = 300)
ggplot() +
  geom_segment(data = Size,
               aes(x = chromosome, xend = chromosome, y = 0, yend = size/10^6),
               lineend = "round", color = "grey70", size = 3) +
  ylab("Genome position (Mbp)") + xlab("") + theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  geom_segment(data = all,
               aes(x = as.numeric(as.factor(CHROM)) - 0.2, 
                   xend = as.numeric(as.factor(CHROM)) + 0.2, 
                   y = POS/10^6, yend = POS/10^6,colour=Markerploidy),
                size = 0.4) +
  scale_color_manual(values = c("red","blue"),name = "Marker ploidy")+
  coord_flip() +
  scale_y_continuous(labels = scales::comma)
dev.off()


### Marker distribution of 3059 markers
Size <- fread("sizes.genome.txt",header = F)
Size <- Size[c(1:12),]
names(Size) <- c("chromosome","size")
Size$chromosome <- sort(unique(merge_all$CHROM))
Size$chromosome<- factor(Size$chromosome,levels = Size$chromosome)

d_all <- read.csv("/Users/sc3339-admin/Desktop/Cranberry/20250214/Raw file/DCran23-8178_MADC.csv")
chr <- stringr::str_split_fixed(unique(d_all$CloneID),pattern = "\\_",n=2)[,1]
pos <- stringr::str_split_fixed(unique(d_all$CloneID),pattern = "\\_",n=2)[,2]
merge_all <- data.frame(SNP=unique(d_all$CloneID), CHROM=tolower(chr),POS=as.numeric(pos))
merge_all$CHROM <- factor(merge_all$CHROM,levels = sort(unique(merge_all$CHROM)))

png("3059SNP_distribution2.png",height = 2300,width = 3500,res = 300)
ggplot() +
  geom_segment(data = Size,
               aes(x = chromosome, xend = chromosome, y = 0, yend = size/10^6),
               lineend = "round", color = "grey70", size = 3) +
  ylab("Genome position (Mbp)") + xlab("") + theme_bw() +
  theme(text = element_text(size = 15), legend.position = "none") +
  geom_segment(data = merge_all,
               aes(x = as.numeric(as.factor(CHROM)) - 0.2, 
                   xend = as.numeric(as.factor(CHROM)) + 0.2, 
                   y = POS/10^6, yend = POS/10^6,color="red"),
               size = 0.4) +
  #scale_color_manual(values = c("red","blue"))+
  coord_flip() +
  scale_y_continuous(labels = scales::comma)
dev.off()
