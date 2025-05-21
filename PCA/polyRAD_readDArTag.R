library(polyRAD)
library(qqman)
library(dplyr)
library(data.table)

setwd("~/Desktop/Cranberry/20250214/PCA")

# Importing the data
d <- read.csv("DCran23-8178_MADC_rmDupTags_snpID_rename_targetSNP_missrate_delete2F1mixedploidy.csv")
botloci <- unique(d$CloneID)

set.seed(2025)

data<-readDArTag("DCran23-8178_MADC_rmDupTags_snpID_rename_targetSNP_missrate_delete2F1mixedploidy.csv",
           botloci=botloci,
           sample.name.row = 1,
           possiblePloidies = list(2),
           taxaPloidy = 2L,
           n.header.rows = 0)

# Look at the RADdata object
data

# View the imported taxa names
samples <- GetTaxa(data)

# PCA plot to check sample clustering
#data_pca <- AddPCA(data)
source("/Users/sc3339-admin/Desktop/Cranberry/20250214/PCA/AddPCA function.R")
data_pca <- AddPCA.RADdata(data)



library(tibble)
pca <- as.data.frame(data_pca$PCA)
pca <- tibble::rownames_to_column(pca, "Sample_ID")

passport <- read.csv("Passport information.csv")
names(passport)[2] <- "Sample_ID"

pca_df <- merge(pca,passport,by="Sample_ID")

rownames(pca_df) <- pca_df[,1]
pca_df<- pca_df[,-1]
unique(pca_df$Population)
  
#all pop
pca_df$Population <- factor(pca_df$Population, levels = c("F1: CNJ97 (diploid)", "F1: CNJ99 (diploid)","F1: parent (diploid)",
                                                          "Diverse (diploid)","Diverse (autotetraploid)","Interspecific hybrid (allotetraploid)","Blueberry (autotetraploid)"))


library(ggrepel)
library(factoextra)
pca_df$label <- NA
which(grepl(pattern = "CNJ97_105_4",rownames(pca_df)))
# Then 'relabel' the points of interest
pca_df[which(grepl(pattern = "CNJ97_105_4",rownames(pca_df))),]$label <- "CNJ97_105_4"
pca_df[which(grepl(pattern = "NJS98_18",rownames(pca_df))),]$label <- "NJS98_18"
pca_df[which(grepl(pattern = "CNJ99_9_96",rownames(pca_df))),]$label <- "CNJ99_9_96"
pca_df[which(grepl(pattern = "Blueberry_01",rownames(pca_df))),]$label <- "Blueberry_01"
pca_df[which(grepl(pattern = "Blueberry_02",rownames(pca_df))),]$label <- "Blueberry_02"
pca_df[which(grepl(pattern = "Macro4x_04",rownames(pca_df))),]$label <- "Macro4x_04"
pca_df[which(grepl(pattern = "Macro4x_01",rownames(pca_df))),]$label <- "Macro4x_01"
pca_df[which(grepl(pattern = "Macro4x_02",rownames(pca_df))),]$label <- "Macro4x_02"
pca_df[which(grepl(pattern = "Macro4x_03",rownames(pca_df))),]$label <- "Macro4x_03"
pca_df[which(grepl(pattern = "Oxy4x_11",rownames(pca_df))),]$label <- "Oxy4x_11"
pca_df[which(grepl(pattern = "Oxy4x_12",rownames(pca_df))),]$label <- "Oxy4x_12"
pca_df[which(grepl(pattern = "Oxy4x_13",rownames(pca_df))),]$label <- "Oxy4x_13"
pca_df[which(grepl(pattern = "Oxy4x_14",rownames(pca_df))),]$label <- "Oxy4x_14"

# pca_df[which(grepl(pattern = "CNJ16_50_16",rownames(pca_df))),]$label <- "CNJ16_50_16"
# pca_df[which(grepl(pattern = "CNJ16_50_18",rownames(pca_df))),]$label <- "CNJ16_50_18"
# pca_df[which(grepl(pattern = "CNJ16_41_31",rownames(pca_df))),]$label <- "CNJ16_41_31"
# pca_df[which(grepl(pattern = "CNJ16_45_24",rownames(pca_df))),]$label <- "CNJ16_45_24"
# pca_df[which(grepl(pattern = "CNJ16_45_26",rownames(pca_df))),]$label <- "CNJ16_45_26"


#parents
pca_df$Parent <- NA
pca_df[which(grepl(pattern = "CNJ97_105_4",rownames(pca_df))),]$Parent <- "Parent of F1"
pca_df[which(grepl(pattern = "NJS98_18",rownames(pca_df))),]$Parent <- "Parent of F1"
pca_df[which(grepl(pattern = "CNJ99_9_96",rownames(pca_df))),]$Parent <- "Parent of F1"
pca_df[which(grepl(pattern = "Blueberry_01",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Blueberry_02",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Macro4x_04",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_11",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_12",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_13",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_14",rownames(pca_df))),]$Parent <- "Parent of hybrids"
pca_df$Parent[is.na(pca_df$Parent)] <- "Unknown"

#write.csv(pca_df,"PCA.csv")

png("PCA.png",width  = 3900,height = 2500,res=300)

ggplot(pca_df,aes(pca_df$PC1,pca_df$PC2,color=Population,shape=Parent)) + 
  geom_point(size=3,alpha=0.8) +
  scale_color_manual(values=c("#d7191c","#2c7bb6","#fa9fb5","#bdbdbd","#abd9e9","#fdae61","#756bb1"))+
 scale_shape_manual(values = c(11,12,16))+
  #stat_ellipse(level = 0.95, show.legend = F)+
  #geom_text(aes(label=Population),vjust = "outward") + 
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme_bw() + theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  #labs(x="PC1 (31.56%)",y="PC2 (6.05%)",title = "")+
  labs(x="PC1 (9.6%)",y="PC2 (6.3%)",title="")+
  theme(text =element_text(size = 15))+
  geom_text_repel(aes(label = pca_df$label),box.padding = 0.5,col="black")

dev.off()
