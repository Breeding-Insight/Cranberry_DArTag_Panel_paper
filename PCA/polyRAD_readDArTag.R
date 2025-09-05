library(polyRAD)
library(qqman)
library(dplyr)
library(data.table)

setwd("~/Desktop/Upoad github for cranberry/PCA")

# Importing the data
d <- read.csv("DCran23-8178_MADC_rmDupTags_snpID_rename_targetSNP_missrate_delete2F1mixedploidy.csv")
# botloci <- unique(d$CloneID)
botloci <- fread("Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_f180bp.botloci",header=FALSE)

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
unique(pca_df$Species)
  
#all pop
pca_df$Species <- factor(pca_df$Species, levels = c("V. macrocarpon (F1: CNJ16-45; 2x)","V. macrocarpon (F1: CNJ16-41; 2x)","V. macrocarpon (F1: parent; 2x)","V. macrocarpon (diverse; 2x)","V. macrocarpon (4x)","V. microcarpum (2x)","V. oxycoccos (4x)","V. meridionale × V. macrocarpon (2x)","V. oxycoccos × V. corymbosum (2x)","V. corymbosum (blueberry; 4x)","V. meridionale (blueberry; 4x)"))


library(ggrepel)
library(factoextra)
pca_df$label <- NA
which(grepl(pattern = "Oxy2x_11",rownames(pca_df)))
# Then 'relabel' the points of interest
pca_df[which(grepl(pattern = "CNJ97_105_4",rownames(pca_df))),]$label <- "CNJ97-105-4 (Mullica Queen)"
pca_df[which(grepl(pattern = "NJS98_18",rownames(pca_df))),]$label <- "NJS98-18"
pca_df[which(grepl(pattern = "CNJ99_9_96",rownames(pca_df))),]$label <- "CNJ99-9-96"
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

pca_df[grepl("^Macro_12$", rownames(pca_df)), "label"] <- "Macro_12 (Murphy's green)"
# pca_df[grepl("^Macro_21$", rownames(pca_df)), "label"]  <- "Macro_21"
# pca_df[grepl("^Macro_22$", rownames(pca_df)), "label"]  <- "Macro_22"
# pca_df[grepl("^Macro_23$", rownames(pca_df)), "label"]  <- "Macro_23"
# pca_df[grepl("^Macro_24$", rownames(pca_df)), "label"]  <- "Macro_24"
pca_df[grepl("^Macro_20$", rownames(pca_df)), "label"] <- "Macro_20 (Yellow Bell)"
pca_df[grepl("^Macro_03$", rownames(pca_df)), "label"]  <- "Macro_03 (Ben Lear)"
pca_df[grepl("^US93_71$", rownames(pca_df)), "label"]  <- "US93_71 (Ben Lear)"
pca_df[grepl("^Macro_16$", rownames(pca_df)), "label"]  <- "Macro_16 (Stevens)"
pca_df[grepl("^US88_32$", rownames(pca_df)), "label"]  <- "US88_32 (Stevens)"
pca_df[grepl("^Macro_05$", rownames(pca_df)), "label"]  <- "Macro_05 (Budd's blues)"
pca_df[grepl("^US88_30$", rownames(pca_df)), "label"]  <- "US88_30 (Budd's blues)"
# pca_df[grepl("^Macro_17$", rownames(pca_df)), "label"]  <- "Macro_17 (Sundance)"
# pca_df[grepl("^Macro_10$", rownames(pca_df)), "label"] <- "Macro_10 (HyRed)"
pca_df[grepl("^Macro_06$", rownames(pca_df)), "label"]  <- "Macro_06 (Crimson Queen)"
pca_df[grepl("^NJS98_23$", rownames(pca_df)), "label"]  <- "NJS98_23 (Crimson Queen)"
# pca_df[grepl("^Macro_18$", rownames(pca_df)), "label"]  <- "Macro_18 (Sweetie)"
pca_df[grepl("^Macro_01$", rownames(pca_df)), "label"]  <- "Macro_01 (#35)"
pca_df[grepl("^US88_51$", rownames(pca_df)), "label"]  <- "US88_51 (#35)"
pca_df[grepl("^Oxy2x_11$", rownames(pca_df)), "label"]  <- "Micro2x_11"
pca_df[grepl("^Macro_07$", rownames(pca_df)), "label"]  <- "Macro_07 (Franklin)"
pca_df[grepl("^US88_13$", rownames(pca_df)), "label"]  <- " US88_13 (Franklin)"
pca_df[grepl("^US94_13$", rownames(pca_df)), "label"]  <- " US94_13 (Franklin)"
pca_df[grepl("^Macro_09$", rownames(pca_df)), "label"]  <- "Macro_09 (Howes)"
pca_df[grepl("^US88_65$", rownames(pca_df)), "label"]  <- "US88_65 (Howes)"
pca_df[grepl("^Macro_11$", rownames(pca_df)), "label"]  <- "Macro_11 (Mullica Queen)"
pca_df[grepl("^CNJ97_105_4$", rownames(pca_df)), "label"]  <- "CNJ97_105_4 (Mullica Queen)"
pca_df[grepl("^Macro_13$", rownames(pca_df)), "label"]  <- "Macro_13 (Pilgrim)"
pca_df[grepl("^US88_93$", rownames(pca_df)), "label"]  <- "US88_93 (Pilgrim)"
# pca_df[grepl("^US88_56$", rownames(pca_df)), "label"]  <- "US88_56 (LEMUNYON)"
# pca_df[grepl("^Macro_14$", rownames(pca_df)), "label"]  <- "Macro_14 (Potter's Favorite)"
pca_df[grepl("^US88_18$", rownames(pca_df)), "label"]  <- "US88_18 (Early Black)"
pca_df[grepl("^US93_85_4$", rownames(pca_df)), "label"]  <- "US93_85_4 (McFarlin)"
pca_df[grepl("^Macro_08$", rownames(pca_df)), "label"]  <- "Macro_08 (Grygleski1)"

#parents
pca_df$Pedigree <- NA
# pca_df[grepl("^Macro_12$", rownames(pca_df)), "Pedigree"] <- "Parent of Macro_21-24"
# pca_df[grepl("^Macro_03$", rownames(pca_df)), "Pedigree"]  <- "Parent of Macro_06&10&17&18"
# pca_df[grepl("^Macro_16$", rownames(pca_df)), "Pedigree"]  <- "Parent of Macro_06&10&17&18"
pca_df[which(grepl(pattern = "CNJ97_105_4",rownames(pca_df))),]$Pedigree <- "Parent of F1"
pca_df[which(grepl(pattern = "NJS98_18",rownames(pca_df))),]$Pedigree <- "Parent of F1"
pca_df[which(grepl(pattern = "CNJ99_9_96",rownames(pca_df))),]$Pedigree <- "Parent of F1"
pca_df[which(grepl(pattern = "Blueberry_01",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Blueberry_02",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Macro4x_04",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_11",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_12",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_13",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
pca_df[which(grepl(pattern = "Oxy4x_14",rownames(pca_df))),]$Pedigree <- "Parent of hybrids"
# pca_df[grepl("^US88_51$", rownames(pca_df)), "Pedigree"]  <- "Parent of Mullica Queen"
# pca_df[grepl("^US88_56$", rownames(pca_df)), "Pedigree"]  <- "Parent of Mullica Queen"
pca_df$Pedigree[is.na(pca_df$Pedigree)] <- "Progeny"

#write.csv(pca_df,"PCA.csv")
library(ggsci)
# png("PCA.png",width  = 3900,height = 2500,res=300)
png("PCA1.png",width  = 7500,height = 5500,res=300)
ggplot(pca_df, aes(PC1, PC2, color = Species, shape = Pedigree)) + 
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(
    values = ggsci::pal_d3("category20")(20)[c(1,2,4,11,7,3,5,9,10,6,12)],
    labels = c(
      expression(" "*italic("V. macrocarpon")*" (F"["1"]*": CNJ16-45; 2x)"),
      expression(" "*italic("V. macrocarpon")*" (F"["1"]*": CNJ16-41; 2x)"),
      expression(" "*italic("V. macrocarpon")*" (F"["1"]*": parent; 2x)"),
      expression(" "*italic("V. macrocarpon")*" (diverse; 2x)"),
      expression(" "*italic("V. macrocarpon")*" (colchicine-generated; 4x)"),
      expression(" "*italic("V. microcarpum")*" (selfs of Micro2x_11; 2x)"),
      expression(" "*italic("V. oxycoccos")*" (4x)"),
      expression(" "*italic("V. meridionale × V. macrocarpon")*" (2x)"),
      expression(" "*italic("V. oxycoccos × V. corymbosum")*" (2x)"),
      expression(" "*italic("V. corymbosum")*" (blueberry; 4x)"),
      expression(" "*italic("V. meridionale")*" (blueberry; 4x)")
    )
  ) +
  
  scale_shape_manual(values = c(11,12,16),#c(11,12, 8, 13, 9, 16)
                     labels = c("Parent of F1" = expression("Parent of F"["1"]))) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 25)
  ) +
  
  labs(x = "PC1 (9.6%)", y = "PC2 (6.3%)", title = "") +
  geom_text_repel(aes(label = label), box.padding = 0.5, size=7,col = "black")

dev.off()

library(plotly)
library("RColorBrewer")
p <- plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~pca_df$Species, 
             #text=~pca_df$label,
             text=~rownames(pca_df),
             colors = brewer.pal(n = 11, name = "Paired")) 
