library(tidyverse)
library(data.table)
setwd("~/Desktop/Cranberry/20250214/LinkageMap")
d <- read.csv("recFRE_all_CNJ99.csv",row.names = 1)
#Chrom <- paste0("chr", formatC(1:12, width = 2, flag = "0"))
Chrom <- unique(stringr::str_split_fixed(rownames(d),pattern = "\\_",n=2)[,1])

SNP_nearby <- data.frame(SNP1 = character(), SNP2 = character(), Value = numeric())
Result <- data.frame(SNP1 = character(), SNP2 = character(), Value = numeric())

for (i in 1:length(Chrom)) {
 
  ### find adjacent SNPs leading to high recombination frequency
  
  Chr <- d[grep(pattern = Chrom[i], rownames(d)), grep(pattern = Chrom[i], colnames(d))]
  indices <- which(Chr > 0.4, arr.ind = TRUE)
  # adjacent_SNP <- indices[abs(indices[, 1] - indices[, 2]) == 1, , drop = FALSE] 
  adjacent_SNP <- lapply(1:5, function(n) indices[abs(indices[, 1] - indices[, 2]) == n, , drop = FALSE])
  adjacent_SNP <- Filter(function(x) nrow(x) > 0, adjacent_SNP)
  new_input <- data.frame(SNP1 = character(), SNP2 = character(), RecFre = numeric())
  
  if (length(adjacent_SNP) > 0) {
    adjacent_SNP_combined <- do.call(rbind, adjacent_SNP)
    new_input <- data.frame(
      SNP1 = rownames(Chr)[adjacent_SNP_combined[, 1]], 
      SNP2 = colnames(Chr)[adjacent_SNP_combined[, 2]], 
      RecFre = Chr[adjacent_SNP_combined],
      stringsAsFactors = FALSE
    )
  
  ### according to adjacent SNPs, to find linked SNP in other chromosomes 

    SNP_output <- d[which(rownames(d) %in% unique(new_input[,1])), -grep(pattern = Chrom[i], colnames(d)), drop = FALSE]
    Result_indices <- which(SNP_output < 0.1, arr.ind = TRUE)
  
    if (length(Result_indices) > 0) {
      new_input2 <- data.frame(
        SNP1 = rownames(SNP_output)[Result_indices[, 1]], 
        SNP2 = colnames(SNP_output)[Result_indices[, 2]], 
        RecFre = SNP_output[Result_indices],
        stringsAsFactors = FALSE
      )
      Result <- rbind(Result, new_input2)
    }
  }
  SNP_nearby <- rbind(SNP_nearby, new_input)
  } 

#write.csv(Result,"SNP_mismatchbased-RECFRE_99.csv",row.names = F)
Result1 <- aggregate(Result[,3],by=Result[1],min)
merged_result <- merge(Result1, Result, by.x = c(names(Result1)[1], "x"), by.y = c(names(Result)[1], names(Result)[3]))
sort(unique(merged_result$SNP1))
merged_result$Chrom <- gsub("^chr","",stringr::str_split_fixed(merged_result$SNP2,pattern = "\\_",n=2)[,1])
merged_result1 <- unique(merged_result[,c(1,4)])

#write.table(merged_result1,"Result_linkedSNP_RECmatrix3_5_CNJ99.txt",row.names = F,quote = F,sep = "\t")

d_SNP <- fread("Result_linkedSNP_RECmatrix3_5_CNJ99.txt")
d_97 <- read.csv("Mappoly_CNJ99_unifysharedparentGeno_deleteThreeSamplesFromPCA.csv")
dim(d_97)
names(d_97)
identical(d_97[which(d_97$SNP %in% d_SNP$SNP1),]$SNP,d_SNP$SNP1)
d_97[which(d_97$SNP %in% d_SNP$SNP1),]$chrom <- d_SNP$Chrom
d_97[which(d_97$SNP %in% d_SNP$SNP1),]$genome_position <- d_SNP$Pos
write.csv(d_97,"Mappoly_CNJ99_unifysharedparentGeno_deleteThreeSamplesFromPCA_reposition.csv",row.names = F)
