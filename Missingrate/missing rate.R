library(tidyverse)
setwd("~/Desktop/Cranberry/20250214/Missing rate")

#calculating counts of total reads as size
dall <- read.csv("DCran23-8178_MADC_rmDupTags_snpID_rename_targetSNP.csv")
d_sum <- aggregate(dall[,4:ncol(dall)], by = dall[2], FUN = sum)
dim(d_sum)

#exclude the first colunm "number",rownames
rownames(d_sum)<-d_sum[,1]
d_sum<- d_sum[,-1]

passport <- read.csv("Passport information.csv")
names(passport)

#F1_97 <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$PopulationType=="CNJ97 (diploid)")$AccessionID_changed)]
#F1_99 <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$PopulationType=="CNJ99 (diploid)")$AccessionID_changed)]
F1 <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$Population=="F1 (diploid)")$AccessionID_changed)]
div_diplo <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$Population=="Diverse (diploid)")$AccessionID_changed)]
div_auto <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$Population=="Diverse (autotetraploid)")$AccessionID_changed)]
hybrid <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$Population=="Interspecific hybrid (allotetraploid)")$AccessionID_changed)]
blueberry <- d_sum[,which(colnames(d_sum) %in% subset(passport,passport$Population=="Blueberry")$AccessionID_changed)]

d_sum_all <- cbind.data.frame(F1,div_diplo,div_auto,hybrid,blueberry)


##sample-level missing rate

missrate <- c()
for (i in 1:ncol(d_sum_all)) {
  missrate[i] <-length(which(d_sum_all[,i] < 10))/nrow(d_sum_all)
}
boxplot(missrate)

mr <- data.frame("ID"=colnames(d_sum_all), "missrate"= missrate)


mr_all <- rbind.data.frame(mr,mr)
mr_all$Population <- rep(c("Entire population", "F1","Diverse (2x)","Diverse (4x)","Interspecific hybrid","Blueberry"),c(372,136,194,22,18,2))
mr_all$Population <- factor(mr_all$Population,levels = c("Entire population", "F1","Diverse (2x)","Diverse (4x)","Interspecific hybrid","Blueberry"))
#write.csv(mr_all,"sample-based missing rate.csv",row.names = F)

library(ggrepel)
library(factoextra)
mr_all$label <- NA
# Then 'relabel' the points of interest
#pop_rate_line_combined[,which(grepl(pattern = "CNJ16_52_57",rownames(pop_rate_line_combined)))]$label <- "CNJ16_52_57"
mr_all[which(grepl(pattern = "Blueberry_01",mr_all$ID)),]$label <- "Blueberry_01"
mr_all[which(grepl(pattern = "Blueberry_02",mr_all$ID)),]$label <- "Blueberry_02"
mr_all[which(grepl(pattern = "CNJ16_52_57",mr_all$ID)),]$label <- "CNJ16_52_57"
mr_all[which(grepl(pattern = "CNJ16_52_62",mr_all$ID)),]$label <- "CNJ16_52_62"

give.n<- function(y) {
  return( 
    data.frame(
      y =1.1*max(y),  #may need to modify this depending on your data
      label = paste('n = ', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n'
      )
    )
  )
}

png("line_missrate_cranberry2.png",height = 2500,width = 4000,res = 300)
ggplot(mr_all,aes(x=mr_all$Population,y=mr_all$missrate*100,
                  fill=mr_all$Population))+geom_boxplot(width=0.5)+ 
  theme_bw()+ ylim(c(0,120))+
  geom_hline(yintercept=95,linetype = "dashed",color="red")+
  geom_hline(yintercept=80,linetype = "dashed",color="red")+
  geom_text(aes(x=0,y = 80, label = "80", color = "black"), size = 5, hjust = 0,inherit.aes = FALSE)+
  geom_text(aes(x=0,y = 95, label = "95", color = "black"), size = 5, hjust = 0,inherit.aes = FALSE)+
  #ggtitle("Sample-based missing rate")+
  theme(axis.text.x=element_text(angle=0,hjust=0.5),
       # plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        #axis.title.x = element_text(face = "bold"), 
        #axis.title.y = element_text(face = "bold"),
        legend.position = "none",text = element_text(size = 20))+
  xlab("Population")+ ylab("Missing rate (%)")+
  stat_summary(fun.data = give.n, geom = "text", fun = median)+
  geom_text_repel(aes(label = mr_all$label),box.padding = 0.35,col="black")
dev.off()



### marker-level missing rate
d_sum_all_delete <- d_sum_all[,-grep(pattern = "CNJ16_52_57|CNJ16_52_62",colnames(d_sum_all))]
write.csv(d_sum_all_delete,"test.csv")
# all
missrate_site <- c()
for (i in 1:nrow(d_sum_all_delete)) {
  missrate_site[i] <-length(which(d_sum_all_delete[i,] < 10))/ncol(d_sum_all_delete)
}

mr_site <- data.frame(SNP=rownames(d_sum_all_delete),missrate=missrate_site)
write.csv(mr_site,"marker-based missing rate.csv",row.names = F)
length(which(mr_site$missrate >=0.95))
mr_site_95 <- mr_site[mr_site$missrate>=0.95,]
boxplot(missrate_site)

## remove samples and markers with >=95% missing rate
d_final <- dall[-which(dall$CloneID %in% mr_site_95$SNP),-grep(pattern = "CNJ16_52_57|CNJ16_52_62",colnames(dall))]
write.csv(d_final,"DCran23-8178_MADC_rmDupTags_snpID_rename_targetSNP_missrate.csv",row.names = F)
