library(mappoly2)
library(dplyr)
setwd("~/Desktop/Cranberry/20250214/LinkageMap")

## CNJ97 population
cranberry.f1 <- read_geno_csv(file.in = "Mappoly_CNJ97_unifysharedparentGeno_deleteTwoSamplesFromPCA_reposition _rmOneBigGap.csv",
                            ploidy.p1 = 2,
                            ploidy.p2 = 2,
                            name.p1 = "NJS98_18",
                            name.p2 = "CNJ97_105_4")

plot(cranberry.f1)
cranberry.f2 <- filter_data(cranberry.f1, mrk.thresh = 0.1, ind.thresh = 0.1) 

cranberry.f3  <- filter_individuals(cranberry.f2,type="Gmat")

#plot(cranberry.f3, type = "density")
cranberry.all <- pairwise_rf(cranberry.f2, mrk.scope = "all",ncpus = 8)
plot(cranberry.all)

cranberry.f3.all <- rf_filter(cranberry.all, thresh.LOD.ph = 5, thresh.LOD.rf = 5, thresh.rf = 0.15, probs = c(0.025, 0.975))
plot(cranberry.f3.all)
#write.csv(cranberry.f3.all$pairwise.rf$rec.mat,file = "recFRE_all_CNJ97.csv")

g <- group(x = cranberry.f3.all, expected.groups = 12, comp.mat = T, inter = F)
g

#write.csv(data.frame(g$groups.snp),"SNP_GROUP_CNJ97.csv")

#Group all
s <- make_sequence(g, 
                   lg = list(c(1,2), c(3,4), c(4,6), 7, c(8,9), c(8,10), 11, c(4,6),12,c(5,9),3,c(2,8)),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom 

#Choose this
s <- make_sequence(g, 
                   lg = list(c(1,2), 4, 6, 7,c(8,9), c(8,10), 11, 4, 12,c(5,9),3,c(2,8)),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom 


#Group after delete one markers due to BIGgap
s <- make_sequence(g, 
                   lg = list(c(1,2), c(3,4), c(4,6), 6, c(7,8), c(7,9), 10, c(6,11),12,c(5,8),3,c(2,7)),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom 

#Choose this
s <- make_sequence(g, 
                   lg = list(c(1,2), 4, 6, 6, c(7,8), c(7,9), 10, 11,12,c(5,8),3,c(2,7)),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom 


s <- order_sequence(s, type = "mds")
s <- order_sequence(s, type = "genome")

#plot_rf_matrix(s, type = "mds", fact = 5)
#plot_rf_matrix(s, type = "genome", fact = 5)

#plot_mds_vs_genome(s)

s <- pairwise_phasing(s, type = "mds",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 10, 
                      max.search.expansion.p2 = 10)
#print(s, type = "mds")

s <- pairwise_phasing(s, 
                      type = "genome",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 10, 
                      max.search.expansion.p2 = 10)

#print(s, type = "genome")


s <- mapping(s, type = "mds", parent = "p1", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p1", ncpus = 8)

#print(s, type = "mds")
#print(s, type = "genome")
#plot_map(s, lg = 1, type = "mds", parent = "p1")
#plot_map(s, lg = 1, type = "genome", parent = "p1")

s <- mapping(s, type = "mds", parent = "p2", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p2", ncpus = 8)

s <- merge_single_parent_maps(s, type = "genome", ncpus = 8, error = 0.05)

s <- merge_single_parent_maps(s, type = "mds", ncpus = 8, error = 0.05)

s <- calc_haplotypes(s, type = "mds", ncpus = 8)
s <- calc_haplotypes(s, type = "genome", ncpus = 8)

plot_map_list(s, type = "mds", parent = "p1p2",col = mp_pal(12))
plot_genome_vs_map(s, type="mds",same.ch.lg = T)

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/LinkageMap_CNJ97.png",height = 2000,width = 3000,res = 300)
plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(12))
dev.off()

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/MapVsPhysical_CNJ97.png",height = 2500,width = 3000,res = 300)
plot_genome_vs_map(s, type="genome",same.ch.lg = T)
dev.off()

map_summary(s,type = "genome")
h <- map_summary(s,type = "genome")
write.csv(h,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/map_summary_CNJ97.csv")

t <- plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(12))
#plot_map(s, lg = 1, type = "mds", parent = "p1p2")
write.csv(t,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/CNJ97_map.csv")

NJS98_18xCNJ97_105_4_map <- s
save(NJS98_18xCNJ97_105_4_map,file = "~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ97_105_4_map.rda")



## CNJ99 population
cranberry.f1 <- read_geno_csv(file.in = "Mappoly_CNJ99_unifysharedparentGeno_deleteThreeSamplesFromPCA_reposition.csv",
                              ploidy.p1 = 2,
                              ploidy.p2 = 2,
                              name.p1 = "NJS98_18",
                              name.p2 = "CNJ99_9_96")

plot(cranberry.f1)
cranberry.f2 <- filter_data(cranberry.f1, mrk.thresh = 0.1, ind.thresh = 0.1) 


cranberry.f3  <- filter_individuals(cranberry.f2,type="Gmat")

#plot(cranberry.f3, type = "density")
cranberry.all <- pairwise_rf(cranberry.f2, mrk.scope = "all", ncpus = 8)
#plot(cranberry.all)

cranberry.f3.all <- rf_filter(cranberry.all, thresh.LOD.ph = 5, thresh.LOD.rf = 5, thresh.rf = 0.15, probs = c(0.025, 0.975))
#write.csv(cranberry.f3.all$pairwise.rf$rec.mat,file = "recFRE_all_CNJ99.csv")


g <- group(x = cranberry.f3.all, expected.groups = 12, comp.mat = T, inter = F)
g

#write.csv(data.frame(g$groups.snp),"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/SNP_GROUP_CNJ99.csv")
#plot(g)

##Group all

s <- make_sequence(g, 
                   lg = list(c(1,2,3), c(5,6,1), c(2,4,1), c(1,8), 9, c(4,8), c(10,11), c(4,11), 9, c(7,12), c(3,4,1), c(5,12)),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom


#s <- make_sequence(g, 
#                   lg = list(c(1,2), c(5,6), c(2,4), 1, 9, c(4,8), c(10,11), 11, 9, c(7,12), c(3,4), c(5,12)),#LG
#                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) #Chrom




s <- order_sequence(s, type = "mds")
s <- order_sequence(s, type = "genome")

#plot_rf_matrix(s, type = "mds", fact = 5)
#plot_rf_matrix(s, type = "genome", fact = 5)

#plot_mds_vs_genome(s)

s <- pairwise_phasing(s, type = "mds",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 10, 
                      max.search.expansion.p2 = 10)
#print(s, type = "mds")

s <- pairwise_phasing(s, 
                      type = "genome",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 10, 
                      max.search.expansion.p2 = 10)

#print(s, type = "genome")


s <- mapping(s, type = "mds", parent = "p1", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p1", ncpus = 8)
#print(s, type = "mds")
#print(s, type = "genome")
#plot_map(s, lg = 1, type = "mds", parent = "p1")
#plot_map(s, lg = 1, type = "genome", parent = "p1")

s <- mapping(s, type = "mds", parent = "p2", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p2", ncpus = 8)


s <- merge_single_parent_maps(s, type = "genome", ncpus = 8, error = 0.05)
#> Multi-locus map estimation
s <- merge_single_parent_maps(s, type = "mds", ncpus = 8, error = 0.05)
#> Multi-locus map estimation
s <- calc_haplotypes(s, type = "mds", ncpus = 8)
s <- calc_haplotypes(s, type = "genome", ncpus = 8)
#plot_haplotypes(s, lg = 1, ind = "CNJ16_41_1")

plot_map_list(s, type = "mds", parent = "p1p2",col = mp_pal(12))
png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/LinkageMap_CNJ99.png",height = 2000,width = 3000,res = 300)
plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(12))
dev.off()
png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/MapVsPhysical_CNJ99.png",height = 2500,width = 3000,res = 300)
plot_genome_vs_map(s, type="genome",same.ch.lg = T)
dev.off()
plot_genome_vs_map(s, type="mds",same.ch.lg = T)
map_summary(s,type = "genome")
h <- map_summary(s,type = "genome")
#write.csv(h,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/map_summary_CNJ99.csv")

t <- plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(12))
#plot_map(s, lg = 1, type = "mds", parent = "p1p2")

#write.csv(t,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/CNJ99_map.csv",row.names = F)

NJS98_18xCNJ99_9_96_map <- s

save(NJS98_18xCNJ99_9_96_map,file = "~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ99_9_96_map.rda")


## consensus.map

load("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ97_105_4_map.rda")
load("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ99_9_96_map.rda")

#combined map

MAPs <- list("NJS98_18&CNJ97_105_4" =NJS98_18xCNJ97_105_4_map, 
             "NJS98_18&CNJ99_9_96" = NJS98_18xCNJ99_9_96_map)
plot_multi_map(MAPs)

prep.maps <- prepare_to_integrate(MAPs,type = "genome")
plot(prep.maps)
plot_shared_markers(prep.maps)

consensus.map <- estimate_consensus_map(prep.maps, ncpus = 8, err = 0.05)
save(consensus.map,file = "~/Desktop/Cranberry/20250214/LinkageMap/Result_map/consensus.map.rda")

plot(consensus.map)

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/combined map.png",height = 2000,width = 3000,res = 300)
plot(consensus.map, only.consensus = TRUE, col = mp_pal(12),horiz=F)
#plot_shared_markers(prep.maps)
dev.off()
h <- plot(consensus.map)
#write.csv(h$data,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/consensus_map_information.csv",row.names = F)


h_map_pos <- plot(consensus.map)
h_map_pos_marker <- h_map_pos$plot_env$maps2

length(which(h_map_pos$plot_env$marker_counts$n=="2"))#shared markers
share_markers <- h_map_pos$plot_env$marker_counts[which(h_map_pos$plot_env$marker_counts$n=="2"),]$mrk.names#shared markers

h_97_marker <- h_map_pos$plot_env$maps1[grep(pattern = "NJS98_18xCNJ97_105_4",h_map_pos$plot_env$maps1$POP),]
h_99_marker <- h_map_pos$plot_env$maps1[grep(pattern = "NJS98_18xCNJ99_9_96",h_map_pos$plot_env$maps1$POP),]


library(ggVennDiagram)
colnames(h_97_marker)[3] <- "CNJ97"
colnames(h_99_marker)[3] <- "CNJ99"
x <- (list("CNJ97"=h_97_marker$CNJ97,"CNJ99"= h_99_marker$CNJ99))

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/marker_venn.png",height = 2000,width = 2200,res = 300)
ggVennDiagram(x,label = "count",edge_size = 0,set_size=5,label_size=5) + 
  #scale_fill_gradient2(low = "green", mid="red", high = "grey",midpoint = 500) +
  scale_fill_distiller(palette = "RdBu")+
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .15))+
  coord_flip()
dev.off()


# 

library(ggpubfigs) 
library(tidyverse)

h_share_markers <- h_map_pos_marker %>% filter(h_map_pos_marker$mrk.names %in% share_markers)
h_97_filter <- h_97_marker %>% filter(!h_97_marker$CNJ97 %in% share_markers)
h_97_filter_2 <- h_map_pos_marker %>% filter(h_map_pos_marker$mrk.names %in% h_97_filter$CNJ97)

h_99_filter <- h_99_marker %>% filter(!h_99_marker$CNJ99 %in% share_markers)
h_99_filter_2 <- h_map_pos_marker %>% filter(h_map_pos_marker$mrk.names %in% h_99_filter$CNJ99)


h_map_physi <- cbind.data.frame(MarkerGroup=c(rep("Shared marker",nrow(h_share_markers)),rep("Unique marker in CNJ97",nrow(h_97_filter_2)),rep("Unique marker in CNJ99",nrow(h_99_filter_2))),
                                           LG=c(h_share_markers$LG,h_97_filter_2$LG,h_99_filter_2$LG),
                                           SNP=c(h_share_markers$mrk.names,h_97_filter_2$mrk.names,h_99_filter_2$mrk.names),
                                           Map_pos=c(h_share_markers$pos,h_97_filter_2$pos,h_99_filter_2$pos)
                                           )

h_map_physi$physical_map <- as.numeric(stringr::str_split_fixed(h_map_physi$SNP,pattern = "_",n=2)[,2])
h_map_physi$Chrom <- NA
h_map_physi$reposition <-NA
d_CNJ97_99 <- read.csv("Result_linkedSNP_RECmatrix3_5_CNJ97&99.csv")
d_CNJ97_99_overlap <- d_CNJ97_99[which(d_CNJ97_99$SNP %in% h_map_physi$SNP),]
identical(h_map_physi[which(h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP),]$SNP,d_CNJ97_99_overlap$SNP)

selected_SNPs <- h_map_physi$SNP[h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP]
d_CNJ97_99_overlap_ordered <- d_CNJ97_99_overlap[match(selected_SNPs, d_CNJ97_99_overlap$SNP), ]
identical(h_map_physi[which(h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP),]$SNP,d_CNJ97_99_overlap_ordered$SNP)

h_map_physi[which(h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP),]$physical_map <- d_CNJ97_99_overlap_ordered$Pos
h_map_physi[which(h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP),]$Chrom <- d_CNJ97_99_overlap_ordered$Chrom
h_map_physi[which(h_map_physi$SNP %in% d_CNJ97_99_overlap$SNP),]$reposition <- d_CNJ97_99_overlap_ordered$Pos
write.csv(h_map_physi,"~/Desktop/Cranberry/20250214/LinkageMap/Result_map/consensus_map_information_reposition.csv",row.names = F)

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/combined_genomeVsmap2.png",height = 2500, width = 3000,res = 300)
ggplot(h_map_physi,aes((physical_map)/1e6,Map_pos,col=MarkerGroup))+
  geom_point(alpha = 0.7,size=3,pch=20)+
  scale_color_manual(values = friendly_pal("contrast_three"))+
  #scale_fill_distiller(palette = "RdBu")+
  labs(subtitle = "", x = "Genome position (Mbp)", y = "Map position (cM)") +
  facet_wrap(~factor(h_map_physi$LG,levels = paste0("lg",c(1:12))))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15), axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),axis.title.x = element_text(size = 15),
                 legend.position = "bottom", plot.subtitle = element_text(hjust = 0.5)) 

dev.off()
h_map_physi$chr <- stringr::str_split_fixed(h_map_physi$SNP,pattern = "\\_",n=2)[,1]
h_map_physi$LG <- factor(h_map_physi$LG,levels = paste0("lg",c(1:12)))
png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/combined_genomeVsmap3.png",height = 3500, width = 3500,res = 300)
ggplot(h_map_physi,aes((physical_map)/1e6,Map_pos,col=MarkerGroup))+
  geom_point(alpha = 0.7,size=3,pch=20)+
  scale_color_manual(values = friendly_pal("contrast_three"))+
  #scale_fill_distiller(palette = "RdBu")+
  labs(subtitle = "", x = "Genome position (Mbp)", y = "Map position (cM)") +
  facet_grid(LG ~ chr)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),
        legend.position = "bottom", plot.subtitle = element_text(hjust = 0.5)) 
dev.off()



## individual map

load("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ97_105_4_map.rda")
load("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/NJS98_18xCNJ99_9_96_map.rda")
d_97 <- plot_map_list(NJS98_18xCNJ97_105_4_map, type = "genome", parent = "p1p2",col = mp_pal(12))
d_97$chr <- stringr::str_split_fixed(d_97$mrk,pattern = "\\_",n=2)[,1]
d_97$physical_map <- as.numeric(stringr::str_split_fixed(d_97$mrk,pattern = "\\_",n=2)[,2])
d_99 <- plot_map_list(NJS98_18xCNJ99_9_96_map, type = "genome", parent = "p1p2",col = mp_pal(12))
d_99$chr <- stringr::str_split_fixed(d_99$mrk,pattern = "\\_",n=2)[,1]
d_99$physical_map <- as.numeric(stringr::str_split_fixed(d_99$mrk,pattern = "\\_",n=2)[,2])
d_99$LG <- factor(d_99$LG,levels = paste0("lg",c(1:12)))

png("~/Desktop/Cranberry/20250214/LinkageMap/Result_map/MapVsPhysical_CNJ99_2.png",height = 3500, width = 3800,res = 300)
ggplot(d_99,aes((physical_map)/1e6,pos,col=LG))+
  geom_point(alpha = 0.7,size=3,pch=20)+
 # scale_color_manual(values = friendly_pal("contrast_three"))+
  #scale_fill_distiller(palette = "RdBu")+
  labs(subtitle = "", x = "Genome position (Mbp)", y = "Map position (cM)") +
  facet_grid(LG ~ chr)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),
        legend.position = "none", plot.subtitle = element_text(hjust = 0.5)) 
dev.off()

