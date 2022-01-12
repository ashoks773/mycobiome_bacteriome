#--- This Scritp is to repeat Analysis after removal of CaptiveChimps-Hodonin
# - Generate All plots again
setwd("~/Box/Gomez_Lab/Project_16S_ITS/Manuscript_Final/NPJ_Biofilms_Microbiomes/NPJ_Submitted/Reviewer3")

library (vegan)
library (ape)
library (pgirmess)
library (labdsv)
library (ggplot2)
library (psych)
library (randomForest)
library (cluster)
library (ade4)
library (pROC)
library(tidyverse)
library(gapminder)
library (phyloseq)

###########################################################################
#---- Figure 1: Alpha-diversity Analysis and PCoA on ITS and 16S seperately
###########################################################################

#------Read in Fungal ASV Table table
asv_table_in <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
#asv_table_in <- read.csv("../ITS_analysis/L6_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))

#-- Additional steps to filter OTUs
ITSasv.sums <- colSums(asv_table_in_t)
ITSasv_filtered <- asv_table_in_t[ , which(ITSasv.sums > 10)] #---Check ASV sum should be more than 10
ITSasv_filtered_filtered <- dropspc(ITSasv_filtered, 5) #--- ASV should be present in more than 5 samples
ITSasv_table <- ITSasv_filtered_filtered[order(rownames(ITSasv_filtered_filtered)),]

#---- New metadata - elephants removed (total 18), less than 1k depht removed (total 18) and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]
ITSmeta <- subset(ITSmeta, Group != "Captive Chimps-Hodonin") # Remove Three Captive Chimps Hodonin Samples: X1FECOB, X2FecJUDY, X4FecTEA

ITS.Obs <- rowSums(ITSasv_table > 0)
ITS.shannon <- diversity(ITSasv_table)
ITS.InvSimpson <- diversity(ITSasv_table, "invsimpson")
ITS.Obs.rare <- rarefy(ITSasv_table, min(rowSums(ITSasv_table)))
ITS_total_diversity <- data.frame(ITS.Obs, ITS.Obs.rare, ITS.shannon, ITS.InvSimpson)

ITS_total_diversity_meta <- merge(ITS_total_diversity, ITSmeta, by=0, all=F)
rownames(ITS_total_diversity_meta) <- ITS_total_diversity_meta$Row.names; ITS_total_diversity_meta$Row.names <- NULL

its_stats <- kruskalmc(ITS_total_diversity_meta$ITS.Obs ~ ITS_total_diversity_meta$Group)
#write.table (its_stats, file = "ITS_Krukshal.txt", sep = "\t")

its_shannon_stats <- kruskalmc(ITS_total_diversity_meta$ITS.shannon ~ ITS_total_diversity_meta$Group)
#write.table (its_shannon_stats, file = "ITS_Krukshal_Shannon.txt", sep = "\t")

#Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred") #- Remove Hodonin Color

#####################################
#----Figure 1A and S1A: Plot Fungal diversity
#####################################
ITS_Observed <- ITS_total_diversity_meta[,c(1,9)]
ITS_Observed_melted <- melt(ITS_Observed, id.vars = "Group")
jpeg("Figure1A.jpg", height = 4, width = 4, units = 'in', res = 600)
#ggplot(Observed_at1k, aes(x=reorder(Group, Observed, FUN=median), y=Observed, color=Group, fill = Group), alpha = 0.1) + geom_boxplot() + ggtitle("Observed") + labs(x="",y="Alpha Diversity Measure") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
ggplot(data = ITS_Observed_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Observed ASVs (Fungi)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
ITS_Shannon <- ITS_total_diversity_meta[,c(3,9)]
ITS_Shannon_melted <- melt(ITS_Shannon, id.vars = "Group")
jpeg("FigureS1.jpg", height = 4, width = 4, units = 'in', res = 600)
#ggplot(Observed_at1k, aes(x=reorder(Group, Observed, FUN=median), y=Observed, color=Group, fill = Group), alpha = 0.1) + geom_boxplot() + ggtitle("Observed") + labs(x="",y="Alpha Diversity Measure") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
ggplot(data = ITS_Shannon_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Shannon (Fungi)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()

############################################
#---Figure 1B: Plot Beta diversity PCoA plot
############################################
ITS_ASV_table_meta <- merge(ITSasv_table, ITSmeta, by=0, all=F)
rownames(ITS_ASV_table_meta) <- ITS_ASV_table_meta$Row.names; ITS_ASV_table_meta$Row.names <- NULL

ITS_species <- ITS_ASV_table_meta[,1:582]
ITS_meta <- ITS_ASV_table_meta[,583:588]

ITS_species_relab <- decostand(ITS_species, method = "total")*100
#write.table (ITS_species_relab, file = "ITS_ASV_filtered_Rel_Abun.txt", sep = "\t")
ITS_bray_dist <- vegdist(ITS_species_relab, method = "bray", binary = TRUE)
ITS_bray_pcoa <- pcoa (ITS_bray_dist)
ITS_bray_pcoa$values[1:2,]
ITS_pc1 <- round(ITS_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
ITS_pc2 <- round(ITS_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
ITS_Bray_PCoA_MATRIX <- ITS_bray_pcoa$vectors[,1:2]
ITS_Bray_PCoA_MATRIX <- data.frame(ITS_Bray_PCoA_MATRIX)
ITS_Bray_PCoA_MATRIX_New <- cbind(ITS_Bray_PCoA_MATRIX, ITS_meta)

#Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

#Shapes <- c(16, 17, 18, 15)

#jpeg("For_colors.jpg", height = 4, width = 6.5, units = 'in', res = 600)
jpeg("Figure1B.jpg", height = 4, width = 6.5, units = 'in', res = 600)
#ggplot(ITS_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", ITS_pc1, "%", sep="")) + ylab(paste("PCo2 - ", ITS_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
ggplot(ITS_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values=c(1,2,4,5,6,7,8,9,10)) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", ITS_pc1, "%", sep="")) + ylab(paste("PCo2 - ", ITS_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(ITS_bray_dist ~ ITS_meta$Group)
#ITS_meta$Group   9    28.518  3.1687   12.59 0.43697  0.001 ***

#############################################################################
#-----Figure 1B Additional: Phylogenetic tree using cumulative ITS ASV abundance
#############################################################################
ITS_relab_meta <- cbind (ITS_meta, ITS_species_relab)
#write.table(ITS_relab_meta, file= "ITS_relab_meta.txt", sep= "\t")

Average_ASV_abundance <- read.csv(file="~/Box/Gomez_Lab/Project_16S_ITS/Combined/Cumulative_ITS_rel_abun.txt", sep="\t", row.names = 1, header=T)
row_names_df_to_remove<-"Captive Chimps-Hodonin" #-- Remove CaptiveChimps Hodonin
Average_ASV_abundance <- Average_ASV_abundance[!(row.names(Average_ASV_abundance) %in% row_names_df_to_remove),]


bray_dist <-vegdist(Average_ASV_abundance, "bray")
hc <- hclust(bray_dist, method="average")
png(filename="Figure1B_histogram.png", height = 2, width = 2.5, units = 'in', res = 300)
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)
#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()

#----- To check Inter-Individual variations
#ITS_bray_dist <- vegdist(ITS_species_relab, method = "bray", binary = TRUE)
#groups <- ITS_meta$Group
#mod3 <- betadisper(ITS_bray_dist, groups, type = "centroid")
#plot(mod3)
#boxplot(mod3)

#---- Unsupervised Clustering
dist_bray <- as.matrix (ITS_bray_dist)
jpeg("FigureS2a.jpg", height = 3.5, width = 3.5, units = 'in', res = 600)
fviz_nbclust(dist_bray, pam, method = "wss") + geom_vline(xintercept = 4, linetype = 2)
dev.off ()

pam <- pam(dist_bray, k=4)
clusters <- pam$clustering
clusters <- data.frame(clusters)
Bray_PCoA_MATRIX_Cluster <- cbind(ITS_Bray_PCoA_MATRIX, clusters)
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "1"] <- "Cluster1"
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "2"] <- "Cluster2"
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "3"] <- "Cluster3"
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "4"] <- "Cluster4"
#Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "5"] <- "Cluster5"

jpeg("FigureS2b.jpg", height = 3.5, width = 3.5, units = 'in', res = 600)
s.class(Bray_PCoA_MATRIX_Cluster[,1:2], fac= as.factor(Bray_PCoA_MATRIX_Cluster$clusters), col = c("darkcyan", "darkorchid1", "red", "darkgreen", "darkgray"), label = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5"))
dev.off ()
Bray_PCoA_MATRIX_Cluster_Meta <- cbind (Bray_PCoA_MATRIX_Cluster, ITS_meta)
write.table (Bray_PCoA_MATRIX_Cluster_Meta, file ="Bray_Clusters_Group.txt", sep = "\t")

############################################
#---Figure 1C: Fungi Bubble plot
############################################
bubble_ITS_reviewer.R

############################################
#---Figure S2: Fungi Taxonomy plot
############################################
Fungi_Taxonomy_plots.R

#---------------------------------- 16S analysis
############################################
#---- For Figure S7A, B, and C
############################################
#--- Read in Bacterial ASV Table table
asv_table_in <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
#asv_table_in <- read.csv("../16S_analysis/L6_feature.table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))
#ab <- table(unlist(asv_table_in_t))

BACasv.sums <- colSums(asv_table_in_t)
BACasv_filtered <- asv_table_in_t[ , which(BACasv.sums > 10)] #---Check ASV sum should be more than 10
BACasv_filtered_filtered <- dropspc(BACasv_filtered, 5) #--- ASV should be present in more than 5 samples
BACasv_table <- BACasv_filtered_filtered[order(rownames(BACasv_filtered_filtered)),]

BACmeta <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/Metadata_Filtered.txt", sep = "\t", row.names = 1)
BACmeta_t <- data.frame(t(BACmeta))
BACmeta <- data.frame(t(BACmeta_t))
BACmeta <- BACmeta[order(rownames(BACmeta)),]
BACmeta <- subset(BACmeta, Group != "Captive Chimps-Hodonin") # Remove Three Captive Chimps Hodonin Samples: X1FECOB, X2FecJUDY, X4FecTEA

BAC.Obs <- rowSums(BACasv_table > 0)
BAC.shannon <- diversity(BACasv_table)
BAC.InvSimpson <- diversity(BACasv_table, "invsimpson")
BAC.Obs.rare <- rarefy(BACasv_table, min(rowSums(BACasv_table)))
BAC_total_diversity <- data.frame(BAC.Obs, BAC.Obs.rare, BAC.shannon, BAC.InvSimpson)

BAC_total_diversity_meta <- merge(BAC_total_diversity, BACmeta, by=0, all=F)
rownames(BAC_total_diversity_meta) <- BAC_total_diversity_meta$Row.names; BAC_total_diversity_meta$Row.names <- NULL

bac_stats <- kruskalmc(BAC_total_diversity_meta$BAC.Obs ~ BAC_total_diversity_meta$Group)
#write.table (bac_stats, file = "Bac_Krukshal.txt", sep = "\t")

bac_shannon_stats <- kruskalmc(BAC_total_diversity_meta$BAC.shannon ~ BAC_total_diversity_meta$Group)
#write.table (bac_shannon_stats, file = "Bac_Krukshal_Shannon.txt", sep = "\t")

Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred") #- Remove Hodonin Color

#--- PLOT Bacterial diversity
BAC_Observed <- BAC_total_diversity_meta[,c(1,9)]
BAC_Observed_melted <- melt(BAC_Observed, id.vars = "Group")
jpeg("FigureS7_obs.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = BAC_Observed_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Observed ASVs (Bacteria)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
BAC_Shannon <- BAC_total_diversity_meta[,c(3,9)]
BAC_Shannon_melted <- melt(BAC_Shannon, id.vars = "Group")
jpeg("FigureS7_Shannon.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = BAC_Shannon_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Shannon (Bacteria)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()

#--- Plot Beta diversity PCoA plot
BAC_ASV_table_meta <- merge(BACasv_table, BACmeta, by=0, all=F)
rownames(BAC_ASV_table_meta) <- BAC_ASV_table_meta$Row.names; BAC_ASV_table_meta$Row.names <- NULL

BAC_species <- BAC_ASV_table_meta[,1:2756]
BAC_meta <- BAC_ASV_table_meta[,2757:2762]

BAC_species_relab <- decostand(BAC_species, method = "total")*100
#write.table (BAC_species_relab, file = "16S_ASV_filtered_Rel_Abun.txt", sep = "\t")
BAC_bray_dist <- vegdist(BAC_species_relab, method = "bray", binary = TRUE)
BAC_bray_pcoa <- pcoa (BAC_bray_dist)
BAC_bray_pcoa$values[1:2,]
BAC_pc1 <- round(BAC_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
BAC_pc2 <- round(BAC_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
BAC_Bray_PCoA_MATRIX <- BAC_bray_pcoa$vectors[,1:2]
BAC_Bray_PCoA_MATRIX <- data.frame(BAC_Bray_PCoA_MATRIX)
BAC_Bray_PCoA_MATRIX_New <- cbind(BAC_Bray_PCoA_MATRIX, BAC_meta)

Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

jpeg("FigureS7C.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(BAC_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values=c(1,2,4,5,6,7,8,9,10)) + 
  scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", BAC_pc1, "%", sep="")) + ylab(paste("PCo2 - ", BAC_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + ylim(-0.50, 0.30) +
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()


adonis2(BAC_bray_dist ~ BAC_meta$Group)
#BAC_meta$Group   9    33.361  3.7068  29.567 0.64572  0.001 ***

#---- Combined Alpha diversity
Total_Fungal_Bacterial_Diversity <- data.frame(ITS_total_diversity_meta, BAC_total_diversity_meta)
#write.table (Total_Fungal_Bacterial_Diversity, file = "Total_Diversity.txt", sep = "\t")


###############################################
#---- Figure 2A: Diversity correlation plot
###############################################
#------- 
jpeg("FigureS8B.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(data = Total_Fungal_Bacterial_Diversity, mapping = aes(x = ITS.Obs, y = BAC.Obs)) + 
  geom_point(mapping = aes(color = Group, shape = Group)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values=c(1,2,4,5,6,7,8,9,10)) + theme_bw() +
  xlab("Observed ASVs (Fungi)") + ylab("Observed ASVs (Bacteria)") +
  geom_smooth()
dev.off ()
Correlation <- corr.test(Total_Fungal_Bacterial_Diversity$BAC.Obs, Total_Fungal_Bacterial_Diversity$ITS.Obs, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=0.36, p=4.775164e-06, and se=0.07


#######################################################
#---Figure 2B: Qiime Procrustes Plot generated in Qiime 
#######################################################
#qiime diversity procrustes-analysis --i-reference 16s/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --i-other ITS/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --p-dimensions 2 --output-dir Procrustes_results
#qiime emperor procrustes-plot --i-reference-pcoa 16s/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --i-other-pcoa ITS/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --m-metadata-file Metadata_filtered.tsv --o-visualization procrustes.qzv --p-ignore-missing-samples

#------#------ Procustes analysis on rarefied data to check the scores and significance
setwd("~/Work/Project_16S_ITS/Combined/Procrustes")
ITS_dist_bray_qiime <- read.csv("ITS_Bray_distance-matrix.tsv", sep = "\t", row.names = 1)
BAC_dist_bray_qiime <- read.csv("16S_Bray_distance-matrix.tsv", sep = "\t", row.names = 1)

procrustes_results <- procrustes(BAC_dist_bray_qiime, ITS_dist_bray_qiime, symmetric = TRUE)
#Procrustes sum of squares: 0.4056 
protest (BAC_dist_bray_qiime, ITS_dist_bray_qiime, permutations = 999)
#Procrustes Sum of Squares (m12 squared):        0.388 
#Correlation in a symmetric Procrustes rotation: 0.7823  
#Significance:  0.001 
mantel(BAC_dist_bray_qiime, ITS_dist_bray_qiime, permutations = 999)
#Mantel statistic r: 0.6266 
#Significance: 0.001

#https://rdrr.io/cran/Evomorph/man/ShapeDist.html
#A matrix containing procrustes distance between shapes. Procrustes distance is the square root of the sum of squared differences in the posititions of the landmarks in two shapes
#--- To calculate Distances between Bacteria and Fungal Composition
library("Evomorph")
Distances <- ShapeDist(shapes = ITS_dist_bray_qiime, reference = BAC_dist_bray_qiime)
Distances <- data.frame(Distances)
row <- row.names(BAC_dist_bray_qiime)
row <- data.frame(row)
Distances_new <- cbind (row, Distances)
row.names(Distances_new) <- Distances_new$row

ITSmeta <- read.csv("../../../Metadata_Filtered.txt", sep = "\t", row.names = 1)

Distances_new_meta <- merge(Distances_new, ITSmeta, by=0, all=F)
rownames(Distances_new_meta) <- Distances_new_meta$Row.names; Distances_new_meta$Row.names <- NULL

##############################################################
# Figure 2C ------- Fungi PCo1 vs Bacteria PCo1
##############################################################
ITS_pc1_BAC_pc1 <- merge(ITS_Bray_PCoA_MATRIX, BAC_Bray_PCoA_MATRIX, by=0, all=F)
rownames(ITS_pc1_BAC_pc1) <- ITS_pc1_BAC_pc1$Row.names; ITS_pc1_BAC_pc1$Row.names <- NULL

ITS_pc1_BAC_pc1_Meta <- merge(ITS_pc1_BAC_pc1, ITS_meta, by=0, all=F)
rownames(ITS_pc1_BAC_pc1_Meta) <- ITS_pc1_BAC_pc1_Meta$Row.names; ITS_pc1_BAC_pc1_Meta$Row.names <- NULL

jpeg("FigureS8A.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(data = ITS_pc1_BAC_pc1_Meta, mapping = aes(x = Axis.1.x, y = Axis.1.y)) + 
  geom_point(mapping = aes(color = Group, shape = Group)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values=c(1,2,4,5,6,7,8,9,10)) + theme_bw() +
  xlab("PCo1- Fungi") + ylab("PCo1- Bacteria") +
  geom_smooth()
dev.off ()

Correlation <- corr.test(ITS_pc1_BAC_pc1_Meta$Axis.1.x, ITS_pc1_BAC_pc1_Meta$Axis.1.y, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")

##############################################################
# Figure 2D ------- Merged to ASV tables and use for PCoA Plot
##############################################################
Table_ITS <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
Table_ITS <- data.frame(t(Table_ITS))
Table_16S <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
Table_16S <- data.frame(t(Table_16S))

Combined_ASVs <- data.frame(Table_16S, Table_ITS)

#-- Additional steps to filter OTUs
combinedASV.sums <- colSums(Combined_ASVs)
CombinedASV_filtered <- Combined_ASVs[ , which(combinedASV.sums > 10)] #---Check ASV sum should be more than 10
CombinedASV_filtered_filtered <- dropspc(CombinedASV_filtered, 5) #--- ASV should be present in more than 5 samples
Combined_ASV_table <- CombinedASV_filtered_filtered[order(rownames(CombinedASV_filtered_filtered)),]

#---- New metadata - elephants removed and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("~/Box/Gomez_Lab/Project_16S_ITS/Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]
ITSmeta <- subset(ITSmeta, Group != "Captive Chimps-Hodonin") # Remove Three Captive Chimps Hodonin Samples: X1FECOB, X2FecJUDY, X4FecTEA

Combined_ASV_table_meta <- merge(Combined_ASV_table, ITSmeta, by=0, all=F)
rownames(Combined_ASV_table_meta) <- Combined_ASV_table_meta$Row.names; Combined_ASV_table_meta$Row.names <- NULL

species <- Combined_ASV_table_meta[,1:3338]
meta <- Combined_ASV_table_meta[,3339:3344]

#species_relab_new <- cbind(BAC_species_relab, ITS_species_relab) #Alternate if want to use
species_relab <- decostand(species, method = "total")*100
bray_dist <- vegdist(species_relab, method = "bray", binary = TRUE)
bray_pcoa <- pcoa (bray_dist)
bray_pcoa$values[1:2,]
pc1 <- round(bray_pcoa$values$Rel_corr_eig[1]*100, 2)
pc2 <- round(bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, meta)

#Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

jpeg("Figure3A.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values=c(1,2,4,5,6,7,8,9,10)) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", pc1, "%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(bray_dist ~ meta$Group)
#meta$Group   9    33.277  3.6974  27.559 0.62947  0.001 ***

#############################################################################
#-----Figure 2D Additional: Phylogenetic tree using cumulative ASV abundance
#############################################################################
#species_relab_meta <- cbind (meta, species_relab)
#write.table(species_relab_meta, file= "species_relative.txt", sep= "\t")

Cumulative_ITS_16S_ASV_abundance <- read.csv(file="~/Box/Gomez_Lab/Project_16S_ITS/Combined/Cumulative_combined_relabun.txt", sep="\t", row.names = 1, header=T)

row_names_df_to_remove<-"Captive Chimps-Hodonin" #-- Remove CaptiveChimps Hodonin
Cumulative_ITS_16S_ASV_abundance <- Cumulative_ITS_16S_ASV_abundance[!(row.names(Cumulative_ITS_16S_ASV_abundance) %in% row_names_df_to_remove),]

bray_dist <-vegdist(Cumulative_ITS_16S_ASV_abundance, "bray")
hc <- hclust(bray_dist)
png(filename="Figure3A_histogram.png", height = 2, width = 2.5, units = 'in', res = 300)
#clus3 = cutree(hc, 4)
#plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE, tip.color = colors[clus3])
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)
#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()

#########################################################
#---Figure S4--- : Bacteria Bubble and taxonomy plot
#########################################################
Bubble_16S.R
Bcateria_Taxonomy_plots.R


##############********************************
#--------#---- Figure 2 - To make Taxa heatmap
###############################################
#-- To make heatmap - Short Genus names top30 - formatted Ordered - by phylogeny as ITS
top_30genus_reorder <- read.csv(file="Sel_30_genus_rel_foramtted_Reordered.txt", sep = "\t", row.names = 1, header =T)
top_30genus_reorder <- subset(top_30genus_reorder, Group != "Captive Chimps-Hodonin") #- Remove captive chimps Hodonin
categories <- top_30genus_reorder[,1]
top_30genus_reorder_pred <- top_30genus_reorder[,2:31]

#Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

aheatmap(sqrt(sqrt(top_30genus_reorder_pred)), color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "spearman", hclustfun = "complete", fontsize=10,  filename="Figure2_Top30_Genus_Reordered.png", scale = "row", Rowv = NA, Colv = NA, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annRow = categories)
