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

#######################################################################################
#---Gorilla Repetative samples were removed to check the patterns as reviewer Asked for
#######################################################################################
#Jab (one extra sample) - we removed Jab2
#Sam2 (one extra sample) - we removed Samson
#SCH2 (two extra samples) we removed Scho2
#-- These Four samples were removed from the Metadata File

#------Read in Fungal ASV Table table
asv_table_in <- read.csv("~/Work/Project_16S_ITS/ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
#asv_table_in <- read.csv("../ITS_analysis/L6_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))


#-- Additional steps to filter OTUs
ITSasv.sums <- colSums(asv_table_in_t)
ITSasv_filtered <- asv_table_in_t[ , which(ITSasv.sums > 10)] #---Check ASV sum should be more than 10
ITSasv_filtered_filtered <- dropspc(ITSasv_filtered, 5) #--- ASV should be present in more than 5 samples
ITSasv_table <- ITSasv_filtered_filtered[order(rownames(ITSasv_filtered_filtered)),]

#---- New metadata - elephants removed (total 18), less than 1k depht removed (total 18) and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("Metadata_Gorilla_removed.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]


Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

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

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
#Shapes <- c(16, 17, 18, 15)

#jpeg("For_colors.jpg", height = 4, width = 6.5, units = 'in', res = 600)
jpeg("ITS_Bray_PCoA_CaptiveWLG_remvoed.jpg", height = 4, width = 6.5, units = 'in', res = 600)
#ggplot(ITS_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", ITS_pc1, "%", sep="")) + ylab(paste("PCo2 - ", ITS_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
ggplot(ITS_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", ITS_pc1, "%", sep="")) + ylab(paste("PCo2 - ", ITS_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis(ITS_bray_dist ~ ITS_meta$Group)
#ITS_meta$Group   9    28.518  3.1687   12.59 0.43697  0.001 ***

#############################################################################
#-----Figure 1B Additional: Phylogenetic tree using cumulative ITS ASV abundance
#############################################################################
ITS_relab_meta <- cbind (ITS_meta, ITS_species_relab)
#write.table(ITS_relab_meta, file= "ITS_relab_meta.txt", sep= "\t")

Average_ASV_abundance <- read.csv(file="Cumulative_ITS_rel_abun.txt", sep="\t", row.names = 1, header=T)
bray_dist <-vegdist(Average_ASV_abundance, "bray")
hc <- hclust(bray_dist, method="average")
png(filename="ITS_Dendogram_withOutColor_Try.png", height = 2, width = 2.5, units = 'in', res = 300)
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)
#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()


