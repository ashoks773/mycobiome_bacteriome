library (vegan)
library (ape)
library (pgirmess)
library (labdsv)
library (ggplot2)
library (psych)
library (randomForest)


###########################################################################
#---- Figure 1: Alpha-diversity Analysis and PCoA on ITS and 16S seperately
###########################################################################

#------Read in Fungal ASV Table table
asv_table_in <- read.csv("../ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
#asv_table_in <- read.csv("../ITS_analysis/L6_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))

#-- Additional steps to filter OTUs
ITSasv.sums <- colSums(asv_table_in_t)
ITSasv_filtered <- asv_table_in_t[ , which(ITSasv.sums > 10)] #---Check ASV sum should be more than 10
ITSasv_filtered_filtered <- dropspc(ITSasv_filtered, 5) #--- ASV should be present in more than 5 samples
ITSasv_table <- ITSasv_filtered_filtered[order(rownames(ITSasv_filtered_filtered)),]

#---- New metadata - elephants removed (total 18), less than 1k depht removed (total 18) and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("../Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]

ITS.Obs <- rowSums(ITSasv_table > 0)
ITS.shannon <- diversity(ITSasv_table)
ITS.InvSimpson <- diversity(ITSasv_table, "invsimpson")
ITS.Obs.rare <- rarefy(ITSasv_table, min(rowSums(ITSasv_table)))
ITS_total_diversity <- data.frame(ITS.Obs, ITS.Obs.rare, ITS.shannon, ITS.InvSimpson)

ITS_total_diversity_meta <- merge(ITS_total_diversity, ITSmeta, by=0, all=F)
rownames(ITS_total_diversity_meta) <- ITS_total_diversity_meta$Row.names; ITS_total_diversity_meta$Row.names <- NULL

its_stats <- kruskalmc(ITS_total_diversity_meta$ITS.Obs ~ ITS_total_diversity_meta$Group)
write.table (its_stats, file = "ITS_Krukshal.txt", sep = "\t")

its_shannon_stats <- kruskalmc(ITS_total_diversity_meta$ITS.shannon ~ ITS_total_diversity_meta$Group)
write.table (its_shannon_stats, file = "ITS_Krukshal_Shannon.txt", sep = "\t")

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

#####################################
#----Figure 1A and S1A: Plot Fungal diversity
#####################################
ITS_Observed <- ITS_total_diversity_meta[,c(1,9)]
ITS_Observed_melted <- melt(ITS_Observed, id.vars = "Group")
jpeg("ITS_Observed.jpg", height = 4, width = 4, units = 'in', res = 600)
#ggplot(Observed_at1k, aes(x=reorder(Group, Observed, FUN=median), y=Observed, color=Group, fill = Group), alpha = 0.1) + geom_boxplot() + ggtitle("Observed") + labs(x="",y="Alpha Diversity Measure") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
ggplot(data = ITS_Observed_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Observed ASVs (Fungi)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
ITS_Shannon <- ITS_total_diversity_meta[,c(3,9)]
ITS_Shannon_melted <- melt(ITS_Shannon, id.vars = "Group")
jpeg("ITS_Shannon.jpg", height = 4, width = 4, units = 'in', res = 600)
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
write.table (ITS_species_relab, file = "ITS_ASV_filtered_Rel_Abun.txt", sep = "\t")
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

jpeg("ITS_Bray_PCoA_Unweighted.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(ITS_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", ITS_pc1, "%", sep="")) + ylab(paste("PCo2 - ", ITS_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis(ITS_bray_dist ~ ITS_meta$Group)
#ITS_meta$Group   9    28.518  3.1687   12.59 0.43697  0.001 ***

#############################################################################
#-----Figure 1B Additional: Phylogenetic tree using cumulative ITS ASV abundance
#############################################################################
ITS_relab_meta <- cbind (ITS_meta, ITS_species_relab)
write.table(ITS_relab_meta, file= "ITS_relab_meta.txt", sep= "\t")

Average_ASV_abundance <- read.csv(file="Cumulative_ITS_rel_abun.txt", sep="\t", row.names = 1, header=T)
bray_dist <-vegdist(Average_ASV_abundance, "bray")
hc <- hclust(bray_dist)
png(filename="ITS_Dendogram_withOutColor.png", height = 2, width = 2.5, units = 'in', res = 300)
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)
#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()

bray_dist <-vegdist(Average_ASV_abundance, "mountford")
hc <- hclust(bray_dist)
png(filename="Check_ITS_Dendogram.png", height = 2, width = 1.5, units = 'in', res = 300)
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)
#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()

############################################
#---Figure 1C: Fungi Bubble plot
############################################
Bubble_ITS.R

#################################################
#-- To find out significantly discriminating Taxa
#################################################
ITS_Family_proportions <- read.csv(file="ITS_family_proportions.txt",sep="\t",header=TRUE,row.names=1)
ITS_Family_proportions <- data.frame(t(ITS_Family_proportions))

ITS_Family_proportions_meta <- merge(ITS_Family_proportions, ITSmeta, by=0, all=F)
rownames(ITS_Family_proportions_meta) <- ITS_Family_proportions_meta$Row.names; ITS_Family_proportions_meta$Row.names <- NULL

ITS_Family_proportions_meta_Group <- ITS_Family_proportions_meta[c(1:323,328)]
#--- Indval Analysis
library (labdsv)
ITS_Family_proportions_meta_Group_filtered <- ITS_Family_proportions_meta_Group[ , which(!apply(ITS_Family_proportions_meta_Group==0,2,all))]
ITS_Family_proportions_meta_Group_filtered[is.na(ITS_Family_proportions_meta_Group_filtered)] <- 0

iva <- indval(ITS_Family_proportions_meta_Group_filtered[,1:319], ITS_Family_proportions_meta_Group_filtered$Group)

gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(ITS_Family_proportions_meta_Group_filtered[,1:319]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="Family_proportions_indvalsummary.txt", sep = "\t")

#################################################################################
#---Supplementary Analysis: PD Whole tree, UniFrac for Fungal, and Taxa Summary
#################################################################################
#-- Make a phyloseq object after filtering of ASVs
# Read in OTU table
#------Read in Fungal ASV Table table
asv_table_in <- read.csv("../ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))

#-- Additional steps to filter OTUs
ITSasv.sums <- colSums(asv_table_in_t)
ITSasv_filtered <- asv_table_in_t[ , which(ITSasv.sums > 10)] #---Check ASV sum should be more than 10
ITSasv_filtered_filtered <- dropspc(ITSasv_filtered, 5) #--- ASV should be present in more than 5 samples
ITSasv_filtered_filtered_t <- data.frame(t(ITSasv_filtered_filtered))
asv_table_in <- as.matrix(ITSasv_filtered_filtered_t)

#---- New metadata - elephants removed (total 18), less than 1k depht removed (total 18) and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
# -- Sample H100 and H9 were not used as they were giving NA values while calculating UniFrac
ITSmeta <- read.table("../Metadata_Filtered_Only_for_Phyloseq_function.txt", sep = "\t", row.names = 1, header=T)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("../ITS_analysis/ITS_taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy_t <- data.frame(t(taxonomy))
taxonomy <- data.frame(t(taxonomy_t))
taxonomy <- as.matrix(taxonomy)

# Read in tree
phy_tree <- read_tree("../ITS_analysis/ITS_tree.nwk")

# Import all as phyloseq objects
ASV <- otu_table(asv_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(ITSmeta)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(ASV)
taxa_names(phy_tree)
# Same sample names
sample_names(ASV)
sample_names(META)

# Finally merge to create Phyloseq object!
ITS_ps <- phyloseq(ASV, TAX, META, phy_tree)

##########################
#----FigureS - Not Included
###########################

source('~/Work/Project_16S_ITS/Combined/estimate_pd.R')
Faith_pd <- estimate_pd(ITS_ps)
Group <- ps@sam_data[,5:6]
Group <- data.frame(Group)
Faith_pd_Group <- cbind (Faith_pd, Group)
Faith_pd_Group <- na.omit(Faith_pd_Group) # Two samples does not have Calculated PD values
write.table (Faith_pd_Group, file="ITS_faith_pd_metadata.txt", sep = "\t")

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

ITS_faithPD <- Faith_pd_Group[,c(1,3)]
ITS_faithPD_melted <- melt(ITS_faithPD, id.vars = "Group")
jpeg("ITS_FaithPD.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = ITS_faithPD_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Faith PD (Fungi)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
kruskalmc (Faith_pd_Group$PD ~ Faith_pd_Group$Group)

####################################
#Figure S1B: Beta diversity analysis 
####################################
#A) Create Phyloseq object on filtered OTU
ITS_ps_normalized  = transform_sample_counts(ITS_ps, function(x) x / sum(x) )
ordUni = ordinate(ITS_ps_normalized, "PCoA", "unifrac", weighted=FALSE)

pc <- round (ordUni$values$Rel_corr_eig, 4)*100
ITS_Uni_PCoA_MATRIX <- ordUni$vectors[,1:2]
ITS_Uni_PCoA_MATRIX <- data.frame(ITS_Uni_PCoA_MATRIX)

ITS_Uni_PCoA_MATRIX_New <- merge(ITS_Uni_PCoA_MATRIX, ITSmeta, by=0, all=F)
rownames(ITS_Uni_PCoA_MATRIX_New) <- ITS_Uni_PCoA_MATRIX_New$Row.names; ITS_Uni_PCoA_MATRIX_New$Row.names <- NULL

jpeg("ITS_UniFrac_PCoA.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(ITS_Uni_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", pc[1], "%", sep="")) + ylab(paste("PCo2 - ", pc[2], "%", sep="")) + ggtitle("UniFrac") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

pcoa1 <- ITS_Uni_PCoA_MATRIX_New[,c(1,7)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group")
jpeg("ITS_Uni_PCoA1_Distances.jpg", height = 2.5, width = 4, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa1_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCo1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) + coord_flip() + theme(legend.position='none')
dev.off ()
pcoa2 <- ITS_Uni_PCoA_MATRIX_New[,c(2,7)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group")
jpeg("ITS_Uni_PCoA2_Distances.jpg", height = 4, width = 2.5, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa2_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

ITS_UniFrac_distances <- UniFrac(ITS_ps_normalized, weighted=FALSE)
adonis(ITS_UniFrac_distances ~ ITS_ps_normalized@sam_data$Group)

############################################
#---Figure S2: Fungi Taxonomy plot
############################################
Fungi_Taxonomy_plots.R

#---------------------------------- 16S analysis
############################################
#---- For Figure 2: Analysis of 16S Data
############################################
#--- Read in Bacterial ASV Table table
asv_table_in <- read.csv("../16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
#asv_table_in <- read.csv("../16S_analysis/L6_feature.table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))
#ab <- table(unlist(asv_table_in_t))

BACasv.sums <- colSums(asv_table_in_t)
BACasv_filtered <- asv_table_in_t[ , which(BACasv.sums > 10)] #---Check ASV sum should be more than 10
BACasv_filtered_filtered <- dropspc(BACasv_filtered, 5) #--- ASV should be present in more than 5 samples
BACasv_table <- BACasv_filtered_filtered[order(rownames(BACasv_filtered_filtered)),]

BACmeta <- read.csv("../Metadata_Filtered.txt", sep = "\t", row.names = 1)
BACmeta_t <- data.frame(t(BACmeta))
BACmeta <- data.frame(t(BACmeta_t))
BACmeta <- BACmeta[order(rownames(BACmeta)),]

BAC.Obs <- rowSums(BACasv_table > 0)
BAC.shannon <- diversity(BACasv_table)
BAC.InvSimpson <- diversity(BACasv_table, "invsimpson")
BAC.Obs.rare <- rarefy(BACasv_table, min(rowSums(BACasv_table)))
BAC_total_diversity <- data.frame(BAC.Obs, BAC.Obs.rare, BAC.shannon, BAC.InvSimpson)

BAC_total_diversity_meta <- merge(BAC_total_diversity, BACmeta, by=0, all=F)
rownames(BAC_total_diversity_meta) <- BAC_total_diversity_meta$Row.names; BAC_total_diversity_meta$Row.names <- NULL

bac_stats <- kruskalmc(BAC_total_diversity_meta$BAC.Obs ~ BAC_total_diversity_meta$Group)
write.table (bac_stats, file = "Bac_Krukshal.txt", sep = "\t")

bac_shannon_stats <- kruskalmc(BAC_total_diversity_meta$BAC.shannon ~ BAC_total_diversity_meta$Group)
write.table (bac_shannon_stats, file = "Bac_Krukshal_Shannon.txt", sep = "\t")

#--- PLOT Bacterial diversity
BAC_Observed <- BAC_total_diversity_meta[,c(1,9)]
BAC_Observed_melted <- melt(BAC_Observed, id.vars = "Group")
jpeg("BAC_Observed.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = BAC_Observed_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Observed ASVs (Bacteria)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
BAC_Shannon <- BAC_total_diversity_meta[,c(3,9)]
BAC_Shannon_melted <- melt(BAC_Shannon, id.vars = "Group")
jpeg("BAC_Shannon.jpg", height = 4, width = 4, units = 'in', res = 600)
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

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

jpeg("BAC_Bray_PCoA_Unweighted.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(BAC_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", BAC_pc1, "%", sep="")) + ylab(paste("PCo2 - ", BAC_pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()


adonis(BAC_bray_dist ~ BAC_meta$Group)
#BAC_meta$Group   9    33.361  3.7068  29.567 0.64572  0.001 ***

#---- Combined Alpha diversity
Total_Fungal_Bacterial_Diversity <- data.frame(ITS_total_diversity_meta, BAC_total_diversity_meta)
write.table (Total_Fungal_Bacterial_Diversity, file = "Total_Diversity.txt", sep = "\t")

###############################################
#---- Figure 2A: Diversity correlation plot
###############################################
#------- 
jpeg("Bacteria_Fungi_Observed.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(data = Total_Fungal_Bacterial_Diversity, mapping = aes(x = ITS.Obs, y = BAC.Obs)) + 
  geom_point(mapping = aes(color = Group, shape = Group)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values =1:20) + theme_bw() +
  xlab("Observed ASVs (Fungi)") + ylab("Observed ASVs (Bacteria)") +
  geom_smooth()
dev.off ()
Correlation <- corr.test(Total_Fungal_Bacterial_Diversity$BAC.Obs, Total_Fungal_Bacterial_Diversity$ITS.Obs, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=0.36, p=4.775164e-06, and se=0.07

jpeg("Bacteria_Fungi_Shannon.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(data = Total_Fungal_Bacterial_Diversity, mapping = aes(x = ITS.shannon, y = BAC.shannon)) +
  geom_point(mapping = aes(color = Group, shape = Group)) + scale_color_manual(values=Colors) + 
  scale_y_continuous(limits = c(3, 5.5)) + #This will remove 4 samples having very low Shannon (Bacteria)
  scale_shape_manual(values =1:20) + theme_bw() +
  xlab("Shannon (Fungi)") + ylab("Shannon (Bacteria)") +
  geom_smooth()
dev.off ()
Correlation_shann <- corr.test(Total_Fungal_Bacterial_Diversity$BAC.shannon, Total_Fungal_Bacterial_Diversity$ITS.shannon, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=0.30, p=8.761658e-05, and se=0.07


#######################################################
#---Figure 2B: Qiime Procrustes Plot generated in Qiime 
#######################################################
#qiime diversity procrustes-analysis --i-reference 16s/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --i-other ITS/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --p-dimensions 2 --output-dir Procrustes_results
#qiime emperor procrustes-plot --i-reference-pcoa 16s/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --i-other-pcoa ITS/core-metrics-results_rare1k/bray_curtis_pcoa_results.qza --m-metadata-file Metadata_filtered.tsv --o-visualization procrustes.qzv --p-ignore-missing-samples

#------#------ Procustes on rarefied data to check the scores and significance
#--- Bray-curtis distance matrices for this analysis
setwd("~/Work/Project_16S_ITS/Combined/Rarefied_1k_for_alpha_beta-diversity/Procrustes")
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


##############################################################
# Figure 2C ------- Fungi PCo1 vs Bacteria PCo1
##############################################################
ITS_pc1_BAC_pc1 <- merge(ITS_Bray_PCoA_MATRIX, BAC_Bray_PCoA_MATRIX, by=0, all=F)
rownames(ITS_pc1_BAC_pc1) <- ITS_pc1_BAC_pc1$Row.names; ITS_pc1_BAC_pc1$Row.names <- NULL

ITS_pc1_BAC_pc1_Meta <- merge(ITS_pc1_BAC_pc1, ITS_meta, by=0, all=F)
rownames(ITS_pc1_BAC_pc1_Meta) <- ITS_pc1_BAC_pc1_Meta$Row.names; ITS_pc1_BAC_pc1_Meta$Row.names <- NULL

jpeg("FungiPC1_vs_BacteriaPC1.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(data = ITS_pc1_BAC_pc1_Meta, mapping = aes(x = Axis.1.x, y = Axis.1.y)) + 
  geom_point(mapping = aes(color = Group, shape = Group)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values =1:20) + theme_bw() +
  xlab("PCo1- Fungi") + ylab("PCo1- Bacteria") +
  geom_smooth()
dev.off ()

Correlation <- corr.test(ITS_pc1_BAC_pc1_Meta$Axis.1.x, ITS_pc1_BAC_pc1_Meta$Axis.1.y, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=0.5229396, p=4.804338e-12

##############################################################
# Figure 2D ------- Merged to ASV tables and use for PCoA Plot
##############################################################
Table_ITS <- read.csv("../ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
Table_ITS <- data.frame(t(Table_ITS))
Table_16S <- read.csv("../16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
Table_16S <- data.frame(t(Table_16S))

Combined_ASVs <- data.frame(Table_16S, Table_ITS)

#-- Additional steps to filter OTUs
combinedASV.sums <- colSums(Combined_ASVs)
CombinedASV_filtered <- Combined_ASVs[ , which(combinedASV.sums > 10)] #---Check ASV sum should be more than 10
CombinedASV_filtered_filtered <- dropspc(CombinedASV_filtered, 5) #--- ASV should be present in more than 5 samples
Combined_ASV_table <- CombinedASV_filtered_filtered[order(rownames(CombinedASV_filtered_filtered)),]

#---- New metadata - elephants removed and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("../Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]

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

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

jpeg("Combined_Bray_PCoA_Unweighted.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", pc1, "%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis(bray_dist ~ meta$Group)
#meta$Group   9    33.277  3.6974  27.559 0.62947  0.001 ***

#############################################################################
#-----Figure 2D Additional: Phylogenetic tree using cumulative ASV abundance
#############################################################################
species_relab_meta <- cbind (meta, species_relab)
write.table(species_relab_meta, file= "species_relative.txt", sep= "\t")

Cumulative_ITS_16S_ASV_abundance <- read.csv(file="Cumulative_combined_relabun.txt", sep="\t", row.names = 1, header=T)
bray_dist <-vegdist(Cumulative_ITS_16S_ASV_abundance, "bray")
hc <- hclust(bray_dist)
png(filename="Cumulative_ITS_16S_Dendogram_withOutColor.png", height = 2, width = 2.5, units = 'in', res = 300)
#colors = c("gray35", "purple3", "darkmagenta", "cadetblue2")
#colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
#clus3 = cutree(hc, 4)
#plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE, tip.color = colors[clus3])
plot(as.phylo(hc), type = "phylogram", cex = 0.5, no.margin = TRUE)

#plot(as.phylo(hc), type = "unrooted", cex = 0.5, no.margin = TRUE)
dev.off ()


###################################################################
#---Supplementary Analysis: PD Whole tree and UniFrac for Bacteria
##################################################################
#-- Make a phyloseq object after filtering of ASVs
# Read in OTU table
#------Read in Bacterial ASV Table table
asv_table_in <- read.csv("../16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))

BACasv.sums <- colSums(asv_table_in_t)
BACasv_filtered <- asv_table_in_t[ , which(BACasv.sums > 10)] #---Check ASV sum should be more than 10
BACasv_filtered_filtered <- dropspc(BACasv_filtered, 5) #--- ASV should be present in more than 5 samples
BACasv_filtered_filtered_t <- data.frame(t(BACasv_filtered_filtered))
asv_table_in <- as.matrix(BACasv_filtered_filtered_t)

BACmeta <- read.csv("../Metadata_Filtered.txt", sep = "\t", row.names = 1)
BACmeta_t <- data.frame(t(BACmeta))
BACmeta <- data.frame(t(BACmeta_t))

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv(file="../16S_analysis/16S_taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy_t <- data.frame(t(taxonomy))
taxonomy <- data.frame(t(taxonomy_t))
taxonomy <- as.matrix(taxonomy)

# Read in tree
phy_tree <- read_tree("../16S_analysis/16S_tree.nwk")

# Import all as phyloseq objects
ASV <- otu_table(asv_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(BACmeta)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(ASV)
taxa_names(phy_tree)
# Same sample names
sample_names(ASV)
sample_names(META)

# Finally merge to create Phyloseq object!
ps <- phyloseq(ASV, TAX, META, phy_tree)

#############################
#----FigureS - Not inlcuded
############################

source('~/Work/Project_16S_ITS/Combined/estimate_pd.R')
Faith_pd <- estimate_pd(ps)
Group <- ps@sam_data[,5:6]
Group <- data.frame(Group)
Faith_pd_Group <- cbind (Faith_pd, Group)
Faith_pd_Group <- na.omit(Faith_pd_Group) # Two samples does not have Calculated PD values
write.table (Faith_pd_Group, file="Bac_faith_pd_metadata.txt", sep = "\t")

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

Bac_faithPD <- Faith_pd_Group[,c(1,3)]
Bac_faithPD_melted <- melt(Bac_faithPD, id.vars = "Group")
jpeg("BAC_FaithPD.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Bac_faithPD_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Faith PD (Bacteria)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
kruskalmc (Faith_pd_Group$PD ~ Faith_pd_Group$Group)

####################################
#Figure S3B: Beta diversity analysis 
####################################
#A) Create Phyloseq object on filtered OTU
ps_normalized  = transform_sample_counts(ps, function(x) x / sum(x) )
ordUni = ordinate(ps_normalized, "PCoA", "unifrac", weighted=FALSE)

pc <- round (ordUni$values$Rel_corr_eig, 4)*100
BAC_Uni_PCoA_MATRIX <- ordUni$vectors[,1:2]
BAC_Uni_PCoA_MATRIX <- data.frame(BAC_Uni_PCoA_MATRIX)

BAC_Uni_PCoA_MATRIX_New <- merge(BAC_Uni_PCoA_MATRIX, BACmeta, by=0, all=F)
rownames(BAC_Uni_PCoA_MATRIX_New) <- BAC_Uni_PCoA_MATRIX_New$Row.names; BAC_Uni_PCoA_MATRIX_New$Row.names <- NULL

jpeg("BAC_UniFrac_PCoA.jpg", height = 4, width = 6.5, units = 'in', res = 600)
ggplot(BAC_Uni_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Group, colour=Group)) + geom_point(size=2) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", pc[1], "%", sep="")) + ylab(paste("PCo2 - ", pc[2], "%", sep="")) + ggtitle("UniFrac") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

pcoa1 <- BAC_Uni_PCoA_MATRIX_New[,c(1,7)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group")
jpeg("BAC_Uni_PCoA1_Distances.jpg", height = 2.5, width = 4, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa1_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCo1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) + coord_flip() + theme(legend.position='none')
dev.off ()
pcoa2 <- BAC_Uni_PCoA_MATRIX_New[,c(2,7)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group")
jpeg("BAC_Uni_PCoA2_Distances.jpg", height = 4, width = 2.5, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa2_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

UniFrac_distances <- UniFrac(ps_normalized, weighted=FALSE)
adonis(UniFrac_distances ~ ps_normalized@sam_data$Group)

#########################################################
#---Figure S4--- : Bacteria Bubble and taxonomy plot
#########################################################
Bubble_16S.R
Bcateria_Taxonomy_plots.R
