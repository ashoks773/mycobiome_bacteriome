# Cohesion.R script

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

#---- Load Mountain Gorilla Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/Moutain_Check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_mountainG <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_mountainG) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load BaAka Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/BaAka_check_New.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_BaAka <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_BaAka) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load Bantu Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/Bantu_check_New.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_Bantu <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_Bantu) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load Mangabey Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/Mangabey_check_New.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_mangabey <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_mangabey) <- c("Positive Cohesion", "Neagtive Cohesion")
#write.table(Cohesion, file="Network_Stats/For_Cohesion/MoutainG_Cohesion.txt", sep="\t")

#---- Load Chimps Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/Chimps_Check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_Chimps <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_Chimps) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load CaptiveChimps Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/CaptiveOstrava_check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_ChimpsOstrava <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_ChimpsOstrava) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load Captive WLG Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/CaptiveWLG_Check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_captiveWLG <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_captiveWLG) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load USA Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/USA_Check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_USA <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_USA) <- c("Positive Cohesion", "Neagtive Cohesion")

#---- Load WLG Network data
load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/WLG_Check.RData")
# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(ccrepe_bac_its$sim.score, 2, pos.mean)
connectedness.neg <- apply(ccrepe_bac_its$sim.score, 2, neg.mean)
# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- (as.matrix(ITSasv_table_relab) %*% connectedness.pos)/100
cohesion.neg <- as.matrix(ITSasv_table_relab) %*% connectedness.neg/100
cohesion.pos <- data.frame(cohesion.pos)
cohesion.neg <- data.frame(cohesion.neg)
cohesin_WLG <- cbind (cohesion.pos, cohesion.neg)
colnames(cohesin_WLG) <- c("Positive Cohesion", "Neagtive Cohesion")

Cohesin_All <- rbind(cohesin_BaAka, cohesin_Bantu, cohesin_ChimpsOstrava, cohesin_captiveWLG, cohesin_Chimps, cohesin_mangabey, cohesin_mountainG, cohesin_USA, cohesin_WLG)


#--- Merge Metadata
meta <- read.csv("~/Work/Project_16S_ITS/Metadata_Filtered.txt", sep = "\t", row.names = 1)
meta <- data.frame(t(meta))
meta <- data.frame(t(meta))
meta <- subset(meta, Group != "Captive Chimps-Hodonin")

Cohesin_All_meta <- merge(Cohesin_All, meta, by=0, all=F)
rownames(Cohesin_All_meta) <- Cohesin_All_meta$Row.names; Cohesin_All_meta$Row.names <- NULL
write.table (Cohesin_All_meta, file="Network_Stats/For_Cohesion/Cohesion_all_meta_Final.txt", sep="\t")

boxplot(Cohesin_All_meta$`Positive Cohesion` ~ Cohesin_All_meta$Group)
boxplot(Cohesin_All_meta$`Neagtive Cohesion` ~ Cohesin_All_meta$Group)

Colors = c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
#--- Positive Cohesion
PC <- Cohesin_All_meta[,c(1,7)]
PC_melted <- melt(PC, id.vars = "Group")
jpeg("Network_Stats/For_Cohesion/PositiveCohesion.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = PC_melted, aes(x=reorder(Group, value, FUN=mean), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("") + labs(x="",y="Positive Cohesion") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
kruskalmc(PC$`Positive Cohesion` ~ PC$Group)

#--- Negative Cohesion
NC <- Cohesin_All_meta[,c(2,7)]
NC_melted <- melt(NC, id.vars = "Group")
jpeg("Network_Stats/For_Cohesion/NegativeCohesion.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = NC_melted, aes(x=reorder(Group, value, FUN=mean), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("") + 
  labs(x="",y="Negative Cohesion") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + 
  ylim(0.0, -0.25) +
  theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
kruskalmc(NC$`Neagtive Cohesion` ~ NC$Group)

#--- Not using this
#--- Ratio of Negative to Positive Cohesion
ratio_neg_vs_pos <- data.frame (Cohesin_All_meta$`Positive Cohesion`/Cohesin_All_meta$`Neagtive Cohesion`)
row.names(ratio_neg_vs_pos) <- row.names(Cohesin_All_meta)
ratio_neg_vs_pos_meta <- merge(ratio_neg_vs_pos, meta, by=0, all=F)
rownames(ratio_neg_vs_pos_meta) <- ratio_neg_vs_pos_meta$Row.names; ratio_neg_vs_pos_meta$Row.names <- NULL
