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
library(Evomorph)
library (reshape2)

#---- Using this Tutorial
#http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-Section-01-Main-Lab.html

#-- For All (These distances both ITS_bray_dist and BAC_bray_dist are coming from Figure1-2_Alpha-Beta R Script)
ITS_bray_dist <- as.matrix(ITS_bray_dist)
BAC_bray_dist <- as.matrix(BAC_bray_dist)
ITS_bray_dist <- data.frame(ITS_bray_dist)
BAC_bray_dist <- data.frame(BAC_bray_dist)

#https://rdrr.io/cran/Evomorph/man/ShapeDist.html
#Procrustes distance provides a measure of coincidence of two point sets xi and yi, i=1..N. For this purpose the variance of point deviations is calculated at the optimal superposition of the sets. It allows to characterize the shape proximity of a given simplex to shape of a reference one.
total_Distances <- ShapeDist(shapes = ITS_bray_dist, reference = BAC_bray_dist)

total_Distances <- data.frame(total_Distances)
row <- row.names(ITS_bray_dist)
row <- data.frame(row)
Distances_new <- cbind (row, total_Distances)
row.names(Distances_new) <- Distances_new$row

Distances_new_meta <- merge(Distances_new, ITSmeta, by=0, all=F)
rownames(Distances_new_meta) <- Distances_new_meta$Row.names; Distances_new_meta$Row.names <- NULL

#write.table(Distances_new_meta, file="Procrustes_Distances.txt", sep = "\t")
Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

procrustes_Dist <- Distances_new_meta[,c(2,7)]
procrustes_Dist_melted <- melt(procrustes_Dist, id.vars = "Group")

jpeg("Procrustes_Dist.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(data = procrustes_Dist_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Procrustes distance") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=12)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
dev.off ()

######################### ********************************
#--- For Reviewer 3 Remove Captive Hodonin and Plot this GAIN
Distances_new_meta <- read.csv(file="Procrustes_Distances.txt", sep = "\t", row.names = 1, header = T)
Distances_new_meta <- subset(Distances_new_meta, Group != "Captive Chimps-Hodonin")
#Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
Colors <- c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

procrustes_Dist <- Distances_new_meta[,c(2,7)]
procrustes_Dist_per = data.frame (procrustes_Dist$total_Distances/100) #-- Intially Distance scale was from 1-100/ I converted that from 0-1
procrustes_Dist_per_meta <- cbind (procrustes_Dist_per, procrustes_Dist$Group)
colnames(procrustes_Dist_per_meta) <- c("Procrustes distance", "Group")
procrustes_Dist_melted <- melt(procrustes_Dist_per_meta, id.vars = "Group")

jpeg("Figure3D.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(data = procrustes_Dist_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("") + labs(x="",y="Procrustes distance") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=12)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
dev.off ()

########
kruskalmc(Distances_new_meta$total_Distances ~ Distances_new_meta$Group)

#------ Trying Procrustes for Each Group
# https://stackoverflow.com/questions/30325739/ggplot2-for-procrustes-rotation-in-vegan
#---------- Procrustes- Analysis in Each Group
ITSmeta <- read.csv("../Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))

#--- For BaAka (Use Filtered ITS relative abunance table, and Bac Relative abundance Table)
baAka <- subset(ITSmeta, Group == "BaAka-Human")

baAka_ITS_meta <- merge(ITS_species_relab, baAka, by=0, all=F)
rownames(baAka_ITS_meta) <- baAka_ITS_meta$Row.names; baAka_ITS_meta$Row.names <- NULL
baAka_ITS_bray_dist <- vegdist(baAka_ITS_meta[,1:582], method = "bray", binary = TRUE)

baAka_Bac_meta <- merge(BAC_species_relab, baAka, by=0, all=F)
rownames(baAka_Bac_meta) <- baAka_Bac_meta$Row.names; baAka_Bac_meta$Row.names <- NULL
baAka_Bac_bray_dist <- vegdist(baAka_Bac_meta[,1:2756], method = "bray", binary = TRUE)

baAka_procrustes <- procrustes(baAka_ITS_bray_dist, baAka_Bac_bray_dist, symmetric = TRUE)

library(ggplot2)
library(grid)

ctest <- data.frame(rda1=baAka_procrustes$Yrot[,1],
                    rda2=baAka_procrustes$Yrot[,2],xrda1=baAka_procrustes$X[,1],
                    xrda2=baAka_procrustes$X[,2],Group=rep(c("BaAka-Human"),each=27))

Colors <- "forestgreen"
png(filename="BaAka.png", height = 1.8, width = 3.2, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()
protest (baAka_Bac_bray_dist, baAka_ITS_bray_dist, permutations = 999)
mantel(baAka_ITS_bray_dist, baAka_Bac_bray_dist, permutations = 999)
#Mantel statistic r: 0.2963, Significance: 0.037
# baAka_ITS_bray_dist <- as.matrix(baAka_ITS_bray_dist)
# baAka_Bac_bray_dist <- as.matrix(baAka_Bac_bray_dist)
# baAka_ITS_bray_dist <- data.frame(baAka_ITS_bray_dist)
# baAka_Bac_bray_dist <- data.frame(baAka_Bac_bray_dist)
# BaAka_Distances <- ShapeDist(shapes = baAka_Bac_bray_dist, reference = baAka_ITS_bray_dist)

#---- For Bantu
bantu <- subset(ITSmeta, Group == "Bantu-Human")

bantu_ITS_meta <- merge(ITS_species_relab, bantu, by=0, all=F)
rownames(bantu_ITS_meta) <- bantu_ITS_meta$Row.names; bantu_ITS_meta$Row.names <- NULL
bantu_ITS_bray_dist <- vegdist(bantu_ITS_meta[,1:582], method = "bray", binary = TRUE)

bantu_Bac_meta <- merge(BAC_species_relab, bantu, by=0, all=F)
rownames(bantu_Bac_meta) <- bantu_Bac_meta$Row.names; bantu_Bac_meta$Row.names <- NULL
bantu_Bac_bray_dist <- vegdist(bantu_Bac_meta[,1:2756], method = "bray", binary = TRUE)

bantu_procrustes <- procrustes(bantu_ITS_bray_dist, bantu_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=bantu_procrustes$Yrot[,1],
                    rda2=bantu_procrustes$Yrot[,2],xrda1=bantu_procrustes$X[,1],
                    xrda2=bantu_procrustes$X[,2],Group=rep(c("Bantu-Human"),each=13))

Colors <- "darkgray"
png(filename="Bantu.png", height = 1.8, width = 3, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group),size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group),size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (bantu_Bac_bray_dist, bantu_ITS_bray_dist, permutations = 999)
mantel(bantu_ITS_bray_dist, bantu_Bac_bray_dist, permutations = 999)

# bantu_ITS_bray_dist <- as.matrix(bantu_ITS_bray_dist)
# bantu_Bac_bray_dist <- as.matrix(bantu_Bac_bray_dist)
# bantu_ITS_bray_dist <- data.frame(bantu_ITS_bray_dist)
# bantu_Bac_bray_dist <- data.frame(bantu_Bac_bray_dist)
# bantu_Distances <- ShapeDist(shapes = bantu_Bac_bray_dist, reference = bantu_ITS_bray_dist)

#---- For CCH
CCH <- subset(ITSmeta, Group == "Captive Chimps-Hodonin")

CCH_ITS_meta <- merge(ITS_species_relab, CCH, by=0, all=F)
rownames(CCH_ITS_meta) <- CCH_ITS_meta$Row.names; CCH_ITS_meta$Row.names <- NULL
CCH_ITS_bray_dist <- vegdist(CCH_ITS_meta[,1:582], method = "bray", binary = TRUE)

CCH_Bac_meta <- merge(BAC_species_relab, CCH, by=0, all=F)
rownames(CCH_Bac_meta) <- CCH_Bac_meta$Row.names; CCH_Bac_meta$Row.names <- NULL
CCH_Bac_bray_dist <- vegdist(CCH_Bac_meta[,1:2756], method = "bray", binary = TRUE)

CCH_procrustes <- procrustes(CCH_ITS_bray_dist, CCH_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=CCH_procrustes$Yrot[,1],
                    rda2=CCH_procrustes$Yrot[,2],xrda1=CCH_procrustes$X[,1],
                    xrda2=CCH_procrustes$X[,2],Group=rep(c("Captive Chimps-Hodonin"),each=3))

Colors <- "darkgoldenrod3"
png(filename="CCH.png", height = 1.6, width = 3.4, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group),size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (CCH_Bac_bray_dist, CCH_ITS_bray_dist, permutations = 999)
mantel(CCH_ITS_bray_dist, CCH_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# CCH_ITS_bray_dist <- as.matrix(CCH_ITS_bray_dist)
# CCH_Bac_bray_dist <- as.matrix(CCH_Bac_bray_dist)
# CCH_ITS_bray_dist <- data.frame(CCH_ITS_bray_dist)
# CCH_Bac_bray_dist <- data.frame(CCH_Bac_bray_dist)
# CCH_Distances <- ShapeDist(shapes = CCH_Bac_bray_dist, reference = CCH_ITS_bray_dist)

#---- For CCO
CCO <- subset(ITSmeta, Group == "Captive Chimps-Ostrava")

CCO_ITS_meta <- merge(ITS_species_relab, CCO, by=0, all=F)
rownames(CCO_ITS_meta) <- CCO_ITS_meta$Row.names; CCO_ITS_meta$Row.names <- NULL
CCO_ITS_bray_dist <- vegdist(CCO_ITS_meta[,1:582], method = "bray", binary = TRUE)

CCO_Bac_meta <- merge(BAC_species_relab, CCO, by=0, all=F)
rownames(CCO_Bac_meta) <- CCO_Bac_meta$Row.names; CCO_Bac_meta$Row.names <- NULL
CCO_Bac_bray_dist <- vegdist(CCO_Bac_meta[,1:2756], method = "bray", binary = TRUE)

CCO_procrustes <- procrustes(CCO_ITS_bray_dist, CCO_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=CCO_procrustes$Yrot[,1],
                    rda2=CCO_procrustes$Yrot[,2],xrda1=CCO_procrustes$X[,1],
                    xrda2=CCO_procrustes$X[,2],Group=rep(c("Captive Chimps-Ostrava"),each=12))

Colors <- "chartreuse1"
png(filename="CCO.png", height = 1.6, width = 3.4, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (CCO_Bac_bray_dist, CCO_ITS_bray_dist, permutations = 999)
mantel(CCO_ITS_bray_dist, CCO_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# CCO_ITS_bray_dist <- as.matrix(CCO_ITS_bray_dist)
# CCO_Bac_bray_dist <- as.matrix(CCO_Bac_bray_dist)
# CCO_ITS_bray_dist <- data.frame(CCO_ITS_bray_dist)
# CCO_Bac_bray_dist <- data.frame(CCO_Bac_bray_dist)
# CCO_Distances <- ShapeDist(shapes = CCO_Bac_bray_dist, reference = CCO_ITS_bray_dist)

#---- For CWLG
CWLG <- subset(ITSmeta, Group == "Captive Western Lowland Gorilla")

CWLG_ITS_meta <- merge(ITS_species_relab, CWLG, by=0, all=F)
rownames(CWLG_ITS_meta) <- CWLG_ITS_meta$Row.names; CWLG_ITS_meta$Row.names <- NULL
CWLG_ITS_bray_dist <- vegdist(CWLG_ITS_meta[,1:582], method = "bray", binary = TRUE)

CWLG_Bac_meta <- merge(BAC_species_relab, CWLG, by=0, all=F)
rownames(CWLG_Bac_meta) <- CWLG_Bac_meta$Row.names; CWLG_Bac_meta$Row.names <- NULL
CWLG_Bac_bray_dist <- vegdist(CWLG_Bac_meta[,1:2756], method = "bray", binary = TRUE)

CWLG_procrustes <- procrustes(CWLG_ITS_bray_dist, CWLG_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=CWLG_procrustes$Yrot[,1],
                    rda2=CWLG_procrustes$Yrot[,2],xrda1=CWLG_procrustes$X[,1],
                    xrda2=CWLG_procrustes$X[,2],Group=rep(c("Captive Western Lowland Gorilla"),each=18))

Colors <- "purple3"
png(filename="CWLG.png", height = 1.8, width = 4, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (CWLG_ITS_bray_dist, CWLG_Bac_bray_dist, permutations = 999)
protest (CWLG_Bac_bray_dist, CWLG_ITS_bray_dist, permutations = 999)
mantel(CWLG_ITS_bray_dist, CWLG_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# CWLG_ITS_bray_dist <- as.matrix(CWLG_ITS_bray_dist)
# CWLG_Bac_bray_dist <- as.matrix(CWLG_Bac_bray_dist)
# CWLG_ITS_bray_dist <- data.frame(CWLG_ITS_bray_dist)
# CWLG_Bac_bray_dist <- data.frame(CWLG_Bac_bray_dist)
# CWLG_Distances <- ShapeDist(shapes = CWLG_Bac_bray_dist, reference = CWLG_ITS_bray_dist)


#---- For Chimps
Chimps <- subset(ITSmeta, Group == "Mangabey")

Chimps_ITS_meta <- merge(ITS_species_relab, Chimps, by=0, all=F)
rownames(Chimps_ITS_meta) <- Chimps_ITS_meta$Row.names; Chimps_ITS_meta$Row.names <- NULL
Chimps_ITS_bray_dist <- vegdist(Chimps_ITS_meta[,1:582], method = "bray", binary = TRUE)

Chimps_Bac_meta <- merge(BAC_species_relab, Chimps, by=0, all=F)
rownames(Chimps_Bac_meta) <- Chimps_Bac_meta$Row.names; Chimps_Bac_meta$Row.names <- NULL
Chimps_Bac_bray_dist <- vegdist(Chimps_Bac_meta[,1:2756], method = "bray", binary = TRUE)

Chimps_procrustes <- procrustes(Chimps_ITS_bray_dist, Chimps_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=Chimps_procrustes$Yrot[,1],
                    rda2=Chimps_procrustes$Yrot[,2],xrda1=Chimps_procrustes$X[,1],
                    xrda2=Chimps_procrustes$X[,2],Group=rep(c("Chimps"),each=11))

Colors <- "darkmagenta"
png(filename="Chimps.png", height = 1.8, width = 2.5, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (Chimps_Bac_bray_dist, Chimps_ITS_bray_dist, permutations = 999)
mantel(Chimps_ITS_bray_dist, Chimps_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# Chimps_ITS_bray_dist <- as.matrix(Chimps_ITS_bray_dist)
# Chimps_Bac_bray_dist <- as.matrix(Chimps_Bac_bray_dist)
# Chimps_ITS_bray_dist <- data.frame(Chimps_ITS_bray_dist)
# Chimps_Bac_bray_dist <- data.frame(Chimps_Bac_bray_dist)
# Chimps_Distances <- ShapeDist(shapes = Chimps_Bac_bray_dist, reference = Chimps_ITS_bray_dist)


#---- For Mangabey
Mangabey <- subset(ITSmeta, Group == "Mangabey")

Mangabey_ITS_meta <- merge(ITS_species_relab, Mangabey, by=0, all=F)
rownames(Mangabey_ITS_meta) <- Mangabey_ITS_meta$Row.names; Mangabey_ITS_meta$Row.names <- NULL
Mangabey_ITS_bray_dist <- vegdist(Mangabey_ITS_meta[,1:582], method = "bray", binary = TRUE)

Mangabey_Bac_meta <- merge(BAC_species_relab, Mangabey, by=0, all=F)
rownames(Mangabey_Bac_meta) <- Mangabey_Bac_meta$Row.names; Mangabey_Bac_meta$Row.names <- NULL
Mangabey_Bac_bray_dist <- vegdist(Mangabey_Bac_meta[,1:2756], method = "bray", binary = TRUE)

Mangabey_procrustes <- procrustes(Mangabey_ITS_bray_dist, Mangabey_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=Mangabey_procrustes$Yrot[,1],
                    rda2=Mangabey_procrustes$Yrot[,2],xrda1=Mangabey_procrustes$X[,1],
                    xrda2=Mangabey_procrustes$X[,2],Group=rep(c("Mangabey"),each=11))

Colors <- "blue"
png(filename="Mangabey.png", height = 1.8, width = 2.7, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (Mangabey_Bac_bray_dist, Mangabey_ITS_bray_dist, permutations = 999)
mantel(Mangabey_ITS_bray_dist, Mangabey_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# Mangabey_ITS_bray_dist <- as.matrix(Mangabey_ITS_bray_dist)
# Mangabey_Bac_bray_dist <- as.matrix(Mangabey_Bac_bray_dist)
# Mangabey_ITS_bray_dist <- data.frame(Mangabey_ITS_bray_dist)
# Mangabey_Bac_bray_dist <- data.frame(Mangabey_Bac_bray_dist)
# Mangabey_Distances <- ShapeDist(shapes = Mangabey_Bac_bray_dist, reference = Mangabey_ITS_bray_dist)


#---- For MG
MG <- subset(ITSmeta, Group == "Mountain Gorilla")

MG_ITS_meta <- merge(ITS_species_relab, MG, by=0, all=F)
rownames(MG_ITS_meta) <- MG_ITS_meta$Row.names; MG_ITS_meta$Row.names <- NULL
MG_ITS_bray_dist <- vegdist(MG_ITS_meta[,1:582], method = "bray", binary = TRUE)

MG_Bac_meta <- merge(BAC_species_relab, MG, by=0, all=F)
rownames(MG_Bac_meta) <- MG_Bac_meta$Row.names; MG_Bac_meta$Row.names <- NULL
MG_Bac_bray_dist <- vegdist(MG_Bac_meta[,1:2756], method = "bray", binary = TRUE)

MG_procrustes <- procrustes(MG_ITS_bray_dist, MG_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=MG_procrustes$Yrot[,1],
                    rda2=MG_procrustes$Yrot[,2],xrda1=MG_procrustes$X[,1],
                    xrda2=MG_procrustes$X[,2],Group=rep(c("Mountain Gorilla"),each=26))

Colors <- "darkcyan"
png(filename="MG.png", height = 2, width = 3.2, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (MG_Bac_bray_dist, MG_ITS_bray_dist, permutations = 999)
mantel(MG_ITS_bray_dist, MG_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# MG_ITS_bray_dist <- as.matrix(MG_ITS_bray_dist)
# MG_Bac_bray_dist <- as.matrix(MG_Bac_bray_dist)
# MG_ITS_bray_dist <- data.frame(MG_ITS_bray_dist)
# MG_Bac_bray_dist <- data.frame(MG_Bac_bray_dist)
# MG_Distances <- ShapeDist(shapes = MG_Bac_bray_dist, reference = MG_ITS_bray_dist)


#---- For USA
USA <- subset(ITSmeta, Group == "USA-Human")

USA_ITS_meta <- merge(ITS_species_relab, USA, by=0, all=F)
rownames(USA_ITS_meta) <- USA_ITS_meta$Row.names; USA_ITS_meta$Row.names <- NULL
USA_ITS_bray_dist <- vegdist(USA_ITS_meta[,1:582], method = "bray", binary = TRUE)

USA_Bac_meta <- merge(BAC_species_relab, USA, by=0, all=F)
rownames(USA_Bac_meta) <- USA_Bac_meta$Row.names; USA_Bac_meta$Row.names <- NULL
USA_Bac_bray_dist <- vegdist(USA_Bac_meta[,1:2756], method = "bray", binary = TRUE)

USA_procrustes <- procrustes(USA_ITS_bray_dist, USA_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=USA_procrustes$Yrot[,1],
                    rda2=USA_procrustes$Yrot[,2],xrda1=USA_procrustes$X[,1],
                    xrda2=USA_procrustes$X[,2],Group=rep(c("USA-Human"),each=12))

Colors <- "darkorchid1"
png(filename="USA.png", height = 2, width = 3, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (USA_Bac_bray_dist, USA_ITS_bray_dist, permutations = 999)
mantel(USA_ITS_bray_dist, USA_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# USA_ITS_bray_dist <- as.matrix(USA_ITS_bray_dist)
# USA_Bac_bray_dist <- as.matrix(USA_Bac_bray_dist)
# USA_ITS_bray_dist <- data.frame(USA_ITS_bray_dist)
# USA_Bac_bray_dist <- data.frame(USA_Bac_bray_dist)
# USA_Distances <- ShapeDist(shapes = USA_Bac_bray_dist, reference = USA_ITS_bray_dist)

#---- For WLG
WLG <- subset(ITSmeta, Group == "Western Lowland Gorilla")

WLG_ITS_meta <- merge(ITS_species_relab, WLG, by=0, all=F)
rownames(WLG_ITS_meta) <- WLG_ITS_meta$Row.names; WLG_ITS_meta$Row.names <- NULL
WLG_ITS_bray_dist <- vegdist(WLG_ITS_meta[,1:582], method = "bray", binary = TRUE)

WLG_Bac_meta <- merge(BAC_species_relab, WLG, by=0, all=F)
rownames(WLG_Bac_meta) <- WLG_Bac_meta$Row.names; WLG_Bac_meta$Row.names <- NULL
WLG_Bac_bray_dist <- vegdist(WLG_Bac_meta[,1:2756], method = "bray", binary = TRUE)

WLG_procrustes <- procrustes(WLG_ITS_bray_dist, WLG_Bac_bray_dist, symmetric = TRUE)

ctest <- data.frame(rda1=WLG_procrustes$Yrot[,1],
                    rda2=WLG_procrustes$Yrot[,2],xrda1=WLG_procrustes$X[,1],
                    xrda2=WLG_procrustes$X[,2],Group=rep(c("Western Lowland Gorilla"),each=19))

Colors <- "darkred"
png(filename="WLG.png", height = 2, width = 3.8, units = 'in', res = 300)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, colour=Group), col= "#666666", shape=21, size=1.2) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Group), col= "#33CCFF", shape=22, size=1.2) +
  #geom_point(aes(x=rda1, y=rda2, colour=Group), size=1.2) +
  #geom_point(aes(x=xrda1, y=xrda2, colour=Group), size=1.2) +
  scale_color_manual(values=Colors) + theme_classic() +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Group),arrow=arrow(length=unit(0.15,"cm")))
dev.off ()

protest (WLG_Bac_bray_dist, WLG_ITS_bray_dist, permutations = 999)
mantel(WLG_ITS_bray_dist, WLG_Bac_bray_dist, permutations = 999)
#Mantel statistic r: -0.4475, Significance: 0.66667
# WLG_ITS_bray_dist <- as.matrix(WLG_ITS_bray_dist)
# WLG_Bac_bray_dist <- as.matrix(WLG_Bac_bray_dist)
# WLG_ITS_bray_dist <- data.frame(WLG_ITS_bray_dist)
# WLG_Bac_bray_dist <- data.frame(WLG_Bac_bray_dist)
# WLG_Distances <- ShapeDist(shapes = WLG_Bac_bray_dist, reference = WLG_ITS_bray_dist)
