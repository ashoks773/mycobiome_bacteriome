#---------------------------------------------------------------------------
#---- FOr Figure3 - Inter-Individual variations in MycoBiome and Bacteriome
#---------------------------------------------------------------------------
BAC_bray_dist <- vegdist(BAC_species_relab, method = "bray", binary = TRUE)
ITS_bray_dist <- vegdist(ITS_species_relab, method = "bray", binary = TRUE)

Bac_dist <- betadisper(BAC_bray_dist, BAC_meta$Group, type = 'centroid')
Bac_dist_new <- Bac_dist$distances
Bac_dist_new <- data.frame(Bac_dist_new)

ITS_dist <- betadisper(ITS_bray_dist, ITS_meta$Group, type = 'centroid')
ITS_dist_new <- ITS_dist$distances
ITS_dist_new <- data.frame(ITS_dist_new)

Bac_ITS_dist <- merge(Bac_dist_new, ITS_dist_new, by=0, all=F)
rownames(Bac_ITS_dist) <- Bac_ITS_dist$Row.names; Bac_ITS_dist$Row.names <- NULL

Bac_ITS_dist_Group <- merge(Bac_ITS_dist, ITS_meta, by=0, all=F)
rownames(Bac_ITS_dist_Group) <- Bac_ITS_dist_Group$Row.names; Bac_ITS_dist_Group$Row.names <- NULL
#write.table (Bac_ITS_dist_Group, file="Bac_ITS_dist_Group.txt", sep="\t")

Dist <- read.csv("Bac_ITS_dist_Group_formatted.txt", sep = "\t", header=T)

Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")
ggplot(Dist, aes(x=Group, y=Dist, color = Group, shape=Mico)) +
  geom_boxplot() + scale_color_manual(values=Colors) + geom_jitter(position=position_jitter(0.2)) +
  scale_shape_manual(values=c(21,23))

#-----
BaAka <- subset (Dist, Group == "BaAka-Human")
jpeg("BaAka.jpg", height = 3, width = 2, units = 'in', res = 600)
BaAka %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='forestgreen') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(BaAka$Dist ~ BaAka$Mico)

Bantu <- subset (Dist, Group == "Bantu-Human")
jpeg("Bantu.jpg", height = 3, width = 2, units = 'in', res = 600)
Bantu %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkgray') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(Bantu$Dist ~ Bantu$Mico)


CCH <- subset (Dist, Group == "Captive Chimps-Hodonin")
jpeg("Hodonin.jpg", height = 3, width = 2, units = 'in', res = 600)
CCH %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkgoldenrod3') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(CCH$Dist ~ CCH$Mico)

CCO <- subset (Dist, Group == "Captive Chimps-Ostrava")
jpeg("Ostrava.jpg", height = 3, width = 2, units = 'in', res = 600)
CCO %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='chartreuse1') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(CCO$Dist ~ CCO$Mico)

CWLG <- subset (Dist, Group == "Captive Western Lowland Gorilla")
jpeg("CaptiveWLG.jpg", height = 3, width = 2, units = 'in', res = 600)
CWLG %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='purple3') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(CWLG$Dist ~ CWLG$Mico)

Chimps <- subset (Dist, Group == "Chimps")
jpeg("Chimps.jpg", height = 3, width = 2, units = 'in', res = 600)
Chimps %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkmagenta') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(Chimps$Dist ~ Chimps$Mico)

Mangabey <- subset (Dist, Group == "Mangabey")
jpeg("Mangabey.jpg", height = 3, width = 2, units = 'in', res = 600)
Mangabey %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='blue') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(Mangabey$Dist ~ Mangabey$Mico)

MG <- subset (Dist, Group == "Mountain Gorilla")
jpeg("Mountain.jpg", height = 3, width = 2, units = 'in', res = 600)
MG %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkcyan') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(MG$Dist ~ MG$Mico)

USA <- subset (Dist, Group == "USA-Human")
jpeg("USA-Human.jpg", height = 3, width = 2, units = 'in', res = 600)
USA %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkorchid1') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(USA$Dist ~ USA$Mico)

WLG <- subset (Dist, Group == "Western Lowland Gorilla")
jpeg("WLG.jpg", height = 3, width = 2, units = 'in', res = 600)
WLG %>%
  # ggplot(aes(Mico,Dist, fill=Mico)) +
  ggplot(aes(Mico, Dist, shape=Mico)) +
  geom_boxplot(fill='darkred') + scale_shape_manual(values=c(21,23)) +
  #scale_color_manual(values=Colors) +
  geom_line(aes(group=Paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Group,group=Paired), size=3, position = position_dodge(0.2)) +
  theme(legend.position = "none") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(colour="black", size=18)) + 
  theme(legend.title = element_blank()) + theme(legend.position='none') + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), axis.text.y = element_text(size = 18, colour = "black"))
dev.off ()
wilcox.test(WLG$Dist ~ WLG$Mico)


#----- To calculate Fold changes
Distances <- read.csv(file="Bac_ITS_dist_Group.txt", sep="\t", row.names = 1, header = T)

MeanBac_Distances <- aggregate(Distances$Bac_dist_new, by = list(Group =Distances$Group), FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
MeanFun_Distances <- aggregate(Distances$ITS_dist_new, by = list(Group =Distances$Group), FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

FC <- MeanFun_Distances$x/MeanBac_Distances$x
FC <- data.frame(FC)
row.names(FC) <- MeanBac_Distances$Group
