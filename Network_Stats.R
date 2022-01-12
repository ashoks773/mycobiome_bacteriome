##############################################################################
#------------- Statistical Analysis on Network statistics- For Figure3B
##############################################################################
Colors = c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

setwd("~/Work/Project_16S_ITS/Manuscript_Final/NPJ_Biofilms_Microbiomes/NPJ_Submitted/Network_ana/Network_Stats")
Cytoscape_stats <- read.csv(file="All_directed_Network.txt", sep="\t",header = T, row.names = 1)

library(ggplot2)

##############################################
#########Let's make figure for Network Stress 
##############################################
NC <- Cytoscape_stats[,c(10,17)]
NC_melted <- melt(NC, id.vars = "Group")
jpeg("Network_Stats/NeighborhoodConnectivity.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = NC_melted, aes(x=reorder(Group, value, FUN=mean), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("") + labs(x="",y="Neighborhood Connectivity") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
Kruskalmc(NC$NeighborhoodConnectivity ~ NC$Group)

degree <- Cytoscape_stats[,c(11,17)]
degree_melted <- melt(degree, id.vars = "Group")
jpeg("Network_Stats/Degree.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = degree_melted, aes(x=reorder(Group, value, FUN=mean), y=value, fill=Group)) + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("") + labs(x="",y="Degree") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black"))
dev.off ()
kruskalmc(degree$Outdegree ~ degree$Group)

#---- Load Modularity Data - Calcualted in Network Attribute Steps
Colors = c("forestgreen", "darkgray", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

setwd("~/Work/Project_16S_ITS/Manuscript_Final/NPJ_Biofilms_Microbiomes/NPJ_Submitted/Network_ana/Network_Stats")
modularity <- read.csv(file="~/Work/Project_16S_ITS/Manuscript_Final/NPJ_Biofilms_Microbiomes/NPJ_Submitted/Network_ana/Modularity.txt", sep="\t",header = T, row.names = 1)

jpeg("Network_Stats/modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(modularity$Modularity, col=Colors, pch=1, lwd = 5)#, xlim=c(-0.5, 12), ylim = c(0.00, 1))
#text (modularity$Modularity, labels = rownames(modularity), pos = 3, cex = 0.7)
dev.off ()

#------------------------------------------------
#---- Number of Correlations vs Number of Samples
Numbers <- read.csv(file="~/Work/Project_16S_ITS/Manuscript_Final/NPJ_Biofilms_Microbiomes/NPJ_Submitted/Network_ana/Numberof_corr_numberof_Samples.txt", sep="\t",header = T, row.names = 1)
jpeg("Network_Stats/NumberofCorrelations_no_of_Samples.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(Numbers$Total.correaltions ~ Numbers$Number.of.Samples, col=Colors, pch=1, lwd = 5)
#model<-glm (Numbers$Total.correaltions ~ Numbers$Number.of.Samples)
cor.test(x=Numbers$Total.correaltions, y=Numbers$Number.of.Samples, method="spearman")
#abline(model)
dev.off ()

