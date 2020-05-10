##############################################################################
#------------- Statistical Analysis on Network statistics- For Figure3B
##############################################################################
Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

setwd("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Correlations/CytoScapeStatistics")
Cytoscape_stats <- read.csv(file="CytoScapeStats-For_AllGroups.txt", sep="\t",header = T)

library(ggplot2)

##############################################
#########Let's make figure for Network Stress 
##############################################
stress <- aggregate(sqrt(Cytoscape_stats$Stress), by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_stress <- do.call(data.frame, stress)
myData_stress$se <- myData_stress$x.sd / sqrt(myData_stress$x.n)
colnames(myData_stress) <- c("Group", "median", "sd", "n", "se")
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = myData_stress$median + myData_stress$se, ymin = myData_stress$median - myData_stress$se)

jpeg("Network_Stress_ordered.jpg", height = 4.8, width = 7, units = 'in', res = 600)
p <- ggplot(data = myData_stress, aes(x=reorder(Group, median, FUN=mean), y = median, fill= Group, colour=Group)) + 
#p <- ggplot(data = myData, aes(x = Group, y = mean, fill= Group, colour=Group)) + 
 scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) +
  labs(x="",y="median (Network Stress)")
p + geom_bar(stat = "identity", position = dodge) +
geom_errorbar(limits, position = dodge, width = 0.25, colour="black") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
axis.title.x=element_blank()) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
stress_stats <- kruskalmc(Cytoscape_stats$Stress ~ Cytoscape_stats$Group)
write.table(stress_stats, file= "stress_Stats.txt", sep="\t")


#################################################
###Let's make figure for Number of direct edges 
#################################################
edges <- aggregate(Cytoscape_stats$NumberOfDirectedEdges, by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_edges <- do.call(data.frame, edges)
myData_edges$se <- myData_edges$x.sd / sqrt(myData_edges$x.n)
colnames(myData_edges) <- c("Group", "median", "sd", "n", "se")
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = myData_edges$median + myData_edges$se, ymin = myData_edges$median - myData_edges$se)

jpeg("NumberofDirectEdges_Odered.jpg", height = 4.8, width = 7, units = 'in', res = 600)
p <- ggplot(data = myData_edges, aes(x=reorder(Group, median, FUN=mean), y = median, fill= Group, colour=Group)) + 
  scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) +
  labs(x="",y="median (Number of Directed Edges)")
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25, colour="black") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
directEdges_stats <- kruskalmc(Cytoscape_stats$NumberOfDirectedEdges ~ Cytoscape_stats$Group)
write.table(directEdges_stats, file= "directEdges_Stats.txt", sep="\t")

#################################################
###Let's make figure for NeighborhoodConnectivity 
#################################################
connectivity <- aggregate(Cytoscape_stats$NeighborhoodConnectivity, by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_connectivity <- do.call(data.frame, connectivity)
myData_connectivity$se <- myData_connectivity$x.sd / sqrt(myData_connectivity$x.n)
colnames(myData_connectivity) <- c("Group", "median", "sd", "n", "se")
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = myData_connectivity$median + myData_connectivity$se, ymin = myData_connectivity$median - myData_connectivity$se)

jpeg("NeighborhoodConnectivity_Ordered.jpg", height = 4.8, width = 7, units = 'in', res = 600)
p <- ggplot(data = myData_connectivity, aes(x=reorder(Group, median, FUN=mean), y = median, fill= Group, colour=Group)) + 
  scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) +
  labs(x="",y="median (Neighborhood Connectivity)")
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25, colour="black") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
NeighborhoodConnectivity <- kruskalmc(Cytoscape_stats$NeighborhoodConnectivity ~ Cytoscape_stats$Group)
write.table(NeighborhoodConnectivity, file= "NeighborhoodConnectivity_Stats.txt", sep="\t")

#--- Note --- If reordered plot shoud be made
p <- ggplot(data = myData, aes(x =reorder(Group, mean, FUN=median), y = mean, fill= Group, colour=Group)) + 
  scale_color_manual(values=Colors) + scale_fill_manual(values=Colors)
