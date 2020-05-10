setwd("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Correlations/CytoScapeStatistics")
Cytoscape_stats <- read.csv(file="CytoScapeStats-For_AllGroups.txt", sep="\t",header = T)

library(ggplot2)
########################################################################
######### Figure for Network Stress vs Observed ASVs Fungi and Bacteria
########################################################################
Colors <- c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred")

stress <- aggregate(sqrt(Cytoscape_stats$Stress), by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_stress <- do.call(data.frame, stress)
myData_stress$se <- myData_stress$x.sd / sqrt(myData_stress$x.n)
colnames(myData_stress) <- c("Group", "median", "sd", "n", "se")

Total_diversity <- read.csv(file="../../Total_Diversity.txt", sep="\t",header = T)
ITS_obs <- aggregate(Total_diversity$ITS.Obs, by = list(Group =Total_diversity$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_ITS_obs <- do.call(data.frame, ITS_obs)
myData_ITS_obs$se <- myData_ITS_obs$x.sd / sqrt(myData_ITS_obs$x.n)
colnames(myData_ITS_obs) <- c("Group", "median", "sd", "n", "se")

pch_site<-c(1:10)[factor(myData_ITS_obs$Group)]
jpeg("FungalDiv_Stress_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_ITS_obs$median, myData_stress$median, main="Fungal diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Fungi))", ylab="Median (Network Stress)", col=Colors, pch=pch_site)
dev.off ()

BAC_obs <- aggregate(Total_diversity$BAC.Obs, by = list(Group =Total_diversity$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_BAC_obs <- do.call(data.frame, BAC_obs)
myData_BAC_obs$se <- myData_BAC_obs$x.sd / sqrt(myData_BAC_obs$x.n)
colnames(myData_BAC_obs) <- c("Group", "median", "sd", "n", "se")

pch_site<-c(1:10)[factor(myData_BAC_obs$Group)]
jpeg("BacDiv_Stress_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_BAC_obs$median, myData_stress$median, main="Bacterial diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Bacteria))", ylab="Median (Network Stress)", col=Colors, pch=pch_site)
dev.off ()


#######################################################
#########Let's make figure for Neiborhood connectivity 
#######################################################
connectivity <- aggregate(Cytoscape_stats$NeighborhoodConnectivity, by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_connectivity <- do.call(data.frame, connectivity)
myData_connectivity$se <- myData_connectivity$x.sd / sqrt(myData_connectivity$x.n)
colnames(myData_connectivity) <- c("Group", "median", "sd", "n", "se")

pch_site<-c(1:10)[factor(myData_ITS_obs$Group)]
jpeg("FungalDiv_Connectivity_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_ITS_obs$median, myData_connectivity$median, main="Fungal diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Fungi))", ylab="Median (Neighborhood Connectivity)", col=Colors, pch=pch_site)
dev.off ()

pch_site<-c(1:10)[factor(myData_BAC_obs$Group)]
jpeg("BacDiv_Connectivity_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_BAC_obs$median, myData_connectivity$median, main="Bacterial diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Bacteria))", ylab="Median (Neighborhood Connectivity)", col=Colors, pch=pch_site)
dev.off ()


#######################################################
#########Let's make figure for Number of direct Edges 
#######################################################
edges <- aggregate(Cytoscape_stats$NumberOfDirectedEdges, by = list(Group =Cytoscape_stats$Group), FUN = function(x) c(median = median(x), sd = sd(x), n = length(x)))
myData_edges <- do.call(data.frame, edges)
myData_edges$se <- myData_edges$x.sd / sqrt(myData_edges$x.n)
colnames(myData_edges) <- c("Group", "median", "sd", "n", "se")

pch_site<-c(1:10)[factor(myData_ITS_obs$Group)]
jpeg("FungalDiv_Edges_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_ITS_obs$median, myData_edges$median, main="Fungal diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Fungi))", ylab="Median (Directed Edges)", col=Colors, pch=pch_site)
dev.off ()

pch_site<-c(1:10)[factor(myData_BAC_obs$Group)]
jpeg("BacDiv_Edges_New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(myData_BAC_obs$median, myData_edges$median, main="Bacterial diversity vs Bacteria-Fungi connections", xlab="Median (Observed ASVs (Bacteria))", ylab="Median (Directed Edges)", col=Colors, pch=pch_site)
dev.off ()
