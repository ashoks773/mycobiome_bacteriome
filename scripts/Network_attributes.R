#-------- BaAka load both Negative and Positive correlations
my_adj_list <- read.csv(file="BaAka_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/BaAka_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/BaAka_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.9163022
l=layout_with_fr(net_pos)
jpeg("Network_Stats/BaAka_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#-- Cohesion
mwBlocks <- cohesive_blocks(net)
mwBlocks
blocks(mwBlocks)
cohesion(mwBlocks)

#----- Bantu
#-------- Bantu load both Negative and Positive correlations
my_adj_list <- read.csv(file="Bantu_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/Bantu_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Bantu_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()


#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.8479379
l=layout_with_fr(net_pos)
jpeg("Network_Stats/Bantu_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#-- Cohesion
mwBlocks <- cohesive_blocks(net)
mwBlocks
blocks(mwBlocks)
cohesion(mwBlocks)

g <- induced_subgraph(net, subcomponent(net, 1))
cohesion(g)

#---Chimps
#-------- Chimps load both Negative and Positive correlations
my_adj_list <- read.csv(file="Chimp_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/Chimps_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Chimps_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.6114407
l=layout_with_fr(net_pos)
jpeg("Network_Stats/Chimps_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()


#---- Mangabey
#-------- Mangabey load both Negative and Positive correlations
my_adj_list <- read.csv(file="Mangabey_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/Mangabey_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Mangabeys_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.6010617
l=layout_with_fr(net_pos)
jpeg("Network_Stats/Mangabey_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#----- For Mountain
#-------- Mountain load both Negative and Positive correlations
my_adj_list <- read.csv(file="Moutain_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/Mountain_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Mountain_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.3038948
l=layout_with_fr(net_pos)
jpeg("Network_Stats/Mountain_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#------- For USA
#-------- USA load both Negative and Positive correlations
my_adj_list <- read.csv(file="USA_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/USA_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/USA_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.0
l=layout_with_fr(net_pos)
jpeg("Network_Stats/USA_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()


#----- For Captive Western Lowland Gorilla
my_adj_list <- read.csv(file="Captive_WLC_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/Captive_WLG_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Captive_WLG_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.5403546
l=layout_with_fr(net_pos)
jpeg("Network_Stats/Captive_WLG_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#----- For Captive Chimps Ostrava
my_adj_list <- read.csv(file="CaptiveOstrava_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/CaptiveOstrava_hubScore_Keynote.txt", sep="\t")
as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/Captive_Ostrava_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.7624629
l=layout_with_fr(net_pos)
jpeg("Network_Stats/CaptiveOstrava_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()

#----- For Western Lowland Gorilla
my_adj_list <- read.csv(file="WLG_neg_pos_r06.txt", sep = "\t", header = T)
my_adj_list <- my_adj_list[,c(1,2,4)]
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
#write.table (hs_new, file="Network_Stats/WLG_hubScore_Keynote.txt", sep="\t")

as <- authority_score(net, weights=NA)$vector
jpeg("Network_Stats/WLG_Hubs_Authorities.jpg", height = 4, width = 4, units = 'in', res = 600)
par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs", vertex.label="")
plot(net, vertex.size=as*30, main="Authorities", vertex.label="")
dev.off ()

#--- Modularity
my_adj_list$weight <- abs(my_adj_list$weight) #-- Convert negative to Postive
my_adj_list_pos <- subset (my_adj_list, weight > 0.6)
net_pos <- graph.data.frame(my_adj_list_pos, directed = FALSE)
ceb <- cluster_edge_betweenness(net_pos) 
modularity(ceb)
# [1] 0.8240825
l=layout_with_fr(net_pos)
jpeg("Network_Stats/WLG_modularity.jpg", height = 4, width = 4, units = 'in', res = 600)
plot(ceb, net_pos, layout=l, vertex.label="")
dev.off ()


###########################################################################################
#----- To Calculate Cohesion We need to use the Relative abundance ASV table of 16S and ITS2
###########################################################################################
Table_ITS <- read.csv("~/Work/Project_16S_ITS/ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
Table_ITS <- data.frame(t(Table_ITS))
Table_16S <- read.csv("~/Work/Project_16S_ITS/16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1)
Table_16S <- data.frame(t(Table_16S))

Combined_ASVs <- data.frame(Table_16S, Table_ITS)

#-- Additional steps to filter OTUs
combinedASV.sums <- colSums(Combined_ASVs)
CombinedASV_filtered <- Combined_ASVs[ , which(combinedASV.sums > 10)] #---Check ASV sum should be more than 10
CombinedASV_filtered_filtered <- dropspc(CombinedASV_filtered, 5) #--- ASV should be present in more than 5 samples
Combined_ASV_table <- CombinedASV_filtered_filtered[order(rownames(CombinedASV_filtered_filtered)),]

#---- New metadata - elephants removed and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
ITSmeta <- read.csv("~/Work/Project_16S_ITS/Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]

Combined_ASV_table_meta <- merge(Combined_ASV_table, ITSmeta, by=0, all=F)
rownames(Combined_ASV_table_meta) <- Combined_ASV_table_meta$Row.names; Combined_ASV_table_meta$Row.names <- NULL

species <- Combined_ASV_table_meta[,1:3338]
meta <- Combined_ASV_table_meta[,3339:3344]

#species_relab_new <- cbind(BAC_species_relab, ITS_species_relab) #Alternate if want to use
species_relab <- decostand(species, method = "total")*100
species_relab_meta <- cbind (species_relab, meta)

############################
#--- For BaAka- 
############################
BaAka_species_relab <- subset(species_relab_meta, Group=="BaAka-Human")
BaAka_species_relab <- BaAka_species_relab[,1:3338]
#-- Additional steps to filter OTUs
baAkaASV.sums <- colSums(BaAka_species_relab)
BaAka_species_relab_filtered <- BaAka_species_relab[ , which(baAkaASV.sums > 0.2)] #Here We are talking about relative abudances
BaAka_species_relab_filtered <- dropspc(BaAka_species_relab_filtered, 5) #--- ASV should be present in more than 5 samples (20% of samples)
#write.table (BaAka_species_relab_filtered, file="Network_Stats/For_Cohesion/BaAka_species_relab_filtered.txt", sep = "\t")

Bantu_species_relab <- subset(species_relab_meta, Group=="Bantu-Human")
Bantu_species_relab <- Bantu_species_relab[,1:3338]
#-- Additional steps to filter OTUs
BantuASV.sums <- colSums(Bantu_species_relab)
Bantu_species_relab_filtered <- Bantu_species_relab[ , which(BantuASV.sums > 0.2)] #Here We are talking about relative abudances
Bantu_species_relab_filtered <- dropspc(Bantu_species_relab_filtered, 3) #--- ASV should be present in more than 3 samples (20% of samples)
write.table (Bantu_species_relab_filtered, file="Network_Stats/For_Cohesion/Bantu_species_relab_filtered.txt", sep = "\t")

CaptiveOstrava_species_relab <- subset(species_relab_meta, Group=="Captive Chimps-Ostrava")
CaptiveOstrava_species_relab <- CaptiveOstrava_species_relab[,1:3338]
#-- Additional steps to filter OTUs
CaptiveOstravaASV.sums <- colSums(CaptiveOstrava_species_relab)
CaptiveOstrava_species_relab_filtered <- CaptiveOstrava_species_relab[ , which(CaptiveOstravaASV.sums > 0.2)] #Here We are talking about relative abudances
CaptiveOstrava_species_relab_filtered <- dropspc(CaptiveOstrava_species_relab_filtered, 3) #--- ASV should be present in more than 3 samples (20% of samples)
write.table (CaptiveOstrava_species_relab_filtered, file="Network_Stats/For_Cohesion/CaptiveOstrava_species_relab_filtered.txt", sep = "\t")

MountainG_species_relab <- subset(species_relab_meta, Group=="Mountain Gorilla")
MountainG_species_relab <- MountainG_species_relab[,1:3338]
#-- Additional steps to filter OTUs
MountainGASV.sums <- colSums(MountainG_species_relab)
MountainG_species_relab_filtered <- MountainG_species_relab[ , which(MountainGASV.sums > 0.2)] #Here We are talking about relative abudances
MountainG_species_relab_filtered <- dropspc(MountainG_species_relab_filtered, 5) #--- ASV should be present in more than 6 samples (20% of samples)
write.table (MountainG_species_relab_filtered, file="Network_Stats/For_Cohesion/MountainG_species_relab_filtered.txt", sep = "\t")
