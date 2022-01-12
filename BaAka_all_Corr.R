#############################################################
#----- For BaAka - Correlations
#############################################################
#---------- First load the complete ASV tables FOR ITS
ITS_asv_table_in <- read.csv("../../ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
ITS_asv_table_in_t <- data.frame(t(ITS_asv_table_in))

meta <- read.csv("BaAka.txt", sep = "\t", row.names = 1)
meta_t <- data.frame(t(meta))
meta <- data.frame(t(meta_t))

ITS_asv_meta <- merge(ITS_asv_table_in_t, meta, by=0, all=F)
rownames(ITS_asv_meta) <- ITS_asv_meta$Row.names; ITS_asv_meta$Row.names <- NULL

ITS_asv_BaAka <- ITS_asv_meta[,1:13249]

ITSasv.sums <- colSums(ITS_asv_BaAka)
ITSasv_filtered <- ITS_asv_BaAka[ , which(ITSasv.sums > 10)] #---Check ASV sum should be more than 10
ITSasv_filtered_filtered <- dropspc(ITSasv_filtered, 3) #--- ASV should be present in atleast 10% of samples
ITSasv_table <- ITSasv_filtered_filtered[order(rownames(ITSasv_filtered_filtered)),]

ITSasv_table_relab <- decostand(ITSasv_table, method = "total")*100


#############################################################
#---------- First load the complete ASV tables - FOR BACTERIA
BAC_asv_table_in <- read.csv(file="../../16S_analysis/16S_feature-table.txt", sep = "\t", row.names = 1, header = T)
BAC_asv_table_in_t <- data.frame(t(BAC_asv_table_in))

BAC_asv_meta <- merge(BAC_asv_table_in_t, meta, by=0, all=F)
rownames(BAC_asv_meta) <- BAC_asv_meta$Row.names; BAC_asv_meta$Row.names <- NULL

BAC_asv_BaAka <- BAC_asv_meta[,1:9272]

BACasv.sums <- colSums(BAC_asv_BaAka)
BACasv_filtered <- BAC_asv_BaAka[ , which(BACasv.sums > 10)] #---Check ASV sum should be more than 10
BACasv_filtered_filtered <- dropspc(BACasv_filtered, 3) #--- ASV should be present in atleast 10% of samples
BACasv_table <- BACasv_filtered_filtered[order(rownames(BACasv_filtered_filtered)),]

BACasv_table_relab <- decostand(BACasv_table, method = "total")*100
write.table (BACasv_table_relab, file="BACasv_table_relab.txt", sep ="\t")
#------ CCREPE Correlations on Filtered Bacterial and ITS relative abundances

library (ccrepe) #minsubjects should be atleast 50%
ccrepe_bac_its <- ccrepe(x = BACasv_table_relab, y = ITSasv_table_relab, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"), iterations = 1000, min.subj = 14)

r_ccrepe_bac_its <- ccrepe_bac_its$sim.score
q_ccrepe_bac_its <- ccrepe_bac_its$p.values

r_ccrepe_bac_its <- melt(r_ccrepe_bac_its)
q_ccrepe_bac_its <- melt (q_ccrepe_bac_its)

q_r_ccrepe_bac_its <- cbind(r_ccrepe_bac_its,q_ccrepe_bac_its)
q_r_ccrepe_bac_its_Final <- q_r_ccrepe_bac_its[,c(1:3,6)]


colnames(q_r_ccrepe_bac_its_Final) <- c("Bacteria", "ITS", "r_value", "p_value")
#q_r_ccrepe_bac_its_Final_formatted <- q_r_ccrepe_bac_its_Final %>% filter(!is.na(q_value))
write.table (q_r_ccrepe_bac_its_Final, file="Bac_ITS_BaAka_corr_Check_Again.txt", sep = "\t")


#------ Igraph to check the KeyNote taxa and other Parameters
# https://kateto.net/netscix2016.html

#------ Load BaAka data
#load ("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/BaAka_check.RData")


#-------- BaAka load both Negative and Positive correlations

my_cor_df <- melt(ccrepe_bac_its$sim.score)
#library("tidyr")
#data4 <- my_cor_df %>% drop_na() 
my_adj_list <- subset(my_cor_df, value > 0.6)
#my_adj_list <- my_cor_df %>% filter(abs(value) > 0.6)
#my_adj_list <- my_cor_df %>% filter(value > 0.6)
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

#write.table (my_adj_list, file = "BaAka_my_adj_list.txt", sep ="\t")
#-- Change IDs F1, F2 for fungi and B1, B2 for bacteria
my_adj_list <- read.csv(file="BaAka_my_adj_list.txt", sep="\t", header = T, row.names = 1)
my_adj_list <- my_adj_list [,c(2,4,5)]

net <- graph.data.frame(my_adj_list, directed = FALSE)
net

#To plot
plot(net, edge.arrow.size=.4,vertex.label=NA)
plot(net, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$media, vertex.label.color="black",
     vertex.label.cex=.7) 

#--- Parameters
edge_density(net, loops=F)

reciprocity(net)
dyad_census(net)
2*dyad_census(net)$mut/ecount(net)

#--- Hub Score
hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
write.table (hs_new, file="BaAka_hubScore_Keynote.txt", sep="\t")

#--- cohesion
##https://rdrr.io/cran/igraph/man/cohesive_blocks.html
mwBlocks <- cohesive_blocks(net)
mwBlocks
blocks(mwBlocks)
cohesion(mwBlocks)
export_pajek(mwBlocks, net, file="BaAka_mwBlocks.paj")
plot(mwBlocks, net)
plot(mwBlocks, net, vertex.label=V(net)$name, margin=-0.2,
     vertex.shape="rectangle", vertex.size=24, vertex.size2=8,
     mark.border=1, colbar=c(NA, NA,"cyan","orange") )

#-- Modularity (This will work only on strictly positives)
#my_adj_list <- my_cor_df %>% filter(value > 0.6)
#net <- graph.data.frame(my_adj_list, directed = FALSE)
ceb <- cluster_edge_betweenness(net) 
modularity(ceb)

###############################
#--- Trying for cohesion Again
install.packages("statnet", dependencies=TRUE)
library(statnet)
g <- graph.data.frame(my_adj_list, directed = TRUE)              # Load with igraph
net1<-network(my_adj_list, directed=TRUE, matrix.type="edgelist")

coreness(g)
kcore <- coreness(g)    # Extract k-cores as a data object.
V(g)$core <- kcore      # Add the cores as a vertex attribute
plot.igraph(g, vertex.color=V(g)$core) # plot

set.seed(18675309)
kco <- kcores(net1) 

gplot(net1, 
      vertex.col=kco, 
      usearrows=FALSE,
      displaylabels=TRUE)

#Girvan-Newman (Edge Betweenness)
g <- as.undirected(g, mode="collapse")
eb <- cluster_edge_betweenness(g)
membership(eb)
V(g)$group <- membership(eb)
sizes(eb)

plot(eb, g,
     vertex.label = V(g)$id,  # Uses index number instead of name
     layout = layout_with_kk, 
     main="Max Modularity Solution")

##-----
V(g)$group <- membership(eb)
# Change the number of communities in the solution.
eb$membership <- cut_at(eb, 2) # 2 communities

plot(eb, g,
     vertex.label = V(g)$id,  # Uses index number instead of name
     layout = layout_with_kk, 
     main="Two-Community Solution")
###-----
gBlocks <- cohesive_blocks(g)
blocks(gBlocks)
plot(gBlocks, g,
     layout = layout_with_kk)



#--- Other parameters
#--- Transivity
transitivity(net, type="global")
transitivity(as.undirected(net, mode="collapse"))
transitivity(net, type="local")

cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))

#V(net)$community <- cfg$membership
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
#plot(net, vertex.color=colrs[V(net)$community])


# Assortability and Homophily
assortativity_nominal(net, V(net), directed=F)
assortativity_degree(net, directed=F)


#=-------
BaAka_corr <- read.csv("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Corr_New_Check/Selected_corrs/BaAka_corr_r06_p05.txt", sep = "\t", header = T)
names(BaAka_corr)[names(BaAka_corr)=="r_value"] <- "r"
names(BaAka_corr)[names(BaAka_corr)=="p_value"] <- "p"

graph_BaAka_corr <- BaAka_corr %>%
  filter(abs(r) > .6) %>%
  filter(abs(p) < 0.05) %>%
  graph_from_data_frame(directed = FALSE)



#------
