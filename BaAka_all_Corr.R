#############################################################
#----- For Mountain Gorilla - Correlations
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

#------ CCREPE Correlations on Filtered Bacterial and ITS relative abundances

library (ccrepe) #minsubjects should be atleast 50%
#ccrepe_bac_its <- ccrepe(x = BACasv_table_relab, y = ITSasv_table_relab, sim.score = nc.score, iterations = 1000, min.subj = 14)
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

my_cor_df <- melt(ccrepe_bac_its$sim.score)
my_adj_list <- my_cor_df %>% filter(abs(value) > 0.6)
#my_adj_list <- my_cor_df %>% filter(value > 0.6)
names(my_adj_list) <- c('from', 'to', 'weight')
class(my_adj_list)
dim(my_adj_list)

net <- graph.data.frame(my_adj_list, directed = FALSE)

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

hs <- hub_score(net, weights=NA)$vector
hs_new <- data.frame (hs)
write.table (hs_new, file="BaAka_hubScore_Keynote.txt", sep="\t")

#--- Transivity
transitivity(net, type="global")
transitivity(as.undirected(net, mode="collapse"))
transitivity(net, type="local")

#-- Modularity (This will work only on strictly positives)
my_adj_list <- my_cor_df %>% filter(value > 0.6)
net <- graph.data.frame(my_adj_list, directed = FALSE)
ceb <- cluster_edge_betweenness(net) 
modularity(ceb)

cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))

#V(net)$community <- cfg$membership
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
#plot(net, vertex.color=colrs[V(net)$community])


# Assortability and Homophily
assortativity_nominal(net, V(net), directed=F)
assortativity_degree(net, directed=F)


#=-------
BaAka_corr <- read.csv("BaAka_corr_r06_p05.txt", sep = "\t", header = T)
names(BaAka_corr)[names(BaAka_corr)=="r_value"] <- "r"
names(BaAka_corr)[names(BaAka_corr)=="p_value"] <- "p"

graph_BaAka_corr <- BaAka_corr %>%
  filter(abs(r) > .6) %>%
  filter(abs(p) < 0.05) %>%
  graph_from_data_frame(directed = FALSE)
