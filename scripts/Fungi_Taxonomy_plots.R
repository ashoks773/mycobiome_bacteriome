#Step 1: ---------- Phyla plot for Each sample in the Group

asv_table_in <- read.csv("../ITS_analysis/ITS_feature-table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))

#-- Additional steps to filter OTUs
asv.sums <- colSums(asv_table_in_t)
asv_table_filtered <- asv_table_in_t[ , which(asv.sums > 10)] #---Check ASV sum should be more than 50
#asv_table_filtered_filtered <- dropspc(asv_table_filtered, 5) #--- ASV should be present in more than 5 samples
asv_table_filtered_t <- data.frame(t(asv_table_filtered))

asv_table_in <- as.matrix(asv_table_filtered_t)

# Read in taxonomy
# Fomatted Separated by kingdom, phylum, class, order, family, genus, species, 
# and No Phyla was aslo renamed as unIdentified
#-- kingdom Chromista, Metazoa, Plantae, Rhizaria, unassigned were removed
taxonomy <- read.csv("../ITS_analysis/ITS_taxonomy_formatted.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("../Metadata_Filtered.txt", sep="\t", row.names = 1, header=T)
metadata_t <- data.frame(t(metadata))
metadata <- data.frame(t(metadata_t))

# Read in tree
phy_tree <- read_tree("../ITS_analysis/ITS_tree.nwk")
# Import all as phyloseq objects
ASV <- otu_table(asv_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(ASV)
taxa_names(phy_tree)
# Same sample names
sample_names(ASV)
sample_names(META)

# Finally merge to create Phyloseq object!
ps <- phyloseq(ASV, TAX, META, phy_tree)
ps_normalized  = transform_sample_counts(ps, function(x) x / sum(x) )

MyPalette <- c("forestgreen", "#5DD0B9", "darkgoldenrod3", "chartreuse1", 
               "purple3", "darkmagenta", "blue", "#1f78b4", "darkcyan", "darkorchid1", "darkred", "#ffff99",
               "#de77ae", "#6a3d9a", "darkgray", "gray32")

jpeg("ITS_Phylum.jpg", height = 10, width = 12, units = 'in', res = 600)
p <- plot_bar(ps_normalized, fill = "Phylum", facet_grid=~Family) + scale_fill_manual(values=MyPalette)
## add facets
p <- p + facet_wrap(~Group, scales="free_x", nrow = 2)
plot(p)
dev.off ()


#Step 2: ---------- Phyla plot for comined Groups (This file has been manually created)
#------- Combined Plot ----- phyla plots only for Groups
asv_table_in <- read.csv("../ITS_analysis/featureTable_forCumulativePhylaPlot.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))
#ab <- table(unlist(asv_table_in_t))

#-- Additional steps to filter OTUs
asv.sums <- colSums(asv_table_in_t)
asv_table_filtered <- asv_table_in_t[ , which(asv.sums > 10)] #---Check ASV sum should be more than 50
#asv_table_filtered_filtered <- dropspc(asv_table_filtered, 5) #--- ASV should be present in more than 5 samples
asv_table_filtered_t <- data.frame(t(asv_table_filtered))

asv_table_in <- as.matrix(asv_table_filtered_t)

# Read in taxonomy
# Fomatted Separated by kingdom, phylum, class, order, family, genus, species, 
# and No Phyla was aslo renamed as unIdentified
#-- kingdom Chromista, Metazoa, Plantae, Rhizaria, unassigned were removed
taxonomy <- read.csv("../ITS_analysis/ITS_taxonomy_formatted.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("../meta_CumulativePhylaPlot.txt", sep="\t", row.names = 1, header=T)
metadata_t <- data.frame(t(metadata))
metadata <- data.frame(t(metadata_t))

# Read in tree
phy_tree <- read_tree("../ITS_analysis//ITS_tree.nwk")
# Import all as phyloseq objects
ASV <- otu_table(asv_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(ASV)
taxa_names(phy_tree)
# Same sample names
sample_names(ASV)
sample_names(META)

# Finally merge to create Phyloseq object!
ps <- phyloseq(ASV, TAX, META, phy_tree)
ps_normalized  = transform_sample_counts(ps, function(x) x / sum(x) )

MyPalette <- c("forestgreen", "#5DD0B9", "darkgoldenrod3", "chartreuse1", 
               "purple3", "darkmagenta", "blue", "#1f78b4", "darkcyan", "darkorchid1", "darkred", "#ffff99",
               "#de77ae", "#6a3d9a", "darkgray", "gray32")

jpeg("ITS_Phylum_Cumulative.jpg", height = 8, width = 5, units = 'in', res = 600)
p <- plot_bar(ps_normalized, fill = "Phylum", facet_grid=~Family) + scale_fill_manual(values=MyPalette)
## add facets
p <- p + facet_wrap(~Group, scales="free_x", nrow = 1)
plot(p)
dev.off ()

#----- Statistical significance of Each Phyla plotted in Step 1 and Step 2


#-- Step 3: Phyla and Family Bubble Plot
source('~/Work/Project_16S_ITS/Combined/Bubble_ITS.R')


#----- Statistical significance of Each Family plotted in Step 3