library("reshape2")
library("stringr")
library("ggplot2")
library("RColorBrewer")

# Path to input files
otu_tab_file <- "OTU_table"
tax_file <- "taxonomy"

#replicates_list <- c("X33_Didon", "X34_Boto", "X35_Ngombi", "X36_Molongo", "X37_Ngonduma", "X38_Ito", "X39_Iyiki", "X40_Epessekele", "X41_Francois", "X42_Ndidje", "X43_Wala", "X44_Zulia", "X45_Maki", "X46_Singa", "X47_Mbusa", "X49_Bonango", "X50_Bajungo", "X51_Malambo", "X52_Lindaki", "X53_Yambi", "X54_Manoi", "X55_Luma", "X56_Adamou", "X57_Yama", "X58_Abuya", "X59_Abeni", "X60_Omer", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25", "X29", "X30", "X31", "X32", "X35", "D-H-S", "Day", "HA", "MK", "Nne", "PA", "T1", "X1FECOB", "X2FecJUDY", "X4FecTEA", "CH1Hope", "CH2Bambari", "CH3Maja", "CH4Zira", "CH10", "CH11", "CH12", "CH13", "CH14", "CH15", "CH16", "CH17", "CH6", "CH7", "CH8", "CH9", "X1", "X2", "X3", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "C6", "X102M", "X249M", "X251M", "X253M", "X87M", "M157", "M158", "M159", "M255", "M45", "M49", "B1", "B2", "B3", "B4", "H1", "K1", "K2", "K3", "K4", "K5", "K6", "M1", "M2", "M3", "M4", "M5", "M7", "MH10", "MH2", "MH3", "MH4", "MH5", "MH6", "MH7", "MH8", "MH9", "X1220", "X11A2", "X121B", "H10", "H100", "H12", "H2", "H3", "H4", "H5", "H6", "H7", "H9", "HB", "X157", "X1_Mata", "X10_Sopo", "X12_Makumba", "X16Susa", "X17_Mayele", "X18Lungu", "X21Mapoki", "X22_Susa", "X23Mongali7516", "X24_Mongali", "X25_Mayele", "X28Wiya", "X3_Mama", "X30_Duma", "X4_Mio", "X5Malui", "X6_samba", "X7_Tembo", "G161", "G6", "Jab", "Jab2", "Sam2", "Samson", "SCH2", "Scho2", "Shorder")
#replicates_groups <- c("BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "USA-Human", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Liberec", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla")

#--- Remove Captive Chimps Hodonin
replicates_list <- c("X33_Didon", "X34_Boto", "X35_Ngombi", "X36_Molongo", "X37_Ngonduma", "X38_Ito", "X39_Iyiki", "X40_Epessekele", "X41_Francois", "X42_Ndidje", "X43_Wala", "X44_Zulia", "X45_Maki", "X46_Singa", "X47_Mbusa", "X49_Bonango", "X50_Bajungo", "X51_Malambo", "X52_Lindaki", "X53_Yambi", "X54_Manoi", "X55_Luma", "X56_Adamou", "X57_Yama", "X58_Abuya", "X59_Abeni", "X60_Omer", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25", "X29", "X30", "X31", "X32", "X35", "D-H-S", "Day", "HA", "MK", "Nne", "PA", "T1", "CH10", "CH11", "CH12", "CH13", "CH14", "CH15", "CH16", "CH17", "CH6", "CH7", "CH8", "CH9", "X1", "X2", "X3", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "C6", "X102M", "X249M", "X251M", "X253M", "X87M", "M157", "M158", "M159", "M255", "M45", "M49", "B1", "B2", "B3", "B4", "H1", "K1", "K2", "K3", "K4", "K5", "K6", "M1", "M2", "M3", "M4", "M5", "M7", "MH10", "MH2", "MH3", "MH4", "MH5", "MH6", "MH7", "MH8", "MH9", "X1220", "X11A2", "X121B", "H10", "H100", "H12", "H2", "H3", "H4", "H5", "H6", "H7", "H9", "HB", "X157", "X1_Mata", "X10_Sopo", "X12_Makumba", "X16Susa", "X17_Mayele", "X18Lungu", "X21Mapoki", "X22_Susa", "X23Mongali7516", "X24_Mongali", "X25_Mayele", "X28Wiya", "X3_Mama", "X30_Duma", "X4_Mio", "X5Malui", "X6_samba", "X7_Tembo", "G161", "G6", "Jab", "Jab2", "Sam2", "Samson", "SCH2", "Scho2", "Shorder")
replicates_groups <- c("BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "BaAka-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Bantu-Human", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "USA-Human", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Captive Chimps-Ostrava", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Chimps", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mangabey", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Mountain Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "USA-Human", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla", "Captive Western Lowland Gorilla")

length(replicates_list) == length(replicates_groups) 

tax_aggr <- "Family"
tax_number <- 20

tax_col <- "Phylum"

file_name <- "Figure1C.svg"
plot_dim <- c(4,8)

otu_tab <- read.table("~/Box/Gomez_Lab/Project_16S_ITS/ITS_analysis/ITS_feature-table.txt", header = TRUE, comment.char = "", sep = "\t")
names(otu_tab)[1] <- "OTU"

tax_tab <- read.table("~/Box/Gomez_Lab/Project_16S_ITS/ITS_analysis/ITS_taxonomy_formatted.tsv", header = TRUE, comment.char = "", sep = "\t", fill = TRUE)
names(tax_tab)[1] <- "OTU"


#-------------- This code removes useless annotations from the taxonomic annotation file, such as 'uncultured bacteria', or 'unkwown family', and fills the blanks by copying the previous taxonomic levels (while adding an 'unassigned' mention).
for (col in 2:ncol(tax_tab)) {
  for (row in 1:nrow(tax_tab)) {
    if (grepl("uncultured",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
    if (grepl("unknown",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
  }
}

# Replace empty cells by 'NA'
tax_tab2 <- as.data.frame(apply(tax_tab, 2, function(x) gsub("^$|^ $", NA, x)))

# Remove columns containing only 'NA'
col_to_remove <- c()

for (col in 2:ncol(tax_tab2)) {
  x <- sum(is.na(tax_tab2[,col]))/nrow(tax_tab2)
  if (x == 1) {
    col_to_remove <- c(col_to_remove, col)
  }
}

if (length(col_to_remove) > 0) {
  tax_tab3 <- tax_tab2[,-col_to_remove]
} else {
  tax_tab3 <- tax_tab2
}

for (col in 2:ncol(tax_tab3)) {
  tax_tab3[,col] <- as.character(tax_tab3[,col])
}

# Fill all NAs

for (col in 2:ncol(tax_tab3)) {
  for (row in 1:nrow(tax_tab3)) {
    if (is.na(tax_tab3[row,col])) {
      if (!grepl("OTU", tax_tab3[row,col-1]) & !grepl("unassigned", tax_tab3[row,col-1])) {
        tax_tab3[row,col] <- paste0("unassigned ", tax_tab3[row,col-1])
      } else {
        tax_tab3[row,col] <- tax_tab3[row,col-1]
      }
    }
  }
}

#---------Compute the relative abundance of OTUs for each sample
otu_counts <- colSums(otu_tab[,-1])
otu_tab2 <- otu_tab
otu_tab2[,-1] <- sweep(otu_tab[,-1], 2, otu_counts, `/`)
otu_tab2[is.na(otu_tab2)] <- 0

# Check that the sum of relative abundances for each sample is 1.
colSums(otu_tab2[,-1])

#----------
m <- merge(otu_tab2, tax_tab3)
dim(m)

taxonomy <- c()
for (row in 1:nrow(m)) {
  taxonomy <- c(taxonomy, paste0(m[row,names(m)==tax_col], ";", m[row,names(m)==tax_aggr]))
}

# Subset from the merged table the selected samples only
m2 <- m[,names(m) %in% replicates_list]

# Aggregate 'm2' based on the selected taxonomic level
m3 <- aggregate(m2, by=list(taxonomy), FUN=sum)

dim(m3)

#-- Export family table to do the statistical test
#write.table (m3, file ="ITS_family_proportions.txt", sep="\t")

m3[1:5,1:4]


#---------
if (tax_number > nrow(m3)) {
  tax_number <- nrow(m3)
}

m3$average <- rowMeans(m3[,-1])
m3.sorted <- m3[order(-m3$average),]

# Aggregate the smaller taxonomic bins together
m3.sorted$selection <- rep("discarded", nrow(m3.sorted))
m3.sorted$selection[1:tax_number] <- "retained"
m3.sorted$Group.1[m3.sorted$selection == "discarded"] <- "Other;Other"
m3.sorted$average <- NULL
m3.sorted$selection <- NULL
m4 <- aggregate(m3.sorted[,-1], by=list(taxonomy=m3.sorted$Group.1), FUN=sum)

m4[m4$taxonomy == "Other;Other", -1]
mean(as.numeric(m4[m4$taxonomy == "Other;Other", -1]))

#------ Transpose m4
n <- m4$taxonomy
m4.t <- as.data.frame(t(m4[,-1]))
colnames(m4.t) <- n
m4.t$sample <- rownames(m4.t)
rownames(m4.t) <- NULL

#---- Calculate mean and standard deviation for each samples
m4.t$replicate <- rep(NA, nrow(m4.t))
for (line in 1:(nrow(m4.t))){
  m4.t$replicate[line] <- replicates_groups[m4.t$sample[line] == replicates_list]
}

# Compute the mean
m4.t.mean <- aggregate(m4.t[,1:(ncol(m4.t)-2)],
                       by = list(m4.t$replicate),
                       FUN = "mean")
names(m4.t.mean)[1] <- "sample"

dim(m4.t.mean)                              
## [1]  6 42

# Compute the standard deviation                             
m4.t.sd <- aggregate(m4.t[,1:(ncol(m4.t)-2)],
                     by = list(m4.t$replicate),
                     FUN = "sd")
names(m4.t.sd)[1] <- "sample"

dim(m4.t.sd) 

#----------
molten.mean <- melt(m4.t.mean, id.vars = "sample")
molten.mean$id <- paste0(molten.mean$sample, "-", molten.mean$variable)

molten.sd <- melt(m4.t.sd, id.vars = "sample")
molten.sd$id <- paste0(molten.sd$sample, "-", molten.sd$variable)

# Merge the dataframes
molten <- merge(molten.mean, molten.sd, by.x = "id", by.y = "id")


#--------- Final arrangements of the dataframe

molten$id <- NULL
molten$sample.y <- NULL
molten$variable.y <- NULL
names(molten) <- c("sample", "taxonomy", "mean", "sd")

molten$tax_col <- str_split_fixed(molten$taxonomy, ";", 2)[,1]
molten$tax_bin <- str_split_fixed(molten$taxonomy, ";", 2)[,2]

# Reorder the taxonomic annotation for the plot
molten <- molten[order(molten$tax_col),]
tax_levels <- as.character(molten$tax_bin[!duplicated(molten$tax_bin)])
tax_levels <- tax_levels[tax_levels != "Other"]
tax_levels <- c(tax_levels, "Other")
molten$tax_bin <- factor(molten$tax_bin, levels = rev(tax_levels))

# Reorder the samples for the plot
#molten$sample <- factor(molten$sample, levels = replicates_groups)

# Remove null values
molten2 <- molten[molten$mean > 0,]

MyPalette <- c("forestgreen", "#5DD0B9", "blue", "darkgoldenrod3", "chartreuse1", 
               "purple3", "darkmagenta", "darkorchid1", "darkred", "#ffff99",
               "#de77ae", "#6a3d9a", "darkgray", "lightgreen", "brown", "magenta", "pink")

#------- Bubble plot
bubble_plot <- ggplot(molten2,aes(sample,tax_bin)) +
  #geom_point(aes(size=mean+sd), shape=16, color = "red") + 
  geom_point(aes(size=mean, fill=tax_col),shape=21,color="black") +  
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_blank()) +
  ylab("Taxonomy (Family)") +
  xlab("Samples") +
  scale_fill_manual(values=MyPalette, name="Taxonomic\nclade") +
 # scale_fill_discrete(name="Taxonomic\nclade") +
  #scale_fill_manual(values= c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3")), name="Taxonomic\nclade") +
  scale_size(name = "Relative\nabundance")

#svg(file="bubble_plot.svg", width = plot_dim[1], height = plot_dim[2])
png(file="Figure1C.png", height = 6, width = 6, units = 'in', res = 600)
bubble_plot
dev.off()
