#---------- Correlations between selected ITS + selected 16S
bac_sel_asv <- read.csv(file="../Bac_Taxa_Analysis/Sel_16S_ASV_filtered_Rel_Abun.txt", sep = "\t", row.names = 1, header = T)
bac_sel_asv_t <- data.frame(t(bac_sel_asv))

its_sel_asv <- read.csv(file="../ITS_Taxa_Analysis/Sel_ITS_ASV_filtered_Rel_Abun.txt", sep = "\t", row.names = 1, header = T)
its_sel_asv_t <- data.frame(t(its_sel_asv))

library (ccrepe)
ccrepe_bac_its <- ccrepe(x = bac_sel_asv_t, y = its_sel_asv_t, sim.score = nc.score, iterations = 50, min.subj = 5)
r_ccrepe_bac_its <- ccrepe_bac_its$sim.score
q_ccrepe_bac_its <- ccrepe_bac_its$q.values

r_ccrepe_bac_its <- melt(r_ccrepe_bac_its)
q_ccrepe_bac_its <- melt (q_ccrepe_bac_its)

q_r_ccrepe_bac_its <- cbind(r_ccrepe_bac_its,q_ccrepe_bac_its)
q_r_ccrepe_bac_its_Final <- q_r_ccrepe_bac_its[,c(1:3,6)]


colnames(q_r_ccrepe_bac_its_Final) <- c("Bacteria", "ITS", "r_value", "q_value")
#q_r_ccrepe_bac_its_Final_formatted <- q_r_ccrepe_bac_its_Final %>% filter(!is.na(q_value))
write.table (q_r_ccrepe_bac_its_Final, file="Bac_ITS_corr.txt", sep = "\t")

# Formt Bac_ITS_corr.txt in vi editor using :%s/^X//g , and :%s/\tX/\t/g commands
#perl pick_details.pl ../Bac_Taxa_Analysis/indvalsummary.txt.taxonomy.txt ../ITS_Taxa_Analysis/indvalsummary.txt.taxonomy.txt Bac_ITS_corr.txt


###############################################################################
#--------- CytoScape Network plots for Figure 3
##############################################################################

##############################################################################
#------------- Top correlated Bacteria and Fungi in Every Group - For Figure5
##############################################################################
setwd("~/Work/Project_16S_ITS/Combined/BAC_ITS_Taxa_Correlations")

Bact_Fungi <- cbind(bac_sel_asv_t, its_sel_asv_t)

#---- New metadata - elephants removed (total 18), less than 1k depht removed (total 18) and 1220, 11A2, and 121B (Unkonwns used as Captive Apes)
meta <- read.csv("../../Metadata_Filtered.txt", sep = "\t", row.names = 1)
meta_t <- data.frame(t(meta))
meta <- data.frame(t(meta_t))
meta <- meta[order(rownames(meta)),]

Bact_Fungi_meta <- merge(Bact_Fungi, meta, by=0, all=F)
rownames(Bact_Fungi_meta) <- Bact_Fungi_meta$Row.names; Bact_Fungi_meta$Row.names <- NULL

Bact_Fungi_meta <- cbind(Bact_Fungi, meta)

jpeg("BaAka-Fecali_Saccharomycetales.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X5edd946e0a26ae1c85ed8c498c328387", y = "X87f1f89372c1542c10504ef404c7ba31",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Faecalibacterium prausnitzii", ylab = "Saccharomycetales")
dev.off ()

jpeg("Bantu-Clostridium_Saccharomycetales.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "c5f7927fc747803122ff7f5a60e8a521", y = "X6aed27fe0a1493fe286a0d9dbd827347",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Clostridium butyricum", ylab = "Saccharomycetales")
dev.off ()

jpeg("Liberac-Paraprevotellaceae_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X4db6fcc5aaf1645e73249672c5a4956e", y = "X125cc88523ee662f0b44d8882a02fefb",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Paraprevotellaceae", ylab = "Unidentified")
dev.off ()
# 100% qcov, 100% seq identity, Sequence ID: KY932438.1 

jpeg("Ostrava-Lachnospiraceae_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X1f2c48155561e7d6b74b67ca7f068922", y = "X799cc660a43d6ee89a7ad4907111c073",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Lachnospiraceae", ylab = "Unidentified")
dev.off ()

jpeg("CaptiveWLG-Coriobacteriaceae_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X10b0ccda77f21a69527f7596b15c1c23", y = "eae7022c2729b4b1ea7a62abe468e799",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Coriobacteriaceae", ylab = "Unidentified")
dev.off ()

jpeg("CaptiveWLG-Coriobacteriaceae_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X10b0ccda77f21a69527f7596b15c1c23", y = "eae7022c2729b4b1ea7a62abe468e799",
          shape = "Group", size = 1,
          cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Coriobacteriaceae", ylab = "Unidentified Nohit")
dev.off ()

jpeg("Chimps-Prevotella_Xylariales.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "ba00d2e0c74c53a4f32a5b6c48859a9c", y = "X3a439a0e92578bb8f02ffdeeafdef664",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Prevotella", ylab = "Xylariales")
dev.off ()

jpeg("Mangabey-Clostridiales_Phoma.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "c11b67d8c052cb2b8eaae7e9bf30d1f9", y = "f43b8a3da88ed03175b1e9ea64eb9bc0",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Clostridiales", ylab = "Phoma")
dev.off ()

jpeg("MountainG-Clostridiales_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X3c9cfbf9c8e038663ccfa9ef241898ca", y = "X473285b616b9b3435278f9a53a1fbf8a",
          shape = "Group", size = 1,
         # cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Clostridiales", ylab = "Unidentified")
dev.off ()
#Blast hit
#Murshidia linstowi isolate DKAT-18 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and large subunit ribosomal RNA gene, partial sequence	
#qcov-99%	evalue-0.0	%iden-97.18%,	Sequence ID: MK968095.1

jpeg("USA-Coprococcus_Saccharomyces.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "X973ca62579b1246f54829278d51475ea", y = "X1c9ff9e2aeeaebc0df06f1e182022bb7",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Coprococcus", ylab = "Saccharomyces")
dev.off ()

jpeg("WLG-Treponema_Unidentified.jpg", height = 3, width = 3, units = 'in', res = 600)
ggscatter(Bact_Fungi_meta, x = "fb620f679ca2c3918529dc0cd826a317", y = "X54b32922901ebad101bbaec01bbccf8f",
          shape = "Group", size = 1,
          #cor.coef = TRUE, cor.method = "spearman", color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          color = "Group",    palette = c("forestgreen", "darkgray", "darkgoldenrod3", "chartreuse1", "purple3", "darkmagenta", "blue", "darkcyan", "darkorchid1", "darkred"),
          ellipse = TRUE, mean.point = TRUE, star.plot = TRUE,
          xlab = "Treponema", ylab = "Unidentified")
dev.off ()
#Blasthit
#Whitfieldia elongata 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence
#qcov-97%	evalue-0.0	%iden-99.18%,	Sequence ID: EU528910.1