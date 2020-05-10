#-- Bacteria relative abundance already calculated
BAC_rel_abun <- read.csv(file = "../16S_ASV_filtered_Rel_Abun.txt", sep = "\t", row.names = 1, header=T)
 
BACmeta <- read.csv("../../Metadata_Filtered.txt", sep = "\t", row.names = 1)
BACmeta_t <- data.frame(t(BACmeta))
BACmeta <- data.frame(t(BACmeta_t))
BACmeta <- BACmeta[order(rownames(BACmeta)),]

BAC_rel_abun_meta <- merge(BAC_rel_abun, BACmeta, by=0, all=F)
rownames(BAC_rel_abun_meta) <- BAC_rel_abun_meta$Row.names; BAC_rel_abun_meta$Row.names <- NULL

library (labdsv)
BAC_rel_abun_meta_filtered <- BAC_rel_abun_meta[ , which(!apply(BAC_rel_abun_meta==0,2,all))]
BAC_rel_abun_meta_filtered[is.na(BAC_rel_abun_meta_filtered)] <- 0
iva <- indval(BAC_rel_abun_meta_filtered[,1:2703], BAC_rel_abun_meta_filtered$Group)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(BAC_rel_abun_meta_filtered[,1:2703]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="indvalsummary.txt", sep = "\t")

#---
Selected_ITS <- read.csv(file="Selected_ITS_ASV_rel_abun.txt", sep="\t", row.names = 1, header = T)
