#-- ITS relative abundance already calculated
ITS_rel_abun <- read.csv(file = "../ITS_ASV_filtered_Rel_Abun.txt", sep = "\t", row.names = 1, header=T)
 
ITSmeta <- read.csv("../../Metadata_Filtered.txt", sep = "\t", row.names = 1)
ITSmeta_t <- data.frame(t(ITSmeta))
ITSmeta <- data.frame(t(ITSmeta_t))
ITSmeta <- ITSmeta[order(rownames(ITSmeta)),]

ITS_rel_abun_meta <- merge(ITS_rel_abun, ITSmeta, by=0, all=F)
rownames(ITS_rel_abun_meta) <- ITS_rel_abun_meta$Row.names; ITS_rel_abun_meta$Row.names <- NULL

library (labdsv)
ITS_rel_abun_meta_filtered <- ITS_rel_abun_meta[ , which(!apply(ITS_rel_abun_meta==0,2,all))]
ITS_rel_abun_meta_filtered[is.na(ITS_rel_abun_meta_filtered)] <- 0
iva <- indval(ITS_rel_abun_meta_filtered[,1:576], ITS_rel_abun_meta_filtered$Group)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(ITS_rel_abun_meta_filtered[,1:576]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="indvalsummary.txt", sep = "\t")
