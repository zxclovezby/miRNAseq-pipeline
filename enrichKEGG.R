library(clusterProfiler)
library(DOSE)

#$comparison $annotation_lib $organism $pq_cutoff
args <-commandArgs(trailingOnly=TRUE)
sample <-args[1]
db<-args[2]
org<-args[3]
significant<-as.numeric(args[4])

#retrieve genes name from cuffdiff gene_exp.diff file
cuffdiff<-read.csv(file=paste0("./",sample),header=T,stringsAsFactors=F)
#cuffdiff[,3]<-sapply(cuffdiff$gene,function(x){strsplit(x,',')[[1]][1]})
cuffdiff<-cuffdiff[-c(which(cuffdiff$geneNames==".")),]
cuffdiff<-cuffdiff[order(cuffdiff$pval),]
duplicated<-duplicated(cuffdiff$geneNames)
cuffdiff<-cuffdiff[!duplicated,]

library(paste0("org.",db,".eg.db"),character.only = TRUE)
egSYMBOL2EG<-toTable(get(paste0("org.",db,".egSYMBOL2EG")))

match_ez<-match(cuffdiff$geneNames,egSYMBOL2EG$symbol)
summary(is.na(match_ez))
cuffdiff$entrezID<-egSYMBOL2EG$gene_id[match_ez]
na_row<-which(is.na(match_ez))
cuffdiff<-cuffdiff[-c(na_row),]

sig_genes<-unlist(subset(cuffdiff,pval <= 0.05,select = 'entrezID'))
kegg<-enrichKEGG(sig_genes, organism = org, pvalueCutoff = significant, pAdjustMethod = "BH", qvalueCutoff = significant, use_internal_data = FALSE)
#kegg2<-setReadable(kegg, OrgDb = paste0("org.",db,".eg.db"))
write.table(summary(kegg) , file =paste0(sample,"_KEGG.tsv") , sep = "\t" , row.names = F )
