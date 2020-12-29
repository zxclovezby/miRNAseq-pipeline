library(topGO)

args <-commandArgs(trailingOnly=TRUE)
sample <-args[1]
db<-args[2]
significant<-as.numeric(args[3])

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
all_genes<-cuffdiff$pval
names(all_genes)<-cuffdiff$geneNames

GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p <significant, description = sample, annot = annFUN.org, mapping = paste0("org.",db,".eg.db"), ID = "symbol")
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p <significant, description = sample, annot = annFUN.org, mapping = paste0("org.",db,".eg.db"), ID = "symbol")
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(p) p <significant, description = sample, annot = annFUN.org, mapping = paste0("org.",db,".eg.db"), ID = "symbol")
resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")

BP_term<-usedGO(GOdata_BP)
MF_term<-usedGO(GOdata_MF)
CC_term<-usedGO(GOdata_CC)
BP_all<-GenTable(GOdata_BP, classicFisher = resultFisher_BP, topNodes = length(BP_term))
MF_all<-GenTable(GOdata_MF, classicFisher = resultFisher_MF, topNodes = length(MF_term))
CC_all<-GenTable(GOdata_CC, classicFisher = resultFisher_CC, topNodes = length(CC_term))
BP_all[,7]<-p.adjust(BP_all$classicFisher, method ='BH')
MF_all[,7]<-p.adjust(MF_all$classicFisher, method ='BH')
CC_all[,7]<-p.adjust(CC_all$classicFisher, method ='BH')
BP_ann<-genesInTerm(GOdata_BP, BP_term)
MF_ann<-genesInTerm(GOdata_MF, MF_term)
CC_ann<-genesInTerm(GOdata_CC, CC_term)
BP_all[,8]<-BP_all$Significant/BP_all$Expected
MF_all[,8]<-MF_all$Significant/MF_all$Expected
CC_all[,8]<-CC_all$Significant/CC_all$Expected

sigGene_BP<-sigGenes(GOdata_BP)
sigGene_MF<-sigGenes(GOdata_MF)
sigGene_CC<-sigGenes(GOdata_CC)
paste0(sample,"_BPsigGene :",length(sigGene_BP))
paste0(sample,"_MFsigGene :",length(sigGene_MF))
paste0(sample,"_CCsigGene :",length(sigGene_CC))

test_BP<-sapply(BP_ann,function(x){intersect(x,sigGene_BP)})
test2_BP=sapply(test_BP,function(x){paste(x,collapse=',')})
test3_BP=data.frame(GO=names(test2_BP),GeneID=test2_BP,stringsAsFactors =F)

test_MF<-sapply(MF_ann,function(x){intersect(x,sigGene_MF)})
test2_MF=sapply(test_MF,function(x){paste(x,collapse=',')})
test3_MF=data.frame(GO=names(test2_MF),GeneID=test2_MF,stringsAsFactors =F)

test_CC<-sapply(CC_ann,function(x){intersect(x,sigGene_CC)})
test2_CC=sapply(test_CC,function(x){paste(x,collapse=',')})
test3_CC=data.frame(GO=names(test2_CC),GeneID=test2_CC,stringsAsFactors =F)

m_BP<-match(BP_all$GO.ID,test3_BP$GO)
m_MF<-match(MF_all$GO.ID,test3_MF$GO)
m_CC<-match(CC_all$GO.ID,test3_CC$GO)
BP_all[,9]<-test3_BP$GeneID[m_BP]
MF_all[,9]<-test3_MF$GeneID[m_MF]
CC_all[,9]<-test3_CC$GeneID[m_CC]
colnames(BP_all)<-c("GO.ID","Term","Annotated","Significant","Expected","p-value","FDR-value","Enrichment","Genes")
colnames(MF_all)<-c("GO.ID","Term","Annotated","Significant","Expected","p-value","FDR-value","Enrichment","Genes")
colnames(CC_all)<-c("GO.ID","Term","Annotated","Significant","Expected","p-value","FDR-value","Enrichment","Genes")
BP_all_sig<-BP_all[which(BP_all[,7]<significant),]
MF_all_sig<-MF_all[which(MF_all[,7]<significant),]
CC_all_sig<-CC_all[which(CC_all[,7]<significant),]
write.table( BP_all_sig , file =paste0("./BP/",sample,"_BP.tsv") , sep = "\t" , row.names = F )
write.table( MF_all_sig , file =paste0("./MF/",sample,"_MF.tsv") , sep = "\t" , row.names = F )
write.table( CC_all_sig , file =paste0("./CC/",sample,"_CC.tsv") , sep = "\t" , row.names = F )
