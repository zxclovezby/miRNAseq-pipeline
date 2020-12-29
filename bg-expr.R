library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
options(stringsAsFactors = F)

args <-commandArgs(trailingOnly=TRUE)

g1 <- args[1]
g2 <- args[2]
rda <- args[3]
num <- args[4]
Org <- args[5]
KEGG.org <- args[6]
html_tail <- args[7]
html_header <- args[8]
map <- args[9]

#load(rda)

abs_FC.f<-function(x){
  x<-x[1]
  if(x>=1){
  k<-x
  }else{
  k<- -(1/x)
  }
  k
  }

logFC.sta<-function(x){
  x=unlist(x)
  a<-mean(x[phenotype_table$group == Groups_un[2]])
  b<-mean(x[phenotype_table$group == Groups_un[1]])
  logFC<-log2(a+1)-log2(b+1)
  logFC
  }
	  
# pData(bg) = read.table("pData.txt", sep = '\t', header = T)
# transcript_fpkm = gexpr(bg)
# whole_tx_table = texpr(bg, 'all')[,-1]
# anno_data = whole_tx_table[,c(1:9)]
# idx = rowSums(transcript_fpkm > 0) >= num
# transcript_fpkm_n = transcript_fpkm[idx,]
# whole_tx_table_n = whole_tx_table[idx,]
# phenotype_table = pData(bg)
r = strsplit(rda, "\\.")
r_n = r[[1]][1]
# dir.create(r_n)
# write.csv(transcript_fpkm_n, paste0(r_n,"/",r_n,"_GeneExpression_Matrix.csv"), row.names=T, quote = F)
# write.csv(whole_tx_table_n, paste0(r_n,"/",r_n,"_TranscriptExpression_Matrix.csv"),quote = F, row.names = T)


# bg = exprfilter(bg, 0, meas = "FPKM")
# bg_C0vsC8 = subset(bg, "group == g1 | group == g2", genomesubset = F)
# stat_results = stattest(bg_C0vsC8, feature='gene', meas='FPKM', covariate='group', getFC=T)
# stat_results$log2FC = log2((stat_results$fc)+1)
# stat_results = stat_results[,c(-1,-3)]
# indices = match(stat_results$id, texpr(bg_C0vsC8, 'all')$gene_id)
# gene_names_for_result = texpr(bg, 'all')$gene_name[indices]
# results_genes = na.omit(data.frame(geneNames=gene_names_for_result,stat_results))
# write.csv(results_genes, paste0(r_n,"/",r_n,"_",g2,"vs",g1,".csv"), row.names = F, quote =F)


# qsig=which(results_genes$pval<0.05&results_genes$log2FC>=1)
# gene_sig = results_genes[qsig,]
# id = gene_sig$id
# write.csv(gene_sig, paste0(r_n,"/",r_n,"_",g2,"vs",g1,"_sig.csv"), row.names = F, quote =F)
# write.table(id, paste0(r_n,"/ID_",r_n,"_",g2,"vs",g1,"_sig.csv"), row.names = F, quote =F, col.names = F, sep = '\t')

#################################### ballgwon QC ploting ##########################################
###### All the Samples ######
# dir.create(paste0(r_n,'/Sample_QC'))
# setwd(paste0(r_n,'/Sample_QC'))
# gene_data = data.frame(gexpr(bg))
# sel<-apply(gene_data,1,max)
# gene_data = gene_data[sel>=1,]
# colnames(gene_data)<-gsub("FPKM[.]","",colnames(gene_data))
# box_data<-log2(gene_data)

# pdf("Boxplot.pdf")
# boxplot(box_data,col="skyblue",las=2)
# dev.off()
# png("Boxplot.png")
# boxplot(box_data,col="skyblue",las=2)
# dev.off()

# #######2.1 PCA分析

# PCA <- prcomp(t(gene_data))
# temp <- predict(PCA)
# n <- rownames(temp)
# n<-gsub("X","",n)
# temp <- data.frame(Sample=n,temp)
# PCA2 <- temp[,1:3] #二维
# PCA2 <- merge(PCA2,phenotype_table,by.x="Sample",by.y="id")
# PCA2$group=as.character(PCA2$group)
# library(ggplot2)
# library(Cairo)
# c <- ggplot()
# p <- c+
#      theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.border = element_rect(fill=NA))+
#      geom_point(data=PCA2,aes(x=PC1,y=PC2,colour = group))+
#      geom_text(data=PCA2,aes(x=PC1,y=PC2,label=Sample),hjust=0.5,vjust=-0.5)+
#      labs(x="PCA1",y="PCA2",title="PCA")
#      pdf("PCA.pdf") #PDF
#      plot(p)
#      dev.off()
#      CairoPNG("PCA.png") #png
#      plot(p)
#      dev.off()

# #####2.2 样品层次聚类

# library(cluster)
# dd <- dist(t(gene_data))
# hc <- hclust(dd)
# pdf("Cluster.pdf")
# plot(hc,hang = -1)
# dev.off()
# png("Cluster.png")
# plot(hc,hang = -1)
# dev.off()

# ####2.3 样品间相关性分布图

# library(lattice)
# cor_data <- as.matrix(cor(gene_data),method="person")
# col<- colorRampPalette(c('white','blue'))
# grid <- expand.grid(x=rownames(cor_data), y=colnames(cor_data))
# grid$z<-as.vector(cor_data)
# pdf("cor_plot.pdf")
# levelplot(cor_data,grid,scales=list(x=list(rot=90)),main="Pearson correlation between samples",col.regions=(col),aspect='fill',panel=function(...) {
# arg <- list(...)
# panel.levelplot(...)
# panel.text(arg$x, arg$y, round(arg$z,3))})
# dev.off()
# png("cor_plot.png")
# levelplot(cor_data,grid,scales=list(x=list(rot=90)),main="Pearson correlation between samples",col.regions=(col),
# aspect='fill',panel=function(...) {
# arg <- list(...)
# panel.levelplot(...)
# panel.text(arg$x, arg$y, round(arg$z,3))})
# dev.off()
# setwd("../../")

#################################### edgeR分析差异表达基因 ##########################################
# whole_gx_table = transcript_fpkm
# whole_gx_table = data.frame(gene_id=rownames(whole_gx_table),whole_gx_table)
# whole_gx_table = merge(anno_data,whole_gx_table,by="gene_id")
# whole_gx_table = whole_gx_table[!duplicated(whole_gx_table$gene_id),]

# dir.create(paste0(r_n,'/Diff_exp'))
# setwd(paste0(r_n,'/Diff_exp'))
# library(edgeR)
# library(Cairo)
VS=paste0(g2,"vs",g1)
#read.delim(paste0("../",r_n,"_gene_count_matrix.csv"),head=T,sep=",",row.names=1)->countTable
# geneHeat1=paste("Gene_Heatmap.",g2,"vs",g1,".pdf",sep="")
# geneHeat2=paste("Gene_Heatmap.",g2,"vs",g1,".png",sep="")
# geneVol1=paste("Gene_Volcanoplot.",g2,"vs",g1,".pdf",sep="")
# geneVol2=paste("Gene_Volcanoplot.",g2,"vs",g1,".png",sep="")
genefile1=paste("Gene_All.",g2,"vs",g1,".xls",sep="")
genefile2=paste("Gene_diff.",g2,"vs",g1,".xls",sep="")
# ##filter data##
# countTable<-countTable[,phenotype_table$id]
# keep <- rowSums(cpm(countTable)> 0.3) >= num
# countTable <- countTable[keep,]

# #######replicated####
# d <- DGEList(counts = countTable, group = phenotype_table$group)
# d <- calcNormFactors(d)
# d <- estimateCommonDisp(d)
# edgeR.de <- exactTest(d)
# ######################
# de <- decideTestsDGE(edgeR.de, adjust.method="fdr",p.value=0.05, lfc=1)

# cpm(d, log=TRUE, prior.count=2)->cpm
# DE=topTags(edgeR.de,adjust.method="fdr",n=400000)$table
# DE=data.frame(Gene_ID=rownames(DE),DE[,-2])
# All<-merge(DE,whole_gx_table,by.x="Gene_ID",by.y="gene_id")
# Groups_un<-unique(phenotype_table$group)

# All$logFC=apply(All[,c((ncol(All)-nrow(phenotype_table)+1):ncol(All))],1,logFC.sta)


# write.table(All,genefile1,sep="\t",quote = F,row.names = F)
# All$Diff<-"notsignificant"
# All$Diff[All$logFC >= 1 & All$FDR <=0.05]<-"UP"
# All$Diff[All$logFC <= -1 & All$FDR <=0.05]<-"Down"
# diff<-All[-which(All$Diff=="notsignificant"),]
# summary_gene <-data.frame(Down = length(All$Diff=="Down"),
#     	     Up = length(All$Diff=="UP"),
#       	     All= nrow(diff))
# write.table(diff,genefile2,sep="\t",quote = F,row.names = F)

# ######### no replicate ##############
# library(NOISeq)
# mydata <- readData(data = countTable, length = NULL, gc = NULL, biotype = NULL, chromosome = NULL, factors = phenotype_table)
# myresults <- noiseq(mydata, factor = "group", condition = c(g1, g2), k = NULL, norm = "n", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
# noiseq = myresults@results[[1]]
# noiseq_n = na.omit(noiseq)
# noiseq_n=noiseq_n[order(noiseq_n$prob, decreasing = T),]
# noiseq_n=data.frame(Gene_ID=rownames(noiseq_n),noiseq_n[,1:5])
# All_n<-merge(noiseq_n,whole_gx_table,by.x="Gene_ID",by.y="gene_id")
# colnames(All_n)[4] = 'logFC'
# All_n$FDR = 1-All_n$prob
# All = All_n[,c(1,4:5,ncol(All_n),2:3,6:7,11,8:10,12:"(ncol(All_n)-1))]


# write.table(All,genefile1,sep="\t",quote = F,row.names = F)
# All$Diff<-"notsignificant"
# All$Diff[All$logFC >= 2 & All$prob >=0.95]<-"UP"
# All$Diff[All$logFC <= -2 & All$prob >=0.95]<-"Down"
# diff<-All[-which(All$Diff=="notsignificant"),]
# summary_gene <-data.frame(Down = sum(All$Diff=="Down"),
#     	     Up = sum(All$Diff=="UP"),
#       	     All= nrow(diff))
# write.table(diff,genefile2,sep="\t",quote = F,row.names = F)


# ################################
# library(pheatmap)
# k<-ncol(diff)-1
# pdf(geneHeat1)
# heatmap_data<-diff[,c(paste0('FPKM.',phenotype_table$id[phenotype_table$group==g1]), paste0('FPKM.',phenotype_table$id[phenotype_table$group==g2]))]
# heatmap_data = heatmap_data[apply(heatmap_data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
# pheatmap(heatmap_data,color = colorRampPalette(c("red","white","blue"))(100),
#          fontsize=5,fontsize_row = 1,scale="row",cluster_cols=F)
# dev.off()
# CairoPNG(geneHeat2)
# pheatmap(heatmap_data,color = colorRampPalette(c("red","white","blue"))(100),fontsize=5,fontsize_row = 1,scale="row",cluster_cols=F)
# dev.off()
# ##############################
# k<-ncol(diff)
# data_plot <- All[,c(2,4,k)]
# data_plot$FDR<--log10(data_plot$FDR)
# colnames(data_plot)<-c("logFC","log10qval","Change")

# library(ggplot2)
# p<-ggplot()+
#    geom_point(data=data_plot,aes(x =logFC,y =log10qval,colour =Change))+#每个组的点单独画
#    #geom_text(data=data_1,aes(x =data_1$Group1,y =data_1$Group2,label=data_1$text), size=3,hjust = -0.1,vjust=1.1)+
#    #geom_text(data=data_2,aes(x =data_2$Group1,y =data_2$Group2,label=data_2$text), size=3,hjust = 1.1,vjust=-0.1)+#加文字注释
#    #geom_smooth(data=data_plot,aes(x =data_plot$Group1,y =data_plot$Group2),method = "lm", se=FALSE, color="black", formula = y ~ x)+#回归线
#    scale_colour_manual(values = c("green","grey","red"))+#颜色是所有颜色一起
#    labs(x="log2FC",y="-log10.Q_value")+ #x,y轴名字
#    geom_hline(yintercept = 1.3)+#加竖线
#    geom_vline(xintercept = 1)+
#    geom_vline(xintercept = -1)+   #加横线
#    theme(panel.grid.major=element_blank(),panel.background=element_blank())+#背景去边框，颜色
#    theme(axis.line = element_line(colour = "black"))+
#    scale_x_continuous(limits = c(-10,10))

# pdf(geneVol1,width = 12,height = 8)
# print(p)
# dev.off()
# CairoPNG(geneVol2,width = 1200,height = 800)
# print(p)
# dev.off()

# # #################################### edgeR分析差异表达 Trans ##########################################

# library(edgeR)
# library(Cairo)
# read.delim(paste0("../",r_n,"_transcript_count_matrix.csv"),head=T,sep=",",row.names=1)->countTable
# transHeat1=paste("transcript_Heatmap.",VS,".pdf",sep="")
# transHeat2=paste("transcript_Heatmap.",VS,".png",sep="")
# transVol1=paste("transcript_Volcanoplot.",VS,".pdf",sep="")
# transVol2=paste("transcript_Volcanoplot.",VS,".png",sep="")
transfile1=paste("transcript_All.",VS,".xls",sep="")
transfile2=paste("transcript_diff.",VS,".xls",sep="")

# ##filter data##
# countTable<-countTable[,phenotype_table$id]
# keep <- rowSums(cpm(countTable)> 0.3) >= num
# countTable <- countTable[keep,]

# #######replicated####

# d <- DGEList(counts = countTable, group = phenotype_table$group)
# d <- calcNormFactors(d)
# d <- estimateCommonDisp(d)
# edgeR.de <- exactTest(d)

# ######################

# de <- decideTestsDGE(edgeR.de, adjust.method="fdr",p.value=0.05, lfc=1)
# cpm(d, log=TRUE, prior.count=2)->cpm
# DE=topTags(edgeR.de,adjust.method="fdr",n=400000)$table
# DE=data.frame(transcript_ID=rownames(DE),DE[,-2])
# All<-merge(DE,whole_tx_table,by.x="transcript_ID",by.y="t_name")
# Groups_un<-unique(phenotype_table$group)
# All$logFC=apply(All[,seq(from = 14, to = ncol(All), by = 2)],1,logFC.sta)
# write.table(All,transfile1,sep="\t",quote = F,row.names = F)
# All$Diff<-"notsignificant"
# All$Diff[All$logFC >= 1 & All$FDR <=0.05]<-"UP"
# All$Diff[All$logFC <= -1 & All$FDR <=0.05]<-"Down"
# diff<-All[-which(All$Diff=="notsignificant"),]
# summary_transcript <-data.frame(Down = length(All$Diff=="Down"),Up = length(All$Diff=="UP"),All= nrow(diff))
# write.table(diff,transfile2,sep="\t",quote = F,row.names = F)
# summary_all <- rbind(summary_gene,summary_transcript)
# summary_all$type<-c('Gene','Transcript')
# summary_all<-summary_all[,c(4,1,2,3)]
# write.table(summary_all,paste0('../Summary_Diff.',VS,'.xls'),sep='\t',row.names = F)


# ####### no replicate #########
# mydata <- readData(data = countTable, length = NULL, gc = NULL, biotype = NULL, chromosome = NULL, factors = phenotype_table)
# myresults <- noiseq(mydata, factor = "group", conditions = c(g1, g2), k = NULL, norm = "n", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
# noiseq = myresults@results[[1]]
# noiseq_n = na.omit(noiseq)
# noiseq_n=noiseq_n[order(noiseq_n$prob, decreasing = T),]
# noiseq_n=data.frame(transcript_ID=rownames(noiseq_n),noiseq_n[,1:5])
# All_n<-merge(noiseq_n,whole_tx_table,by.x="transcript_ID",by.y="t_name")
# colnames(All_n)[4] = 'logFC'
# All_n$FDR = 1-All_n$prob
# All = All_n[,c(1,4:5,ncol(All_n),2:3,6:7,11,8:10,12:"(ncol(All_n)-1))]

# write.table(All,transfile1,sep="\t",quote = F,row.names = F)
# All$Diff<-"notsignificant"
# All$Diff[All$logFC >= 2 & All$prob >=0.95]<-"UP"
# All$Diff[All$logFC <= -2 & All$prob >=0.95]<-"Down"
# diff<-All[-which(All$Diff=="notsignificant"),]
# summary_transcript <-data.frame(Down = sum(All$Diff=="Down"),
#     	     Up = sum(All$Diff=="UP"),
#       	     All= nrow(diff))
# write.table(diff,transfile2,sep="\t",quote = F,row.names = F)
# summary_all <- rbind(summary_gene,summary_transcript)
# summary_all$type<-c('Gene','Transcript')
# summary_all<-summary_all[,c(4,1,2,3)]
# write.table(summary_all,paste0('../Summary_Diff.',VS,'.xls'),sep='\t',row.names = F)

# ################################
# library(pheatmap)
# k<-ncol(diff)-1
# pdf(transHeat1)
# heatmap_data<-diff[,c(paste0('FPKM.',phenotype_table$id[phenotype_table$group==g1]), paste0('FPKM.',phenotype_table$id[phenotype_table$group==g2]))]
# heatmap_data = heatmap_data[apply(heatmap_data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
# pheatmap(heatmap_data,color = colorRampPalette(c("red","white","blue"))(100),fontsize=5,fontsize_row = 1,scale="row",cluster_cols=F)
# dev.off()
# CairoPNG(transHeat2)
# pheatmap(heatmap_data,color = colorRampPalette(c("red","white","blue"))(100),fontsize=5,fontsize_row = 1,scale="row",cluster_cols=F)
# dev.off()
# ##############################
# k<-ncol(diff)
# data_plot <- All[,c(2,4,k)]
# data_plot$FDR<--log10(data_plot$FDR)
# colnames(data_plot)<-c("logFC","log10qval","Change")

# library(ggplot2)
# p<-ggplot()+
#   geom_point(data=data_plot,aes(x =logFC,y =log10qval,colour =Change))+#每个组的点单独画
#   #geom_text(data=data_1,aes(x =data_1$Group1,y =data_1$Group2,label=data_1$text), size=3,hjust = -0.1,vjust=1.1)+
#   #geom_text(data=data_2,aes(x =data_2$Group1,y =data_2$Group2,label=data_2$text), size=3,hjust = 1.1,vjust=-0.1)+#加文字注释
#   #geom_smooth(data=data_plot,aes(x =data_plot$Group1,y =data_plot$Group2),method = "lm", se=FALSE, color="black", formula = y ~ x)+#回归线
#   scale_colour_manual(values = c("green","grey","red"))+#颜色是所有颜色一起
#   labs(x="log2FC",y="-log10.Q_value")+ #x,y轴名字
#   geom_hline(yintercept = 1.3)+#加竖线
#   geom_vline(xintercept = 1)+
#   geom_vline(xintercept = -1)+   #加横线
#   theme(panel.grid.major=element_blank(),panel.background=element_blank())+#背景去边框，颜色
#   theme(axis.line = element_line(colour = "black"))+
#   scale_x_continuous(limits = c(-10,10))

# pdf(transVol1,width = 12,height = 8)
# print(p)
# dev.off()
# CairoPNG(transVol2,width = 1200,height = 800)
# print(p)
# dev.off()
# setwd("../../")
############################################# ENRICHMENT Gene ###################################################

dir.create(paste0(r_n,'/Enrichment'))
setwd(paste0(r_n,'/Enrichment'))
library('clusterProfiler')
library('annotate')
library('ggplot2')
library('topGO')
library('Rgraphviz')
library('AnnotationHub')
library(KEGG.db)
library(paste0('org.',Org,'.eg.db'),character.only = TRUE)
org <- paste0('org.',Org,'.eg.db')
data <-read.delim(paste0('../Diff_exp/',genefile2),sep='\t',stringsAsFactors = F, na.strings = ".")
#data.1<- na.omit(data)
data.1 = data[data$Gene_ID %like% "ENSG",]
ids <- bitr(data.1$Gene_ID, fromType="ENSEMBL", toType= "ENTREZID", OrgDb=org)
#ids <- bitr(data.1$Gene_ID, fromType="SYMBOL", toType= "ENTREZID", OrgDb=org)
enrichGO.MF <- enrichGO(gene = ids$ENTREZID, OrgDb=org,ont = "MF", qvalueCutoff=0.05, readable =T)
enrichGO.BP <- enrichGO(gene = ids$ENTREZID, OrgDb=org,ont = "BP", qvalueCutoff = 0.05, readable =T)
enrichGO.CC <- enrichGO(gene = ids$ENTREZID, OrgDb=org,ont = "CC", qvalueCutoff = 0.05, readable =T)
dir.create(paste('Gene_diff',VS,sep = '.'))
setwd(paste('Gene_diff',VS,sep = '.'))
dir.create('GO')
setwd('GO')

#####修改doplot函数，使其更美观
dotplot_new <- function(object, x="geneRatio", colorBy="qvalue", showCategory=10, split=NULL, font.size=11, title="") {
	    colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
	    if (x == "geneRatio" || x == "GeneRatio") {
	    x <- "GeneRatio"
	    size <- "Count"
	    } else if (x == "count" || x == "Count") {
            x <- "Count"
	    size <- "GeneRatio"
	    } else {
	    stop("x should be geneRatio or count...")
	    }
	    df <- fortify(object, showCategory = showCategory, split=split)	## already parsed in fortify
	    ## df$GeneRatio <- parse_ratio(df$GeneRatio)
	    idx <- order(df$GeneRatio, decreasing = FALSE)
	    df$Description <- factor(df$Description, levels=df$Description[idx])
	    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
	       geom_point() +
	       scale_color_gradient2(low = "red", high = "green",mid = "black", midpoint = (min(df$qvalue)+max(df$qvalue))/2, space = "rgb", guide = "colourbar", breaks=seq(min(df$qvalue),max(df$qvalue),(max(df$qvalue)-min(df$qvalue))/4)) +
	       ylab("") +
	       ggtitle(title) +
	       theme_dose(font.size)
	    }

pdf(file='dotplot_MF_of_GO_result.pdf',width = 12,height = 14)
dotplot_new(enrichGO.MF,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_MF_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.MF,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
pdf(file='MF_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.MF)
dev.off()
png(file='MF_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.MF)
dev.off()

####输出BP水平的气泡图（如果报错请查看BP水平富集分析是否不存在显著结果，下同）
pdf(file='dotplot_BP_of_GO_result.pdf',width = 15,height = 15)
dotplot_new(enrichGO.BP,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_BP_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.BP,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()

#####输出BP水平的有向无环图
pdf(file='BP_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.BP)
dev.off()
png(file='BP_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.BP)
dev.off()

#####输出CC水平的气泡图
pdf(file='dotplot_CC_of_GO_result.pdf',width = 15,height = 15)
dotplot_new(enrichGO.CC,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_CC_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.CC,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()

#####输出CC水平的有向无环图
pdf(file='CC_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.CC,firstSigNodes=10)
dev.off()
png(file='CC_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.CC,firstSigNodes=10)
dev.off()

######准备做GO分析条形图的数据框
bar.cc<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "CC", level = 2,readable = TRUE)
bar.BP<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "BP", level = 2,readable = TRUE)
bar.MF<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "MF", level = 2,readable = TRUE)
GO.histo<-rbind(bar.cc@result,bar.BP@result,bar.MF@result)
type<-c(rep('Cellular component',nrow(bar.cc@result)),rep('Biological process',nrow(bar.BP@result)),rep('Molecular function',nrow(bar.MF@result)))
GO.histo$type<-type
GO.histo.1<-GO.histo[nrow(GO.histo):1,]
GO.histo.1<-GO.histo.1[which(GO.histo.1$Count!=0),]

######输出GO分析条形图
DES_level <-factor((1:nrow(GO.histo.1)),labels=GO.histo.1$Description)
pdf('hist_GO_result.pdf',width = 14,height = 12)
ggplot(data=GO.histo.1,aes(x=DES_level,y=GO.histo.1$Count))+geom_bar(stat="identity",aes(fill=type))+
	labs(x='GO term',y='Number of genes',title = "GO Terms")+
	coord_flip()+
	theme_bw()+
	theme(panel.border=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	plot.title=element_text(hjust = 0.5))
dev.off()
png('hist_GO_result.png',width = 2400,height = 1280)

ggplot(data=GO.histo.1,aes(x=DES_level,y=GO.histo.1$Count))+geom_bar(stat="identity",aes(fill=type))+
	labs(x='GO term',y='Number of genes',title = "GO Terms")+
	coord_flip()+
	theme_bw()+
	theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))
dev.off()
write.table(GO.histo.1,'go_level2_annotation.xls',sep = '\t',row.names = F)

######################
GO.result<-rbind(enrichGO.BP@result,enrichGO.CC@result,enrichGO.MF@result)
GO.result$ontology<-c(rep('Biological process',nrow(enrichGO.BP@result)),rep('Cellular component',nrow(enrichGO.CC@result)),rep('Molecular function',nrow(enrichGO.MF@result)))

#######################################################################
##########  此部分代码为替换富集结果中的gene id为gene symbol ##########
#######################################################################

#id.anno<-merge(ids,data.1,by.x='ENSEMBL',by.y='gene_id')
#for(i in 1:nrow(id.anno)){
#  GO.result$geneID<-gsub(id.anno$ENTREZID[i],id.anno$gene_name[i],GO.result$geneID)
#}
write.table(GO.result,'go_enrich_result.xls',sep = '\t',row.names = F)

########################################################################
############################  KEGG  ####################################
########################################################################

dir.create(paste0(r_n,'/KEGG'))
setwd(paste0(r_n, '/KEGG'))

ids <- bitr(data.1$Gene_ID, fromType="SYMBOL", toType= "UNIPROT", OrgDb=org)
KEGG <- enrichKEGG(gene = ids$UNIPROT,organism = KEGG.org,keyType = 'uniprot', qvalueCutoff=0.05)

#####输出KEGG的气泡图（如果报错请查看BP水平富集分析是否不存在显著结果，下同）
pdf(file='dotplot_of_KEGG_result.pdf',width = 15,height = 15)
dotplot_new(KEGG,colorBy = "qvalue",title='Statistics of KEGG enrichment')
dev.off()
png(file='dotplot_of_KEGG_result.png',width = 980,height = 1460)
dotplot_new(KEGG,colorBy = "qvalue",title='Statistics of KEGG enrichment')
dev.off()

############### level2 水平柱状图
KEGG.hist <- enrichKEGG(gene = ids$UNIPROT,organism = KEGG.org,keyType = 'uniprot', qvalueCutoff=1,pvalueCutoff = 1)
KEGG.histo.1<-KEGG.hist@result
anno<-read.delim(map,stringsAsFactors = F)
see<-merge(KEGG.histo.1,anno,by='ID',all=F)
KEGG.histo.2<-data.frame(table(see$B))
for(i in 1:nrow(KEGG.histo.2)){
      temp<-strsplit(paste(see[which(see$B==KEGG.histo.2$Var1[i]),]$geneID,collapse = '/'),'/')[[1]]
      temp<-unique(temp)
      KEGG.histo.2$Freq[i]<-length(temp)
      KEGG.histo.2$Protein_uniprot_id[i]<-paste(temp,collapse = '/')
      }

DES_level <-factor((1:nrow(KEGG.histo.2)),labels=KEGG.histo.2$Var1)
pdf('hist_KEGG_result.pdf',width = 10,height = 10)
ggplot(data=KEGG.histo.2,aes(x=DES_level,y=KEGG.histo.2$Freq))+geom_bar(stat="identity")+
	labs(x='KEGG term',y='Number of genes',title = "The KEGG Terms")+
	coord_flip()+
	theme_bw()+
	theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))
dev.off()

png('hist_KEGG_result.png',width = 980,height = 1460)
ggplot(data=KEGG.histo.2,aes(x=DES_level,y=KEGG.histo.2$Freq))+geom_bar(stat="identity")+
  labs(x='KEGG term',y='Number of genes',title = "The KEGG Terms")+
  coord_flip()+
  theme_bw()+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))
dev.off()

colnames(KEGG.histo.2)[1:2]<-c('KEGG Term','Protein numbers')
write.table(KEGG.histo.2,'KEGG_level2_annotation.xls',sep = '\t',row.names = F)

########### 生成网络图片
result<-KEGG@result
result$map<-apply(result,1,function(x){
	up<-unlist(strsplit(x[8],'/'))[which(unlist(strsplit(x[8],'/'))%in%ids$UNIPROT[which(ids$ENSEMBL%in%data.1$Gene_ID[data.1$Diff=='UP'])])]
	up<-as.character(bitr(up, fromType="UNIPROT", toType= "ENTREZID", OrgDb=org)$ENTREZID)
	down<-unlist(strsplit(x[8],'/'))[which(unlist(strsplit(x[8],'/'))%in%ids$UNIPROT[which(ids$SYMBOL%in%data.1$gene_name[data.1$Diff=='Down'])])]
	down<-as.character(bitr(down, fromType="UNIPROT", toType= "ENTREZID", OrgDb=org)$ENTREZID)
	gene1 = paste0(up,'%09red',collapse = '/')
	gene2 = paste0(down,'%09green',collapse = '/')
	return(
	paste0("http://www.kegg.jp/kegg-bin/show_pathway?@",x[1],"/reference%3dwhite/",gene1,'/',gene2))								   })

#######################################################################
##########  此部分代码为替换富集结果中的gene id为gene symbol ##########
#######################################################################

id.anno<-merge(ids,data.1,by.x='ENSEMBL',by.y='gene_id')
#for(i in 1:nrow(id.anno)){
#  result$geneID<-gsub(id.anno$UNIPROT[i],id.anno$gene_name[i],result$geneID)
#}
#######################################################################

write.table(result,'KEGG_enrich_result.xls',sep = '\t',row.names = F)

##############生成网页报告
header<-readLines(html_header)
header<-paste(header,"\n",sep = "",collapse = "")
header.1 = "<h1> Pathway Mapping</h1>\n   <h2>result\n</h2>\n<p class='fst'> <p class='snd'> </p>\n"
header.2 = "<p><p></p><figcaption>Table KEGG pathway enrichment.</figcaption>\n"
header.3 = paste0('<div class="tableHeader">\n<table class="table"><thead><tr>',paste("<th>",colnames(result)[-10],"</th>",sep="",collapse = ""),"</tr></thead>\n<tbody>\n")
make_report = function(x){
  temp = paste0("<a href=",x[10]," target='new'>",x[1],"</a> ")
  temp2 = paste("<td>",c(temp,x[2:9]),"</td>",sep = "",collapse = "")
  return(paste("<tr>",temp2,"</tr>\n",sep = ""))
  }

get_report<-apply(result, 1, make_report)
tail.1 = readLines(html_tail)
tail.1<-paste(tail.1,collapse = '')
report = paste0(header,header.1,header.2,header.3,paste(get_report[1:length(get_report)],sep='',collapse = ""),tail.1)
write.table(report,"kegg.html",row.names = F,col.names = F,sep = "",quote=F)

############################################# ENRICHMENT Trans ###################################################

dir.create(paste('Transcript_diff',VS,sep = '.'))
setwd(paste('Transcript_diff',VS,sep = '.'))
data <-read.delim(paste0('../../Diff_exp/',transfile2),sep='\t',stringsAsFactors = F)
data.1 = data[data$transcript_ID %like% "ENST",]
ids <- bitr(data.1$transcript_ID, fromType="ENSEMBLTRANS", toType= "ENTREZID", OrgDb=org)
enrichGO.MF <- enrichGO(gene = ids$ENTREZID,OrgDb=org,ont = "MF", qvalueCutoff=0.05, readable =T)
enrichGO.BP <- enrichGO(gene = ids$ENTREZID,  OrgDb=org,ont = "BP", qvalueCutoff = 0.05, readable =T)
enrichGO.CC <- enrichGO(gene = ids$ENTREZID,  OrgDb=org,ont = "CC", qvalueCutoff = 0.05, readable =T)
dir.create('GO')
setwd('GO')
#修改doplot函数，使其更美观
dotplot_new <- function(object, x="geneRatio", colorBy="qvalue", showCategory=10, split=NULL, font.size=11, title="") {
	    colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
	    if (x == "geneRatio" || x == "GeneRatio") {
            x <- "GeneRatio"
	    size <- "Count"
	    } else if (x == "count" || x == "Count") {
            x <- "Count"
	    size <- "GeneRatio"
	    } else {
	    stop("x should be geneRatio or count...")
	    }
	    df <- fortify(object, showCategory = showCategory, split=split)  ## already parsed in fortify
	    ## df$GeneRatio <- parse_ratio(df$GeneRatio)
	    idx <- order(df$GeneRatio, decreasing = FALSE)
	    df$Description <- factor(df$Description, levels=df$Description[idx])
	    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
            	       geom_point() + scale_color_gradient2(low = "red", high = "green",mid = "black", midpoint = (min(df$qvalue)+max(df$qvalue))/2, space = "rgb", guide = "colourbar", breaks=seq(min(df$qvalue),max(df$qvalue),(max(df$qvalue)-min(df$qvalue))/4)) +
		       ylab("") + ggtitle(title) + theme_dose(font.size)
	}

pdf(file='dotplot_MF_of_GO_result.pdf',width = 12,height = 14)
dotplot_new(enrichGO.MF,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_MF_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.MF,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
pdf(file='MF_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.MF)
dev.off()
png(file='MF_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.MF)
dev.off()

####输出BP水平的气泡图（如果报错请查看BP水平富集分析是否不存在显著结果，下同）
pdf(file='dotplot_BP_of_GO_result.pdf',width = 15,height = 15)
dotplot_new(enrichGO.BP,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_BP_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.BP,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()

#####输出BP水平的有向无环图
pdf(file='BP_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.BP)
dev.off()
png(file='BP_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.BP)
dev.off()

#####输出CC水平的气泡图
pdf(file='dotplot_CC_of_GO_result.pdf',width = 15,height = 15)
dotplot_new(enrichGO.CC,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
png(file='dotplot_CC_of_GO_result.png',width = 980,height = 1460)
dotplot_new(enrichGO.CC,colorBy = "qvalue",title='Statistics of GO enrichment')
dev.off()
######输出CC水平的有向无环图
pdf(file='CC_DAG_plot_of_GO_result.pdf',width = 15,height = 15)
plotGOgraph(enrichGO.CC,firstSigNodes=10)
dev.off()
png(file='CC_DAG_plot_of_GO_result.png',width = 1280,height = 1280)
plotGOgraph(enrichGO.CC,firstSigNodes=10)
dev.off()
######准备做GO分析条形图的数据框
bar.cc<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "CC", level = 2,readable = T)
bar.BP<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "BP", level = 2,readable = TRUE)
bar.MF<-groupGO(as.character(ids$ENTREZID), OrgDb=org, ont = "MF", level = 2,readable = TRUE)
GO.histo<-rbind(bar.cc@result,bar.BP@result,bar.MF@result)
type<-c(rep('Cellular component',nrow(bar.cc@result)),rep('Biological process',nrow(bar.BP@result)),rep('Molecular function',nrow(bar.MF@result)))
GO.histo$type<-type
GO.histo.1<-GO.histo[nrow(GO.histo):1,]
GO.histo.1<-GO.histo.1[which(GO.histo.1$Count!=0),]

p#####输出GO分析条形图
DES_level <-factor((1:nrow(GO.histo.1)),labels=GO.histo.1$Description)
pdf('hist_GO_result.pdf',width = 14,height = 12)
ggplot(data=GO.histo.1,aes(x=DES_level,y=GO.histo.1$Count))+geom_bar(stat="identity",aes(fill=type))+
	labs(x='GO term',y='Number of genes',title = "GO Terms")+
	coord_flip()+theme_bw()+theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))

dev.off()
png('hist_GO_result.png',width = 2400,height = 1280)
ggplot(data=GO.histo.1,aes(x=DES_level,y=GO.histo.1$Count))+geom_bar(stat="identity",aes(fill=type))+
	labs(x='GO term',y='Number of genes',title = "GO Terms")+
	coord_flip()+theme_bw()+theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))

dev.off()
write.table(GO.histo.1,'go_level2_annotation.xls',sep = '\t',row.names = F)

######################
GO.result<-rbind(enrichGO.BP@result,enrichGO.CC@result,enrichGO.MF@result)
GO.result$ontology<-c(rep('Biological process',nrow(enrichGO.BP@result)),rep('Cellular component',nrow(enrichGO.CC@result)),rep('Molecular function',nrow(enrichGO.MF@result)))

#######################################################################
##########  此部分代码为替换富集结果中的gene id为gene symbol ##########
#######################################################################
#id.anno<-merge(ids,data.1,by.x='ENSEMBL',by.y='gene_id')
#for(i in 1:nrow(id.anno)){
#  GO.result$geneID<-gsub(id.anno$ENTREZID[i],id.anno$gene_name[i],GO.result$geneID)
#}
write.table(GO.result,'go_enrich_result.xls',sep = '\t',row.names = F)
###################################################################################
##################   KEGG
###################################################################################

dir.create('../KEGG')
setwd('../KEGG')
ids <- bitr(data.1$transcript_ID, fromType="ENSEMBLTRANS", toType= "UNIPROT", OrgDb=org)
KEGG <- enrichKEGG(gene = ids$UNIPROT,organism = KEGG.org,keyType = 'uniprot', qvalueCutoff=0.05)

####输出KEGG的气泡图（如果报错请查看BP水平富集分析是否不存在显著结果，下同）
pdf(file='dotplot_of_KEGG_result.pdf',width = 15,height = 15)
dotplot_new(KEGG,colorBy = "qvalue",title='Statistics of KEGG enrichment')
dev.off()
png(file='dotplot_of_KEGG_result.png',width = 980,height = 1460)
dotplot_new(KEGG,colorBy = "qvalue",title='Statistics of KEGG enrichment')
dev.off()

############### level2 水平柱状图
KEGG.hist <- enrichKEGG(gene = ids$UNIPROT,organism = KEGG.org,keyType ='uniprot', qvalueCutoff=1,pvalueCutoff = 1)
KEGG.histo.1<-KEGG.hist@result
anno<-read.delim(map,stringsAsFactors = F)
see<-merge(KEGG.histo.1,anno,by='ID',all=F)
KEGG.histo.2<-data.frame(table(see$B))
for(i in 1:nrow(KEGG.histo.2)){
  temp<-strsplit(paste(see[which(see$B==KEGG.histo.2$Var1[i]),]$geneID,collapse = '/'),'/')[[1]]
  temp<-unique(temp)
  KEGG.histo.2$Freq[i]<-length(temp)
  KEGG.histo.2$Protein_uniprot_id[i]<-paste(temp,collapse = '/')
  }

DES_level <-factor((1:nrow(KEGG.histo.2)),labels=KEGG.histo.2$Var1)
pdf('hist_KEGG_result.pdf',width = 10,height = 10)
ggplot(data=KEGG.histo.2,aes(x=DES_level,y=KEGG.histo.2$Freq))+geom_bar(stat="identity")+
	labs(x='KEGG term',y='Number of genes',title = "The KEGG Terms")+
	coord_flip()+theme_bw()+theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))

dev.off()
png('hist_KEGG_result.png',width = 980,height = 1460)
ggplot(data=KEGG.histo.2,aes(x=DES_level,y=KEGG.histo.2$Freq))+geom_bar(stat="identity")+
	labs(x='KEGG term',y='Number of genes',title = "The KEGG Terms")+
	coord_flip()+theme_bw()+theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), plot.title=element_text(hjust = 0.5))

dev.off()
colnames(KEGG.histo.2)[1:2]<-c('KEGG Term','Protein numbers')
write.table(KEGG.histo.2,'KEGG_level2_annotation.xls',sep = '\t',row.names = F)

# ########### 生成网络图片
# result<-KEGG@result
# result$map<-apply(result,1,function(x){
# 	up<-unlist(strsplit(x[8],'/'))[which(unlist(strsplit(x[8],'/'))%in%ids$UNIPROT[which(ids$ENSEMBLTRANS%in%data.1$transcript_ID[which(data.1$Diff=='UP')])])]
# 	up<-as.character(bitr(up, fromType="UNIPROT", toType= "ENTREZID", OrgDb=org)$ENTREZID)
# 	down<-unlist(strsplit(x[8],'/'))[which(unlist(strsplit(x[8],'/'))%in%ids$UNIPROT[which(ids$ENSEMBLTRANS%in%data.1$transcript_ID[which(data.1$Diff=='Down')])])]
# 	down<-as.character(bitr(down, fromType="UNIPROT", toType= "ENTREZID", OrgDb=org)$ENTREZID)
# 	gene1 = paste0(up,'%09red',collapse = '/')
# 	gene2 = paste0(down,'%09green',collapse = '/')
# 	return(
# 	paste0("http://www.kegg.jp/kegg-bin/show_pathway?@",x[1],"/reference%3dwhite/",gene1,'/',gene2)
# 	)
# })

#######################################################################
##########  此部分代码为替换富集结果中的gene id为gene symbol ##########
#######################################################################
#id.anno<-merge(ids,data.1,by.x='ENSEMBL',by.y='gene_id')
#for(i in 1:nrow(id.anno)){
#  result$geneID<-gsub(id.anno$UNIPROT[i],id.anno$gene_name[i],result$geneID)
#}
#######################################################################
write.table(result,'KEGG_enrich_result.xls',sep = '\t',row.names = F)
