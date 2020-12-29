library(cummeRbund)

cuff<-readCufflinks()
g1<-samples(genes(cuff))[1]
g2<-samples(genes(cuff))[2]

#Density Plot
png(filename ="images/genes_dens.png", width = 768, height = 768)
csDensity(genes(cuff),replicates=T)
dev.off()

#Box Plot
png(filename ="images/genes_bplot.png", width = 768, height = 768)
csBoxplot(genes(cuff),replicates=T)
dev.off()

#Volcano Plot
png(filename ="images/genes_volc.png", width = 768, height = 768)
csVolcano(genes(cuff),g1,g2,alpha=0.05, showSignificant=TRUE)
dev.off()

#Scatter Plot
png(filename ="images/genes_scatp.png", width = 768, height = 768)
csScatter(genes(cuff),g1,g2,smooth=T)
dev.off()

#Heatmap
genes_diff<- diffData(genes(cuff))
sig_genes<- subset(genes_diff, (significant ==  'yes'))
top_50<-sig_genes[order(sig_genes$q_value,decreasing=T)[1:50],]
myGenes <- getGenes(cuff, top_50$gene_id)
png(filename ="images/genes_heatm.png", width = 768, height = 768)
csHeatmap(myGenes,cluster='both',replicates=T)
dev.off()
