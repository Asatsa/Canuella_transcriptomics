## DEG top 1000 genes from the edger results
library(DESeq2)
closeAllConnections()

y = DGEList(counts=data, group=Source)
order(y$tagwise.dispersion, decreasing=TRUE)[1:1000]

countData <- read.csv("~/Thesis/Canuella_transcriptomics/results/results/Trinity.isoform.counts.matrix.tsv", row.names=1,header = TRUE, sep = "\t")
##usingnormalised data
#countData <- read.csv("~/Thesis/Canuella_transcriptomics/results/results/Trinity.isoform.TMM.EXPR.matrix.tsv", row.names=1,header = TRUE, sep = "\t")
metaData <- read.table("~/Thesis/Canuella_transcriptomics/results/results/sample1.txt", header = TRUE, sep = "")
#Construct DESEQDataSet Object

dds <- DESeqDataSetFromMatrix(countData=round(countData), 
                              colData=metaData, 
                              design=~temperature)
dds
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)



resLFC3 <- results(dds, lfcThreshold=2)
table(resLFC3$padj < 0.05)
summary(resLFC3)
##We subset the results table to these genes and then sort it by the log2 fold change 
#estimate to get the significant genes with the strongest down-regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])


#and with the strongest up-regulation:
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

options(max.print = 99999)
#define file name
sink("diffexpressed_dat3.txt")

#define data frame to write to file
idx <- subset(resLFC3,resLFC3$padj < 0.05)##earlier i made a mistake and put 0.05 for the final results
as.data.frame(rownames(idx))

#close the external connection
sink()

png(file="MAplot_deseq_2706.png")
plotMA(res, ylim=c(-2,2))
dev.off()





par(mfrow=c(2,3))
png(file="TRINITY_DN3190_c0_g1_i1.png")
plotCounts(dds, gene="TRINITY_DN3190_c0_g1_i1", intgroup="temperature")
dev.off()
png(file="TRINITY_DN8535_c0_g1_i5.png")
plotCounts(dds, gene="TRINITY_DN8535_c0_g1_i5", intgroup="temperature")
dev.off()
png(file="TRINITY_DN247_c0_g1_i4.png")
plotCounts(dds, gene="TRINITY_DN247_c0_g1_i4", intgroup="temperature")
dev.off()
png(file="TRINITY_DN28064_c2_g1_i1.png")
plotCounts(dds, gene="TRINITY_DN28064_c2_g1_i1", intgroup="temperature")
dev.off()
png(file="TRINITY_DN230014_c0_g1_i2.png")
plotCounts(dds, gene="TRINITY_DN230014_c0_g1_i2", intgroup="temperature")
dev.off()
png(file="TRINITY_DN5097_c2_g1_i4.png")
plotCounts(dds, gene="TRINITY_DN5097_c2_g1_i4", intgroup="temperature")
dev.off()
## transcripp id of interest
DN3190_c0_g1_i1
DN8535_c0_g1_i5
DN247_c0_g1_i4
DN28064_c2_g1_i1
DN230014_c0_g1_i2
DN5097_c2_g1_i4
DN5510_c0_g2_i3
DN71883_c0_g1_i1
DN5107_c0_g1_i10
DN239852_c0_g1_i1
DN18599_c0_g1_i5
DN294403_c0_g1_i1
DN53887_c0_g1_i1
DN6271_c0_g1_i7
DN210730_c0_g1_i1
DN818_c0_g1_i6
DN161529_c0_g1_i1
DN5460_c4_g3_i2
DN295646_c0_g1_i1
DN39_c0_g1_i41
DN39_c0_g1_i73
DN43878_c0_g1_i1
DN4757_c0_g1_i2
DN14867_c0_g1_i4
DN104329_c0_g1_i4
DN4321_c0_g1_i5
DN247866_c0_g2_i1
DN103288_c0_g1_i12
DN310_c0_g1_i7
DN5688_c0_g1_i3
DN97530_c0_g1_i1
DN13495_c0_g1_i1
DN2076_c0_g2_i1
DN2176_c0_g2_i14
DN4740_c0_g1_i3
DN114064_c0_g1_i1
DN28264_c0_g2_i1
### id of LC-PUFA transcrippts
par(mfrow=c(2,3))
png(file="TRINITY_DN37924_c0_g1_i1.png")
plotCounts(idx, gene="TRINITY_DN37924_c0_g1_i1", intgroup="temperature")
dev.off()

plotCounts(res05, gene="TRINITY_DN10382_c0_g1_i2", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN14119_c0_g2_i3", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN1389_c0_g1_i1", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN6216_c0_g3_i1", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN5426_c0_g1_i1", intgroup="temperature")
par(mfrow=c(2,3))
plotCounts(dds, gene="TRINITY_DN120834_c0_g1_i11", intgroup="temperature")###this gene is not there
plotCounts(dds, gene="TRINITY_DN58850_c0_g1_i1", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN9614_c0_g2_i10", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN29323_c0_g1_i1", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN8101_c0_g1_i3", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN2851_c0_g1_i1", intgroup="temperature")
plotCounts(dds, gene="TRINITY_DN2010_c0_g1_i7", intgroup="temperature")

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
png(file="volcanoplot_DESEQ.png")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
png(file="PCA_dESEQ2706.png")
plotPCA(vsdata, intgroup="temperature") #using the DESEQ2 plotPCA fxn we can

dev.off()
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="temperature", 
                returnData=TRUE)
library("ggplot2")
png(file="plot_last.png")
ggplot(d, aes(x=temperature, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
dev.off()



vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")

##plottingthe heatmap the code that actually worked
library('org.Hs.eg.db')
library(ComplexHeatmap)
idx.df <- as.data.frame(idx)
idx.df <- idx.df[(idx.df$baseMean>100)&abs(idx.df$log2FoldChange>2.0),]
#idx.df$symbol <- mapIds(org.Hs.eg.db,keys = row.names(idx.df), keytype = "REFSEQ", column = "SYMBOL")
mat <- counts(dds,normalized=T)[row.names(idx.df),]
mat1 <- t(apply(mat, 1 ,scale))
colnames(mat1) <- (metaData$id)
h <- Heatmap(mat1,cluster_rows = T, cluster_columns = T,column_labels = colnames(mat1),name = 'Z score')#,row_labels = [idx.df(rownames(mat1),]&symbol)
png('heatmapdeseqlog2fc.png',res=250,width = 2000,height = 4000)
print(h)
dev.off()
resultsNames(dds)

#getting upregulated and downregulated genes from edgeR results
library(edgeR)
edgeRfor <- read.table("~/Thesis/Canuella_transcriptomics/results/edgeR.622143.dir/Trinity.isoform.counts.matrix.T20_vs_T23.edgeR.DE_results", row.names=1,header = TRUE, sep = "\t")
fc_res <- edgeRfor[, -c(1:2)] 

out <- topTags(fc_res, n = Inf, p = 0.05)$fc_res
fit <- glmFit(fc_res, metaData)
DE.up <- fc_res[fc_res$logFC>2, ]
DE.up1 <- DE.up[DE.up$FDR<0.05,]
DE.down <- fc_res[fc_res$logFC<2, ]
DE.down1 <- DE.up[DE.down$FDR<0.05,]

