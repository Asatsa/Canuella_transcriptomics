#MIXOMICS ANALYSIS
library(tidyverse)
library(GGally)
library(corrr)
library(reshape2)
library(mixOmics)

## setting up data
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
#y <- DGEList("~/Thesis/Canuella_transcriptomics/results/counts.matrix",header = TRUE)
seqdata <- read.delim("~/Thesis/Canuella_transcriptomics/results/edgeR.622143.dir/Trinity.isoform.counts.matrix.T20_vs_T23.edgeR.count_matrix.txt", row.names=1, stringsAsFactors = FALSE)
logcounts <- cpm(seqdata,log=TRUE)
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
select_var1 <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)
highly_variable_lcpm1 <- logcounts[select_var1,]
dim(highly_variable_lcpm1)
write.table(highly_variable_lcpm1,"top1001.txt",sep = "\t")

gene <- read.csv("~/Thesis/Canuella_transcriptomics/results/edgeR.622143.dir/Trinity.isoform.counts.matrix.T20_vs_T23.edgeR.count_matrix.txt", row.names = 1,header = TRUE,sep = "\t")
#gene <- gene %>% relocate(6:10)
gene <- highly_variable_lcpm##top 500
gene <- highly_variable_lcpm1##top100
gene <- t(gene)
gene <- as.matrix(gene)

#gene <- read.csv("~/Thesis/Canuella_transcriptomics/results/edgeR.622143.dir/diffExpr.P1e-3_C2.matrix.log2.centered.dat.csv", row.names = 1,header = TRUE,sep = "\t")

FA <- read.csv("~/Documents/Thesis/masshunter/FA_RESULTS/FA_relative.csv", row.names = 1, header = TRUE, sep = ",")
FA <- FA[c(11:20),]
rownames(FA) <- rownames(gene)
factor <- as.data.frame(FA[,c(30)])
rownames(factor) <- rownames(gene)
colnames(factor)[1] <- "temperature"
factor$temperature <- as.factor(factor$temperature)
str(factor)
FA <- FA[,-c(11, 13, 16, 20, 23, 25, 30)] # remove FA that are 0, 19:0 and temperature factor
FA <- as.matrix(FA)

dim(gene) # check the dimensions of the X dataframe
dim(FA) # check the dimensions of the Y dataframe
png(file="samplecorrelation100C.png")
im1 <- imgCor(gene, FA, sideColors = c("purple", "green"))
dev.off()

grid1 <- seq(0.001, 0.4, length = 30) 
grid2 <- seq(0.001, 0.4, length = 30)

tune.rcc <- tune.rcc(gene, FA, grid1 = grid1, grid2 = grid2, validation = "loo")
tune.rcc
opt.l1 <- tune.rcc$opt.lambda1 # extract the optimal lambda values
opt.l2 <- tune.rcc$opt.lambda2

rcc.cv <- rcc(FA,gene, method = "ridge", lambda1 = opt.l1, lambda2 = opt.l2)

rcc.shrink <- rcc(FA,gene, method = "shrinkage")
rcc.shrink$lambda

png(file="crossvalidation100C.png")
plot(rcc.cv, type = "barplot", main = "Cross Validation")
dev.off()

png(file="shrinkage100C.png")
plot(rcc.shrink, type = "barplot", main = "Cross Validation")
dev.off()

png(file="Plotofindividuals100C.png")
plotIndiv(rcc.cv, comp = 1:2, ind.names = rownames(factor), group = factor$temperature, rep.space = "XY-variate", 
          legend = TRUE, title = "Plot of Individuals")
dev.off()

plotIndiv(rcc.shrink, comp = 1:2, ind.names = rownames(factor), group = factor$temperature, rep.space = "XY-variate", 
          legend = TRUE, title = "Plot of Individuals")
png(file="Arrowplot100C.png")
plotArrow(rcc.cv, group = factor$temperature, col.per.group = color.mixo(1:2), title = 'Arrow Plot')
dev.off()
png(file="varplot100C1.png",res=250,width = 3000,height = 3000)
plotVar(rcc.cv, var.names = c(TRUE, TRUE), cex = c(4, 4), cutoff = 0.5, title = 'Plot of Variables(CV), cutoff = 0.5')
dev.off()
png(file="varplot200C1.png",res=250,width = 3000,height = 3000)
plotVar(rcc.shrink, var.names = c(TRUE, TRUE), cex = c(4, 4), cutoff = 0.5, title = 'Plot of Variables(shrink), cutoff = 0.5')
dev.off()
png(file="network100C1.png",res=250,width = 2000,height = 3000)
rcc_network <- network(rcc.cv, comp = 1:2, interactive = F, lwd.edge = 3, cutoff = 0.634)
dev.off()
png(file="network2001.png")
rcc_network_full <- network(rcc.cv, comp = 1:2, interactive = FALSE, lwd.edge = 3)
dev.off()

cancor <- data.frame(rcc_network_full$M)
cancor <- data.frame(t(cancor))
write.csv(cancor, file = "./output/Asatsa_mixomics_top100C_matrix1.csv")
png(file="corelationmatrixlast100C1.png",res=250,width = 2000,height = 2000)
cim(rcc.cv, comp = 1:2, xlab = "genes", ylab = "fatty acids")
dev.off()