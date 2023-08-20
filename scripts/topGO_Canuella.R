##        Gene ontology enrichment analysis of RNA-seq data          ## 
## from the ACCLIMATE experiment using Bioconductor, edgeR and topGO ##
## ----------------------------------------------------------------- ##


## Date GO analysis:    June 2023

## Assembly server:     we11sv01.ugent.be -> http://www.kennybogaert.eu/we11sv01-usertips/


## R script source:     https://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

## ---------------------------------------------------------------- ##

## 1. Setting up
## -------------

BiocManager::install("Rgraphviz")
library(topGO)
library(Rgraphviz)

## 2. Input data
## -------------

## You need: annotations.txt interestinggenes.txt output_file


#geneID2GO = readMappings(file = "~/Downloads/annotations.txt")
geneID2GO = readMappings(file = "~/Transcriptomics/go_annotations1.txt")
#genesOfInterest = read.csv("~/Transcriptomics/allgenes.txt", header = F)

genesOfInterest = read.csv("~/Transcriptomics/diffexpressedgenetrue.txt", header = F)#edgeRdata
#genesOfInterest = read.csv("~/Transcriptomics/diffexpressed_dat.txt", header = F)#deseq data
#genesOfInterest = read.csv("~/Transcriptomics/diffexpressed_dat1.txt", header = F)#deseq data with log2 fold set at 2 not TMM data
genesOfInterest = read.csv("~/Transcriptomics/diffexpressed_dat2.txt", header = F)# all parameters, TMM data,log2fold at 1
# set the output file
sink(file = "./output_topgolo2foldat1.txt")

# read in the 'gene universe' file
geneUniverse <- names(geneID2GO)
#geneUniverse <- gsub("\\\"", "", geneUniverse)

# read in the genes of interest
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse%in% genesOfInterest))
names(geneList) <- geneUniverse


#myInterestingGenes <- sample(geneID2GO, length(geneID2GO) / 10)
#geneList <- factor(as.integer(geneID2GO %in% myInterestingGenes))
#names(geneList) <- geneID2GO
#str(geneList)

# build the GOdata object in topGO
myGOdata <- new("topGOdata", description="ACCLIMATE", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
#saveRDS(myGOdata, file = "GOdata_deseq_withadjustedlcfold.rds")#non tTMM log2fold data
saveRDS(myGOdata, file = "GOdata_deseq_withadjustedlcfold.rds")# TMM data at 1
#myGOdata <- readRDS("GOdata_deseq.rds")
# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.01:
mysummary <- summary(attributes(resultTopgo)$score <= 0.01)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.01

# print out the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif,numChar=1000)
allRes
write.table(allRes, file = "./output/topGOdeseqTMM.txt")
# print a graph (to a pdf file) with the top 'numsignif' results:
#printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = "./output/output_TopGOdeseq_", useInfo = "all", pdfSW = TRUE)
#dev.off()

# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink() 

###generating plots

require(ggplot2)
library(scales)
goEnrichment <- allRes
goEnrichment$topgoFisher <- as.numeric(allRes$topgoFisher)
goEnrichment <- goEnrichment[allRes$topgoFisher < 0.05,] # filter terms for KS p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","topgoFisher")]
goEnrichment
ntop <- 30
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(unique(ggdata$Term))) # fixes orderrev(unique(dataset$variable))
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(topgoFisher), size = -log10(topgoFisher), fill = -log10(topgoFisher))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO processes(topgo data)',
    subtitle = 'Top 30 terms ordered by topgoFisher p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()
png("goseqdeseqTMMfold130.png",res=250,width = 5500,height = 4000)
print(gg1)
dev.off()