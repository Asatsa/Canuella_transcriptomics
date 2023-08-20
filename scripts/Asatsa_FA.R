# Fatty acid analysis plots

library(tidyverse)
library(GGally)
library(corrr)
library(reshape2)
library(mixOmics)
library(ggplot2)


#gene <- read.csv("./data/diffExpr.P1e-3_C2.matrix.log2.centered.dat.csv", row.names = 1,header = TRUE,sep = "\t")
gene <- read.csv("~/Thesis/Canuella_transcriptomics/results/edgeR.622143.dir/diffExpr.P1e-3_C2.matrix.log2.centered.dat.txt", row.names = 1,header = TRUE,sep = "\t")
gene <- gene %>% relocate(6:10)
gene <- t(gene)
gene <- as.matrix(gene)
data1 <- read.table("~/Documents/Thesis/masshunter/FA_RESULTS/FA_raw1.csv",row.names = 1, header = T, sep = ",")
FA <- read.csv("~/Documents/Thesis/masshunter/FA_RESULTS/FA_relative.csv", row.names = 1, header = TRUE, sep = ",")
FA <- FA[c(11:20),]#lab copepods
FA1 <- FA[c(6:10),]#algae
FA <- FA[,-c(20)] # remove the lab standard 19:0
FA2 <- FA[c(1:5,11:20),]#all copepods
FA3 <- FA[c(1:5),]#field copepods
rownames(FA) <- rownames(gene)
FA <- FA[,-c(11, 13, 16, 20, 23, 25)]# remove FA that are 0, 19:0 and temperature factor

FA3 <- FA3[,-c(20)]# remove FA that are 0, 19:0 and temperature factor

FA_long <- gather(data = FA, key = FA, value = percentage, 'X12.0':'X22.6.n.3', factor_key = TRUE)
FA_LONGALGAE <- gather(data = FA1, key = FA1, value = percentage, 'X12.0':'X22.6.n.3', factor_key = TRUE)
FA_copepodsf <- gather(data = FA3, key = FA3, value = percentage, 'X12.0':'X22.6.n.3', factor_key = TRUE)#filed copepods
FA_copepods <- gather(data = FA2, key = FA2, value = percentage, 'X12.0':'X22.6.n.3', factor_key = TRUE)#allcopepods
write.table(FA_long,"RIGHTDATAFORSTATS.txt",sep = '\t')
FAlist <- unique(FA_long$FA)

FA_subset <- list()

for (i in 1:length(FAlist)){
  FA_subset[[i]] <- subset(FA_long, FA==FAlist[i])
}

for (i in 1:length(FAlist)){
    if(i == "17") next       #you add this bit of code when you want to exclude certain FAs from the analysis (e.g. C17:)               
  ttest <- t.test(percentage ~ Temperature, data = FA_subset[[i]])
  print(paste(FAlist[i]))
  print(ttest)
}

for (i in 1:length(FAlist)){
    if(i == 17) next      # you add this bit of code when you want to exclude certain FAs from the analysis (e.g. C17:)               
  wilcoxon <- wilcox.test(percentage ~ Temperature, data = FA_subset[[i]])
  print(paste(FAlist[i]))
  print(wilcoxon)
}


FAplot <- function(dat, name){
  
  p = ggplot(dat, aes(x= Temperature, y= percentage, fill = Temperature)) + 
    geom_boxplot()+
    ggtitle(paste0(name))+
    theme_bw()

 ggsave(filename = paste0("./output/", name, ".png"), width=4, height=3, unit='in', dpi=300)
}  

for (i in 1:length(FAlist)){
  print(i)
  print(FAplot(dat = FA_subset[[i]], name = FAlist[i]))
}

ggplot(FA_LONGALGAE, aes(x = FA, y=percentage, fill = factor(Temperature))) + 
  geom_col() +
  geom_errorbar(aes(ymin = percentage - se, ymax = percentage + se), width =.2) + 
  geom_label(aes(label = paste(percentage, "\ub1", se)), nudge_y = 5, size = 2,
             label.size = 0, label.r = unit(0, "pt"), fill = "#ebebeb") +
  geom_point(size = 1) +
  theme(legend.position = "right", axis.text.x = element_text(face = "bold", size = 10, colour = "black"),
  #axis.text.x = element_text(colour = "black", face = "bold", size = 7), 
  legend.text = element_text(size = 7, face ="bold", colour ="black"),
  axis.title.y = element_text(face = "bold", size = 10),
  legend.title = element_text(size = 14, colour = "black", face = "bold"),
  legend.key=element_blank(),
  labs(x = "FA", y = "percentage"),
  scale_colour_manual(values = c("#009E73", "#E69F00")))



p1<- ggplot(FA_long, aes(x = FA, y=percentage, fill = factor(Temperature))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin= percentage - se, ymax = percentage + se), width=.2,
                position=position_dodge(.9))+
theme(axis.text.x = element_text(face = "bold", size = 5, colour = "black"))
print(p1)
# Finished bar plot
p1+labs(title="Relative abundance per FA", x="FA", y = "% concentration")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))

se <- sd(FA_long$percentage)
#t.test(FA_long$percentage ~ FA_long$Temperature, data = FA_long$FA)


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
##plot for field copepods

df4 <-  data_summary(FA_copepodsf, varname="percentage", 
                     groupnames=c("FA3", "Temperature"))
png("relative abundancefieldcopepodsfixedaxis1.png",res=250,width = 2000,height = 2000)
z<- ggplot(df4, aes(x = FA3, y=percentage, fill = factor(FA3))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin= percentage - sd, ymax = percentage + sd), width=.01,
                position=position_dodge(.9)) +
  labs(title="Relative concentration per FA for field samples", x="FA", y ="% concentration")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.text.x = element_text(angle = 45, size = 10, face = 'bold', hjust = 1.10),
        legend.position="none")

z
dev.off()
##plot for lab copepods only
df1 <- data_summary(FA_long, varname="percentage", 
                    groupnames=c("FA", "Temperature"))
png("relative abundancelabcopepodsfixedaxis1.png", res=250,width = 2000,height = 2000)
z<- ggplot(df1, aes(x = FA, y=percentage, fill = factor(Temperature))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin= percentage - sd, ymax = percentage + sd), width=.2,
                position=position_dodge(.9)) +
  labs(title="Relative abundance per FA for lab samples", x="FA", y ="% concentration")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.text.x = element_text(angle = 45, size = 10, face = 'bold', hjust = 1.10),
        legend.position="none")
 
z
dev.off()
# Finished bar plot

z+labs(title="Relative abundance per FA", x="FA", y ="% concentration")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))
dev.off()

df3 <- data_summary(FA_LONGALGAE, varname="percentage", 
                    groupnames=c("FA1", "Temperature"))
png("relative abundancealgaefixedaxis.png")
AL<- ggplot(df3, aes(x = FA1, y=percentage,fill = factor(FA1), show.legend = F))+
  geom_bar(stat="identity",color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin= percentage - sd, ymax = percentage + sd), width=.3,
                position=position_dodge(.9)) +
  labs(title="Relative abundance per FA in algae samples", x="FA1", y ="% concentration")+
  theme(axis.text.x = element_text(angle = 45, size = 12, face = 'bold', hjust = 1.10),
        legend.position="none")
  
  
 
AL
dev.off()


write.csv(df2,"table with sd.csv")
sd <- sd(FA_copepodsf$percentage)
sd <- sd(FA_long$percentage)
sd <- sd(FA_LONGALGAE$percentage)

sd(FA$X22.6.n.3)
dha <- (FA$X22.6.n.3)
dha1 <-  c(0.4218547,0.4072803,0.3877807,0.4400876,0.4565700)
dha2 <- c( 0.2375991, 0.3176449, 0.2579014, 0.4247258, 0.3554029)
sd(dha2)
#t.test between field samples and labsamples
t.test(FA_long$percentage, FA_copepodsf$percentage)


theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
  
  show.legend = F
