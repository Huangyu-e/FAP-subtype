##Fig.4A
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
sample=as.vector(colnames(mRNA_expression))
express=as.vector(mRNA_expression[which(rownames(mRNA_expression)=="CD274"),])
PDL1_exp <- data.frame(cbind(sample,express))
PDL1_exp$sample_1 <- gsub("\\.", "-", PDL1_exp$sample)
merge_score <- merge(cluster,PDL1_exp,by.x="sample",by.y="sample_1")
merge_score$express <- as.numeric(merge_score$express)
merge_score$express <- log10(merge_score$express+1)
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
library(ggpubr)
compare_means(express  ~ clus,data = merge_score)
	  
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
ggbetweenstats(
  data = merge_score,
  x = clus,
  y = express,centrality.type="nonparametric"
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()

##Fig.4B
cbio <- read.table("TCGA_data/combined_study_clinical_data.tsv",header = T,sep = "\t",stringsAsFactors = F)
cbio_1 <- cbio[,c(2,29,30,31)]
cbio_1$TMB <- cbio_1$Mutation.Count/38
merge_score <- merge(cluster,cbio_1,by.x="sample",by.y="Patient.ID")
merge_score$TMB <- log10(merge_score$TMB+1)
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
compare_means(TMB  ~ clus,data = merge_score)
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
ggbetweenstats(
  data = merge_score,
  x = clus,
  y = TMB,centrality.type="nonparametric"
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()

##Fig.4C
merge_score$MSI.MANTIS.Score1 <- log10(merge_score$MSI.MANTIS.Score+1)
compare_means(MSI.MANTIS.Score1 ~ clus,data = merge_score)
ggbetweenstats(
  data = merge_score,
  x = clus,
  y = MSI.MANTIS.Score1,centrality.type="nonparametric"
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()

##Fig.4D
Thorsson <- read.table("TCGA_data/Thorsson.txt",header = T,sep = "\t",stringsAsFactors = F)
BRCA_Thorsson <- Thorsson[which(Thorsson[,2]=="BRCA"),]
BRCA_Thorsson_1 <- merge(cluster,BRCA_Thorsson,by.x="sample",by.y="TCGA.Participant.Barcode")
BRCA_Thorsson_1$clus <- factor(BRCA_Thorsson_1$clus ,levels=c("C1","C2","C3","C4"))
BRCA_Thorsson_1$SNV.Neoantigens <- log10(BRCA_Thorsson_1$SNV.Neoantigens+1)
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =SNV.Neoantigens,centrality.type="nonparametric"
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(SNV.Neoantigens ~ clus,data = BRCA_Thorsson_1)

##Fig.4E
BRCA_Thorsson_1$Indel.Neoantigens <- log10(BRCA_Thorsson_1$Indel.Neoantigens+1)
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =Indel.Neoantigens,centrality.type="nonparametric",
 )+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(Indel.Neoantigens ~ clus,data = BRCA_Thorsson_1)

##Fig.4F
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =Leukocyte.Fraction,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(Leukocyte.Fraction ~ clus,data = BRCA_Thorsson_1)

##Fig.4G
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =Lymphocytes,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(Lymphocytes ~ clus,data = BRCA_Thorsson_1)

##Fig.4H
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =TIL.Regional.Fraction,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(TIL.Regional.Fraction ~ clus,data = BRCA_Thorsson_1)

##Fig.4I
BRCA_Thorsson_1$TCR <- log10(BRCA_Thorsson_1$TCR.Richness+1)
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =TCR,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)  
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(TCR ~ clus,data = BRCA_Thorsson_1)

##Fig.4J
BRCA_Thorsson_1$BCR <- log10(BRCA_Thorsson_1$BCR.Richness+1)
ggbetweenstats(
  data = BRCA_Thorsson_1,
  x = clus,
  y =BCR,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(BCR ~ clus,data = BRCA_Thorsson_1)

##Fig.4K
library(IOBR)
estimate <- deconvo_tme(eset = mRNA_expression, method = "estimate")
estimate <- as.data.frame(estimate)
estimate$sample <- gsub("\\.", "-", estimate[,1])
merge_score <- merge(cluster,estimate,by="sample")
merge_score$clus <- factor(merge_score$clus ,levels=c("C1","C2","C3","C4"))
ggbetweenstats(
  data = merge_score,
  x = clus,
  y =ImmuneScore_estimate,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(ImmuneScore_estimate ~ clus,data = merge_score)

##Fig.4L
library(GSVA)
xcell <- deconvo_tme(eset = mRNA_expression, method = "xcell", arrays = TRUE)
xcell <- as.data.frame(xcell)
xcell$sample <- gsub("\\.", "-", xcell[,1])
merge_score <- merge(cluster,xcell,by="sample")
merge_score$clus <- factor(merge_score$clus ,levels=c("C1","C2","C3","C4"))
library(ggpubr)
ggbetweenstats(
  data = merge_score,
  x = clus,
  y =ImmuneScore_xCell,centrality.type="nonparametric",
 point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 2, stroke = 0)  
)+ scale_color_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
compare_means(ImmuneScore_xCell ~ clus,data = merge_score)

#Fig.4M
mcpcounter <- deconvo_tme(eset = mRNA_expression, method = "mcpcounter")
mcpcounter <- as.data.frame(mcpcounter)
mcpcounter$sample <- gsub("\\.", "-", mcpcounter[,1])
merge_score <- merge(cluster,mcpcounter,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
library(ggpubr)
final_map <- c()
for(i in 4:11){
    each_map <- merge_score[,c(i,2)]
	each_map$cell_type <- colnames(merge_score)[i]
	colnames(each_map)[1] <- "score"
	final_map <- rbind(final_map,each_map)
}
final_map_1 <- final_map[which(!is.na(final_map[,1])),]
final_map <- final_map_1
final_map <- as.data.frame(final_map)
final_map[,1] <- as.numeric(final_map[,1])
final_map$clus <- factor(final_map$clus,levels=c("C1","C2","C3","C4"))
final_map[,1] <- log10(final_map[,1]+1)
ggplot(final_map, aes(x=cell_type, y=score, fill=clus)) + 
    geom_boxplot()+theme_classic()+scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

library(ggpubr)
stat <- c()
for(i in 4:13){
    each_hallmark <- merge_score[,c(2,i)]
	colnames(each_hallmark)[2] <- "score"
	each_P <- compare_means(score~clus,data=each_hallmark,method="kruskal.test")$p.adj
	each_stat <- c(colnames(merge_score)[i],each_P)
	stat <- rbind(stat,each_stat)
}
stat <- as.data.frame(stat)
stat[,2] <- as.numeric(stat[,2])
stat[which(stat[,2]>=0.05),]

#Fig.4N
Heetal <- read.table("pathway/Heetal.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
genesetList <- list()
for(i in 1:nrow(Heetal)){
    each_signature <- Heetal[i,]
	genesetList[[Heetal[i,1]]]= as.vector(as.character(unlist(strsplit(each_signature[,2],split=", "))))
	
}
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
mRNA_expression <- as.matrix(mRNA_expression)
library(GSVA)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", rownames(GSVA_score_1))
merge_score <- merge(cluster,GSVA_score_1,by="sample")
library(ggpubr)
stat <- c()
for(i in 3:ncol(merge_score)){
    each_hallmark <- merge_score[,c(2,i)]
	colnames(each_hallmark)[2] <- "score"
	each_P <- compare_means(score~clus,data=each_hallmark,method="kruskal.test")$p
	each_stat <- c(colnames(merge_score)[i],each_P)
	stat <- rbind(stat,each_stat)
}
stat <- as.data.frame(stat)
stat[,2] <- as.numeric(stat[,2])
stat[which(stat[,2]>=0.05),]

final_map <- c()
for(i in 3:ncol(merge_score)){
    each_map <- merge_score[,c(i,2)]
	each_map$cell_type <- colnames(merge_score)[i]
	colnames(each_map)[1] <- "score"
	final_map <- rbind(final_map,each_map)
}
final_map_1 <- final_map[which(!is.na(final_map[,1])),]
final_map <- final_map_1
final_map <- as.data.frame(final_map)
final_map[,1] <- as.numeric(final_map[,1])
final_map$clus <- factor(final_map$clus,levels=c("C1","C2","C3","C4"))
ggplot(final_map, aes(x=cell_type, y=score, fill=clus)) + 
    geom_boxplot()+theme_classic()+scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

