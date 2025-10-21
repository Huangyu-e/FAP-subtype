#Fig.7A
library(clusterProfiler)
ann_colors = list(Type = c(C1 = "#E59F01", C2 = "#0073B3",C3 = "#CC79A7",C4 = "#019E73"))
DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
DEG_1 <- DEG[which(DEG[,6]<0.01&abs(DEG[,2])>1),]
up <- DEG_1 [which(DEG_1 [,2]>0),]
eg <- bitr(rownames(up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway1 <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
								
DEG <- read.table("C2DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
DEG_1 <- DEG[which(DEG[,6]<0.01&abs(DEG[,2])>1),]
up <- DEG_1 [which(DEG_1 [,2]>0),]
eg <- bitr(rownames(up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway2 <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
																
DEG <- read.table("C3DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
DEG_1 <- DEG[which(DEG[,6]<0.01&abs(DEG[,2])>1),]
up <- DEG_1 [which(DEG_1 [,2]>0),]
eg <- bitr(rownames(up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway3 <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
  
DEG <- read.table("C4DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
DEG_1 <- DEG[which(DEG[,6]<0.01&abs(DEG[,2])>1),]
up <- DEG_1 [which(DEG_1 [,2]>0),]
eg <- bitr(rownames(up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway4 <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)

union_pathway <- unique(c(KEGG_pathway1[1:5,4],KEGG_pathway2[1:5,4],KEGG_pathway3[1:5,4],KEGG_pathway4[1:5,4]))
C1path <- KEGG_pathway1[which(!is.na(match(KEGG_pathway1[,4],union_pathway))),]
C1path$type <- "C1"
C2path <- KEGG_pathway2[which(!is.na(match(KEGG_pathway2[,4],union_pathway))),]
C2path$type <- "C2"
C3path <- KEGG_pathway3[which(!is.na(match(KEGG_pathway3[,4],union_pathway))),]
C3path$type <- "C3"
C4path <- KEGG_pathway4[which(!is.na(match(KEGG_pathway4[,4],union_pathway))),]
C4path$type <- "C4"
pathways <- rbind(C1path,C2path,C3path,C4path)
pathways_1 <- pathways[,c(4,10,15)]
pathways_1$log10P <- (-log10(pathways_1$pvalue))
pathways_1$Description <- factor(pathways_1$Description,levels=rev(unique(c(KEGG_pathway1[1:5,4],KEGG_pathway2[1:5,4],KEGG_pathway3[1:5,4],KEGG_pathway4[1:5,4]))))

library(ggplot2)
ggplot(pathways_1, aes(x=type,y=Description)) +
  geom_point(aes(size=log10P),color="#F24F4F") +
  scale_size_continuous(range = c(1,7)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "lines"))

##Fig.7B
library(clusterProfiler)
DEG <- read.table("GEO/GSE22220/C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
DEG_1 <- DEG[which(abs(DEG[,1])>1&DEG[,5]<0.01),]
up <- DEG_1 [which(DEG_1 [,2]>0),]
eg <- bitr(rownames(up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway1 <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
KEGG_pathway1_1 <- KEGG_pathway1[1:10,]
KEGG_pathway1_1$log10P <- (-log10(KEGG_pathway1_1$pvalue))
KEGG_pathway1_1$Description <- factor(KEGG_pathway1_1$Description,levels=rev(KEGG_pathway1_1$Description))							
ggplot(KEGG_pathway1_1 , aes(x = log10P, y = Description)) +
  geom_bar(stat = "identity",fill= "#F24F4F") +
  theme_minimal() + theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black") 
  )+
  labs(title = "C1", x = "-log10P", y = "Description") 
  
##Fig.7C
nonsilent_mutation <- read.table("GDCdata/TCGA-BRCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/nonsilent_mutation.txt",header = F,sep = "\t",stringsAsFactors = F)
sample_trans <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	return(paste(split_name[1:3],collapse = "-"))
}
nonsilent_mutation$sample <- unlist(lapply(nonsilent_mutation[,16],sample_trans))
nonsilent_mutation_1 <- nonsilent_mutation[,c(5,141)]
cluster <- read.table("cluster.txt",header = T,sep = "\t",stringsAsFactors = F)
nonsilent_mutation_2 <- merge(nonsilent_mutation_1,cluster,by="sample")
mutation_matrix <- matrix(rep(0,4*length(unique(nonsilent_mutation_2[,2]))),nrow=4)
rownames(mutation_matrix) <- c("C1","C2","C3","C4")
colnames(mutation_matrix) <- unique(nonsilent_mutation_2[,2])
for(i in rownames(mutation_matrix)){
    each_i <- nonsilent_mutation_2[which(nonsilent_mutation_2[,3]==i),]
	for(j in unique(each_i[,2])){
	    each_chr <- each_i[which(each_i[,2]==j),]
		mutation_matrix[which(rownames(mutation_matrix)==i),which(colnames(mutation_matrix)==j)] <- nrow(each_chr)
	}
}
my.data <- mutation_matrix[,c(3,11,14,19,13,12,10,21,7,6,18,15,16,2,8,17,9,23,1,22,20,4,5,24)]
grid.col <- c("#E59F01","#0073B3","#CC79A7","#019E73",rep("gray",24))
library(circlize) 
circos.par(gap.degree = c(rep(2, nrow(my.data)-1), 10, rep(2, ncol(my.data)-1), 10),
           start.degree = 180)
chordDiagram(my.data,
             directional = TRUE,
             diffHeight = 0.06,
             grid.col = grid.col,
             transparency = 0.5,reduce = 0)
circos.clear()

##Fig.7D
library(ggpubr)
sample_stat <- data.frame(sample=names(table(nonsilent_mutation_2[,1])),count=as.vector(table(nonsilent_mutation_2[,1])))
merge_score <- merge(sample_stat,cluster,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
merge_score$count <- log10(merge_score$count+1)
ggplot(merge_score, aes(y = count, x = clus, fill = clus)) +
  geom_violin() +scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
  labs(x = "cluster",
      y = "log10 (Mutation Load+1)")
	  
##Fig.7E
sample_cnv <- read.table("TCGA-BRCA.masked_cnv_DNAcopy.tsv/TCGA-BRCA.masked_cnv_DNAcopy.tsv",header = T,sep = "\t",stringsAsFactors = F)
type_trans <- function(x){return(unlist(strsplit(x,split="-"))[4])}
sample_type <- unlist(lapply(sample_cnv$sample,type_trans))
sample_cnv_1 <- sample_cnv[which(sample_type=="01A"),]
sample_trans <- function(x){
  split_name <- unlist(strsplit(x,split="-"))
  return(paste(split_name[1:3],collapse = "-"))
}
sample_cnv_1$sample_1 <- unlist(lapply(sample_cnv_1$sample,sample_trans)) 
sample_cnv_2 <- sample_cnv_1[,c(6,2,3,4,5)]
C1 <- sample_cnv_2[which(!is.na(match(sample_cnv_2[,1],cluster[which(cluster[,2]=="C1"),1]))),]
C1 <- C1[,-1]
C2 <- sample_cnv_2[which(!is.na(match(sample_cnv_2[,1],cluster[which(cluster[,2]=="C2"),1]))),]
C2 <- C2[,-1]
C3 <- sample_cnv_2[which(!is.na(match(sample_cnv_2[,1],cluster[which(cluster[,2]=="C3"),1]))),]
C3 <- C3[,-1]
C4 <- sample_cnv_2[which(!is.na(match(sample_cnv_2[,1],cluster[which(cluster[,2]=="C4"),1]))),]
C4 <- C4[,-1]
library(maftools)
gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_95.BRCA.txt", 
                     gisticAmpGenesFile="amp_genes.conf_95.BRCA.txt", 
                     gisticDelGenesFile="del_genes.conf_95.BRCA.txt", 
                     gisticScoresFile="scores.BRCA.gistic", isTCGA=TRUE)
gisticChromPlot(gistic=gistic, markBands="all")
gisticOncoPlot(gistic = gistic, 
               sortByAnnotation = TRUE, top = 10)
library(OmicCircos)
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="")
circos(R=390, cir="hg19", W=4,   type="chr", 
       print.chr.lab=TRUE, scale=TRUE);
circos(R=340, cir="hg19", W=40,  mapping=C1, col="#faeccc",     
       col.v=4,   type="ml3", B=T, lwd=1, cutoff=0);
circos(R=300, cir="hg19", W=40,  mapping=C2,  col="#cce3f0",     
       col.v=5,   type="ml3", B=T, lwd=1, cutoff=0);
circos(R=260, cir="hg19", W=40,  mapping=C3, col="#f5e4ed",      
       col.v=5,   type="ml3", B=T, lwd=1, cutoff=0);
circos(R=220, cir="hg19", W=40,  mapping=C4,  col="#ccece3",     
       col.v=5,   type="ml3", B=T, lwd=1, cutoff=0)

##Fig.7F
library(ggpubr)
sample_stat <- data.frame(sample=names(table(sample_cnv_2[,1])),count=as.vector(table(sample_cnv_2[,1])))
merge_score <- merge(sample_stat,cluster,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
merge_score$count <- log10(merge_score$count+1)
ggplot(merge_score, aes(y = count, x = clus, fill = clus)) +
  geom_violin() +scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
  labs(x = "cluster",
      y = "log10 (CNV Load+1)")

