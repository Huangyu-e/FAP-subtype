##Fig.3A
pathwaytype <- read.table("Ferroptosispathway_type.txt",header = T,sep = "\t",stringsAsFactors = F)
pathway <- read.table("Ferroptosispathway.txt",header = T,sep = "\t",stringsAsFactors = F)
Promotion_pathway <- list()
Suppression_pathway <- list()
for(i in unique(pathwaytype[,1])){
    each_pathway <- pathway[which(pathway[,2]==i),]
	if(pathwaytype[which(pathwaytype[,1]==i),2]=="Promotion"){
	    Promotion_pathway[[i]] <- each_pathway[,1]
	}
	else if(pathwaytype[which(pathwaytype[,1]==i),2]=="Suppression"){
	    Suppression_pathway[[i]] <- each_pathway[,1]
	}
}
Promotion_pathway <- list(pathway[which(!is.na(match(pathway[,2],pathwaytype[which(pathwaytype[,2]=="Promotion"),1]))),1])	 
##Ferrdb
Ferroptosis_driver <- read.table("ferroptosis_driver.txt",header = T,sep = "\t",stringsAsFactors = F,quote="") 
Ferroptosis_driver_1 <- Ferroptosis_driver[which(Ferroptosis_driver$testin=="Human"|Ferroptosis_driver$testin=="Human, mice"|Ferroptosis_driver$testin=="Human, rat"|Ferroptosis_driver$testin=="Human, porcine"),]
##KEGG
library(clusterProfiler)
library(org.Hs.eg.db)  
library(KEGGREST)
kegg_code <- "hsa04216"
pathway_info <- keggGet(kegg_code)
gene_entries <- pathway_info[[1]]$GENE
entrez_ids <- gene_entries[seq(1, length(gene_entries), 2)]
gene_symbols <- bitr(entrez_ids, fromType = "ENTREZID",
                     toType = "SYMBOL", OrgDb = org.Hs.eg.db)

driver_1 <- list(Reduce(intersect,list(driver=Ferroptosis_driver_1[,3],Promotion=Promotion_pathway[[1]],IOBR_ferro=gene_symbols[,2])))

library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
ssgsea <- ssgseaParam(mRNA_expression,driver_1)
signature_score <- gsva(ssgsea)
GSVA_score_1 <- as.data.frame(t(signature_score))
GSVA_score_1$ID <- rownames(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", GSVA_score_1$ID)
merge_score <- merge(cluster,GSVA_score_1,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))

ggplot(merge_score, aes(y = V1, x = clus, fill = clus)) +
  geom_violin() +scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
  labs(x = "cluster",
      y = "Ferroptosis")#
compare_means(V1 ~ clus ,data = merge_score)	

##Fig.3B
Thorsson <- read.table("TCGA_data/Thorsson.txt",header = T,sep = "\t",stringsAsFactors = F)
BRCA_Thorsson <- Thorsson[which(Thorsson[,2]=="BRCA"),]
BRCA_Thorsson_1 <- merge(cluster,BRCA_Thorsson,by.x="sample",by.y="TCGA.Participant.Barcode")
BRCA_Thorsson_1$clus <- factor(BRCA_Thorsson_1$clus,levels=c("C1","C2","C3","C4"))
library(ggpubr)

ggplot(BRCA_Thorsson_1, aes(y = Intratumor.Heterogeneity, x = clus, fill = clus)) +
  geom_violin() +scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+theme_classic()
  labs(x = "cluster",
      y = "Intratumor.Heterogeneity")  
compare_means(Intratumor.Heterogeneity ~ clus,data = BRCA_Thorsson_1)

##Fig.3C
library(ggridges)
library(clusterProfiler)
hsa_kegg <- clusterProfiler::download_KEGG("hsa")
pathway_gene <- hsa_kegg [[1]]
stem_pathway <- pathway_gene[which(pathway_gene[,1]=="hsa04550"),]
eg <- bitr(stem_pathway[,2],fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")
stem_sig <- list(stemn=eg[,2])
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
mRNA_expression <- as.matrix(mRNA_expression)
library(GSVA)
ssgsea <- ssgseaParam(mRNA_expression,stem_sig)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", rownames(GSVA_score_1))

merge_score <- merge(cluster,GSVA_score_1,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
library(ggpubr)
library(ggridges)
ggplot(merge_score, aes(x = stemn  , y = clus, fill = clus)) +
  geom_density_ridges(alpha = 0.7) +
  theme_ridges() +scale_fill_manual(values=c("#E59F01","#0073B3","#CC79A7","#019E73"))+
  labs(title = "Ridgeline Plot of Sepal Lengths by Species",
       x = "stemn",
      y = "cluster")
compare_means(stemn ~ clus,data = merge_score)##两两之间全都显著

##Fig.3D
library(IOBR)
cluster <- read.table("cluster.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
EMT <- list(WNT_target=signature_collection$WNT_target,Pan_F_TBRs=signature_collection$Pan_F_TBRs,
EMT1=signature_collection$EMT1,EMT2=signature_collection$EMT2,EMT3=signature_collection$EMT3)
mRNA_expression <- as.matrix(mRNA_expression)
library(GSVA)
ssgsea <- ssgseaParam(mRNA_expression,EMT)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", rownames(GSVA_score_1))
merge_score <- merge(cluster,GSVA_score_1,by="sample")
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))

library(ggpubr)
stat <- c()
for(i in 3:ncol(merge_score)){
    each_hallmark <- merge_score[,c(2,i)]
	colnames(each_hallmark)[2] <- "score"
	each_P <- compare_means(score~clus,data=each_hallmark,method="kruskal.test")$p.adj
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


##Fig.3E  
genesetList <- list()
pathway_list <- c("BER.txt","FA.txt","HR.txt","MMR.txt","NER.txt","NHEJ.txt")
for(pathway in pathway_list){
    pathway_path <- paste("pathway/DDR",pathway,sep="/")
	each_pathway <- read.table(pathway_path,header = T,sep = "\t",stringsAsFactors = F)
	genesetList[[pathway]] <- each_pathway[,1]
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
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))

library(ggpubr)
stat <- c()
for(i in 3:ncol(merge_score)){
    each_hallmark <- merge_score[,c(2,i)]
	colnames(each_hallmark)[2] <- "score"
	each_P <- compare_means(score~clus,data=each_hallmark,method="kruskal.test")$p.adj
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

##Fig.3F 	
genesetList <- list()
pathway_list <- list.files("pathway/cancerpath")
for(pathway in pathway_list){
    pathway_path <- paste("pathway/cancerpath",pathway,sep="/")
	each_pathway <- read.table(pathway_path,header = F,sep = "\t",stringsAsFactors = F)
	genesetList[[pathway]] <- each_pathway[,1]
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
merge_score$clus <- factor(merge_score$clus,levels=c("C1","C2","C3","C4"))
library(ggpubr)
ggboxplot(
  merge_score, x = "clus", y = "WNT.txt",fill="clus" ,palette = c("#E59F01","#0073B3","#CC79A7","#019E73"),width=0.5
  )+stat_compare_means()+ylab("WNT.txt")
compare_means(WNT.txt ~ clus,data = merge_score)

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
