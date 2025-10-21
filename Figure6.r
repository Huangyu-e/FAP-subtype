#Fig.6A
C1DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
library(IOBR)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down
BRCA_clinical <- read.table("TCGA_data/clinical.txt",header = T,sep = "\t",stringsAsFactors = F,fill=T)
merge_data <- merge(BRCA_clinical,GSVA_score_1,by.x="bcr_patient_barcode",by.y="sample")
merge_data$type <- "high"
merge_data[which(merge_data$score<=median(merge_data$score)),10] <- "low"
library(survival)
library(survminer)
fit <- survfit(Surv(OS.time, OS) ~ type, data = merge_data)
ggsurvplot(fit,
       pval = TRUE,
       risk.table = TRUE,
       risk.table.col = "strata",
       ggtheme = theme_classic(),
       palette = c("#f4a261","#2a9d8f")
       )

##Fig.6B
mRNA_expression <- read.table("GEO/GSE22220/GSE22220-GPL6098_series_matrix.txt/GSE22220-GPL6098_series_matrix.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNAID_map <- read.table("GEO/GSE22220/mRNA_IDmap.txt",header = T,sep = "\t",stringsAsFactors = F,quote="")
map_symbol <- mRNAID_map[match(mRNA_expression[,1],mRNAID_map[,1]),2]
mRNA_expression[,1] <- map_symbol
miRNA_exp_mean <- aggregate(.~ID_REF, data = mRNA_expression, mean)
rownames(miRNA_exp_mean) <- miRNA_exp_mean[,1]
miRNA_exp_mean <- miRNA_exp_mean[,-1]
mRNA_expression <- miRNA_exp_mean
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
library(IOBR)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down
GSVA_score_1$sample <- rownames(GSVA_score_1)

mRNA_sample <- read.table("GEO/GSE22220/mRNA_sample.txt",header = T,sep = "\t",stringsAsFactors = F)
miRNA_sample <- read.table("GEO/GSE22220/miRNA_sample.txt",header = T,sep = "\t",stringsAsFactors = F)
sample_map <- read.table("GEO/GSE22220/GSE22220_miRNA_mRNA_sample_associations (1).txt/GSE22220_miRNA_mRNA_sample_associations.txt",header = T,sep = "\t",stringsAsFactors = F)
sample_1 <- merge(mRNA_sample,sample_map,by.x="Sample_title",by.y="mRNA.sample.Name")
colnames(sample_1)[2] <- "mRNA_ID"
sample_2 <- merge(miRNA_sample,sample_1,by.x="Sample_title",by.y="miRNA.sample.Name")
colnames(sample_2)[2] <- "miRNA_ID"
merge_data <- merge(GSVA_score_1,sample_2,by.x="sample",by.y="mRNA_ID")
clinical <- read.table("GEO/GSE22220/clinical.txt",header = F,sep = "\t",stringsAsFactors = F,fill=T)##210 sample
clinical_2 <- clinical[,c(16,17,2)]
time_sel <- function(x){
    return(unlist(strsplit(x,split=": "))[2])
}
clinical_2$RFS.time <- unlist(lapply(clinical_2[,2],time_sel))
clinical_2$status <- unlist(lapply(clinical_2[,1],time_sel))
clinical_2$RFS.time <- as.numeric(clinical_2$RFS.time)
clinical_2$status <- as.numeric(clinical_2$status)
clinical_3 <- merge(clinical_2,merge_data,by.x="V2",by.y="miRNA_ID")

clinical_3$type <- "high"
clinical_3[which(clinical_3$score<=median(clinical_3$score)),12] <- "low"

library(survival)
library(survminer)
fit <- survfit(Surv(RFS.time, status) ~ type, data = clinical_3)
ggsurvplot(fit,
       pval = TRUE,
       risk.table = TRUE,
       risk.table.col = "strata",
       ggtheme = theme_classic(),
       palette = c("#f4a261","#2a9d8f")
       )

##Fig.6C
pre_drug <- read.table("TCGA_data/drug/drug.txt",header = T,sep = "\t",stringsAsFactors = F)	
drug_correction <-  read.table("GDCdata/TCGA-BRCA/Transcriptome_Profiling/clinical/drugCorrection.txt",head=TRUE,sep="\t",fill=TRUE,quote="")
drug <- merge(pre_drug,drug_correction,by.x="pharmaceutical_therapy_drug_name",by.y="OldName")
C1DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
library(IOBR)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("\\.", "-", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down
merge_data <- merge(drug,GSVA_score_1,by.x="bcr_patient_barcode",by.y="sample")
response <- c("Complete Response","Partial Response","Stable Disease","Clinical Progressive Disease")
merge_cluster_1 <- merge_data[which(!is.na(match(merge_data[,3],response))),]
select_drug <- names(table(merge_cluster_1[,4]))[which(table(merge_cluster_1[,4])>20)]
merge_cluster_2 <- merge_cluster_1[which(!is.na(match(merge_cluster_1[,4],select_drug))),]
library(dplyr)
merge_cluster_2 <- merge_cluster_2 %>%
  mutate(response = case_when(
    treatment_best_response %in% c("Complete Response", "Partial Response") ~ "0",
    treatment_best_response %in% c("Stable Disease", "Clinical Progressive Disease") ~ "1",
    TRUE ~ treatment_best_response
  ))
merge_cluster_2$response <- as.numeric(merge_cluster_2$response)
i="Fluorouracil" 
each_drug <- merge_cluster_2[which(merge_cluster_2[,4]==i),]
library(pROC)
library(ggplot2)
rocobj <- roc(each_drug[,8],each_drug[,7])
auc <- round(auc(rocobj),4)
auc

ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#e07a5f")+
geom_abline(intercept = 0,
              slope = 1,
              color = "black",
              size = 1,
              linetype = "dashed")+
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
  theme_classic()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(rocobj$auc,2)))

##Fig.6D
i="Epirubicin" 
each_drug <- merge_cluster_2[which(merge_cluster_2[,4]==i),]
library(pROC)
library(ggplot2)
rocobj <- roc(each_drug[,8],each_drug[,7])
auc <- round(auc(rocobj),4)
auc

ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#e07a5f")+
geom_abline(intercept = 0,
              slope = 1,
              color = "black",
              size = 1,
              linetype = "dashed")+
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
  theme_classic()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(rocobj$auc,2)))
  
##Fig.6E
i="Tamoxifen" 
each_drug <- merge_cluster_2[which(merge_cluster_2[,4]==i),]
library(pROC)
library(ggplot2)
rocobj <- roc(each_drug[,8],each_drug[,7])
auc <- round(auc(rocobj),4)
auc

ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#e07a5f")+
geom_abline(intercept = 0,
              slope = 1,
              color = "black",
              size = 1,
              linetype = "dashed")+
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
  theme_classic()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(rocobj$auc,2)))
  
##Fig.6F
i="Anastrozole" 
each_drug <- merge_cluster_2[which(merge_cluster_2[,4]==i),]
library(pROC)
library(ggplot2)
rocobj <- roc(each_drug[,8],each_drug[,7])
auc <- round(auc(rocobj),4)
auc

ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#e07a5f")+
geom_abline(intercept = 0,
              slope = 1,
              color = "black",
              size = 1,
              linetype = "dashed")+
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
  theme_classic()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(rocobj$auc,2)))