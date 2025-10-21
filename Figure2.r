##Fig.2A
candidate <- read.table("candidate.txt",header = F,sep = "\t",stringsAsFactors = F)
RPM <- read.table("TCGA_data/BRCA_miRNA_RPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
type_trans <- function(x){
    split_name <- unlist(strsplit(x,split="\\."))
	return(split_name[4])
}
sample_type <- unlist(lapply(colnames(RPM),type_trans))
RPM_1 <- RPM[,which(sample_type=="01A")]
candidate_RPM <- RPM_1[which(!is.na(match(rownames(RPM_1),candidate[,1]))),]
sample_trans <- function(x){
  split_name <- unlist(strsplit(x,split="\\."))
	merge_1 <- paste(split_name[1],split_name[2],sep="-")
	merge_2 <- paste(merge_1,split_name[3],sep="-")
	return(merge_2)
}
trans_name <- unlist(lapply(colnames(RPM_1),sample_trans))
colnames(candidate_RPM) <- trans_name
candidate_RPM_1 <- t(apply(candidate_RPM,1,scale))
colnames(candidate_RPM_1) <- trans_name

library(ConsensusClusterPlus)
N <- 4 
beta_select <- as.matrix(candidate_RPM_1)
result_cs <- ConsensusClusterPlus(beta_select,
								  maxK = 6,  seed = 123456, plot = 'pdf')
clusterCS <- result_cs[[N]][["consensusClass"]]
cluster <- data.frame(sample=names(clusterCS),clus_1=clusterCS)
cluster[which(cluster[,2]==4),3] <- "C1"
cluster[which(cluster[,2]==1),3] <- "C2"
cluster[which(cluster[,2]==2),3] <- "C3"
cluster[which(cluster[,2]==3),3] <- "C4"
colnames(cluster)[3] <- "clus"
cluster <- cluster[,-2]

##Fig.2B
BRCA_clinical <- read.table("TCGA_data/clinical.txt",header = T,sep = "\t",stringsAsFactors = F,fill=T)
patient_1 <- merge(cluster,BRCA_clinical,by.x="sample",by.y="bcr_patient_barcode")
patient_1$clus <- factor(patient_1$clus,levels=c("C1","C2","C3","C4"))
library(survival)
library(survminer)
fit <- survfit(Surv(OS.time, OS) ~ clus, data = patient_1)
ggsurvplot(fit,
       pval = TRUE,
       risk.table = TRUE,
       risk.table.col = "strata",
       ggtheme = theme_classic(),
       palette = c("#E59F01", "#0073B3","#CC79A7","#019E73")
       )

#Fig.2C
clinical_1 <- read.table("TCGA_data/clinical_1.txt",header = T,sep = "\t",stringsAsFactors = F)
clinical_1[which(clinical_1$ajcc_pathologic_tumor_stage=="Stage X"|clinical_1$ajcc_pathologic_tumor_stage=="[Discrepancy]"|clinical_1$ajcc_pathologic_tumor_stage=="[Not Available]"),9] <- NA
clinical_1$age <- "upper65"
clinical_1[which(clinical_1$age_at_diagnosis<=65),14] <- "lower65"
clinical_1[which(clinical_1$ajcc_pathologic_tumor_stage=="Stage I"),9] <- "StageI"
clinical_1[which(clinical_1$ajcc_pathologic_tumor_stage=="Stage II"),9] <- "StageII"
clinical_1[which(clinical_1$ajcc_pathologic_tumor_stage=="Stage III"),9] <- "StageIII"
clinical_1[which(clinical_1$ajcc_pathologic_tumor_stage=="Stage IV"),9] <- "StageIV"
Thorsson <- read.table("TCGA_data/Thorsson.txt",header = T,sep = "\t",stringsAsFactors = F)
BRCA_Thorsson <- Thorsson[which(Thorsson[,2]=="BRCA"),]
subtype <- BRCA_Thorsson[,c(1,4)]
subtype_1 <- merge(subtype,BRCA_clinical,by.x="TCGA.Participant.Barcode",by.y="bcr_patient_barcode")
merge_1 <- merge(subtype,clinical_1,by.x="TCGA.Participant.Barcode",by.y="bcr_patient_barcode")
merge_2 <- merge(cluster,merge_1,by.x="sample",by.y="TCGA.Participant.Barcode",all.x=T)
merge_3 <- merge_2[match(cluster[,1],merge_2[,1]),]
TCGA.Subtype <- merge_3[,c(2,11)]
TCGA.Subtype_1 <- TCGA.Subtype[which(!is.na(TCGA.Subtype[,2])),]
C1_StageI <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="StageI"),])
C2_StageI <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="StageI"),])
C3_StageI <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="StageI"),])
C4_StageI <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="StageI"),])
C1_StageII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="StageII"),])
C2_StageII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="StageII"),])
C3_StageII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="StageII"),])
C4_StageII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="StageII"),])
C1_StageIII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="StageIII"),])
C2_StageIII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="StageIII"),])
C3_StageIII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="StageIII"),])
C4_StageIII <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="StageIII"),])
C1_StageIV <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="StageIV"),])
C2_StageIV <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="StageIV"),])
C3_StageIV <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="StageIV"),])
C4_StageIV <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="StageIV"),])
chisq_matrix <- matrix(c(C1_StageI,C2_StageI,C3_StageI,C4_StageI,C1_StageII,C2_StageII,C3_StageII,C4_StageII,C1_StageIII,C2_StageIII,C3_StageIII,C4_StageIII,C1_StageIV,C2_StageIV,C3_StageIV,C4_StageIV),ncol=4)
chisq.test(chisq_matrix)#0.0001456
TCGA.Subtype_1 <- merge_3[,c(2,11)]
TCGA.Subtype_1 <- na.omit(TCGA.Subtype_1)
sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"),]
C1_1 <- length(which(sample_1[,2]=="StageI"))/nrow(sample_1)
C1_2 <- length(which(sample_1[,2]=="StageII"))/nrow(sample_1)
C1_3 <- length(which(sample_1[,2]=="StageIII"))/nrow(sample_1)
C1_4 <- length(which(sample_1[,2]=="StageIV"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"),]
C2_1 <- length(which(sample_1[,2]=="StageI"))/nrow(sample_1)
C2_2 <- length(which(sample_1[,2]=="StageII"))/nrow(sample_1)
C2_3 <- length(which(sample_1[,2]=="StageIII"))/nrow(sample_1)
C2_4 <- length(which(sample_1[,2]=="StageIV"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"),]
C3_1 <- length(which(sample_1[,2]=="StageI"))/nrow(sample_1)
C3_2 <- length(which(sample_1[,2]=="StageII"))/nrow(sample_1)
C3_3 <- length(which(sample_1[,2]=="StageIII"))/nrow(sample_1)
C3_4 <- length(which(sample_1[,2]=="StageIV"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"),]
C4_1 <- length(which(sample_1[,2]=="StageI"))/nrow(sample_1)
C4_2 <- length(which(sample_1[,2]=="StageII"))/nrow(sample_1)
C4_3 <- length(which(sample_1[,2]=="StageIII"))/nrow(sample_1)
C4_4 <- length(which(sample_1[,2]=="StageIV"))/nrow(sample_1)

data <- data.frame(risk.type=c(rep("C1",4),rep("C2",4),rep("C3",4),rep("C4",4)),TCGA_subtype=rep(c("StageI","StageII","StageIII","StageIV"),4),
                   ratio=c(C1_1,C1_2,C1_3,C1_4,C2_1,C2_2,C2_3,C2_4,C3_1,C3_2,C3_3,C3_4,C4_1,C4_2,C4_3,C4_4))
data$risk.type <- factor(data$risk.type,levels=c("C1","C2","C3","C4"))
data$TCGA_subtype <- factor(data$TCGA_subtype,levels=c("StageI","StageII","StageIII","StageIV"))
ggplot(data, aes(fill=TCGA_subtype, y=ratio, x=risk.type)) + 
    geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("#BED9E8","#AEBBD9","#7698C3","#477DB4"))+
	theme_classic()

##Fig.2D
TCGA.Subtype <- merge_3[,c(2,16)]
TCGA.Subtype_1 <- TCGA.Subtype[which(!is.na(TCGA.Subtype[,2])),]
C1_upper65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="upper65"),])
C2_upper65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="upper65"),])
C3_upper65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="upper65"),])
C4_upper65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="upper65"),])
C1_lower65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="lower65"),])
C2_lower65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="lower65"),])
C3_lower65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="lower65"),])
C4_lower65 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="lower65"),])
chisq_matrix <- matrix(c(C1_upper65,C2_upper65,C3_upper65,C4_upper65,C1_lower65,C2_lower65,C3_lower65,C4_lower65),ncol=2)
chisq.test(chisq_matrix)#1.789e-05

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"),]
C1_1 <- length(which(sample_1[,2]=="lower65"))/nrow(sample_1)
C1_2 <- length(which(sample_1[,2]=="upper65"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"),]
C2_1 <- length(which(sample_1[,2]=="lower65"))/nrow(sample_1)
C2_2 <- length(which(sample_1[,2]=="upper65"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"),]
C3_1 <- length(which(sample_1[,2]=="lower65"))/nrow(sample_1)
C3_2 <- length(which(sample_1[,2]=="upper65"))/nrow(sample_1)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"),]
C4_1 <- length(which(sample_1[,2]=="lower65"))/nrow(sample_1)
C4_2 <- length(which(sample_1[,2]=="upper65"))/nrow(sample_1)

data <- data.frame(risk.type=c(rep("C1",2),rep("C2",2),rep("C3",2),rep("C4",2)),TCGA_subtype=rep(c("lower65","upper65"),4),
                   ratio=c(C1_1,C1_2,C2_1,C2_2,C3_1,C3_2,C4_1,C4_2))
data$risk.type <- factor(data$risk.type,levels=c("C1","C2","C3","C4"))
data$TCGA_subtype <- factor(data$TCGA_subtype,levels=c("upper65","lower65"))
ggplot(data, aes(fill=TCGA_subtype, y=ratio, x=risk.type)) + 
    geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("#AEBBD9","#BED9E8"))+
	theme_classic()

##Fig.2E
TCGA.Subtype <- merge_3[,c(2,3)]
TCGA.Subtype_1 <- TCGA.Subtype[which(!is.na(TCGA.Subtype[,2])),]
C1_LumA <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="BRCA.LumA"),])
C2_LumA <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="BRCA.LumA"),])
C3_LumA <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="BRCA.LumA"),])
C4_LumA <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="BRCA.LumA"),])
C1_LumB <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="BRCA.LumB"),])
C2_LumB <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="BRCA.LumB"),])
C3_LumB <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="BRCA.LumB"),])
C4_LumB <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="BRCA.LumB"),])
C1_Basal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="BRCA.Basal"),])
C2_Basal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="BRCA.Basal"),])
C3_Basal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="BRCA.Basal"),])
C4_Basal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="BRCA.Basal"),])
C1_Normal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="BRCA.Normal"),])
C2_Normal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="BRCA.Normal"),])
C3_Normal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="BRCA.Normal"),])
C4_Normal <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="BRCA.Normal"),])
C1_Her2 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="BRCA.Her2"),])
C2_Her2 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="BRCA.Her2"),])
C3_Her2 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="BRCA.Her2"),])
C4_Her2 <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="BRCA.Her2"),])

chisq_matrix <- matrix(c(C1_LumA,C2_LumA,C3_LumA,C4_LumA,C1_LumB,C2_LumB,C3_LumB,C4_LumB,C1_Basal,C2_Basal,C3_Basal,C4_Basal,C1_Normal,C2_Normal,C3_Normal,C4_Normal,C1_Her2,C2_Her2,C3_Her2,C4_Her2),ncol=5)
chisq.test(chisq_matrix)#< 2.2e-16
TCGA.Subtype_1 <- merge_3[,c(2,3)]
TCGA.Subtype_1 <- na.omit(TCGA.Subtype_1)
sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"),]
LumA_1 <- length(which(sample_1[,2]=="BRCA.LumA"))/nrow(sample_1)
Her2_1 <- length(which(sample_1[,2]=="BRCA.Her2"))/nrow(sample_1)
LumB_1 <- length(which(sample_1[,2]=="BRCA.LumB"))/nrow(sample_1)
Normal_1 <- length(which(sample_1[,2]=="BRCA.Normal"))/nrow(sample_1)
Basal_1 <- length(which(sample_1[,2]=="BRCA.Basal"))/nrow(sample_1)

sample_2 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"),]
LumA_2 <- length(which(sample_2[,2]=="BRCA.LumA"))/nrow(sample_2)
Her2_2 <- length(which(sample_2[,2]=="BRCA.Her2"))/nrow(sample_2)
LumB_2 <- length(which(sample_2[,2]=="BRCA.LumB"))/nrow(sample_2)
Normal_2 <- length(which(sample_2[,2]=="BRCA.Normal"))/nrow(sample_2)
Basal_2 <- length(which(sample_2[,2]=="BRCA.Basal"))/nrow(sample_2)

sample_3 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"),]
LumA_3 <- length(which(sample_3[,2]=="BRCA.LumA"))/nrow(sample_3)
Her2_3 <- length(which(sample_3[,2]=="BRCA.Her2"))/nrow(sample_3)
LumB_3 <- length(which(sample_3[,2]=="BRCA.LumB"))/nrow(sample_3)
Normal_3 <- length(which(sample_3[,2]=="BRCA.Normal"))/nrow(sample_3)
Basal_3 <- length(which(sample_3[,2]=="BRCA.Basal"))/nrow(sample_3)

sample_4 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"),]
LumA_4 <- length(which(sample_4[,2]=="BRCA.LumA"))/nrow(sample_4)
Her2_4 <- length(which(sample_4[,2]=="BRCA.Her2"))/nrow(sample_4)
LumB_4 <- length(which(sample_4[,2]=="BRCA.LumB"))/nrow(sample_4)
Normal_4 <- length(which(sample_4[,2]=="BRCA.Normal"))/nrow(sample_4)
Basal_4 <- length(which(sample_4[,2]=="BRCA.Basal"))/nrow(sample_4)

data <- data.frame(risk.type=c(rep("C1",5),rep("C2",5),rep("C3",5),rep("C4",5)),TCGA_subtype=rep(c("LumA","Her2+","LumB","Normal-like","Basal-like"),4),
                   ratio=c(LumA_1,Her2_1,LumB_1,Normal_1,Basal_1,LumA_2,Her2_2,LumB_2,Normal_2,Basal_2,LumA_3,Her2_3,LumB_3,Normal_3,Basal_3,LumA_4,Her2_4,LumB_4,Normal_4,Basal_4))
data$risk.type <- factor(data$risk.type,levels=c("C1","C2","C3","C4"))
data$TCGA_subtype <- factor(data$TCGA_subtype,levels=c("Normal-like","Her2+","LumA","LumB","Basal-like"))
ggplot(data, aes(fill=TCGA_subtype, y=ratio, x=risk.type)) + 
    geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("#A4C7E2","#BED9E8","#AEBBD9","#7698C3","#477DB4"))+
	theme_classic()

##Fig.2F
status <- read.table("TCGA_data/status.txt",header = T,sep = "\t",stringsAsFactors = F)	
library(dplyr)
status$subtype <- ifelse(
  status$er_status_by_ihc == "Negative" & status$pr_status_by_ihc == "Negative" & status$her2_status_by_ihc == "Negative", 
  "TNBC",
  "noTNBC"
)
merge_data <- merge(cluster,status,by.x="sample",by.y="bcr_patient_barcode")
TCGA.Subtype_1 <- merge_data[,c(2,6)]
TCGA.Subtype_1 <- na.omit(TCGA.Subtype_1)
C1_TNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="TNBC"),])
C1_noTNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"&TCGA.Subtype_1[,2]=="noTNBC"),])
C2_TNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="TNBC"),])
C2_noTNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"&TCGA.Subtype_1[,2]=="noTNBC"),])
C3_TNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="TNBC"),])
C3_noTNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"&TCGA.Subtype_1[,2]=="noTNBC"),])
C4_TNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="TNBC"),])
C4_noTNBC <- nrow(TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"&TCGA.Subtype_1[,2]=="noTNBC"),])
chisq_matrix <- matrix(c(C1_TNBC,C1_noTNBC,C2_TNBC,C2_noTNBC,C3_TNBC,C3_noTNBC,C4_TNBC,C4_noTNBC),ncol=4)
chisq.test(chisq_matrix)

sample_1 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C1"),]
TNBC_1 <- length(which(sample_1[,2]=="TNBC"))/nrow(sample_1)
noTNBC_1 <- length(which(sample_1[,2]=="noTNBC"))/nrow(sample_1)

sample_2 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C2"),]
TNBC_2 <- length(which(sample_2[,2]=="TNBC"))/nrow(sample_2)
noTNBC_2 <- length(which(sample_2[,2]=="noTNBC"))/nrow(sample_2)

sample_3 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C3"),]
TNBC_3 <- length(which(sample_3[,2]=="TNBC"))/nrow(sample_3)
noTNBC_3 <- length(which(sample_3[,2]=="noTNBC"))/nrow(sample_3)

sample_4 <- TCGA.Subtype_1[which(TCGA.Subtype_1[,1]=="C4"),]
TNBC_4 <- length(which(sample_4[,2]=="TNBC"))/nrow(sample_4)
noTNBC_4 <- length(which(sample_4[,2]=="noTNBC"))/nrow(sample_4)

data <- data.frame(clus=c(rep("C1",2),rep("C2",2),rep("C3",2),rep("C4",2)),type=rep(c("TNBC","noTNBC"),4),
                   ratio=c(TNBC_1,noTNBC_1,TNBC_2,noTNBC_2,TNBC_3,noTNBC_3,TNBC_4,noTNBC_4))

data$clus <- factor(data$clus,levels=c("C1","C2","C3","C4"))
data$type <- factor(data$type,levels=c("noTNBC","TNBC"))
ggplot(data, aes(fill=type, y=ratio, x=clus)) + 
    geom_bar(position="fill", stat="identity")+scale_fill_manual(values=c("#C0D6EA","#ACBBDC"))+
	theme_classic()

##Fig.2G
ID_map <- read.table("GEO/GSE40267/ID_map.txt",header = T,sep = "\t",stringsAsFactors = F)
library(stringr)
strsub <- function(x){
    start_index <- str_locate(x, "MIMAT")
	result <- substr(x, start_index[1,1], start_index[1,2] + 7)
	return(result)
}
ID_map$ID <- unlist(lapply(ID_map$Accession_String,strsub))
miRNA_exp <- read.table("GEO/GSE40267/GSE40267_series_matrix.txt",header = T,sep = "\t",stringsAsFactors = F)
miRNA_exp_1 <- unique(merge(ID_map,miRNA_exp,by.x="miRNA_ID",by.y="ID_REF"))
miRNA_exp_2 <- miRNA_exp_1[,-c(1,2)]
rownames(miRNA_exp_2) <- miRNA_exp_2[,1]
miRNA_exp_2 <- miRNA_exp_2[,-1]
miRNA_exp_3 <- miRNA_exp_2[which(!is.na(match(rownames(miRNA_exp_2),candidate[,1]))),]#62
miRNA_exp_4 <- t(apply(miRNA_exp_3,1,scale))
colnames(miRNA_exp_4) <- colnames(miRNA_exp_3)

library(ConsensusClusterPlus)
N <- 4
beta_select <- as.matrix(miRNA_exp_4)
result_cs <- ConsensusClusterPlus(beta_select,
								  maxK = 6,  seed = 123456, plot = 'pdf')
clusterCS <- result_cs[[N]][["consensusClass"]]
cluster <- data.frame(sample=names(clusterCS),clus=clusterCS)
cluster[which(cluster[,2]==4),3] <- "C4"
cluster[which(cluster[,2]==1),3] <- "C1"
cluster[which(cluster[,2]==2),3] <- "C2"
cluster[which(cluster[,2]==3),3] <- "C3"
cluster$clus <- cluster[,3]
cluster <- cluster[,-3]

clinical <- read.table("GEO/GSE40267/clinical.txt",header = F,sep = "\t",stringsAsFactors = F)##181 sample
clinical_1 <- clinical[,c(21,23,43)]
time_sel <- function(x){
    return(unlist(strsplit(x,split=": "))[2])
}
clinical_1$OS.time <- unlist(lapply(clinical_1[,1],time_sel))
clinical_1$OS <- unlist(lapply(clinical_1[,2],time_sel))
clinical_2 <- clinical_1[which(clinical_1[,4]!="NA"),]
clinical_2$OS_status <- 1
clinical_2[which(clinical_2$OS=="NA"),6] <- 0
clinical_2$OS.time <- as.numeric(clinical_2$OS.time)
clinical_2$OS_status <- as.numeric(clinical_2$OS_status)
patient_1 <- merge(cluster,clinical_2,by.x="sample",by.y="V43")
patient_1$clus <- factor(patient_1$clus,levels=c("C1","C2","C3","C4"))
patient_1$OS.time  <- as.numeric(patient_1$OS.time )
patient_1$OS_status <- as.numeric(patient_1$OS_status)
library(survival)
library(survminer)
fit <- survfit(Surv(OS.time , OS_status) ~ clus, data = patient_1)
ggsurvplot(fit,
       pval = TRUE,
       risk.table = TRUE,
       risk.table.col = "strata",
       ggtheme = theme_classic(),
       palette = c("#E59F01", "#0073B3","#CC79A7","#019E73")
       )

##Fig.2H
setwd("D:/课题相关/Ferroptosis miRNA/Ferroptosis miRNA/Ferroptosis miRNA")
ID_map <- read.table("GEO/GSE22220/ID_map.txt",header = T,sep = "\t",stringsAsFactors = F)
ID_map_1 <- ID_map[,c(1,7)]
ID_map_2 <- c()
for(i in 1:nrow(ID_map_1)){
    if(grepl(",",ID_map_1[i,2])){
	    split_name <- unlist(strsplit(ID_map_1[i,2],split=","))
		each_map <- data.frame(ID=rep(ID_map_1[i,1],length(split_name)),TargetMatureName=split_name)
	}
	else{
	    each_map <- ID_map_1[i,]
	}
	ID_map_2 <- rbind(ID_map_2,each_map)
}

library(miRBaseConverter)
nameto22 <- miRNAVersionConvert(ID_map_2[,2], targetVersion = "v22", exact = TRUE,verbose = TRUE)
ID_map_3 <- merge(ID_map_2,nameto22,by.x="TargetMatureName",by.y="OriginalName")
write.table(ID_map_3,"GEO/GSE22220/ID_map_deal.txt",col.names = T,row.names = T,sep="\t",quote=F)##处理&连接的字符串

ID_map <- read.table("GEO/GSE22220/ID_map_deal.txt",header = T,sep = "\t",stringsAsFactors = F)
miRNA_exp <- read.table("GEO/GSE22220/GSE22220-GPL8178_series_matrix.txt",header = T,sep = "\t",stringsAsFactors = F)
miRNA_exp_1 <- unique(merge(ID_map,miRNA_exp,by.x="ID",by.y="ID_REF"))
miRNA_exp_2 <- miRNA_exp_1[,-c(1,2,3)]
miRNA_exp_3 <- miRNA_exp_2[which(!is.na(miRNA_exp_2[,1])),]
miRNA_exp_4 <- aggregate(. ~ Accession, data = miRNA_exp_3, mean)
rownames(miRNA_exp_4) <- miRNA_exp_4[,1]
miRNA_exp_4 <- miRNA_exp_4[,-1]
miRNA_exp_5 <- miRNA_exp_4[which(!is.na(match(rownames(miRNA_exp_4),candidate[,1]))),]#43
miRNA_exp_6 <- t(apply(miRNA_exp_5,1,scale))
colnames(miRNA_exp_6) <- colnames(miRNA_exp_5)

library(ConsensusClusterPlus)
N <- 2
beta_select <- as.matrix(miRNA_exp_6)
result_cs <- ConsensusClusterPlus(beta_select,
								  maxK = 6,  seed = 123456, plot = 'pdf')
clusterCS <- result_cs[[N]][["consensusClass"]]
cluster <- data.frame(sample=names(clusterCS),clus=clusterCS)
cluster[which(cluster[,2]==1),3] <- "C1"
cluster[which(cluster[,2]==2),3] <- "C2"
cluster$clus <- cluster[,3]
cluster <- cluster[,-3]

clinical <- read.table("GEO/GSE22220/clinical.txt",header = F,sep = "\t",stringsAsFactors = F,fill=T)##210 sample
clinical_2 <- clinical[,c(16,17,2)]
time_sel <- function(x){
    return(unlist(strsplit(x,split=": "))[2])
}
clinical_2$RFS.time <- unlist(lapply(clinical_2[,2],time_sel))
clinical_2$status <- unlist(lapply(clinical_2[,1],time_sel))
clinical_2$RFS.time <- as.numeric(clinical_2$RFS.time)
clinical_2$status <- as.numeric(clinical_2$status)
patient_1 <- merge(cluster,clinical_2,by.x="sample",by.y="V2")
library(survival)
library(survminer)
fit <- survfit(Surv(RFS.time , status) ~ clus, data = patient_1)
ggsurvplot(fit,
       pval = TRUE, #conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       #linetype = "strata", # Change line type by groups
       #surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#E59F01", "#0073B3","#CC79A7","#019E73")
       )
