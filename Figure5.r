##Fig.5A
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
mRNA_expression <- read.table("TPM.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)##181 sample
PDL1_exp <- data.frame(sample=colnames(mRNA_expression),express=as.numeric(mRNA_expression[which(rownames(mRNA_expression)=="CD274"),]))
PDL1_exp$sample_1 <- gsub("\\.", "-", PDL1_exp$sample)
merge_score <- merge(GSVA_score_1,PDL1_exp,by.x="sample",by.y="sample_1")
merge_score$express1 <- log10(merge_score$express+1)
p <- ggplot(merge_score, aes(x=score, y=express1)) +
  geom_point(color = "#69b7c9") +
  geom_smooth(method=lm ,color='black', fill="#e6e6e6", se=TRUE) +
  theme_bw() +
  theme(text = element_text(family = "serif", size = 20),
        panel.grid = element_blank()
        )+
  theme(text = element_text(size = 12)) +
  stat_cor(method = "spearman")
library(ggExtra)
ggMarginal(p, type = "density", fill = "#bcd2ed")

##Fig.5B
cbio <- read.table("TCGA_data/combined_study_clinical_data.tsv",header = T,sep = "\t",stringsAsFactors = F)
cbio_1 <- cbio[,c(2,29,30,31)]
cbio_1$TMB <- cbio_1$Mutation.Count/38
merge_score <- merge(GSVA_score_1,cbio_1,by.x="sample",by.y="Patient.ID")
merge_score$TMB1 <- log10(merge_score$TMB+1)
library(ggplot2)
p <- ggplot(merge_score, aes(x=score, y=TMB1)) +
  geom_point(color = "#69b7c9") +
  geom_smooth(method=lm ,color='black', fill="#e6e6e6", se=TRUE) +
  theme_bw() +
  theme(text = element_text(family = "serif", size = 20),
        panel.grid = element_blank()
        )+
  theme(text = element_text(size = 12)) +
  stat_cor(method = "spearman") 
cor.test(merge_score$score,merge_score$TMB1,method="spearman")
library(ggExtra)
ggMarginal(p, type = "density", fill = "#bcd2ed")

##Fig.5C
C1DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNA_expression <- read.table("GEO/GSE194040/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
mRNA_expression <- na.omit(mRNA_expression)
library(IOBR)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("X", "", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down
clinical <- read.table("GEO/GSE194040/clinical.txt",header = T,sep = "\t",stringsAsFactors = F)
immune_clinical <- clinical[which(clinical[,7]=="arm: Paclitaxel + Pembrolizumab"),]
merge_score <- merge(GSVA_score_1,immune_clinical,by.x="sample",by.y="patient")
library(pROC)
library(ggplot2)
merge_score[,4] <- as.numeric(merge_score[,4])
rocobj1 <- roc(merge_score[,9],merge_score[,4],direction="<")
auc <- round(auc(rocobj1),4)
auc
ggroc(rocobj1,linetype=1,size=2,legacy.axes = T,color="#e07a5f")+
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
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(rocobj1$auc,2)))

##Fig.5D 
C1DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
mRNA_expression <- read.table("GEO/GSE194040/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
mRNA_expression <- as.matrix(mRNA_expression)
mRNA_expression <- na.omit(mRNA_expression)
library(GSVA)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("X", "", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down
clinical <- read.table("GEO/GSE194040/clinical.txt",header = T,sep = "\t",stringsAsFactors = F)
immune_clinical <- clinical[which(clinical[,7]=="arm: Paclitaxel + Pembrolizumab"),]
merge_score <- merge(GSVA_score_1,immune_clinical,by.x="sample",by.y="patient")
library(pROC)
library(ggplot2)
rocobjscore <- roc(merge_score[,9],merge_score[,4],direction="<")
aucscore <- round(auc(rocobjscore),4)

PDL1_exp <- data.frame(sample_1=colnames(mRNA_expression),express=mRNA_expression[which(rownames(mRNA_expression)=="CD274"),])
PDL1_exp$sample <- gsub("X", "", PDL1_exp$sample_1)
merge_score <- merge(immune_clinical,PDL1_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjPDL1 <- roc(merge_score[,6],merge_score[,9],direction="<")
aucPDL1 <- round(auc(rocobjPDL1),4)

PD1_exp <- data.frame(sample_1=colnames(mRNA_expression),express=mRNA_expression[which(rownames(mRNA_expression)=="PDCD1"),])
PD1_exp$sample <- gsub("X", "", PD1_exp$sample_1)
merge_score <- merge(immune_clinical,PD1_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjPD1 <- roc(merge_score[,6],merge_score[,9],direction="<")
aucPD1 <- round(auc(rocobjPD1),4)

CTLA4_exp <- data.frame(sample_1=colnames(mRNA_expression),express=mRNA_expression[which(rownames(mRNA_expression)=="CTLA4"),])
CTLA4_exp$sample <- gsub("X", "", CTLA4_exp$sample_1)
merge_score <- merge(immune_clinical,CTLA4_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjCTLA4 <- roc(merge_score[,6],merge_score[,9],direction="<")
aucCTLA4 <- round(auc(rocobjCTLA4),4)

CD8T <- read.table("signature/CD8_T_cell.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
clinical <- read.table("GEO/GSE194040/clinical.txt",header = T,sep = "\t",stringsAsFactors = F)
immune_clinical <- clinical[which(clinical[,7]=="arm: Paclitaxel + Pembrolizumab"),]
cd8texp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),CD8T[,1]))),]
cd8texp_1 <- colMeans(cd8texp)
cd8t_exp <- data.frame(sample_1=colnames(mRNA_expression),express=cd8texp_1)
cd8t_exp$sample <- gsub("X", "", cd8t_exp$sample_1)
merge_score <- merge(immune_clinical,cd8t_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjcd8t <- roc(merge_score[,6],merge_score[,9],direction="<")
auccd8t <- round(auc(rocobjcd8t),4)

CYT <- read.table("signature/CYT.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
CYTexp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),CYT[,1]))),]
library(psych)
CYT_score <- unlist(lapply(CYTexp,geometric.mean,2))
CYT_exp <- data.frame(sample_1=colnames(mRNA_expression),express=as.numeric(CYT_score))
CYT_exp$sample <- gsub("X", "", CYT_exp$sample_1)
merge_score <- merge(immune_clinical,CYT_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjCYT <- roc(merge_score[,6],merge_score[,9],direction="<")
aucCYT <- round(auc(rocobjCYT),4)

CTL <- read.table("signature/CTL.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
CTLexp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),CTL[,1]))),]
CTLexp_1 <- colMeans(CTLexp)
CTL_exp <- data.frame(sample_1=colnames(mRNA_expression),express=CTLexp_1)
CTL_exp$sample <- gsub("X", "", CTL_exp$sample_1)
merge_score <- merge(immune_clinical,CTL_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjCTL <- roc(merge_score[,6],merge_score[,9],direction="<")
aucCTL <- round(auc(rocobjCTL),4)

IFN <- read.table("signature/IFN.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
IFNexp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),IFN[,1]))),]
IFNexp_1 <- colMeans(IFNexp)
IFN_exp <- data.frame(sample_1=colnames(mRNA_expression),express=IFNexp_1)
IFN_exp$sample <- gsub("X", "", IFN_exp$sample_1)
merge_score <- merge(immune_clinical,IFN_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjIFN <- roc(merge_score[,6],merge_score[,9],direction="<")
aucIFN <- round(auc(rocobjIFN),4)

Tinf <- read.table("signature/T_cell_inflamed.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Tinfexp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),Tinf[,1]))),]
Tinfexp_1 <- colMeans(Tinfexp)
Tinf_exp <- data.frame(sample_1=colnames(mRNA_expression),express=Tinfexp_1)
Tinf_exp$sample <- gsub("X", "", Tinf_exp$sample_1)
merge_score <- merge(immune_clinical,Tinf_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjTinf<- roc(merge_score[,6],merge_score[,9],direction="<")
aucTinf <- round(auc(rocobjTinf),4)

library(sparrow)
LRRC15 <- read.table("signature/LRRC15.CAF",header=F,sep="\t",quote="",stringsAsFactors=F)
LRRC15exp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),LRRC15[,1]))),]
LRRC15exp_1 <- eigenWeightedMean(LRRC15exp)
LRRC15_exp <- data.frame(sample_1=colnames(mRNA_expression),express=as.numeric(LRRC15exp_1$score))
LRRC15_exp$sample <- gsub("X", "", LRRC15_exp$sample_1)
merge_score <- merge(immune_clinical,LRRC15_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjLRRC15<- roc(merge_score[,6],merge_score[,9],direction="<")
aucLRRC15 <- round(auc(rocobjLRRC15),4)

NLRP3 <- read.table("signature/NLRP3.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
library(GSVA)
signature <- list(NLRP3[,1])
ssgsea <- ssgseaParam(mRNA_expression,signature)
GSVA_score <- gsva(ssgsea)
NLRP3_exp <- data.frame(sample_1=colnames(mRNA_expression),express=as.numeric(GSVA_score))
NLRP3_exp$sample <- gsub("X", "", NLRP3_exp$sample_1)
merge_score <- merge(immune_clinical,NLRP3_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjNLRP3<- roc(merge_score[,6],merge_score[,9],direction="<")
aucNLRP3 <- round(auc(rocobjNLRP3),4)

CD8 <- read.table("signature/cd8.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
CD8exp <- mRNA_expression[which(!is.na(match(rownames(mRNA_expression),CD8[,1]))),]
CD8exp_1 <- colSums(CD8exp)
CD8_exp <- data.frame(sample_1=colnames(mRNA_expression),express=CD8exp_1)
CD8_exp$sample <- gsub("X", "", CD8_exp$sample_1)
merge_score <- merge(immune_clinical,CD8_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjCD8<- roc(merge_score[,6],merge_score[,9],direction="<")
aucCD8 <- round(auc(rocobjCD8),4)## 0.6333

library(IOBR)
IPS <- deconvo_ips(mRNA_expression,plot=FALSE)
IPS_exp <- data.frame(sample_1=colnames(mRNA_expression),express=IPS$IPS_IPS)
IPS_exp$sample <- gsub("X", "", IPS_exp$sample_1)
merge_score <- merge(immune_clinical,IPS_exp,by.x="patient",by.y="sample")
library(pROC)
library(ggplot2)
rocobjIPS<- roc(merge_score[,6],merge_score[,9],direction="<")
aucIPS <- round(auc(rocobjIPS),4)

library(pROC)
library(ggplot2)
auc <- data.frame(Name = c("C1Sig","PDL1","PD1","CD8T","CYT","CTL","IFNR","Tinf","NLRP3","CD8","IPS"),
AUC = c(aucscore,aucPDL1,aucPD1,auccd8t,aucCYT,aucCTL,aucIFN,aucTinf,aucNLRP3,aucCD8,aucIPS))
auc <- auc[order(auc$AUC),]
auc$fill <- c(rep("EBC691",10),"e07a5f")
auc$text_size <- 12
auc$text_color <- rep("black",11)

mydata <- auc
rx<-10.5 * 35.43307
ry<-7.5 * 35.43307
rr<-5.00 * 35.43307
circos <- data.frame(paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .2 * rr, "\" style=\"stroke:#b7b7b7; stroke-width:1; stroke-dasharray: 10 5; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .4 * rr, "\" style=\"stroke:#6c6b6b; stroke-width:2; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .6 * rr, "\" style=\"stroke:#b7b7b7; stroke-width:1; stroke-dasharray: 10 5; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .8*rr, "\" style=\"stroke:#6c6b6b; stroke-width:2; fill:none\"/>", sep = "")
)
angle_base <- 360/nrow(mydata)
for (i in 1: nrow(mydata)) {
  mydata[i, 6] <- angle_base * (i-1)
} 
names(mydata)[6] <- "angle"
cx1<-rx - 35/(nrow(mydata)/6) - 20 * mydata$AUC
cy1<-ry - mydata$AUC * 59 / 0.25
cx2<-rx + 35/(nrow(mydata)/6) + 20 * mydata$AUC
cy2<-cy1
mydata$bezier <- paste("<path d=\"M",rx,",",ry," C", cx1, ",", cy1, ",", cx2, ",", cy2, ",",rx,",",ry,"\" transform=\"rotate(", mydata$angle, " ",rx," ",ry,")\"", " fill=\"#", mydata$fill, "\" stroke = \"black\" stroke-width = \"0.5px\" />", sep = "")
tx<-rx - nchar(mydata$Name) * 3
ty<-ry - rr - 6
mydata$text <- paste("<text x=\"", tx, "\" y=\"", ty, "\" font-size=\"", mydata$text_size, "\" transform=\"rotate(", mydata$angle, " ",rx," ",ry,")\"", " fill=\"#", mydata$text_color, "\" >", mydata$Name, "</text>", sep = "")
first_line <- data.frame("<?xml version=\"1.0\" standalone=\"no\"?>",
                         "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                         "",
                         paste("<svg id=\"svg\" width=\"744.0945\" height=\"1052.362\">", "\t")#画布设置为21cm * 29.7cm（宽高）的A4纸
)
write.table(first_line[1, 1], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(first_line[1, 2], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 3], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 4], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 5], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 1], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 2], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 3], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 4], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(mydata$bezier, "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(mydata$text, "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
last_line <- data.frame(paste("</svg>"))
write.table(last_line[1, 1], "radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

#Fig.5E
C1DEG <- read.table("C1DEG.txt",header = T,sep = "\t",stringsAsFactors = F)
mRNA_expression <- read.table("TIGER/expression_bulk/Melanoma-Nathanson_2017.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("X", "", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down

clinical <- read.table("TIGER/cliniclal_bulk/Melanoma-Nathanson_2017.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,fill=T)
clinical <- clinical[which(clinical[,7]!="UNK"),]
clinical[which(clinical[,7]=="R"),18] <- 1
clinical[which(clinical[,7]=="N"),18] <- 0
merge_score <- merge(GSVA_score_1,clinical,by.x="sample",by.y="sample_id")
ggboxplot(
  merge_score, x = "response_NR", y = "score",fill="response_NR" ,palette = c("#DDECC4","#F8B291","red","blue"),width=0.5
  )+stat_compare_means()+ylab("score")
compare_means(score ~ response_NR,data = merge_score)

library(pROC)
library(ggplot2)
merge_score[,21] <- as.numeric(merge_score[,21])
rocobj <- roc(merge_score[,21],merge_score[,4],direction="<")
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

#Fig.5F
mRNA_expression <- read.table("TIGER/expression_bulk/Melanoma-PRJEB23709.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("X", "", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down

clinical <- read.table("TIGER/cliniclal_bulk/Melanoma-PRJEB23709.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,fill=T)
clinical <- clinical[which(clinical[,7]!="UNK"),]
clinical[which(clinical[,7]=="R"),18] <- 1
clinical[which(clinical[,7]=="N"),18] <- 0
merge_score <- merge(GSVA_score_1,clinical,by.x="sample",by.y="sample_id")
ggboxplot(
  merge_score, x = "response_NR", y = "score",fill="response_NR" ,palette = c("#DDECC4","#F8B291","red","blue"),width=0.5)+
 stat_compare_means()+ylab("score")
compare_means(score ~ response_NR,data = merge_score)

library(pROC)
library(ggplot2)
merge_score[,21] <- as.numeric(merge_score[,21])
rocobj <- roc(merge_score[,21],merge_score[,4],direction="<")
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

#Fig.5G
mRNA_expression <- read.table("TIGER/expression_bulk/NSCLC_GSE126044.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,row.names=1)
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- gsub("X", "", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down

clinical <- read.table("TIGER/cliniclal_bulk/NSCLC_GSE126044.Response.tsv",header = T,sep = "\t",stringsAsFactors = F,fill=T)
clinical <- clinical[which(clinical[,7]!="UNK"),]
clinical[which(clinical[,7]=="R"),18] <- 1
clinical[which(clinical[,7]=="N"),18] <- 0
merge_score <- merge(GSVA_score_1,clinical,by.x="sample",by.y="sample_id")

library(pROC)
library(ggplot2)
merge_score[,21] <- as.numeric(merge_score[,21])
rocobj <- roc(merge_score[,21],merge_score[,4],direction="<")
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
  
##Fig.5H
mRNA_expression <- read.table("TIGER/Miao/expression.txt",header = T,sep = "\t",stringsAsFactors = F)
library(clusterProfiler)
eg <- bitr(mRNA_expression[,1],fromType="ENSEMBL",toType="SYMBOL",OrgDb = org.Hs.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("SYMBOL","GENETYPE"),
  keytype = "ENTREZID"
)
protein_coding_genes <- unique(gene_info$SYMBOL[gene_info$GENETYPE == "protein-coding"])
eg_1 <- eg[which(!is.na(match(eg[,2],protein_coding_genes))),]
mRNA_expression_1 <- merge(eg_1,mRNA_expression,by.x="ENSEMBL",by.y="gene_id")
mRNA_expression_1 <- mRNA_expression_1[,-1]
mRNA_expression_mean <- aggregate(.~SYMBOL, data = mRNA_expression_1, mean)
rownames(mRNA_expression_mean) <- mRNA_expression_mean[,1]
mRNA_expression <- mRNA_expression_mean[,-1]
up <- rownames(C1DEG[which(C1DEG[,2]>0),])
down <- rownames(C1DEG[which(C1DEG[,2]<0),])
genesetList <- list(up=up,down=down)
library(GSVA)
mRNA_expression <- as.matrix(mRNA_expression)
ssgsea <- ssgseaParam(mRNA_expression,genesetList)
GSVA_score <- gsva(ssgsea)
GSVA_score_1 <- t(GSVA_score)
GSVA_score_1 <- as.data.frame(GSVA_score_1)
GSVA_score_1$sample <- sub("^(([^_]*_[^_]*))_.*$", "\\1", rownames(GSVA_score_1))
GSVA_score_1$score <- GSVA_score_1$up-GSVA_score_1$down

clinical <- read.table("TIGER/Miao/clinical.txt",header = T,sep = "\t",stringsAsFactors = F,fill=T)##没有响应信息
clinical <- clinical[which(clinical$response_category!="intermediate benefit"),]
clinical[which(clinical[,20]=="clinical benefit"),21] <- 1
clinical[which(clinical[,20]=="no clinical benefit"),21] <- 0
merge_score <- merge(GSVA_score_1,clinical,by.x="sample",by.y="patient_id")

library(pROC)
library(ggplot2)
merge_score[,21] <- as.numeric(merge_score[,21])
rocobj <- roc(merge_score[,21],merge_score[,4],direction="<")
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
