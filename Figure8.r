library(cmapR)
library(Biobase)
library(GeneExpressionSignature) 
celline <- read.table("cellinfo_beta_1.txt",header = T,sep = " ",stringsAsFactors = F,quote="")
breast_cell <- as.vector(celline[which(celline$primary_disease=="breast cancer"),1])
gct_file<-"level5_beta_trt_cp_n720216x12328.gctx"
ds_row <- parse_gctx(gct_file, rid=c(1)) 
pattern <- paste(breast_cell, collapse = "|")
matching_columns <- grep(pattern, ds_row@cid)
colname<-strsplit(ds_row@cid,split = ":") 
small_molecules<-c()
for(i in 1:length(colname)){ 
	if(length(grep("^BRD-",colname[[i]][2]))==1){
		tmp<-strsplit(colname[[i]][2],split = "-")
		small_molecules[i]<-paste(tmp[[1]][1],"-",tmp[[1]][2],sep="")
	}
	else{
		small_molecules[i]<-colname[[i]][2]
	}
}
write.table(ds_row@cid,file="colnames.txt",sep="\t",quote=FALSE)
unique_SM <- unique(small_molecules)
files=file("compoundinfo_beta.txt",open="r")
LINCS_small_molecules<-c()
while ( TRUE ) {
  line = readLines(files, n = 1)
  if ( length(line) == 0 ) {
    break;
  }
  drug<-strsplit(line,split = "\t")
  LINCS_small_molecules<-rbind(LINCS_small_molecules,drug[[1]][1:2])	
}
common_SM<-intersect(LINCS_small_molecules[,1],unique_SM)
ds_col <- parse_gctx(gct_file, cid=c(1))
write.table(ds_col@rid ,file="rownames.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
for(i in 1:length(common_SM)){  
	posi_SM<-which(small_molecules==common_SM[i])
	posi<-intersect(matching_columns,posi_SM)
	posi
	if(length(posi)>0){
		ds <- parse_gctx(gct_file, cid=c(posi))
		exprs<-ds@mat
		exprs_1 <- rowMeans(exprs)
		filename=paste("mean/",common_SM[i],sep="")
		write.table(exprs(MergingSet),file=filename,quote=FALSE,row.names=FALSE,col.names=FALSE)
		print(end-start)
		print(length(posi))
	}
	print(i)
}


library(clusterProfiler)
library(stringr)
library(ggplot2)
library(enrichplot)
library(fgsea) 
C1DEG <- read.table("C1DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C1DEG[which(C1DEG[,2]>1&C1DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="neg",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C1up/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C1up_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C1up/C1_up.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C1DEG <- read.table("C1DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C1DEG[which(C1DEG[,2]<(-1)&C1DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="pos",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C1down/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C1down_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C1down/C1_down.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C2DEG <- read.table("C2DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C2DEG[which(C2DEG[,2]>1&C2DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="neg",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C2up/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C2up_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C2up/C2_up.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C2DEG <- read.table("C2DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C2DEG[which(C2DEG[,2]<(-1)&C2DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="pos",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C2down/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C2down_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C2down/C2_down.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C3DEG <- read.table("C3DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C3DEG[which(C3DEG[,2]>1&C3DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="neg",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C3up/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C3up_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C3up/C3_up.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C3DEG <- read.table("C3DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C3DEG[which(C3DEG[,2]<(-1)&C3DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="pos",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C3down/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C3down_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C3down/C3_down.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C4DEG <- read.table("C4DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C4DEG[which(C4DEG[,2]>1&C4DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="neg",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C4up/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C4up_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C4up/C4_up.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

C4DEG <- read.table("C4DEG.txt",header=TRUE,sep = "\t",stringsAsFactors = F,row.names=1) 
genes<- rownames(C4DEG[which(C4DEG[,2]<(-1)&C4DEG[,6]<0.01),])
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
allfile<-list.files()
allfile_remove<-list.files()
pdffile <- grep("*.pdf", allfile_remove)
file.remove(allfile_remove[pdffile]) 
all_result <- c()
for(i in 1:length(allfile_remove)){
	input_file<-allfile_remove[i]
	drug_response<-read.table(input_file,header=TRUE,row.names=1)
	drug_response_1 <-as.vector(drug_response[,1])
	names(drug_response_1) <- as.character(rownames(drug_response))
	geneList <- drug_response_1
	geneList <- sort(geneList,decreasing=TRUE)
	common_gene<-intersect(rownames(drug_response),gs[,2])
	gsea_sets <- list(genes = common_gene)
	egmt <- fgsea(pathways = gsea_sets, stats = geneList,  nproc = 1,scoreType="pos",)
	title_text<-paste("up_regulated NES: ",signif(egmt$NES,3),"pval: ",signif(egmt$pval,3),sep="")
	pdf(file=paste("~/LINCS/mean/C4down/","_",i,"_",input_file,"_",signif(egmt$NES,3),"pval ",signif(egmt$pval,3),"C4down_regulated.pdf",sep=""),width=9,height=6)
	print(plotEnrichment(gsea_sets[[1]],    geneList,ticksSize = 0.1) + labs(title=title_text))
	dev.off()
	result<-as.data.frame(cbind(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	write.table(result,file="~/LINCS/mean/C4down/C4_down.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
	print(i)
    print (result)	
}

##Fig.8A
up <- read.table(file="mean/C1up/C1_up.txt",sep="\t")
colnames(up)<-c("drug_up","ES_up","NES_up","pvalue_up")
down<-read.table(file="mean/C1down/C1_down.txt",sep="\t")
colnames(down)<-c("drug_down","ES_down","NES_down","pvalue_down")
merge_drug <- merge(up,down,by.x="drug_up",by.y="drug_down")
WTCS<-c()
for(i in 1:nrow(merge_drug)){
	if(sign(merge_drug$NES_up[i]) == sign(merge_drug$NES_down[i])){
		WTCS[i]<-0
	}
	else{
		WTCS[i]<-(merge_drug$NES_up[i]-merge_drug$NES_down[i])/2
	}
}
merge_result <- cbind(merge_drug,WTCS)
merge_result_1 <- merge_result[which(merge_result$pvalue_up<0.05&merge_result$pvalue_down<0.05),]
SM <- read.table(file="LINCS_small_molecules.txt",sep="\t",header=TRUE,quote="")
merge_result_2 <- unique(merge(merge_result_1,SM,by.x="drug_up",by.y="pert_id"))
merge_result <- merge_result_2[order(merge_result_2$WTCS),]
order_WTCS <- merge_result[order(abs(merge_result$WTCS),decreasing=T),]
WTCS_top <- order_WTCS[1:10,]
WTCS_up <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_up,type="up")
WTCS_down <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_down,type="down")
WTCS_data <- rbind(WTCS_up,WTCS_down)
WTCS_data$drug <- factor(WTCS_top$cmap_name,levels=WTCS_top$cmap_name)
ggplot(WTCS_data, aes(x = drug, y = NES, fill = type)) +
  geom_bar(stat = "identity",alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = c("#f2cf80","#E59F01"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = NES),vjust = -0.5)

##Fig.8B
up <- read.table(file="mean/C2up/C2_up.txt",sep="\t")
colnames(up)<-c("drug_up","ES_up","NES_up","pvalue_up")
down<-read.table(file="mean/C2down/C2_down.txt",sep="\t")
colnames(down)<-c("drug_down","ES_down","NES_down","pvalue_down")
down_1 <- down[match(up[,1],down[,1]),]
WTCS<-c()
for(i in 1:nrow(up)){
	if(sign(up$NES[i]) == sign(down_1$NES[i])){
		WTCS[i]<-0
	}
	else{
		WTCS[i]<-(up$NES[i]-down_1$NES[i])/2
	}
}
merge_result<-cbind(up,down_1,WTCS)
posi<-intersect(which(merge_result$pvalue_up<0.05),which(merge_result$pvalue_down<0.05))
merge_result<-merge_result[posi,]
merge_result_1 <- unique(merge(merge_result,SM,by.x="drug_up",by.y="pert_id"))
merge_result <- merge_result_1[order(merge_result_1$WTCS),]
order_WTCS <- merge_result[order(abs(merge_result$WTCS),decreasing=T),]
WTCS_top <- order_WTCS[1:10,]
WTCS_up <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_up,type="up")
WTCS_down <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_down,type="down")
WTCS_data <- rbind(WTCS_up,WTCS_down)
WTCS_data$drug <- factor(WTCS_top$cmap_name,levels=WTCS_top$cmap_name)
ggplot(WTCS_data, aes(x = drug, y = NES, fill = type)) +
  geom_bar(stat = "identity",alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = c("#7fb9d9", "#0073B3"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = NES),vjust = -0.5)

##Fig.8C
up <- read.table(file="mean/C3up/C3_up.txt",sep="\t")
colnames(up)<-c("drug_up","ES_up","NES_up","pvalue_up")
down<-read.table(file="mean/C3down/C3_down.txt",sep="\t")
colnames(down)<-c("drug_down","ES_down","NES_down","pvalue_down")
down_1 <- down[match(up[,1],down[,1]),]
WTCS<-c()
for(i in 1:nrow(up)){
	if(sign(up$NES[i]) == sign(down_1$NES[i])){
		WTCS[i]<-0
	}
	else{
		WTCS[i]<-(up$NES[i]-down_1$NES[i])/2
	}
}
merge_result<-cbind(up,down_1,WTCS)
posi<-intersect(which(merge_result$pvalue_up<0.05),which(merge_result$pvalue_down<0.05))
merge_result<-merge_result[posi,]
merge_result_1 <- unique(merge(merge_result,SM,by.x="drug_up",by.y="pert_id"))
merge_result <- merge_result_1[order(merge_result_1$WTCS),]
order_WTCS <- merge_result[order(abs(merge_result$WTCS),decreasing=T),]
WTCS_top <- order_WTCS[1:10,]
WTCS_up <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_up,type="up")
WTCS_down <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_down,type="down")
WTCS_data <- rbind(WTCS_up,WTCS_down)
WTCS_data$drug <- factor(WTCS_top$cmap_name,levels=WTCS_top$cmap_name)
ggplot(WTCS_data, aes(x = drug, y = NES, fill = type)) +
  geom_bar(stat = "identity",alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = c("#e5bcd3", "#CC79A7"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = NES),vjust = -0.5)
  
##Fig.8D
up <- read.table(file="mean/C4up/C4_up.txt",sep="\t")
colnames(up)<-c("drug_up","ES_up","NES_up","pvalue_up")
down<-read.table(file="mean/C4down/C4_down.txt",sep="\t")
colnames(down)<-c("drug_down","ES_down","NES_down","pvalue_down")
down_1 <- down[match(up[,1],down[,1]),]
WTCS<-c()
for(i in 1:nrow(up)){
	if(sign(up$NES[i]) == sign(down_1$NES[i])){
		WTCS[i]<-0
	}
	else{
		WTCS[i]<-(up$NES[i]-down_1$NES[i])/2
	}
}
merge_result<-cbind(up,down_1,WTCS)
posi<-intersect(which(merge_result$pvalue_up<0.05),which(merge_result$pvalue_down<0.05))
merge_result<-merge_result[posi,]
merge_result_1 <- unique(merge(merge_result,SM,by.x="drug_up",by.y="pert_id"))
merge_result <- merge_result_1[order(merge_result_1$WTCS),]
order_WTCS <- merge_result[order(abs(merge_result$WTCS),decreasing=T),]
WTCS_top <- order_WTCS[1:10,]
WTCS_up <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_up,type="up")
WTCS_down <- data.frame(drug=WTCS_top$cmap_name,NES=WTCS_top$NES_down,type="down")
WTCS_data <- rbind(WTCS_up,WTCS_down)
WTCS_data$drug <- factor(WTCS_top$cmap_name,levels=WTCS_top$cmap_name)
ggplot(WTCS_data, aes(x = drug, y = NES, fill = type)) +
  geom_bar(stat = "identity",alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = c("#80ceb9", "#019E73"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = NES),vjust = -0.5)