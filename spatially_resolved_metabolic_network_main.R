##########################################################
#Load tab-separated expression table, normalize, and filter down to expressed genes. 
prostate<-read.table("./prostate-twelve/P1.2.tsv",sep="\t")
prostate_coord<-colnames(prostate)
x<-gsub("X","",sapply(strsplit(prostate_coord,"x"),function(x) x[1]))
x<-as.numeric(x)
y<-sapply(strsplit(prostate_coord,"x"),function(x) x[2])
y<-as.numeric(y)
df<-data.frame(x=x,y=y,total_counts=colSums(prostate))
write.csv(df,file="prostate_1.2_total.csv")


genes<-sapply(strsplit(rownames(prostate)," "),function(x) x[1])
total_reads<-colSums(prostate)
total_reads<-total_reads/median(total_reads)
prostate_norm<-prostate/matrix(rep(total_reads,nrow(prostate)),nrow=nrow(prostate),ncol=ncol(prostate),byrow = T)

expressed<-rowSums(prostate_norm>1)>=10
prostate_norm<-prostate_norm[expressed,]
prostate<-prostate[expressed,]
genes<-genes[expressed]
write.csv(prostate,file="prostate_1.2_expressed_genes.csv")

#prostate_1.2_expressed_genes.csv is the input to Python script
"prostate_spatialDE.ipynb"

#############################################################
#Metabolic pathway enrichment
source("recon_enrichment.R")
library(pheatmap)

sde_expressed<-read.csv("prostate_1.2_spatialDE_expressed_genes.csv")
sde_expressed<-arrange(sde_expressed,qval,pval)
sde_expressed$g<-as.character(sde_expressed$g)
sde_expressed$g<-sapply(strsplit(sde_expressed$g," "),function(x) x[1])

pathways_12<-recon_enrichment(sde_expressed$g[sde_expressed$qval<0.5],sde_expressed$g,ensg2symbol)

sde_expressed<-read.csv("prostate_3.3_spatialDE_expressed_genes.csv")
sde_expressed<-arrange(sde_expressed,qval,pval)
sde_expressed$g<-as.character(sde_expressed$g)
sde_expressed$g<-sapply(strsplit(sde_expressed$g," "),function(x) x[1])

pathways_33<-recon_enrichment(sde_expressed$g[sde_expressed$qval<0.5],sde_expressed$g,ensg2symbol)


sde_expressed<-read.csv("prostate_2.4_spatialDE_expressed_genes.csv")
sde_expressed<-arrange(sde_expressed,qval,pval)
sde_expressed$g<-as.character(sde_expressed$g)
sde_expressed$g<-sapply(strsplit(sde_expressed$g," "),function(x) x[1])

pathways_24<-recon_enrichment(sde_expressed$g[sde_expressed$qval<0.5],sde_expressed$g,ensg2symbol)

mat<-cbind(pathways_12$pvals,pathways_24$pvals,pathways_33$pvals)
mat<- -log10(mat)
colnames(mat)<- paste("Section",c("1.2","2.4","3.3"))
rownames(mat)<-pathways_12$pathways
mat<-mat[rowSums(mat>2)>=1,]
mat<-mat[c(4,3,1,2),c(3,1,2)]

rownames(mat)[2]<-"Glycolysis"
rownames(mat)[1]<-"OXPHOS"
pdf("./metabolic_pathway_heatmap.pdf",width=3.3,height=2)
pheatmap(mat,cluster_rows = F,cluster_cols = F,fontsize_row = 12,fontsize_col = 12)
dev.off()

###Make spatial expression plot###
single_spatial_plot('ACSL5',prostate_norm,genes,x,y,marker)

##########Re-analysis of laser-capture microdissected RNA-seq data######
setwd("./E-GEOD-6099/")
files<-list.files("./",".txt")
first_file<-read.table(files[1],header=T,sep="\t",stringsAsFactors = F,na.strings = "null")
prostate_cancer<-matrix(0,nrow=nrow(first_file),ncol=length(files))
for (i in seq_along(files)){
  temp<-read.table(files[i],header=T,sep="\t",stringsAsFactors = F,na.strings = "null")
  prostate_cancer[,i]<-temp$VALUE
}

colnames(prostate_cancer)<-sapply(strsplit(files,"_"),function(x) x[1])
rownames(prostate_cancer)<-first_file$Reporter.Identifier

probe_annot<-read.table("./GPL2013.annot",header=T,sep="\t",quote="",comment.char = "",nrows = 20000)
sample_info<-read.table("./E-GEOD-6099.sdrf.txt",header=T,sep="\t",comment.char = "",quote="")
sample_info<-sample_info[-grep("Pool",sample_info$Comment..Sample_characteristics.),]
sample_info$Source.Name<-gsub(" [12]","",sample_info$Source.Name)
ind<-match(colnames(prostate_cancer),sample_info$Source.Name)
sample_info<-sample_info[ind,]
sample_info$condition<-gsub("_[0-9]+[a-z]?$","",sample_info$Comment..Sample_source_name.)
sample_info$condition<-gsub("[0-9]$","",sample_info$condition)
sample_info<-sample_info[,c(1,6,7,13,26)]

save(probe_annot,prostate_cancer,sample_info,file="GSE6099_processed.RData")
load("GSE6099_processed.RData")
bad<-grep("^$",probe_annot$Gene.symbol)
probe_annot<-probe_annot[-bad,]
prostate_cancer<-prostate_cancer[-bad,]

mean_exprs<-rowMeans(prostate_cancer,na.rm = T)
rank_genes<-data.frame(id=probe_annot$ID,genes=probe_annot$Gene.symbol,expression=mean_exprs)
rank_genes<-arrange(rank_genes,genes,desc(expression))
bad<-duplicated(rank_genes$genes)
rank_genes<-rank_genes[!bad,]

good<- probe_annot$ID %in% rank_genes$id
prostate_cancer<-prostate_cancer[good,]
probe_annot<-probe_annot[good,]

sample_info$condition[sample_info$condition=="EPI_ADJ_PCA"]<-"Normal"
good<- sample_info$condition %in% c("Normal","PCA","PIN","MET_HR")
prostate_cancer<-prostate_cancer[,good]
sample_info<-sample_info[good,]
genes2plot<-c("SLC14A1")
genes2plot<-unique(recon3d_genes$symbol)
genes2plot<-genes2plot[genes2plot %in% probe_annot$Gene.symbol]
genes2plot<-as.character(genes2plot)
for (i in seq_along(genes2plot)){
  gene2plot<-genes2plot[i]
  df<-data.frame(Expression=prostate_cancer[probe_annot$Gene.symbol==gene2plot,],Condition=sample_info$condition)
  df$Condition<-factor(df$Condition,levels=c("Normal","PIN","PCA","MET_HR"),ordered=T)
  
  pdf(paste(gene2plot,"_GSE6099.pdf",sep=""),width=4,height=3.3)
  print(ggplot(data=df,aes(x=Condition,y=Expression)) + geom_violin(aes(fill=Condition)) +geom_point(position=position_dodge(width=0.75)) +scale_fill_manual(values=c("green","blue","orange","red"))+ ylab("Expression") + xlab("Condition") + theme(legend.position ="None",plot.title = element_text(hjust = 0.5,face='bold',size=17),axis.text = element_text(size=15),axis.title = element_text(size=17)) + ggtitle(gene2plot) + xlab("Disease condition") + ylab("Expression (relative)"))
  dev.off()
  
}

