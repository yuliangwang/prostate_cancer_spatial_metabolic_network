recon_enrichment<-function(target,background,ensg2symbol){
  pathways<-list.files("~/Dropbox/research_projects/epigenetic_metabolic_interaction/recon204_pathways","txt",full.names = T)
  pathways<-pathways[-grep("Mis|Exchange|Una",pathways)]
  pathway_tissue_mat<-rep(1,length(pathways))
  names(pathway_tissue_mat)<-gsub(".txt","",gsub("recon204_","",sapply(strsplit(pathways,"/"),function(x) x[[8]])))
  n_total<-length(background)
  for (i in seq_along(pathways)){
    pathway_genes<-as.numeric(read.table(pathways[i],header=F,sep="\t",stringsAsFactors = F)$V1)
    ind<-match(pathway_genes,ensg2symbol$entrezgene)
    pathway_genes<-ensg2symbol$hgnc_symbol[ind[!is.na(ind)]]
    if (length(pathway_genes)>=10){
      x<- length(intersect(pathway_genes,target))-1
      m<- length(target)
      n<- n_total - m
      k<- sum(pathway_genes %in% background)
      pathway_tissue_mat[i]<- phyper(x,m,n,k,lower.tail = F)
    }
  }
  df<-data.frame(pathways=names(pathway_tissue_mat),pvals=pathway_tissue_mat,fdr=p.adjust(pathway_tissue_mat,method="BH"))
}