pacman::p_load(SoupX,Seurat,ggplot2,glue,dplyr,DoubletFinder)
# library(DropletUtils)
source('Visualization.R')
source('process.R')
papath <- '/data/matrix_byManual/' # parent path
sam <- c('ZF1_OT','ZF2_OT')
datalist <- list()
for(i in 1:length(sam)){
  path <- paste0(papath,sam[i])
  lib <- dir(path) # different library
  for(j in 1:length(lib)){
    fil <- glue('{path}/{lib[j]}/FilterMatrix')
    raw <- glue('{path}/{lib[j]}/RawMatrix')
    # remove ambient RNA
    tmp <- run_soupx(fil,raw) # a basic filtered and adjusted Seurat object
    tmp$library <- lib[j]
    # basic QC
    tmp <- subset(tmp,nFeature_RNA>400&nFeature_RNA<3000&nCount_RNA > 800)
    # remove doublet
    tmp <- rundf(tmp)
    if(j==1){
      data=tmp
    }else{
      data=merge(data,tmp) # same sample different library merged
    }
  }
  data$sam <- sam[i]
  data$sample <- strsplit(sam[i],split = '_')[[1]][1] # individual
  datalist <- c(datalist,data)
}

# integration from different sample (only two samples here) 
data <- DoIntegration(datalist,method='sct',resolution=0.8,dims=1:15)
saveRDS(data,'zf_OT_snRNA_integrated_filtered_0229.rds') # 28822

pdf('dimplot_integrated_zf_0229.pdf')
DimPlot(data,label=T)
DimPlot(data,group.by = 'sam')
DimPlot(data,group.by = 'sample')
DimPlot(data,group.by = 'library')
dev.off()

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
genes <- intersect(major_markers,rownames(data))
pdf('featureplot_zf_sn_majorMarkers_0229.pdf')
for(i in genes){
  print(FeaturePlot(data,features = i,order = T))
}
dev.off()

Idents(data) <- 'seurat_clusters'
data <- subset(data,idents=0:24) # ncells<100 is filtered, 28685 rest
anno <- data.frame(seurat_clusters=c(1,2,3,13,14,10,22,8,19,24,4,5,23,18,9,21,
                                     11,17,0,6,12,15,16,7,20),
                   celltype=c(rep('EXN',7),rep('AST',5),'BLD','OPC','OLI','VLMC',
                              rep('INN',9)))
data$celltype <- anno[match(data$seurat_clusters,anno$seurat_clusters),'celltype']
pdf('dimplot_zf_OT_snRNA_celltype_0301.pdf')
DimPlot(data,group.by = 'celltype',label = T)
dev.off()
saveRDS(data,'zf_OT_snRNA_integrated_filtered_annotated_0301.rds') # 50111 cells

# manually remove cells expressing multiple markers 
meta <- data@meta.data
df=Seurat::Embeddings(data,reduction = 'umap')
meta <- cbind(meta,df)
library(ggplot2)
library(plotly)
plot_ly(data=meta,x=meta$UMAP_1,y=meta$UMAP_2,color = factor(meta$seurat_clusters))
library(dplyr)
tmp <- subset(meta,between(meta$UMAP_1,-1.6,-0.8) & between(meta$UMAP_2,-7.6,-6.6)) # cells omit
meta <- meta[setdiff(rownames(meta),rownames(tmp)),]   
data <- data[,rownames(meta)] # 28685 cells left


Idents(data) <- 'celltype'
exn <- subset(data,idents='EXN')
exn <- DietSeurat(exn)
saveRDS(exn,'zf_OT_EXN_snRNA_0301.rds') # 9862
inn <- subset(data,idents='INN')
inn <- DietSeurat(inn)
saveRDS(inn,'zf_OT_INN_snRNA_0301.rds') # 11359



