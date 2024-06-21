pacman::p_load(SoupX,Seurat,ggplot2,glue,dplyr,DoubletFinder)
# library(DropletUtils)
source('Visualization.R')
source('process.R')
papath <- '/data/matrix_byManual/' # parent path
sam <- c('PG1_OT_L','PG1_OT_R','PG2_OT_L','PG2_OT_R')
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
saveRDS(data,'pigeon_OT_snRNA_integrated_filtered_0229.rds')

pdf('dimplot_integrated_pigeon_0229.pdf')
DimPlot(data,label=T)
DimPlot(data,group.by = 'sam')
DimPlot(data,group.by = 'sample')
DimPlot(data,group.by = 'library')
dev.off()

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
genes <- intersect(major_markers,rownames(data))
pdf('featureplot_pg_sn_majorMarkers_0229.pdf')
for(i in genes){
  print(FeaturePlot(data,features = i,order = T))
}
dev.off()

qcplot(data,group='sam',show=F)
pdf('featureplot_qc_0301.pdf')
FeaturePlot(data,features = c('nFeature_RNA','nCount_RNA'))
dev.off()

# annotation
deg5=FindMarkers(data,ident.1=5,only.pos=T)
deg19=FindMarkers(data,ident.1=19,only.pos=T)
write.csv(deg19,'deg_cluster19_0304.CSV')
data$hemisphere <- sapply(strsplit(data$sam,'_'), '[', 3) # brain hemisphere info
anno <- data.frame(seurat_clusters=c(5,6,12,22,10,1,2,11,15,16,26,20,3,13,27,28,
                                     21,24,17,0,7,14,18,4,8,9,23,25,19),
                   celltype=c(rep('AST',3),rep('EXN',8),'OPC',rep('OLI',4),'VLMC',
                              rep('INN',9),'BLD','MIC','MIX')) #
# BLD for blood cell  MIX for mixture
deg23=FindMarkers(data,ident.1='23',only.pos=T) # A306-00008843RA -- HBA1 and HBA2
write.csv(deg23,'deg_cluster23_0301.CSV')
data$celltype <- anno[match(data$seurat_clusters,anno$seurat_clusters),'celltype']
pdf('dimplot_pg_OT_snRNA_celltype_0301.pdf')
DimPlot(data,group.by = 'celltype',label = T)
dev.off()
saveRDS(data,'pigeon_OT_snRNA_integrated_filtered_annotated_0301.rds') # 50111 cells


# another part
Idents(data) <- 'celltype'
exn <- subset(data,idents='EXN')
pdf('featureplot_exn_byCounts_0304.pdf')
FeaturePlot(exn,features=c('GAD1','SLC17A6'),slot='counts')
dev.off()
exn$inn_score <- exn@assays$RNA@counts['GAD1',]+exn@assays$RNA@counts['GAD2',]
exn=subset(exn,inn_score<1.0) 
exn <- DietSeurat(exn)
saveRDS(exn,'pigeon_OT_EXN_snRNA_0301.rds') # 14130
inn <- subset(data,idents='INN')
inn <- DietSeurat(inn)
saveRDS(inn,'pigeon_OT_INN_snRNA_0301.rds') # 20998




