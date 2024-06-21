source('Visualization.R')
source('tools.R')

library(Seurat)

######### human and bird integration ###########
# 2024/02/22
# integrate neurons from 4 species
library(Seurat)
library(ggplot2)
library(dplyr)

# gene filter -- only one2one orthologs for 4 species were kept
ort <- data.table::fread("./mart_export.txt")
ort[ort==''|ort=='ortholog_one2many'] <- NA
ort <- na.omit(ort) # any row including NA is removed
ort=ort[!(duplicated(ort$`Gene name`)|duplicated(ort$`Mouse gene name`)|duplicated(ort$`Zebra finch gene name`)),]

dim(ort) # 11216 genes left
write.csv(ort,'./orthologs_3species_one2one_0222.CSV')

##zf
d1=readRDS('/data/zebrafinch/zf_OT_EXN_snRNA_0301.rds')
d2=readRDS('/data/zebrafinch/zf_OT_INN_snRNA_0301.rds')
zf=merge(d1,d2)
DefaultAssay(zf)='RNA'
zf <- DietSeurat(zf,assays = 'RNA',counts = T,data = T) # only RNA assay is required
# change gene names
ort2 <- ort[which(ort$`Zebra finch gene name`%in%rownames(zf)),]
zf <- zf[ort2$`Zebra finch gene name`,]
rownames(zf@assays$RNA@counts) <- ort2$`Gene name`
rownames(zf@assays$RNA@data) <- ort2$`Gene name`
head(row.names(zf),10)
g1 <- row.names(zf)
# split by sample
zf$species='zebrafinch'
zf$sample=zf$sample
zflist <- SplitObject(zf,split.by = 'sample')

## pg
d1=readRDS('/data/pigeon/pigeon_OT_EXN_snRNA_0301.rds')
d2=readRDS('/data/pigeon/pigeon_OT_INN_snRNA_0301.rds')
pg=merge(d1,d2)
DefaultAssay(pg)='RNA'
pg <- DietSeurat(pg,assays = 'RNA',counts = T,data = T) # only RNA assay is required
# change gene names
ort2 <- ort[which(ort$`Zebra finch gene name`%in%rownames(pg)),] #
pg <- pg[ort2$`Zebra finch gene name`,]
rownames(pg@assays$RNA@counts) <- ort2$`Gene name`
rownames(pg@assays$RNA@data) <- ort2$`Gene name`
head(row.names(pg),10)
g2 <- row.names(pg)
# split by sample
pg$species='pigeon'
pg$sample=pg$sample
pglist <- SplitObject(pg,split.by = 'sample')

## hum
d1=readRDS('/data/human_sn/human_SC_snRNA_EXN_0303.rds')
d2=readRDS('/data/human_sn/human_SC_snRNA_INN_0303.rds')
hum=merge(d1,d2)
DefaultAssay(hum)='RNA'
hum <- DietSeurat(hum,assays = 'RNA',counts = T,data = T) # 26744 cells
hum$species <- 'human'
hum$sample <- hum$development_stage
humlist <- SplitObject(hum,split.by = 'sample')
g4 <- rownames(hum)

## integration
gc()
g <- Reduce(intersect,list(g1,g2,g4))
print(length(g))
datalist <- c(zflist,pglist,humlist)
datalist <- lapply(datalist, function(x){
  x <- x[g,]
  return(x)
})
source('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/process.R')
source('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/Visualization.R')
source('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/tools.R')
data <- DoIntegration(datalist,method = 'sct',dims = 1:15,resolution = 0.8)
saveRDS(data,'ot_neuron_4species_integrated_0417.rds')

pdf('dimplot_ot_neuron_4species_integrated_0417.pdf')
DimPlot(data,cols = col37)
DimPlot(data,group.by = 'species')
DimPlot(data,group.by = 'sample')
dev.off()

# feature plot
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
pdf('featureplot_ExnInn_markers_0222.pdf')
FeaturePlot(data,features = c('SLC17A6','SLC17A7','GAD1','GAD2'))
dev.off()



####### futher analysis ##########
data=readRDS('ot_neuron_4species_integrated_0417.rds') # only human and birds
data$mspecies <- ifelse(data$species=='human','human','bird')
# consist proportion
mtr <- consist(data@meta.data,level = 'seurat_clusters',consist = 'mspecies',plot = F)
tmp=melt(mtr,id.vars='seurat_clusters',c('bird_percent_related','human_percent_related')) 
pdf('seurat_proportion_Res0.4_HumBird_0427.pdf',8,4)
ggplot(tmp,aes(x=seurat_clusters,y=value))+geom_bar(aes(fill=variable),stat= 'identity',position='stack')+scale_fill_manual(values=c("#91CF60","#FC8D59"))
dev.off()


# species-related clusters deg analysis
data=FindClusters(data,resolution=0.4)
deg <- read.csv("./deg_allneurons_HumBird_res0.4_nofilter_0418.CSV",row.names = 1)
deg <- subset(deg,deg$p_val_adj<0.05&deg$avg_log2FC>1)
deg1 <- subset(deg,deg$cluster==1)
deg1 <- deg1[order(deg1$avg_log2FC,decreasing = T),]
deg12 <- subset(deg,deg$cluster==12)
deg12 <- deg12[order(deg12$avg_log2FC,decreasing = T),]
deg14 <- subset(deg,deg$cluster==14)
deg14 <- deg14[order(deg14$avg_log2FC,decreasing = T),]

library(clusterProfiler)
library(org.Hs.eg.db)
# compare go
gene <- list(gene1=deg1$gene,gene2=deg12$gene,gene3=deg14$gene)
gene <- lapply(gene, function(x){
  ids=bitr(x, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  x <- ids$ENTREZID
  return(x)
})
sam.go.BP <- compareCluster(gene,
                            fun = "enrichGO",
                            OrgDb = "org.Hs.eg.db",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05
)
df <- sam.go.BP@compareClusterResult
pdf('./compareGo_deg_0418.pdf',6,10)
dotplot(sam.go.BP,font.size=5,title="GO of 3cluster cells",showCategory = 10) 
dev.off()
# deg visualization
require(ComplexHeatmap)
top5 <- unique(c(deg1$gene[1:5],deg12$gene[1:5],deg14$gene[1:5]))
mtr <- scheatmap(data,features = top5,idents = 'seurat_clusters',plot = T)
col <- colorRampPalette(c("#4DBBD5","white","#E64B35"))(20)
pdf('heatmap_top5deg_clu1_12_14_0427.pdf',12,6)
Heatmap(mtr,cluster_columns = F,cluster_rows = F,show_column_names = T,col=col)
dev.off()


pdf('vlnplot_CBLN2_CABP7.pdf',10,10)
VlnPlot(data,features=c('CBLN2','CABP7'),pt.size=0,ncol=1)
dev.off()
pdf('featureplot_CABP7_OTX2_CBLN2_split_0418.pdf',10,6)
FeaturePlot(data,features='CABP7',split.by='mspecies')
FeaturePlot(data,features='OTX2',split.by='mspecies',order=T)
FeaturePlot(data,features='CBLN2',split.by='mspecies',order=T)
dev.off()



########### Focus on the thalamic projection cell group to explore ###########
setwd('clu10')
clu10 <- subset(data,idents=10) #
pdf('dimplot_0418.pdf')
DimPlot(clu10,label=T)
DimPlot(clu10,label=T,group.by='species')
DimPlot(clu10,label=T,group.by='mspecies')
dev.off()
df=Embeddings(clu10,reduction='umap')
clu10$UMAP_1 <- df[,1]
clu10$UMAP_2 <- df[,2]
clu10 <- subset(clu10,UMAP_1 < 0 & UMAP_2 > 0 & UMAP_1 > -4) # 
DefaultAssay(clu10) <- 'RNA';clu10 <- NormalizeData(clu10)
pdf('featureplot_CABP7_CBLN2_0418.pdf',10,5)
FeaturePlot(clu10,features='CABP7',split.by='mspecies')
FeaturePlot(clu10,features='CBLN2',split.by='mspecies',order=T)
dev.off()
Idents(clu10) <- 'mspecies'
deg <- FindMarkers(clu10,ident.1='bird',min.pct=0.01,only.pos=F) # 
deg <- deg[order(deg$avg_log2FC,decreasing = T),]
write.csv(deg,'deg_clu10_bird2hum_0419.CSV')



## 0614 
data=readRDS('./ot_neuron_4species_integrated_0417.rds')
data=FindClusters(data,resolution=0.4)
DefaultAssay(data)='RNA'
pdf('figS3_supplement_dimplot_featureplot_0614.pdf')
DimPlot(data,group.by = 'sample',cols = RColorBrewer::brewer.pal(7,'Accent'))
DimPlot(data,group.by = 'sample',cols = RColorBrewer::brewer.pal(7,'Accent'))+NoLegend()
FeaturePlot(data,features = 'GAD1')
FeaturePlot(data,features = 'SLC17A6')
dev.off()
meta=data@meta.data
meta$ann <- factor(meta$ann,levels = c(paste0('Exn',1:6),paste0('Inn',1:10)))
pdf('figS3_CellCount_0614.pdf',8,8)
p <- ggplot(data=meta,aes(x=ann,y=..count..,fill=species))
p + geom_bar(stat = 'count',position = 'dodge',color='black')+ 
  theme(panel.background=element_rect(fill='transparent',
                                      color ="gray"),
        axis.text.x = element_text(angle = 90, hjust = 0.5,
                                   vjust = 0.5,color = "black",size=9))+
  scale_fill_brewer(palette = "Set3")
dev.off()

