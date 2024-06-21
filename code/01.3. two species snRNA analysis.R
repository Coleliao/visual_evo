pacman::p_load(dplyr,ggplot2,Seurat,ggsci)
source("process.R")
source("tools.R")
source("Visualization.R")
zf <- readRDS('/data/zebrafinch/zf_OT_snRNA_integrated_filtered_annotated_0301.rds')
pg <- readRDS('/data/pigeon/pigeon_OT_snRNA_integrated_filtered_annotated_0301.rds')
zf$species <- 'zebrafinch'
pg$species <- 'pigeon'

# plot and filter
Idents(pg) <- 'celltype'
pg <- subset(pg,idents='MIX',invert=T) # mixture omit
zf$celltype <- ifelse(zf$celltype=='MIC','BLD',zf$celltype) # re-annotated MIC into BLD
zf$celltype <- factor(zf$celltype,levels = c('EXN','INN','AST','OLI','OPC','VLMC','BLD'))
pg$celltype <- factor(pg$celltype,levels = c('EXN','INN','AST','OLI','OPC','VLMC','BLD','MIC'))
saveRDS(zf,'zf_OT_snRNA_integrated_filtered_annotated_0314.rds')
saveRDS(pg,'pg_OT_snRNA_integrated_filtered_annotated_0314.rds')
# dimplot
cols <- pal_npg("nrc", alpha = 0.8)(9)[c(1:7,9)] # alpha 
pdf('dimplot_zf_pg_0314.pdf',9,8)
DimPlot(zf,group.by = 'celltype',label = F,cols = cols,pt.size = 0.5)
DimPlot(zf,group.by = 'celltype',label = T,cols = cols,pt.size = 0.5)
DimPlot(pg,group.by = 'celltype',label = F,cols = cols,pt.size = 0.5)
DimPlot(pg,group.by = 'celltype',label = T,cols = cols,pt.size = 0.5)
dev.off()
pdf('dimplot_zf_pg_groupby_sample_0612.pdf',9,8)
DimPlot(zf,group.by = 'sample',label = F,cols = col37,pt.size = 0.5)
DimPlot(pg,group.by = 'sample',label = F,cols = col37,pt.size = 0.5)
dev.off()
# featureplot
pdf("vlnplot_markergenes_pg_0314.pdf",10,8)
genes <- c('SLC17A6','GAD1','AQP4','PLP1','PDGFRA','COL1A2','A306-00008843RA','PTPRC')
VlnPlot(pg,features = genes,group.by = 'celltype',assay = 'RNA',pt.size = 0,stack = T)+scale_fill_manual(values=cols)
dev.off()
pdf("vlnplot_markergenes_zf_0314.pdf",8,6)
genes <- c('SLC17A6','GAD1','AQP4','PLP1','PDGFRA','COL1A2','HBAD')
VlnPlot(zf,features = genes,group.by = 'celltype',assay = 'RNA',pt.size=0,stack = T)+scale_fill_manual(values=cols)
dev.off()
# correlation
Idents(pg)='celltype'
deg1 <- FindAllMarkers(pg,only.pos = T)
Idents(zf)='celltype'
deg2 <- FindAllMarkers(zf,only.pos = T)
deg1 <- subset(deg1,deg1$p_val_adj<0.05)
deg2 <- subset(deg2,deg2$p_val_adj<0.05)
write.csv(deg1,'deg_pg_majorcelltype_0314.CSV')
write.csv(deg2,'deg_zf_majorcelltype_0314.CSV')
gene1 <- deg1 %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n=50) # top50 deg
gene2 <- deg2 %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n=50) # top50 deg
gene <- unique(c(gene1$gene,gene2$gene)) # 582 genes
source("/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/tools.R")
res <- corSeurat(pg,zf,idents1 = 'celltype',genes = gene,idents2 = 'celltype',plot = F)
pdf('cor_heatmap_2birds_snRNA_0314.pdf')
cols <- colorRampPalette(colors = rev(c("#831A21","#A13D3B","#c16c55","#daa37d",
                                    "#ECD0B4","#F2EBE5","#E8EDF1","#C8D6E7","#9EBCDB"
                                    ,"#7091C5","#4E70AF","#375093")))(100)
print(corheatmap(r=res$r,p = res$p.adj, cluster_row = F, cluster_col = F,color=cols))
dev.off()

## integration
DefaultAssay(zf) <- 'RNA'
DefaultAssay(pg) <- 'RNA'
zf <- DietSeurat(zf,assays = 'RNA')
pg <- DietSeurat(pg,assays = 'RNA')
genes <- intersect(rownames(zf),rownames(pg))
zf <- zf[genes,];pg <- pg[genes,]
zf$species='zebrafinch';pg$species='pigeon'
data <- merge(zf,pg)
data <- NormalizeData(data)
genes <- c('SNAP25','SLC17A6','GAD1','AQP4','PLP1','PDGFRA','PTPRC','COL1A2')
genes <- intersect(genes,rownames(data))
pdf("dotplot_markergenes_0308.pdf")
DotPlot(data,features = genes,group.by = 'celltype',split.by = 'species',scale = F)
dev.off()

# feature plot
# # projection markers from literature
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
## neuron markers 
pdf('featureplot_NeuronMarkers_2birds_0314.pdf')
genes <- intersect(rownames(data),neu_markers)
for(i in genes){print(FeaturePlot(data,features = i,order = T))}
dev.off()
# # projection markers from literature
pdf('featureplot_VisionMarker_2birds_0314.pdf')
genes <- c('NPNT','NTSR1','CBLN2','GRP','TAC1','RORB','CDH7','PITX2','GDA',
           'TMEM132C','RXFP2','PMFBP1')
genes <- intersect(rownames(data),genes)
for(i in genes){print(FeaturePlot(data,features = i,order = T))}
dev.off()
# pigeon
pdf('featureplot_VisionMarker_zf_0314.pdf')
genes <- c('NPNT','NTSR1','CBLN2','GRP','TAC1','RORB','CDH7','PITX2','GDA',
           'TMEM132C','RXFP2','PMFBP1')
genes <- intersect(rownames(zf),genes)
for(i in genes){print(FeaturePlot(zf,features = i,order = T))}
dev.off()


## 0319 neuron only
data <- ScaleData(data)
DefaultAssay(data) <- 'integrated'
data <- FindClusters(data,resolution = 0.2)
pdf('dimplot_0319.pdf')
DimPlot(data,label = T)
dev.off()
DefaultAssay(data) <- 'RNA'
deg=FindAllMarkers(data,min.pct=0.01,only.pos=T,logfc.threshold=0.5)
top5 <- deg %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 5)
pdf('heatmap_topdeg5_neuron_0319.pdf',15,15)
DoHeatmap(data,features = top5$gene,group.by = 'seurat_clusters')
dev.off()

## visualization
data=readRDS('integ_neuron_2birds_0317.rds')
DefaultAssay(data) <- 'integrated'
data <- FindClusters(data,resolution = 0.2)
pdf('dimplot_neurons_2birds_0323.pdf',6,6)
DimPlot(data,label = T,cols = col37,pt.size = 1)
DimPlot(data,label = F,cols = col37,pt.size = 1)
DimPlot(data,label = F,pt.size = 0.1,group.by = 'species')
dev.off()
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
pdf('featureplot_Fgenes_2birds_0323.pdf',6,6)
FeaturePlot(data,features = 'CCK',order = T,pt.size = 0.2)+scale_color_gradientn(colors=c("#F7F7E9","firebrick3"))
FeaturePlot(data,features = 'CABP7',pt.size = 0.2)+scale_color_gradientn(colors=c("#F7F7E9","firebrick3"))
FeaturePlot(data,features = 'OTX2',order = T,pt.size = 0.2)+scale_color_gradientn(colors=c("#F7F7E9","firebrick3"))
FeaturePlot(data,features = 'CBLN2',order = T,pt.size = 0.2)+scale_color_gradientn(colors=c("#F7F7E9","firebrick3"))
FeaturePlot(data,features = 'SLC18A3',order = T,pt.size = 0.2)+scale_color_gradientn(colors=c("#F7F7E9","firebrick3"))
dev.off()



