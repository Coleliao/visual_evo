
##
sc <- readRDS('/data/human_sn/human_midbrain_SC.rds') # cells
lgn <- readRDS('/data/human_sn/human_Thalamus_LG.rds') # 48856 cells
human <- merge(sc,lgn)
# gene name
e2n=read.csv("/jdfssz3/ST_STOMICS/P22Z10200N0664/2.CrossSpecies_Spatial/xiangya/Primary/VisualCtx/SC_feature_names.csv")
all(rownames(human)==e2n$ensembl_ids)
rownames(human@assays$RNA@counts)=e2n$feature_name
rownames(human@assays$RNA@data)=e2n$feature_name
Idents(human)='development_stage'
human@assays$RNA@counts=human@assays$RNA@data
m <- gsub('_',"-",row.names(human@assays$RNA@counts))
row.names(human@assays$RNA@counts) <- m
row.names(human@assays$RNA@data) <- m
human[["RNA"]]@meta.features <- data.frame(row.names = rownames(human[["RNA"]]))
datalist <- SplitObject(human,split.by='development_stage')
human <- DoIntegration(object.list = datalist,method = 'log',resolution = 0.4)
## major type annotation
pdf('FeaturePlot_human_0427.pdf',10,10)
FeaturePlot(human,features = c('SLC17A6','SLC17A7','GAD1','GAD2'),ncol = 2)
dev.off()
m <- table(human$supercluster_term)
m <- names(m[m>100])
Idents(human) <- 'supercluster_term'
human <- subset(human,idents=m)
pdf('DimPlot_human_0429.pdf',8,8)
DimPlot(human,label = T,group.by = 'seurat_clusters')
DimPlot(human,label = T,group.by = 'development_stage')
DimPlot(human,label = T,group.by = 'roi')
DimPlot(human,label = T,group.by = 'cell_type')
DimPlot(human,label = T,group.by = 'supercluster_term') #+NoLegend()
dev.off()
# re-annotate splatter
human$celltype <- ifelse(human$supercluster_term=='Splatter',
                         ifelse(human$seurat_clusters%in%c(4,14,16,23,29),'Splatter_excitatory','Splatter_inhibitory'),
                         human$supercluster_term)
pdf('DimPlot_human_celltype_0429.pdf',12,8)
DimPlot(human,label = T,group.by = 'celltype')
dev.off()
saveRDS(human,'human_LG_SC_reintegrated_0429.rds')

human$major <- as.vector(human$celltype)
human$major <- ifelse(human$major%in%c('Thalamic excitatory','Splatter_excitatory','Miscellaneous'),'Exn',human$major)
human$major <- ifelse(human$major%in%c('Midbrain-derived inhibitory','Splatter_inhibitory',
                                       'MGE interneuron','CGE interneuron','LAMP5-LHX6 and Chandelier'),'Inn',human$major)
human$major <- ifelse(human$major%in%c('Exn','Inn'),human$major,'Non-neuron')
human$region <- sapply(strsplit(human$roi,split = ' '),'[',2)
human$region_major <- paste0(human$region,'_',human$major)

pdf('human_snRNA_SC_LG_CellComposition_0615.pdf')
meta <- human@meta.data
p <- ggplot(data=meta,aes(x=region,y=..count../sum(..count..),fill=celltype))
p + geom_bar(stat = 'count',position = 'fill')+ 
  theme(panel.background=element_rect(fill='transparent',
                                      color ="gray"),
        axis.text.x = element_text(angle = 90, hjust = 0.5,
                                   vjust = 0.5,color = "black",size=9))+
  labs(title = 'cell proportion in human SC and LG', x = 'tissue',y="cell proportion")+
  scale_fill_manual(values = col37)
dev.off()


# deg 
pdf('DimPlot_human_0505.pdf',8,8)
DimPlot(human,label = F,group.by = 'development_stage',cols = c('#F28E2B','#4E79A7','#A0CBE8'),pt.size = 0.1)+NoLegend()
DimPlot(human,label = T,group.by = 'development_stage',cols = c('#F28E2B','#4E79A7','#A0CBE8'),pt.size = 0.1)
DimPlot(human,label = F,group.by = 'roi',cols = c("#32A251","#ACD98D"),pt.size = 0.1,order = T)+NoLegend()
DimPlot(human,label = T,group.by = 'roi',cols = c("#32A251","#ACD98D"),pt.size = 0.1,order = T)
DimPlot(human,group.by = 'celltype',pt.size = 0.1,cols = sample(col37,16))+NoLegend()
DimPlot(human,label = T,group.by = 'celltype',cols = sample(col37,16))
dev.off()


Idents(human) <- 'region_major'
edeg <- FindMarkers(human,ident.1='LG_Exn',ident.2='SC_Exn',min.pct=0.05,only.pos=F,logfc.threshold=0.05)
ideg <- FindMarkers(human,ident.1='LG_Inn',ident.2='SC_Inn',min.pct=0.05,only.pos=F,logfc.threshold=0.05)
edeg$celltype='Exn';edeg$gene=row.names(edeg)
ideg$celltype='Inn';ideg$gene=row.names(ideg)
deg <- rbind(edeg,ideg)
#deg <- subset(deg,deg$p_val_adj<0.05)
write.csv(deg,'deg_human_LG2SC_INNandEXN_nofilter_0429.CSV')

library(scRNAtoolVis)
deg$cluster <- deg$celltype
pdf('deg_volcanoplot_LGvsSC_human_0505.pdf',6,6)
markerVocalno(markers = deg,
              topn = 10,
              labelCol = RColorBrewer::brewer.pal(3,'Accent')[2:3],
              pforce = T,nforce = T)
dev.off()


# GO
cp=deg
deg <- subset(deg,deg$avg_log2FC>1.0)
gene <- deg$gene
ids=bitr(gene, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
erich.go.BP = enrichGO(gene = ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

pacman::p_load('DOSE','ggplot2','aPEAR','clusterProfiler','cols4all')
df <- erich.go.BP@result
df$GeneRatio <- parse_ratio(df$GeneRatio)
df$BgRatio <- parse_ratio(df$BgRatio)
df <- df[order(df$pvalue,decreasing = F),]
df$Description <- factor(df$Description,levels = unique(df$Description))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11),
                 plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) 
p1 <- ggplot(data = df[1:10,], aes(x = GeneRatio, y = Description, fill = -log10(pvalue))) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9,'YlOrRd')) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +
  labs(x = 'GeneRatio', y = '') +
  geom_text(aes(x = 0.03, 
                label = Description),
            size=5,
            hjust = 0)+ 
  theme_classic() + mytheme2
pdf('D:/Documents/BGI/B_OT视顶盖_Part2/results/04.SC-LGd差异/LG_Exn_up_GO_0507.pdf',8,8)
print(p1)
dev.off()
