# 鸟类OT 分别与人的SC和LGN的比较 2024/03/25
# 区分EXN和INN （先不管non-neuron了） 分TF和non-TF做
# ?鸟类需要区分亚群去做一下吗(取决于是想讨论到细胞类型 还是差异基因水平) ? V1需要放进来吗
# 看HVG
library(Seurat)
options(stringsAsFactors=F)
data <- readRDS('./neuron_integration/integ_neuron_2birds_0317.rds')

### Part1 human data processing ####
# ensemble ID to gene name
sc <- readRDS('/data/human_sn/human_midbrain_SC.rds') # cells
lgn <- readRDS('/data/human_sn/human_Thalamus_LG.rds') # 48856 cells
v1 <- readRDS('/data/human_sn/Human_primary_visual_cortex_V1C.rds') # 33484 cells
sc$region='SC';lgn$region='LGd';v1$region='V1'
human <- merge(sc,c(lgn,v1)) # 112247 cells
e2n=read.csv("/jdfssz3/ST_STOMICS/P22Z10200N0664/2.CrossSpecies_Spatial/xiangya/Primary/VisualCtx/SC_feature_names.csv")
all(rownames(human)==e2n$ensembl_ids)
rownames(human@assays$RNA@counts)=e2n$feature_name
rownames(human@assays$RNA@data)=e2n$feature_name
Idents(human) <- 'cell_type'
human <- subset(human,idents='neuron') # 97978 nuclei
# down sample
Idents(human) <- 'region'
human <- subset(human,downsample=20000)
saveRDS(human,'human_snRNA_3region_genenamed_0327.rds') # 
human@assays$RNA@counts <- human@assays$RNA@data


source('process.R')
source('tools.R')

# 
DefaultAssay(data) <- 'RNA'
gene <- intersect(rownames(data),rownames(human)) # 10208
human <- human[gene,]
human <- quickSeurat(human,method='sct',resolution=0.2)
pdf('DimPlot_human_0327.pdf')
DimPlot(human,label = T)
DimPlot(human,label = T,group.by = 'region')
DimPlot(human,label = T,group.by = 'development_stage')
DimPlot(human,label = T,group.by = 'supercluster_term')
dev.off()
pdf('FeaturePlot_human_0327.pdf',10,5)
FeaturePlot(human,features = c('SCL13A6','GAD1'))
dev.off()
# clu0 spatter 
human$se_r0.2 <- human$seurat_clusters
human <- FindClusters(human,resolution = 0.8)
pdf('DimPlot_human_res1.5_0327.pdf')
DimPlot(human,label = T)
dev.off()
meta <- human@meta.data
meta$se_r0.2=as.vector(meta$se_r0.2) # 
meta$selfcluster <- ifelse(meta$seurat_clusters%in%c(1,32)&meta$se_r0.2=='0','0_EX',meta$se_r0.2)
meta$selfcluster <- ifelse(meta$seurat_clusters%in%c(21,6,22,28)&meta$se_r0.2==0,'0_IN',meta$selfcluster)
meta$major <- ifelse(meta$selfcluster%in%c('0','0_IN','15','19','7','8','9','1'),'INN','EXN')
human$major <- meta$major
pdf('FeaturePlot_human_0327.pdf',5,5)
FeaturePlot(human,features = 'SLC16A7')
FeaturePlot(human,features = 'SLC16A6')
FeaturePlot(human,features = 'GAD1')
DimPlot(human,label = T,group.by = 'major')
dev.off()


#### part2 correlation between humans and birds #####
human$region_type <- paste(human$region,human$major,sep='_')
source("/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/tools.R")
human[["RNA"]]@meta.features <- data.frame(row.names = rownames(human[["RNA"]]))
res <- corSeurat(data,human,idents1 = 'celltype',idents2 = 'region_type',plot = F)
pdf('cor_heatmap_OT_humanRegions_snRNA_0328.pdf')
print(corheatmap(r=res$r,p = res$p.adj, cluster_row = F, cluster_col = F))
dev.off()
saveRDS(human,'human_regions_analysed_0328.rds')

#
Idents(human) <- 'region_type'
deg <- FindAllMarkers(human,only.pos = T,logfc.threshold = 0.5)
deg <- subset(deg,deg$p_val_adj<0.05)
write.csv(deg,'human_region_major_deg_0328.csv') # 
tf <- data.table::fread('/data/0.ortholog/TF/TF/Homo_sapiens_TF.txt',header = T)
tf <- tf$Symbol[which(tf$Symbol!='')] # 1637 
gene1 <- intersect(unique(deg$gene),tf) # tf
gene2 <- setdiff(unique(deg$gene),tf) # non-tf
res1 <- corSeurat(data,human,idents1 = 'celltype',idents2 = 'region_type',gene=gene1, plot = F)
res2 <- corSeurat(data,human,idents1 = 'celltype',idents2 = 'region_type',gene=gene2, plot = F)
pdf('cor_heatmap_OT_humanRegions_snRNA_TF_0328.pdf')
print(corheatmap(r=res1$r,p = res1$p.adj, cluster_row = F, cluster_col = F))
print(corheatmap(r=res2$r,p = res2$p.adj, cluster_row = F, cluster_col = F))
dev.off()

#### Part3 deg overlap #####
human <- readRDS('/data/bird_human/human_regions_analysed_0328.rds')
human <- NormalizeData(human)
# deg1 LG_EXN VS. SC_EXN
# deg2 LG_INN VS. SC_INN
Idents(human) <- 'region_type'
deg1 <- FindMarkers(human,ident.1='LGd_EXN',ident.2='SC_EXN',min.pct=0.01) # 
deg2 <- FindMarkers(human,ident.1='LGd_INN',ident.2='SC_INN',min.pct=0.01)

# exn deg1+deg2
bex <- readRDS('./ot_EXN_4species_integrated_0303.rds')
bex$mspecies <- ifelse(bex$species%in%c('pigeon','zebrafinch'),'bird',bex$species)
Idents(bex) <- 'mspecies'
# deg3 OT_EX vs. SC_EX
DefaultAssay(bex) <- 'RNA'
bex <- NormalizeData(bex)
deg3 <- FindMarkers(bex,ident.1='bird',ident.2='human')
# inn
bin <- readRDS('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/ot_part2/04.ot_sn_integ/inn/ot_INN_4species_integrated_0303.rds')
bin$mspecies <- ifelse(bin$species%in%c('pigeon','zebrafinch'),'bird',bin$species)
Idents(bin) <- 'mspecies'
# deg4 OT_IN vs SC_IN
DefaultAssay(bin) <- 'RNA'
bin <- NormalizeData(bin)
deg4 <- FindMarkers(bin,ident.1='bird',ident.2='human')

deg <- list(deg1_lg_ex=deg1,
            deg2_lg_in=deg2,
            deg3_ot_ex=deg3,
            deg4_ot_in=deg4)
saveRDS(deg,'deg_domain_visual_pathway_0408.rds')
library(VennDiagram)
deglist <- lapply(deg, function(x){
  x <- subset(x,x$avg_log2FC > 0.75 & x$p_val_adj < 0.05)
  x <- x[order(x$avg_log2FC,decreasing = T),]
  x <- rownames(x)
  return(x)
})
vplot <- venn.diagram(x=deglist,filename = NULL)
vplot <- venn.diagram(x=deglist[c(1,3)],filename = NULL,fill=RColorBrewer::brewer.pal(3, "Set1")[1:2])
pdf('venn_shared_up_deg_0409.pdf',5,5)
print(grid.draw(vplot))
dev.off()
# EXN
venn_list=deglist[c(1,3)]
venn_list=deglist[c(2,4)]
inter <- get.venn.partitions(venn_list)
#for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
#res=inter[,c(1:3,6,7)]
shared <- unlist(inter[1,'..values..'])
write.table(res,'venn_inter.txt',sep='\t',row.names=FALSE)
# EXN
up_shared <- c('CDKL5',"PPFIBP1", "NR2F1", "GABBR2")
tmp1=subset(human,idents=c('SC_EXN','LGd_EXN'))
tmp2=subset(bex,idents=c('bird','human'))
tmp2$mspecies <- factor(tmp2$mspecies,levels = c('human','bird'))
pdf('vln_exn_shared_updeg_humLG_humSC_0409.pdf',10,4)
VlnPlot(tmp1,features = up_shared,pt.size = F,split.by = 'region_type',stack = T)
VlnPlot(tmp2,features = up_shared,pt.size = F,split.by = 'mspecies',stack = T)
dev.off()

# 0615 visualize DEG from OT_EXN VS. SC_EXN
deg <- readRDS("deg_domain_visual_pathway_0408.rds")
deg=deg$deg3_ot_ex
deg$gene <- rownames(deg)
library(ggplot2)
library(ggrepel)
deg$threshold = factor(ifelse(deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) >= 1.0, ifelse(deg$avg_log2FC>= 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ran<-runif(nrow(deg),min = 25,max =50)
ran<-100^-ran
deg$p_val_adj<-deg$p_val_adj+ran 
pdf('volcano_deg_OtEx_ScEx_0615.pdf',8,6)
ggplot(deg,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#E41A1C","#03396B","#CDDBDB"))+
  geom_text_repel(
    data = deg[deg$p_val<0.05&abs(deg$avg_log2FC)>3.0,],
    aes(label = gene),
    size = 4,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log10 (p_val)')+
  xlab('avg_log2FC')+
  geom_vline(xintercept=c(-1.0,1.0),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
dev.off()

