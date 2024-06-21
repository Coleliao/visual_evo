library(Seurat)
library(dplyr)
setwd("/data/pigeon/deg_vln_OT_LR_new/")
data <- readRDS("/data/pigeon/pigeon_OT_snRNA_integrated_filtered_annotated_0301.rds")
DefaultAssay(data)="RNA"
Idents(data)="celltype"

##
pdf('dimplot_hemisphere_pg_0430.pdf',6,6)
DimPlot(data,group.by = 'hemisphere',pt.size = 0.05,cols = c("#d14a43","#8cd0ab"),order = T)
dev.off()
# to show cell type different
Idents(data) <- 'celltype'
data <- subset(data,idents=c('EXN','INN','OPC','OLI','AST','VLMC'))
genes <- c('EPHA6','ADGRL2','NR2F1','EPHA3','YWHAG','AQP4','SYT1','NEFL')
pdf('vlnplot_deg_0430.pdf',8,6)
VlnPlot(data,features = genes,group.by = 'celltype',
        split.by = 'hemisphere',pt.size = 0,stack = T,cols = c("#d14a43","#8cd0ab"))
dev.off()


#EXN
EXN=subset(data,idents="EXN")
Idents(EXN)="hemisphere"
deg_EXN=FindMarkers(EXN,ident.1="L",min.pct = 0.25)
write.csv(deg_EXN,paste0('deg_FindMarkers_EXN_LR_',Sys.Date(),'.CSV'))
deg_EXN2=deg_EXN %>% subset(deg_EXN$p_val_adj<0.05)
write.csv(deg_EXN2,paste0('deg_FindMarkers_EXN_LR_filter_',Sys.Date(),'.CSV'))

#INN
INN=subset(data,idents="INN")
Idents(INN)="hemisphere"
deg_INN=FindMarkers(INN,ident.1="L",min.pct = 0.25)
write.csv(deg_INN,paste0('deg_FindMarkers_INN_LR_',Sys.Date(),'.CSV'))
deg_INN2=deg_INN %>% subset(deg_INN$p_val_adj<0.05)
write.csv(deg_INN2,paste0('deg_FindMarkers_INN_LR_filter_',Sys.Date(),'.CSV'))


##OLI
OLI=subset(data,idents="OLI")
Idents(OLI)="hemisphere"
deg_OLI=FindMarkers(OLI,ident.1="L",min.pct = 0.25)
write.csv(deg_OLI,paste0('deg_FindMarkers_OLI_LR_',Sys.Date(),'.CSV'))
deg_OLI2=deg_OLI %>% subset(deg_OLI$p_val_adj<0.05)
write.csv(deg_OLI2,paste0('deg_FindMarkers_OLI_LR_filter_',Sys.Date(),'.CSV'))

##AST
AST=subset(data,idents="AST")
Idents(AST)="hemisphere"
deg_AST=FindMarkers(AST,ident.1="L",min.pct = 0.25)
write.csv(deg_AST,paste0('deg_FindMarkers_AST_LR_',Sys.Date(),'.CSV'))
deg_AST2=deg_AST %>% subset(deg_AST$p_val_adj<0.05)
write.csv(deg_AST2,paste0('deg_FindMarkers_AST_LR_filter_',Sys.Date(),'.CSV'))

#OPC  
OPC=subset(data,idents="OPC")
Idents(OPC)="hemisphere"
deg_OPC=FindMarkers(OPC,ident.1="L",min.pct = 0.25)
write.csv(deg_OPC,paste0('deg_FindMarkers_OPC_LR_',Sys.Date(),'.CSV'))
deg_OPC2=deg_OPC%>% subset(deg_OPC$p_val_adj<0.05)
write.csv(deg_OPC2,paste0('deg_FindMarkers_OPC_LR_filter_',Sys.Date(),'.CSV'))


#VLMC
VLMC=subset(data,idents="VLMC")
Idents(VLMC)="hemisphere"
deg_VLMC=FindMarkers(VLMC,ident.1="L",min.pct = 0.25)
write.csv(deg_VLMC,paste0('deg_FindMarkers_VLMC_LR_',Sys.Date(),'.CSV'))
deg_VLMC2=deg_VLMC%>% subset(deg_VLMC$p_val_adj<0.05)
write.csv(deg_VLMC2,paste0('deg_FindMarkers_VLMC_LR_filter_',Sys.Date(),'.CSV'))

library(Seurat)
library(reshape2)
library(RColorBrewer)
EXN<- read.csv('./deg_FindMarkers_EXN_LR_filter_2024-04-28.CSV',row.names = 1)
INN<- read.csv('./deg_FindMarkers_INN_LR_filter_2024-04-28.CSV',row.names = 1)
AST<- read.csv('./deg_FindMarkers_AST_LR_filter_2024-04-28.CSV',row.names = 1)
OLI<- read.csv('./deg_FindMarkers_OLI_LR_filter_2024-04-28.CSV',row.names = 1)
OPC<- read.csv('./deg_FindMarkers_OPC_LR_filter_2024-04-28.CSV',row.names = 1)
VLMC<- read.csv('./deg_FindMarkers_VLMC_LR_filter_2024-04-28.CSV',row.names = 1)

EXN1<- subset(EXN,EXN$avg_log2FC>0.25)
INN1<- subset(INN,INN$avg_log2FC>0.25)
AST1<- subset(AST,AST$avg_log2FC>0.25)
OLI1<- subset(OLI,OLI$avg_log2FC>0.25)
OPC1<- subset(OPC,OPC$avg_log2FC>0.25)
VLMC1<- subset(VLMC,VLMC$avg_log2FC>0.25)

EXN1$degree="L_up"
INN1$degree="L_up"
AST1$degree="L_up"
OLI1$degree="L_up"
OPC1$degree="L_up"
VLMC1$degree="L_up"
data=list(EXN=EXN1,INN=INN1,AST=AST1,OLI=OLI1,OPC=OPC1,VLMC=VLMC1)
table_list <- lapply(data, function(x) {
  table(x[,"degree"])  
})
mydata <- do.call(cbind, table_list)
mydata=as.data.frame(melt(mydata))
mydata$value <- as.factor(mydata$value)
myColors <- brewer.pal(6, "Oranges")
pdf("DEG_num_Lup.pdf")
ggplot(mydata) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_manual(values = myColors) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_fixed()
dev.off()

EXN2<- subset(EXN,EXN$avg_log2FC<(-0.25))
INN2<- subset(INN,INN$avg_log2FC<(-0.25))
AST2<- subset(AST,AST$avg_log2FC<(-0.25))
OLI2<- subset(OLI,OLI$avg_log2FC<(-0.25))
OPC2<- subset(OPC,OPC$avg_log2FC<(-0.25))
VLMC2<- subset(VLMC,VLMC$avg_log2FC<(-0.25))

EXN2$degree="R_up"
INN2$degree="R_up"
AST2$degree="R_up"
OLI2$degree="R_up"
OPC2$degree="R_up"
VLMC2$degree="R_up"
data=list(EXN=EXN2,INN=INN2,AST=AST2,OLI=OLI2,OPC=OPC2,VLMC=VLMC2)
table_list <- lapply(data, function(x) {
  table(x[,"degree"])  
})
mydata <- do.call(cbind, table_list)
mydata=as.data.frame(melt(mydata))
mydata$value <- as.factor(mydata$value)
myColors <- brewer.pal(6, "Blues")
pdf("DEG_num_up.pdf")
ggplot(mydata) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_manual(values = myColors) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_fixed()
dev.off()

#deg GO
Idents(data)="celltype"
sub=subset(data,idents=c("EXN","INN","AST"))
sub=NormalizeData(sub)
Idents(sub)='hemisphere'
deg=FindMarkers(sub,ident.1="L",ident.2="R",min.pct = 0.25)
write.csv(deg,paste0('deg_FindMarkers_EXNINNAST_LR_',Sys.Date(),'.CSV'))

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
data<- read.csv('./deg_FindMarkers_EXNINNAST_LR_2024-06-13.CSV',row.names = 1)
data=subset(data,data$p_val_adj<0.05)
data$gene=rownames(data)
ids=bitr(data$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ids <- na.omit(ids)
data=merge(data,ids,by.x='gene',by.y='SYMBOL')
data2 <- subset(data,abs(data$avg_log2FC)>0.25)
enrich.go.BP = enrichGO(gene = data2$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)
pdf(file="enrichGO_OT_LR_EIA_0.25_new.pdf", width=16, height=9)
p <- dotplot(enrich.go.BP, showCategory=20, font.size=20, title="GO of OT_LR_EIA_0.25") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  geom_point(size=6)
print(p)
dev.off()
