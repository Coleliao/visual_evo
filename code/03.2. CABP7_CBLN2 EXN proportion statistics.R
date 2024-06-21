library(ggplot2)
library(ggpubr)
library(RColorBrewer)

data=readRDS('./ot_EXN_4species_integrated_0303.rds')
meta=data@meta.data
meta$CABP7_count=data@assays$RNA@counts['CABP7',]
meta$CABP7_data=data@assays$RNA@data['CABP7',]
meta$CBLN2_count=data@assays$RNA@counts['CBLN2',]
meta$CBLN2_data=data@assays$RNA@data['CBLN2',]
table(meta$sample)

# average expression level 
m <- meta %>% group_by(sample) %>% 
  summarise(mean1=mean(CABP7_data))
n <- meta %>% group_by(sample) %>% 
  summarise(mean1=mean(CBLN2_data))
m <- tmp$mean_CAPB7_data
n <- tmp$mean_CBLN2_data

# relative cellnumber of CABP7+ cells and CBLN2+ cells
meta1 <- subset(meta,meta$CABP7_count>2)
meta2 <- subset(meta,meta$CBLN2_count>2)
a <- table(meta1$sample)/table(meta$sample)
b <- table(meta2$sample)/table(meta$sample)

res <- list(mean_CABP7_data=m,mean_CBLN2_data=n, # expression
            CABP7_cellnumber=a,CBLN2_cellnumber=b)
saveRDS(res,'4species_CABP7_CBLN2_res_0426.rds')

res <- readRDS('D:/下载/4species_CABP7_CBLN2_res_0426.rds')
res <- data.frame(mean_CABP7_data=res$mean_CABP7_data$mean1,mean_CBLN2_data=res$mean_CBLN2_data$mean1, # expression
                  CABP7_cellnumber=as.vector(res$CABP7_cellnumber),CBLN2_cellnumber=as.vector(res$CBLN2_cellnumber))
res$mspecies <- c(rep('human',3),rep('mouse',5),rep('bird',4))

# visualization and  significance test
# cellnumber proportion (major fig)
res$CABP7_cellnumber <- 100*res$CABP7_cellnumber # 100%
res$CBLN2_cellnumber <- 100*res$CBLN2_cellnumber
pdf('./CABP7_CBLN2_4species_cell_proportion_0426.pdf',4,8)
ggboxplot(data = res, x = "mspecies", y = "CABP7_cellnumber",
          fill = "mspecies", palette =col ,width = 0.5)+ stat_compare_means(label.y = 5) 
ggboxplot(data = res, x = "mspecies", y = "CBLN2_cellnumber",
          fill = "mspecies", palette =col ,width = 0.5)+ stat_compare_means(label.y = 15) 
dev.off()
