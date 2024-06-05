# Figure 3

library(Seurat)
library(Matrix)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(magrittr)
library(dplyr)
library(stringr)
library(scales)
library(ggalluvial)
library(gridExtra)
library(grid)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggalluvial)
library(ggpubr)
library(scCancer)

clu.colors <-c( '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

setwd('./Analysis/Data/')
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

flow.meta <- read.csv('Flow_information.csv', header = TRUE, check.names = F)
cluster.meta <- read.csv('Control_cluster_information.csv', header = TRUE, check.names = F)
cluster.meta$Patient <- substr(cluster.meta$`Control Cluster`, 1, 2)

smp <- readRDS('Combine_PDAC_Seurat_Object.RDS')
table(smp$Clusters)

#----------------------------------------------------------------------------------------#
# Figure 3A, Drug response identification for control clusters. 
options(repr.plot.width=5.5, repr.plot.height=7)
cluster.meta$`Drug response type` <- factor(cluster.meta$`Drug response type` , 
                                            levels = c( 'Absolutely Intrinsic Resistance', 'Relatively Intrinsic Resistance', 'Sensitivity'), 
                                            ordered = T)

p1 <- ggplot(cluster.meta, aes(x=log2(`RDR score`), y=log2(`ADR score`), color=`Drug response type`)) +
  geom_point(size = 4) + 
  scale_color_manual(values  = c('#EE2C2C', "#FFC300", "#1ABC9C"))+
  xlab('log2(RDR score)') +ylab('log2(ADR score)')+labs(colour="Drug response")+
  theme_classic()+
  theme(
    legend.key.height=unit(1.8, "line"),  legend.position = 'bottom', legend.direction = "vertical",
    legend.key = element_rect(fill='transparent', color='transparent'), 
    legend.key.size = unit(0.75, "cm"),
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 22), 
    plot.title = element_text(size = 16, vjust = 0.5), 
    panel.background = element_blank(),
    strip.text = element_text(size=16),
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20, angle = 90, vjust= 1.5), 
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(colour="red", size=.5, linetype=1))+
  geom_text_repel(size = 5,  label= cluster.meta$`Control Cluster`, 
                  box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.2, "lines"), 
                  max.overlaps = 20)+
  geom_vline(aes(xintercept = log2(0.9)), colour = 'gray55', linetype="longdash", alpha = 0.6)+
  geom_hline(aes(yintercept = log2(0.9)), colour = 'gray55', linetype="longdash", alpha = 0.6)

print(p1)

pdf('Cluster_Intrinsic_sensitivity.pdf', width = 5.5, height = 7)
print(p1)
dev.off()


#----------------------------------------------------------------------------------------#
# Figure 3B, Venn diagrams showing the consistently upregulated genes 
# when comparing AIR clusters with RIR and sensitivity (SEN) clusters across three patients. 

smp.sub <- subset(smp, subset = Clusters %in% cluster.meta$`Control Cluster`)
smp.sub$DR_Type <- cluster.meta$`Drug response type`[match(smp.sub$Clusters, cluster.meta$`Control Cluster`)]
table(smp.sub$DR_Type)

smp.sub$DR_Type_short <- as.character(smp.sub$DR_Type)
smp.sub$DR_Type_short[smp.sub$DR_Type_short  == 'Absolutely Intrinsic Resistance'] <- 'AIR'
smp.sub$DR_Type_short[smp.sub$DR_Type_short  == 'Relatively Intrinsic Resistance'] <- 'RIR'
smp.sub$DR_Type_short[smp.sub$DR_Type_short  == 'Sensitivity'] <- 'SEN'

smp.sub$Sample_DR_Type <-  paste0(smp.sub$patient, '_', smp.sub$DR_Type_short)
table(smp.sub$DR_Type_short)
table(smp.sub$Sample_DR_Type)


# Differentially expressed genes of AIR compared with RIR
A_vs_R_results <- data.frame(Gene = c(), wilcox.test = c(), wilcox.test_BH = c(), avg_logFC = c(), pct.1 = c(), pct.2 = c(), Type = c())
for (pa in c('P1', 'P2', 'P4')){
  print(pa)
  air.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Absolutely Intrinsic Resistance' & cluster.meta$Patient == pa)]
  air.clus <- as.character(air.clus)
  print(paste0('Absolute IR: ', air.clus))
  
  rir.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Relatively Intrinsic Resistance' & cluster.meta$Patient == pa)]
  rir.clus <- as.character(rir.clus)
  print(paste0('Relative IR: ', paste(rir.clus, collapse = ', ')))
  
  sensi.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Sensitivity' & cluster.meta$Patient == pa)]
  sensi.clus <- as.character(sensi.clus)
  print(paste0('Sensitivity: ', sensi.clus))
  
  compare.1 <- air.clus
  compare.2 <- rir.clus
  
  expr <- GetAssayData(object = smp.sub, slot = "data")
  cells <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  
  cells.1 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1))]
  cells.2 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.2))]
  
  clusters <- smp.sub$Clusters[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  group <- clusters
  group[which(group %in% compare.1)] <- 'Absolutely Intrinsic Resistance'
  group[which(group %in% compare.2)] <- 'Relatively Intrinsic Resistance'
  print(table(group))
  
  pct.1 <-  apply(expr[, cells.1, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  pct.2 <-  apply(expr[, cells.2, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  genes_use <- names(pct.1)[which(pct.1 > 0.1)]
  
  wicox_test_results <- data.frame(Gene =  genes_use)
  wicox_test_results$wilcox.test <- apply(expr[genes_use, cells], 1,
                                          function(x) unlist(wilcox.test(as.numeric(x) ~ group, exact=FALSE)[3]))
  wicox_test_results$wilcox.test_BH <- p.adjust(wicox_test_results$wilcox.test, method = "BH")
  
  bulk.1 <- expr[genes_use, cells.1]
  bulk.2 <- expr[genes_use, cells.2]
  avg_logFC <- rowMeans(bulk.1) - rowMeans(bulk.2)
  wicox_test_results$avg_logFC <- avg_logFC
  wicox_test_results$pct.1 <- pct.1[wicox_test_results$Gene]
  wicox_test_results$pct.2 <- pct.2[wicox_test_results$Gene]
  rownames(wicox_test_results) <- wicox_test_results$Gene
  
  wicox_test_results$Type <- paste0(pa, '_AIR_vs_RIR')
  A_vs_R_results <- rbind(A_vs_R_results,  wicox_test_results)                                        
}


# Differentially expressed genes of AIR compared with SEN
A_vs_S_results <- data.frame(Gene = c(), wilcox.test = c(), wilcox.test_BH = c(), 
                             avg_logFC = c(), pct.1 = c(), pct.2 = c(), Type = c())
for (pa in c('P1', 'P2', 'P4')){
  print(pa)
  air.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Absolutely Intrinsic Resistance' & cluster.meta$Patient == pa)]
  air.clus <- as.character(air.clus)
  print(paste0('Absolute IR: ', air.clus))
  
  rir.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Relatively Intrinsic Resistance' & cluster.meta$Patient == pa)]
  rir.clus <- as.character(rir.clus)
  print(paste0('Relative IR: ', paste(rir.clus, collapse = ', ')))
  
  sensi.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Sensitivity' & cluster.meta$Patient == pa)]
  sensi.clus <- as.character(sensi.clus)
  print(paste0('Sensitivity: ', sensi.clus))
  
  
  compare.1 <- air.clus
  compare.2 <- sensi.clus
  
  expr <- GetAssayData(object = smp.sub, slot = "data")
  cells <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  
  cells.1 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1))]
  cells.2 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.2))]
  
  clusters <- smp.sub$Clusters[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  group <- clusters
  group[which(group %in% compare.1)] <- 'Absolutely Intrinsic Resistance'
  group[which(group %in% compare.2)] <- 'Sensitivity'
  print(table(group))
  
  pct.1 <-  apply(expr[, cells.1, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  pct.2 <-  apply(expr[, cells.2, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  genes_use <- names(pct.1)[which(pct.1 > 0.1)]
  
  wicox_test_results <- data.frame(Gene =  genes_use)
  wicox_test_results$wilcox.test <- apply(expr[genes_use, cells], 1,
                                          function(x) unlist(wilcox.test(as.numeric(x) ~ group, exact=FALSE)[3]))
  wicox_test_results$wilcox.test_BH <- p.adjust(wicox_test_results$wilcox.test, method = "BH")
  
  bulk.1 <- expr[genes_use, cells.1]
  bulk.2 <- expr[genes_use, cells.2]
  avg_logFC <- rowMeans(bulk.1) - rowMeans(bulk.2)
  wicox_test_results$avg_logFC <- avg_logFC
  wicox_test_results$pct.1 <- pct.1[wicox_test_results$Gene]
  wicox_test_results$pct.2 <- pct.2[wicox_test_results$Gene]
  rownames(wicox_test_results) <- wicox_test_results$Gene
  
  wicox_test_results$Type <- paste0(pa, '_AIR_vs_SEN')
  A_vs_S_results <- rbind(A_vs_S_results,  wicox_test_results)                                        
}


# Differentially expressed genes of RIR compared with SEN
R_vs_S_results <- data.frame(Gene = c(), wilcox.test = c(), wilcox.test_BH = c(), 
                             avg_logFC = c(), pct.1 = c(), pct.2 = c(), Type = c())
for (pa in paste0('P', 1:6)){
  print(pa)
  air.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Absolutely Intrinsic Resistance' & cluster.meta$Patient == pa)]
  air.clus <- as.character(air.clus)
  print(paste0('Absolute IR: ', air.clus))
  
  rir.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Relatively Intrinsic Resistance' & cluster.meta$Patient == pa)]
  rir.clus <- as.character(rir.clus)
  print(paste0('Relative IR: ', paste(rir.clus, collapse = ', ')))
  
  sensi.clus <- cluster.meta$`Control Cluster`[which(cluster.meta$`Drug response type` == 'Sensitivity' & cluster.meta$Patient == pa)]
  sensi.clus <- as.character(sensi.clus)
  print(paste0('Sensitivity: ', sensi.clus))
  
  
  compare.1 <- rir.clus
  compare.2 <- sensi.clus
  
  expr <- GetAssayData(object = smp.sub, slot = "data")
  cells <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  
  cells.1 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.1))]
  cells.2 <- colnames(smp.sub)[which(smp.sub$Clusters %in% c(compare.2))]
  
  clusters <- smp.sub$Clusters[which(smp.sub$Clusters %in% c(compare.1, compare.2))]
  group <- clusters
  group[which(group %in% compare.1)] <- 'Relatively Intrinsic Resistance'
  group[which(group %in% compare.2)] <- 'Sensitivity'
  print(table(group))
  
  pct.1 <-  apply(expr[, cells.1, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  pct.2 <-  apply(expr[, cells.2, drop = F],
                  MARGIN = 1,
                  FUN = function(x) {
                    return(sum(x > 0) / length(x = x))
                  }
  )
  genes_use <- names(pct.1)[which(pct.1 > 0.1)]
  
  wicox_test_results <- data.frame(Gene =  genes_use)
  wicox_test_results$wilcox.test <- apply(expr[genes_use, cells], 1,
                                          function(x) unlist(wilcox.test(as.numeric(x) ~ group, exact=FALSE)[3]))
  wicox_test_results$wilcox.test_BH <- p.adjust(wicox_test_results$wilcox.test, method = "BH")
  
  bulk.1 <- expr[genes_use, cells.1]
  bulk.2 <- expr[genes_use, cells.2]
  avg_logFC <- rowMeans(bulk.1) - rowMeans(bulk.2)
  wicox_test_results$avg_logFC <- avg_logFC
  wicox_test_results$pct.1 <- pct.1[wicox_test_results$Gene]
  wicox_test_results$pct.2 <- pct.2[wicox_test_results$Gene]
  rownames(wicox_test_results) <- wicox_test_results$Gene
  
  wicox_test_results$Type <- paste0(pa, '_RIR_vs_SEN')
  R_vs_S_results <- rbind(R_vs_S_results,  wicox_test_results)                                        
}


markers <- rbind(A_vs_R_results, R_vs_S_results, A_vs_S_results)
table(markers$Type)
#write.csv(markers, 'Intrinsic_resistance_compare_3groups_wilcox_test.csv')


# Venn plots
markers.all <- read.csv( '/home/glab/lxq/PDAC_writing/github/Mateiral/Intrinsic_resistance_compare_3groups_wilcox_test.csv', 
                         row.names = 1, check.names=F)
markers.all$Patient <- substr(markers.all$Type, 1, 2)
markers.all$Group <- substr(markers.all$Type, 4, 15 )


thre <- 0
markers <- markers.all

genes.P1.AvsR <- markers$Gene[which(markers$Patient == 'P1' &  markers$Group == 'AIR_vs_RIR' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P1.AvsS <- markers$Gene[which(markers$Patient == 'P1' &  markers$Group == 'AIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P1.RvsS <- markers$Gene[which(markers$Patient == 'P1' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]

genes.P2.AvsR <- markers$Gene[which(markers$Patient == 'P2' &  markers$Group == 'AIR_vs_RIR' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P2.AvsS <- markers$Gene[which(markers$Patient == 'P2' &  markers$Group == 'AIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P2.RvsS <- markers$Gene[which(markers$Patient == 'P2' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]

genes.P4.AvsR <- markers$Gene[which(markers$Patient == 'P4' &  markers$Group == 'AIR_vs_RIR' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P4.AvsS <- markers$Gene[which(markers$Patient == 'P4' &  markers$Group == 'AIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P4.RvsS <- markers$Gene[which(markers$Patient == 'P4' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]

genes.P3.RvsS <- markers$Gene[which(markers$Patient == 'P3' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P5.RvsS <- markers$Gene[which(markers$Patient == 'P5' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]
genes.P6.RvsS <- markers$Gene[which(markers$Patient == 'P6' &  markers$Group == 'RIR_vs_SEN' &
                                      markers$avg_logFC > thre & markers$wilcox.test_BH < 0.05)]

library(venn) 
library(tidyverse)
library(stringr)

set.seed(2023)

venn.data = 
  list(P1 = genes.P1.AvsR, 
       P2 = genes.P2.AvsR, 
       P4 = genes.P4.AvsR)
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')
options(repr.plot.width=2.5, repr.plot.height=2.5)

pdf('vennplot_AIRvsRIR.pdf', height = 2.5, width = 2.5)
venn.result =
  venn(venn.data, ilabels = TRUE, 
       zcolor =colors_use[c(1, 2, 3)] , opacity =  0.3, size = 30, cexil = 1.2, cexsn = 1.4, box= F, borders = T, sncs = 1.22, ilcs = 1.03);
dev.off()

venn.result =
  venn(venn.data, ilabels = TRUE, 
       zcolor =colors_use[c(1, 2, 3)] , opacity =  0.3, size = 30, cexil = 1.2, cexsn = 1.4, box= F, borders = T, sncs = 1.22, ilcs = 1.03);



venn.data = 
  list(P1 = genes.P1.AvsS, 
       P2 = genes.P2.AvsS, 
       P4 = genes.P4.AvsS)
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

venn.result =
  venn(venn.data, ilabels = TRUE, 
       zcolor =colors_use[c(1, 2, 3)] , opacity =  0.3, size = 30, cexil = 1.2, cexsn = 1.4, box= F, borders = T, sncs = 1.22, ilcs = 1.03);


pdf('vennplot_AIRvsSEN.pdf', height = 2.5, width = 2.5)
venn.result =
  venn(venn.data, ilabels = TRUE, 
       zcolor =colors_use[c(1, 2, 3)] , opacity =  0.3, size = 30, cexil = 1.2, cexsn = 1.4, box= F, borders = T, sncs = 1.22, ilcs = 1.03);
dev.off()


#----------------------------------------------------------------------------------------#
# Figure 3C, Violin plots showing typical upregulated genes associated with stemness (CD44, CTSB, and EPCAM) and ER stress (ERP29, ERGIC3, and PRKCSH).

smp.sub$Sample_DR_Type <- factor(smp.sub$Sample_DR_Type, levels = c(
  'P1_AIR', 'P1_RIR',  'P1_SEN',
  'P2_AIR', 'P2_RIR',  'P2_SEN',
  'P4_AIR', 'P4_RIR',  'P4_SEN',
  'P3_RIR',  'P3_SEN',
  'P5_RIR',  'P5_SEN',
  'P6_RIR',  'P6_SEN'), ordered = T)

groups <- names(table(smp.sub$Sample_DR_Type))

my_comparisons_height <- c(0, 1, 2, 0, 1, 2, 0,  1, 2,  0, 0, 0) 


smp.sub2 <- subset(smp.sub, subset = orig.ident %in% c('P1-Ctrl', 'P2-Ctrl', 'P4-Ctrl'))
my_comparisons <- list(P11 = c(groups[1], groups[2]), 
                       P13 = c(groups[2], groups[3]), 
                       P12 = c(groups[1], groups[3]), 
                       
                       P21 = c(groups[4], groups[5]), 
                       P213 = c(groups[5], groups[6]), 
                       P22 = c(groups[4], groups[6]), 
                       
                       P41 = c(groups[7], groups[8]), 
                       P43 = c(groups[8], groups[9]), 
                       P42 = c(groups[7], groups[9])
)

Genes = c('CD44', 'CTSB', 'EPCAM') # Figure 3C
#Genes = c('ERP29', 'ERGIC3', 'PRKCSH') # Figure 3C
#Genes = c('CD9', 'TMBIM6', 'CLU') # Supp Figure 5B
#Genes = c('PDIA6', 'P4HB', 'GSTM3') # Supp Figure 5B

options(repr.plot.width=5.6, repr.plot.height=6.5)

label_y= 4.2 + (my_comparisons_height) *0.8

p1 <- VlnPlot(smp.sub2, features=Genes, group.by = 'Sample_DR_Type', split.by = 'patient', slot = 'data', 
              cols = colors_use, pt.size = 1, stack=T, flip=T)+

  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test',p.adjust.method="BH",label.x=1,label.y=label_y,  size =4.8)+
  NoLegend()+
  scale_y_continuous(limits = c(-0.1,6.4))+
  theme(
    axis.text.x = element_text(colour = "black",size=14, angle = 30, hjust = 1, vjust = 1), 
    axis.text.y = element_text(colour = "black",size=14), 
    axis.title.y = element_text(colour = "black",size=16), 
    strip.text = element_text(colour = "black",size=16, face = 'plain', hjust = 0.1), 
    axis.line = element_line(colour = "black",size = 0.5)  )+
  xlab('')
p1

pdf('Intrinsic Resistance_vlnplot_genes_1.pdf', height = 6.5, width = 5.6)
print(p1)
dev.off()


#----------------------------------------------------------------------------------------#
# Figure 3D-E, Violin plots showing scores for three signatures in AIR, RIR, and sensitivity clusters (left).
# Kaplan-Meier curves with log-rank tests show the relationship between signature expressions and patient overall survival within the TCGA PAAD cohort. 
AddModuleScore_bulk <- function (object, features, nbin = 24, ctrl = 100, 
                                 k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = 1)
  
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = rownames(x = object))
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ", 
              paste(missing.features, collapse = ", "), 
              ifelse(test = FALSE, yes = ", attempting to find updated synonyms", 
                     no = ", not searching for symbol synonyms"), 
              call. = FALSE, immediate. = TRUE)
    }
    return(intersect(x = x, y = rownames(x = object)))
  })
  
  cluster.length <- length(x = features)
  
  pool <- rownames(x = object)
  data.avg <- apply(object,1, mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- apply(object[features.use,], 2, mean)
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- apply(data.use, 2, mean)
  }
  features.scores.use <- features.scores - ctrl.scores
  
  colnames(features.scores.use) <- colnames(object)
  rownames(features.scores.use) <- name
  return(features.scores.use)
}

# read TCGA data
TCGA_survival.df <- read.csv( 'TCGA_PAAD_Survival_info.csv', header = T, row.names = 1, check.names = F)
TCGA_pheno.df <- read.csv( 'TCGA_PAAD_Phenotype_info.csv', header = T, row.names = 1, check.names = F)
TCGA_expr.mat <- read.csv( 'TCGA_PAAD_FPKM.csv', header = T, row.names = 1, check.names = F)

TCGA_pheno.df <- TCGA_pheno.df[which(TCGA_pheno.df$primary_diagnosis.diagnoses == 'Infiltrating duct carcinoma, NOS' &
                                       TCGA_pheno.df$sample_type.samples == 'Primary Tumor'), ]
TCGA_survival.df <- TCGA_survival.df[TCGA_pheno.df$submitter_id.samples, ]
TCGA_expr.mat <- TCGA_expr.mat[, TCGA_pheno.df$submitter_id.samples]
colnames(TCGA_survival.df) <- c('sample', 'OS_STATUS', 'Patient', 'OS_DAYS')
TCGA_survival.df$OS_MONTHS <- round(TCGA_survival.df$OS_DAYS / 30, 4)


# Figure 3d
pathway.sel <- 'Stemness'

# downloaded form http://dibresources.jcbose.ac.in/ssaha4/bcscdb/search.php
# https://academic.oup.com/database/article/doi/10.1093/database/baac082/6725752
stemness.df <- read.csv('/home/glab/lxq/reference/CSC_Biomarker_2022_All.csv', check.names = F, header = F)


colnames(stemness.df) <- c('Gene', 'TM', 'Direction', 'HGNC_ID', 'Disease', 'Disease_Type', 'Cell_Line', 'Drug', 'Method',
                           'CONFIDENCE_SCORING', 'GLOBAL_SCORING', 'PUBMED_ID')
pdac.stemness.df <- stemness.df[which(stemness.df$Disease == 'Pancreatic Cancer'), ]
dim(pdac.stemness.df)
pdac.stemness.df$Disease_Type[is.na(pdac.stemness.df$Disease_Type)] <- 'Undetected'
table(pdac.stemness.df$Disease_Type)

sel1 <- pdac.stemness.df[which(pdac.stemness.df$Disease == 'Pancreatic Cancer' & 
                                 (pdac.stemness.df$Disease_Type %in% c('Ductal Adenocarcinoma', 'Ductal Carcinoma', 'Carcinoma')) & 
                                 pdac.stemness.df$Cell_Line == 'Primary'), ]

stemness.genes1 <- intersect( unique(sel1$Gene), rownames(smp))
print(paste(stemness.genes1, collapse = ', '))
print(length(stemness.genes1))


sel.df <- pdac.stemness.df[which(pdac.stemness.df$Disease == 'Pancreatic Cancer' & 
                                   (pdac.stemness.df$Disease_Type %in% c('Ductal Adenocarcinoma', 'Ductal Carcinoma', 'Carcinoma')) & 
                                   pdac.stemness.df$Direction == 'Up' &
                                   pdac.stemness.df$Cell_Line != 'Primary'), ]

table(sel.df$Disease_Type)
table(sel.df$Cell_Line)
genes.all <- unique(sel.df$Gene)
stemness.genes2 <- c()
for (gene in genes.all){
  tmp.df <- sel.df[which(sel.df$Gene == gene), ]
  origin <- length(unique(tmp.df$PUBMED_ID))
  cellline <- length(unique(tmp.df$Cell_Line))
  if (origin >= 1 & cellline > 1){
    stemness.genes2 <- c(stemness.genes2, gene)
  }
}

stemness.genes2 <- intersect( stemness.genes2, rownames(smp))
print(paste(stemness.genes2, collapse = ', '))
print(length(stemness.genes2))


print(length(unique(c(stemness.genes1, stemness.genes2))))
print(paste(intersect(rownames(smp), unique(c(stemness.genes1, stemness.genes2))), collapse = ', '))


smp.sub$Sample_DR_Type <- factor(smp.sub$Sample_DR_Type, levels = c(
  'P1_AIR', 'P1_RIR',  'P1_SEN',
  'P2_AIR', 'P2_RIR',  'P2_SEN',
  'P4_AIR', 'P4_RIR',  'P4_SEN',
  'P3_RIR',  'P3_SEN',
  'P5_RIR',  'P5_SEN',
  'P6_RIR',  'P6_SEN'), ordered = T)

groups <- names(table(smp.sub$Sample_DR_Type))
my_comparisons2 <- list(P11 = c(groups[1], groups[2]), 
                        P13 = c(groups[2], groups[3]), 
                        P12 = c(groups[1], groups[3]), 
                        
                        P21 = c(groups[4], groups[5]), 
                        P213 = c(groups[5], groups[6]), 
                        P22 = c(groups[4], groups[6]), 
                        
                        P41 = c(groups[7], groups[8]), 
                        P43 = c(groups[8], groups[9]), 
                        P42 = c(groups[7], groups[9]), 
                        
                        P3 = c(groups[10], groups[11]), 
                        P5 = c(groups[12], groups[13]), 
                        P6 = c(groups[14], groups[15]))


my_comparisons_height2 <- c(0, 1, 2, 0, 1, 2, 0,  1, 2,  0, 0, 0) 



smp.sub <- AddModuleScore(smp.sub, features =list(unique(c(stemness.genes1, stemness.genes2))), name = 'Pathway' )


my_comparisons_height <- c(0, 1, 2, 0, 1, 2, 0,  1, 2,  0, 1, 0) 

aa <- max(smp.sub$Pathway1)
bb <- min(smp.sub$Pathway1)

label_y= aa + (my_comparisons_height) *0.15
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')



options(repr.plot.width=9, repr.plot.height=4.8)

p2 <- VlnPlot(smp.sub, 'Pathway1', group.by = 'Sample_DR_Type', split.by = 'patient', slot = 'data', cols = colors_use[c(1, 2, 4, 3, 5, 6)])+
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "solid",color="red", size=1)+
  NoLegend()+
  stat_compare_means(comparisons = my_comparisons2,method = 'wilcox.test',p.adjust.method="BH",label.x=1,label.y=label_y, size = 5.5)+
  ylim(c(bb, max(label_y)*1.08))+
  theme(legend.position="none", 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18))+
  ggtitle('')+xlab("")+ylab("Signature Score")

print(p2)

# pdf('vlnplot_stemness.pdf', height = 5, width = 9)
# print(p2)
# dev.off()


# Survival analysis

scores <- AddModuleScore_bulk(TCGA_expr.mat, features= list( unique(c(stemness.genes1, stemness.genes2))))
scores <- as.vector(scores)
names(scores) <- colnames(TCGA_expr.mat)

patients <- intersect(rownames(TCGA_survival.df),colnames(TCGA_expr.mat))  
length(patients)
tcga.df <- cbind(TCGA_survival.df[patients, c('OS_STATUS','OS_MONTHS')], scores[patients])
tcga.df$OS_STATUS <- as.character(tcga.df$OS_STATUS) 
tcga.df$OS_STATUS <- as.numeric(substr(tcga.df$OS_STATUS, 1, 1))


colnames(tcga.df)[3] <- 'Value'  
tcga.df$group <- "Median"
tcga.df$group[tcga.df$Value > quantile(tcga.df$Value, 0.5)] <- "High"
tcga.df$group[tcga.df$Value <= quantile(tcga.df$Value, 0.5)] <- "Low"
tcga.df$group <- factor(tcga.df$group, levels = c('High', 'Low'), ordered = T)

surv_object <- Surv(time = tcga.df$OS_MONTHS, event = tcga.df$OS_STATUS)
fit <- survfit(surv_object ~ group, data =  tcga.df)

options(repr.plot.width=3.5, repr.plot.height=4.3)

p.surv <- ggsurvplot(fit, pval = T,
                     legend.labs=c("High" ,  "Low"),
                     palette = c("#DE77AE", "#7FBC41"),
                     # axes.offset = F,
                     #test.for.trend	= T,
                     #legend = "right",
                     ggtheme = theme_bw(),
                     legend.title = "Group", 
                     
                     conf.int = F,
                     # Add risk table
                     risk.table = TRUE,
                     tables.height = 0.25,
                     #tables.y.text = F, 
                     tables.theme = theme_cleantable())
p.surv$plot <- p.surv$plot + 
  xlab('Months')+
  theme(
    #       legend.key.height=unit(1.8, "line"),
    legend.position=c(0.8, 0.78),
    legend.background=element_blank(),
    legend.spacing.y = unit(0.3, "cm"), 
    plot.title = element_text(size = 14, vjust = 0.5), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    strip.text = element_text(size=14),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1))+ guides(fill=guide_legend(
        keywidth=0.5,
        keyheight=0.5,
        default.unit="inch")
      )

print(p.surv)

#----------------------------------------------------------------------------#
# Figure 3e
gmt <- read.gmt("/home/glab/lxq/reference/c5.go.bp.v2022.1.Hs.symbols.gmt")

pathway.sel <- 'GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY'
genes.path <- gmt$gene[which(gmt$term == pathway.sel)]
length(genes.path)

genes_use <- intersect(genes.path, rownames(smp))
length(genes_use)
paste(genes_use, collapse = ', ')

smp.sub <- AddModuleScore(smp.sub, features =list(genes_use), name = 'Pathway' )

aa <- max(smp.sub$Pathway1)
bb <- min(smp.sub$Pathway1)

label_y= aa + (my_comparisons_height) *0.055 -0.02
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

library(ggpubr)

options(repr.plot.width=9, repr.plot.height=5)

p1 <- VlnPlot(smp.sub, 'Pathway1', group.by = 'Sample_DR_Type', split.by = 'patient', slot = 'data', cols = colors_use[c(1, 2, 4, 3, 5, 6)])+
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "solid",color="red", size=1)+
  NoLegend()+
  stat_compare_means(comparisons = my_comparisons2,method = 'wilcox.test',p.adjust.method="BH",label.x=1,label.y=label_y, size = 5)+
  #ylim(c(bb+0.05, max(label_y)*1.055+0.01))+
  scale_y_continuous(limits = c(bb, max(label_y)*1.055+0.01))+
  theme(legend.position="none", 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18))+
  ggtitle('')+xlab("")+ylab("Signature Score")

p1

# survival analysis
scores <- AddModuleScore_bulk(TCGA_expr.mat, features= list( genes_use))
scores <- as.vector(scores)
names(scores) <- colnames(TCGA_expr.mat)

patients <- intersect(rownames(TCGA_survival.df),colnames(TCGA_expr.mat))  
length(patients)
tcga.df <- cbind(TCGA_survival.df[patients, c('OS_STATUS','OS_MONTHS')], scores[patients])
tcga.df$OS_STATUS <- as.character(tcga.df$OS_STATUS) 
tcga.df$OS_STATUS <- as.numeric(substr(tcga.df$OS_STATUS, 1, 1))


colnames(tcga.df)[3] <- 'Value'  
tcga.df$group <- "Median"
tcga.df$group[tcga.df$Value > quantile(tcga.df$Value, 0.5)] <- "High"
tcga.df$group[tcga.df$Value <= quantile(tcga.df$Value, 0.5)] <- "Low"
tcga.df$group <- factor(tcga.df$group, levels = c('High', 'Low'), ordered = T)

surv_object <- Surv(time = tcga.df$OS_MONTHS, event = tcga.df$OS_STATUS)
fit <- survfit(surv_object ~ group, data =  tcga.df)

options(repr.plot.width=3.5, repr.plot.height=4.3)

p.surv <- ggsurvplot(fit, pval = T,
                     legend.labs=c("High" ,  "Low"),
                     palette = c("#DE77AE", "#7FBC41"),
                     # axes.offset = F,
                     #test.for.trend	= T,
                     #legend = "right",
                     ggtheme = theme_bw(),
                     legend.title = "Group", 
                     
                     conf.int = F,
                     # Add risk table
                     risk.table = TRUE,
                     tables.height = 0.25,
                     #tables.y.text = F, 
                     tables.theme = theme_cleantable())
p.surv$plot <- p.surv$plot + 
  xlab('Months')+
  theme(
    #       legend.key.height=unit(1.8, "line"),
    legend.position=c(0.8, 0.78),
    legend.background=element_blank(),
    legend.spacing.y = unit(0.3, "cm"), 
    plot.title = element_text(size = 14, vjust = 0.5), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    strip.text = element_text(size=14),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1))+ guides(fill=guide_legend(
        keywidth=0.5,
        keyheight=0.5,
        default.unit="inch")
      )

print(p.surv)

#----------------------------------------------------------------------------#
# Figure 3f
gmt <- read.gmt("/home/glab/lxq/reference/c5.go.bp.v2022.1.Hs.symbols.gmt")

pathway.sel <-'GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS'
genes.path <- gmt$gene[which(gmt$term == pathway.sel)]
length(genes.path)

genes_use <- intersect(genes.path, rownames(smp))
length(genes_use)
paste(genes_use, collapse = ', ')


smp.sub <- AddModuleScore(smp.sub, features =list(genes_use), name = 'Pathway' )

aa <- max(smp.sub$Pathway1)
bb <- min(smp.sub$Pathway1)

label_y= aa + (my_comparisons_height) *0.058-0.06
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

library(ggpubr)

options(repr.plot.width=9, repr.plot.height= 5)

p1 <- VlnPlot(smp.sub, 'Pathway1', group.by = 'Sample_DR_Type', split.by = 'patient', slot = 'data', cols = colors_use[c(1, 2, 4, 3, 5, 6)])+
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "solid",color="red", size=1)+
  NoLegend()+
  stat_compare_means(comparisons = my_comparisons2,method = 'wilcox',p.adjust.method="BH",label.x=1,label.y=label_y, size = 5.5)+
  scale_y_continuous(limits = c(bb+0.015, max(label_y)*1.08))+
  theme(legend.position="none", 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18))+
  ggtitle('')+xlab("")+ylab("Signature Score")
p1

# pdf('vlnplot_ER_Stress.pdf', height = 4.8, width = 9)
# print(p1)
# dev.off()


# survival analysis
scores <- AddModuleScore_bulk(TCGA_expr.mat, features= list( genes_use))
scores <- as.vector(scores)
names(scores) <- colnames(TCGA_expr.mat)

patients <- intersect(rownames(TCGA_survival.df),colnames(TCGA_expr.mat))  
length(patients)
tcga.df <- cbind(TCGA_survival.df[patients, c('OS_STATUS','OS_MONTHS')], scores[patients])
tcga.df$OS_STATUS <- as.character(tcga.df$OS_STATUS) 
tcga.df$OS_STATUS <- as.numeric(substr(tcga.df$OS_STATUS, 1, 1))


colnames(tcga.df)[3] <- 'Value'  
tcga.df$group <- "Median"
tcga.df$group[tcga.df$Value > quantile(tcga.df$Value, 0.5)] <- "High"
tcga.df$group[tcga.df$Value <= quantile(tcga.df$Value, 0.5)] <- "Low"
tcga.df$group <- factor(tcga.df$group, levels = c('High', 'Low'), ordered = T)

surv_object <- Surv(time = tcga.df$OS_MONTHS, event = tcga.df$OS_STATUS)
fit <- survfit(surv_object ~ group, data =  tcga.df)

options(repr.plot.width=3.5, repr.plot.height=4.3)

p.surv <- ggsurvplot(fit, pval = T,
                     legend.labs=c("High" ,  "Low"),
                     palette = c("#DE77AE", "#7FBC41"),
                     # axes.offset = F,
                     #test.for.trend	= T,
                     #legend = "right",
                     ggtheme = theme_bw(),
                     legend.title = "Group", 
                     
                     conf.int = F,
                     # Add risk table
                     risk.table = TRUE,
                     tables.height = 0.25,
                     #tables.y.text = F, 
                     tables.theme = theme_cleantable())
p.surv$plot <- p.surv$plot + 
  xlab('Months')+
  theme(
    #       legend.key.height=unit(1.8, "line"),
    legend.position=c(0.8, 0.78),
    legend.background=element_blank(),
    legend.spacing.y = unit(0.3, "cm"), 
    plot.title = element_text(size = 14, vjust = 0.5), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    strip.text = element_text(size=14),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1))+ guides(fill=guide_legend(
        keywidth=0.5,
        keyheight=0.5,
        default.unit="inch")
      )

print(p.surv)



#----------------------------------------------------------------------------#
# Supp Figure 5c-e
gmt <- read.gmt("/home/glab/lxq/reference/h.all.v2022.1.Hs.symbols.gmt")

# pathways used to show
pathway.sel <-'HALLMARK_ADIPOGENESIS'
pathway.sel <-'HALLMARK_FATTY_ACID_METABOLISM'
pathway.sel <-'HALLMARK_OXIDATIVE_PHOSPHORYLATION'

genes.path <- gmt$gene[which(gmt$term == pathway.sel)]
genes_use <- intersect(genes.path, rownames(smp))

smp.sub <- AddModuleScore(smp.sub, features =list(genes_use), name = 'Pathway' )

aa <- max(smp.sub$Pathway1)
bb <- min(smp.sub$Pathway1)

label_y= aa + (my_comparisons_height) *0.058-0.06
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

library(ggpubr)

options(repr.plot.width=9, repr.plot.height= 5)

p1 <- VlnPlot(smp.sub, 'Pathway1', group.by = 'Sample_DR_Type', split.by = 'patient', slot = 'data', cols = colors_use[c(1, 2, 4, 3, 5, 6)])+
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "solid",color="red", size=1)+
  NoLegend()+
  stat_compare_means(comparisons = my_comparisons2,method = 'wilcox',p.adjust.method="BH",label.x=1,label.y=label_y, size = 5.5)+
  scale_y_continuous(limits = c(bb+0.015, max(label_y)*1.08))+
  theme(legend.position="none", 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18))+
  ggtitle('')+xlab("")+ylab("Signature Score")
p1


# survival analysis
scores <- AddModuleScore_bulk(TCGA_expr.mat, features= list( genes_use))
scores <- as.vector(scores)
names(scores) <- colnames(TCGA_expr.mat)

patients <- intersect(rownames(TCGA_survival.df),colnames(TCGA_expr.mat))  
length(patients)
tcga.df <- cbind(TCGA_survival.df[patients, c('OS_STATUS','OS_MONTHS')], scores[patients])
tcga.df$OS_STATUS <- as.character(tcga.df$OS_STATUS) 
tcga.df$OS_STATUS <- as.numeric(substr(tcga.df$OS_STATUS, 1, 1))


colnames(tcga.df)[3] <- 'Value'  
tcga.df$group <- "Median"
tcga.df$group[tcga.df$Value > quantile(tcga.df$Value, 0.5)] <- "High"
tcga.df$group[tcga.df$Value <= quantile(tcga.df$Value, 0.5)] <- "Low"
tcga.df$group <- factor(tcga.df$group, levels = c('High', 'Low'), ordered = T)

surv_object <- Surv(time = tcga.df$OS_MONTHS, event = tcga.df$OS_STATUS)
fit <- survfit(surv_object ~ group, data =  tcga.df)

options(repr.plot.width=3.5, repr.plot.height=4.3)

p.surv <- ggsurvplot(fit, pval = T,
                     legend.labs=c("High" ,  "Low"),
                     palette = c("#DE77AE", "#7FBC41"),
                     # axes.offset = F,
                     #test.for.trend	= T,
                     #legend = "right",
                     ggtheme = theme_bw(),
                     legend.title = "Group", 
                     
                     conf.int = F,
                     # Add risk table
                     risk.table = TRUE,
                     tables.height = 0.25,
                     #tables.y.text = F, 
                     tables.theme = theme_cleantable())
p.surv$plot <- p.surv$plot + 
  xlab('Months')+
  theme(
    #       legend.key.height=unit(1.8, "line"),
    legend.position=c(0.8, 0.78),
    legend.background=element_blank(),
    legend.spacing.y = unit(0.3, "cm"), 
    plot.title = element_text(size = 14, vjust = 0.5), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    strip.text = element_text(size=14),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1))+ guides(fill=guide_legend(
        keywidth=0.5,
        keyheight=0.5,
        default.unit="inch")
      )

print(p.surv)

