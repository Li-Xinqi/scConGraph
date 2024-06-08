library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggExtra)
library(gridExtra)
library(grid)
library(magrittr)
library(dplyr)
library(Seurat)
library(scales)
library(ggalluvial)
library(gridExtra)
library(grid)
library(patchwork)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(scCancer)
library(ggalluvial)
library(pheatmap)
library(survival)
library(survminer)


setwd('./Analysis/')
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')

# read info
flow.meta <- read.csv('Flow_information.csv', header = TRUE, check.names = F)
smp <- readRDS('Combine_PDAC_Seurat_Object.RDS')

#---------------------------------------------------------------------------------------------------#
# Figure 4A: Transcriptional change patterns for acquired resistant flows

# Calculate differential genes between aligned treated and control clusters
markers <- data.frame( p_val =c(),  avg_log2FC=c(), pct.1=c(), `pct.2`=c(), p_val_adj=c(), Gene = c(), Cluster = c(), Compare = c())
for (i in 1:nrow(flow.meta)){
  print(i)
  compare.1 <- flow.meta$Source[i]
  compare.2 <- flow.meta$Target[i]
  markers.sub <- FindMarkers(smp, ident.1 = compare.2, ident.2 = compare.1, min.pct = 0.1, logfc.threshold = 0)
  markers.sub$Gene <- rownames(markers.sub)
  markers.sub$Cluster <- compare.2
  markers.sub$Compare <- compare.1
  markers <- rbind(markers, markers.sub)
}

markers$Flow <- paste0(markers$Compare, '->', markers$Cluster)
table(markers$Flow)

markers$Flow <-  paste0(substr(markers$Cluster, 1, 2), ': ', substr(markers$Compare, 4, 5), '->', substr(markers$Cluster, 4, 5))
#write.csv(markers, './Figure4/Flow_Acquired_Resistance_DEGs_LogFC0_Pct10.csv', row.names = F)


# Log2FC value matrix for acquired resistant flows
markers <- markers[which(markers$Flow %in% flow.meta$Flow[flow.meta$Drug_Response_Type == 'Acquired Resistance'] & 
                           markers$p_val_adj < 0.05), ]
genes_use <- unique(markers$Gene)
length(genes_use)

## Use 12165 markers
bulk <- matrix(0, nrow = nrow(smp), ncol = length(unique(markers$Flow)))
rownames(bulk) <- rownames(smp)
colnames(bulk) <- unique(markers$Flow)
for (clus in colnames(bulk)){
  print(clus)
  compare.1 <- paste0(substr(clus, 1, 2),  '_', substr(clus, 5, 6))
  compare.2 <- paste0(substr(clus, 1, 2),  '_', substr(clus, 9, 10))
  bulk.1 <- smp@assays$RNA@data[, which(smp$Clusters %in% compare.1)]
  bulk.2 <- smp@assays$RNA@data[, which(smp$Clusters %in% compare.2)]
  bulk.delta <- log2(rowMeans(expm1(bulk.2))+1) - log2(rowMeans(expm1(bulk.1))+1)
  bulk[names(bulk.delta),clus] <- bulk.delta
}
bulk.use <- bulk[genes_use, ]
dim(bulk.use)


# Create Seurat object for flows of acquired resistance
cluster_bulk <- CreateSeuratObject(counts = bulk.use , project = "cluster", min.cells = 1)
dim(cluster_bulk)

cluster_bulk@assays$RNA@data <- cluster_bulk@assays$RNA@counts
cluster_bulk <- FindVariableFeatures(cluster_bulk, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cluster_bulk)
cluster_bulk <- ScaleData(cluster_bulk, features = all.genes)

cluster_bulk <- RunPCA(cluster_bulk, features = VariableFeatures(object = cluster_bulk), npcs = 30 )
cluster_bulk <- FindNeighbors(cluster_bulk, dims = 1:20, k.param =5)
cluster_bulk <- FindClusters(cluster_bulk, resolution = 0.5)
cluster_bulk <- RunUMAP(cluster_bulk, dims = 1:20)

cluster_bulk$Sample <- colnames(cluster_bulk)
cluster_bulk$Group <- substr( colnames(cluster_bulk), 1 , 2)

cluster_bulk <- FindClusters(cluster_bulk, resolution = 1)
cluster_bulk$seurat_clusters[which(colnames(cluster_bulk) == 'P2: C3->G4')] <- '2' # P2: C3->G4 is between Quiescence and EMT, but it is more like EMT.
table(cluster_bulk$seurat_clusters)
Idents(cluster_bulk) <- cluster_bulk$seurat_clusters

cluster_bulk$Mode <-  as.character(cluster_bulk$seurat_clusters)
cluster_bulk$Mode[which(cluster_bulk$Mode == '0')] <- 'Quiescence'
cluster_bulk$Mode[which(cluster_bulk$Mode == '1')] <- 'Cell Cycle I'
cluster_bulk$Mode[which(cluster_bulk$Mode == '2')] <- 'EMT'
cluster_bulk$Mode[which(cluster_bulk$Mode == '3')] <- 'Cell Cycle II'
table(cluster_bulk$Mode)


# plot Figure 4A, group by patterns
library(ggrepel)
clu.colors <-c( '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

options(repr.plot.width=5.5, repr.plot.height=6.2)
cluster_bulk$UMAP1 <- Embeddings(object = cluster_bulk[["umap"]])[,1]
cluster_bulk$UMAP2 <- Embeddings(object = cluster_bulk[["umap"]])[,2]
cluster_bulk$Flow <- colnames(cluster_bulk)

p <- cluster_bulk@meta.data  %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Mode)) + 
  geom_point( size = 6) + 
  scale_color_manual(values  =  clu.colors[c(2, 4, 3, 1 )])+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill='transparent', color='white'), 
        plot.title = element_text(hjust = 0.5,size=16), 
        legend.title =  element_text(colour = "black",size=20),
        axis.text = element_text(colour = "black",size=14), 
        axis.title = element_text(colour = "black",size=14), 
        axis.line = element_line(colour = "black",size = 0.5), 
        axis.ticks = element_line(colour= "black",size=0.5),
        legend.key = element_rect(fill='transparent', color='transparent'), 
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(colour = "black", face="plain", size=20))+
  labs(color='Pattern')+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "left"))
print(p)


# plot Figure 4A, group by patients

options(repr.plot.width=5.5, repr.plot.height=6.2)
p <- cluster_bulk@meta.data  %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Group)) + 
  geom_point( size = 6) + 
  scale_color_manual(values  =  colors_use[c(1, 2, 4, 3, 5, 6)])+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill='transparent', color='white'), 
        plot.title = element_text(hjust = 0.5,size=16), 
        legend.title =  element_text(colour = "black",size=20), 
        axis.text = element_text(colour = "black",size=14), 
        axis.title = element_text(colour = "black",size=14), 
        axis.line = element_line(colour = "black",size = 0.5), 
        axis.ticks = element_line(colour= "black",size=0.5),
        legend.key = element_rect(fill='transparent', color='transparent'), 
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(colour = "black", face="plain", size=20))+
  labs(color='Patient-derived')+
  guides(color=guide_legend(nrow=2,byrow=T, title.position = "left"))
print(p)

#saveRDS(cluster_bulk, "./Figure4/Acquired_Resistance_Cluster_bulk.RDS")




#---------------------------------------------------------------------------------------------------#
# Figure 4C: Typical genes consistently upregulated for each pattern. 
expr <- GetAssayData(cluster_bulk, slot = 'data')
mode0.cluster <- colnames(cluster_bulk)[which(cluster_bulk$seurat_clusters == 0)]
mode0.mean <- apply(expr[, mode0.cluster], 1, mean)
mode0.percent.up <- apply(expr[, mode0.cluster] > 0, 1, sum)/length(mode0.cluster)
mode0.percent.down <- apply(expr[, mode0.cluster] < 0, 1, sum)/length(mode0.cluster)
mode0.perct <- ifelse(mode0.mean > 0, mode0.percent.up, mode0.percent.down)


mode1.cluster <- colnames(cluster_bulk)[which(cluster_bulk$seurat_clusters == 1)]
mode1.mean <- apply(expr[, mode1.cluster], 1, mean)
mode1.percent.up <- apply(expr[, mode1.cluster] > 0, 1, sum)/length(mode1.cluster)
mode1.percent.down <- apply(expr[, mode1.cluster] < 0, 1, sum)/length(mode1.cluster)
mode1.perct <- ifelse(mode1.mean > 0, mode1.percent.up, mode1.percent.down)


mode2.cluster <- colnames(cluster_bulk)[which(cluster_bulk$seurat_clusters == 2)]
mode2.mean <- apply(expr[, mode2.cluster], 1, mean)
mode2.percent.up <- apply(expr[, mode2.cluster] > 0, 1, sum)/length(mode2.cluster)
mode2.percent.down <- apply(expr[, mode2.cluster] < 0, 1, sum)/length(mode2.cluster)
mode2.perct <- ifelse(mode2.mean > 0, mode2.percent.up, mode2.percent.down)


mode3.cluster <- colnames(cluster_bulk)[which(cluster_bulk$seurat_clusters == 3)]
mode3.mean <- apply(expr[, mode3.cluster], 1, mean)
mode3.percent.up <- apply(expr[, mode3.cluster] > 0, 1, sum)/length(mode3.cluster)
mode3.percent.down <- apply(expr[, mode3.cluster] < 0, 1, sum)/length(mode3.cluster)
mode3.perct <- ifelse(mode3.mean > 0, mode3.percent.up, mode3.percent.down)

genes.dot.df <- data.frame(Gene = rownames(cluster_bulk), 
                           Mode0_Mean = mode0.mean,  
                           Mode1_Mean = mode1.mean, 
                           Mode2_Mean = mode2.mean, 
                           Mode3_Mean = mode3.mean , 
                           Mode0_Percent = mode0.perct,  
                           Mode1_Percent = mode1.perct, 
                           Mode2_Percent = mode2.perct, 
                           Mode3_Percent = mode3.perct  
)
head(genes.dot.df)

#write.csv(genes.dot.df, 'Genes_mean_percent.csv')
genes.dot.df <- read.csv( 'Genes_mean_percent.csv', row.names = 1, header = T, check.names = F)
genes_use1 <- c(
  # Quiescence
  'CDKN1A', 'TP53INP2', 'PLAUR', 'PIM1', 'IRS2', 'ADM', 'VEGFA', 'PDGFA', 'TJP1', 'CLDN7', 'F11R'
)

genes_use2 <- c(
  # Cell Cycle I
  'STMN1', 'MKI67', 'CKS2','CCNB1', 'PCNA', 'RPA3',  'BRCA1', 'NPM1', 'CDK4',  'KPNA2'
)

genes_use3 <- c(
  #EMT
  'MMP14', 'FN1', 'COL6A2','EREG','IGFBP3','TWIST1', 'KRT19', 'IFIT3', 'STAT2', 'CD74', 'IFI44', 'ABCC10', 'ABCB8', 'ABCA12'
)

genes_use4 <- c(
  # Cell Cycle II
  'S100A4', 'S100A11','HMGCS1',  'FDFT1', 'ADIRF','DHCR7', 'EBP', 'SC5D', 'RHEB'
)

plot.df <- melt(genes.dot.df[, c('Gene', 'Mode0_Mean', 'Mode1_Mean', 'Mode2_Mean', 'Mode3_Mean')])
dim(plot.df)

perct.df <- melt(genes.dot.df[, c('Gene', 'Mode0_Percent', 'Mode1_Percent', 'Mode2_Percent', 'Mode3_Percent')])
dim(perct.df)

all(plot.df$Gene == perct.df$Gene)
all(substr(plot.df$variable, 1, 5) == substr(perct.df$variable, 1, 5))

colnames(plot.df) <- c('Gene', 'Mode', 'Mean')
plot.df$Percent <- perct.df$value
plot.df$Mode <- substr(plot.df$Mode, 1, 5)

library(ggthemes)

genes_use <- unique(c(genes_use2, genes_use4, genes_use3, genes_use1))


options(repr.plot.width=25, repr.plot.height=4)

plot.df.sub <- plot.df[which(plot.df$Gene %in% genes_use), ]
plot.df.sub$Gene <- factor(plot.df.sub$Gene, levels = unique(genes_use), ordered = T)
plot.df.sub$Mode <- paste0('Mode' , as.numeric(substr(plot.df.sub$Mode, 5, 5))+1)
plot.df.sub$Mode2 <-  plot.df.sub$Mode

plot.df.sub$Mode2[which(plot.df.sub$Mode == 'Mode1')] <- 'Quiescence'
plot.df.sub$Mode2[which(plot.df.sub$Mode == 'Mode2')] <- 'Cell Cycle I'
plot.df.sub$Mode2[which(plot.df.sub$Mode == 'Mode4')] <- 'Cell Cycle II'
plot.df.sub$Mode2[which(plot.df.sub$Mode == 'Mode3')] <- 'EMT'
plot.df.sub$Mode2 <- factor(plot.df.sub$Mode2, levels = c('Quiescence', 'EMT', 'Cell Cycle II', 'Cell Cycle I'), ordered = T)

plot.df.sub$Mean[which(plot.df.sub$Mean > 1)] <- 1
plot.df.sub$Mean[which(plot.df.sub$Mean < -1)] <- -1


# Figure 4C
p8 <- ggplot(plot.df.sub, aes(x = Gene, y = Mode2)) + 
  geom_point(aes(fill = Mean, size = Percent), alpha = 1, shape=21, colour = 'black') +
  # scale_fill_gradient2(  low = muted("dodgerblue"),
  #   mid = "white",
  #   high = muted("red"),
  #   midpoint = 0)+
  scale_fill_gradient2(  low = muted("#24693D"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0)+
  theme_bw()+
  scale_size(range = c(3, 13))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 18),
        axis.text.y = element_text( size = 20),
        axis.line = element_line(colour = "black"), 
        legend.box = "horizontal", 
        legend.title = element_text( size = 18),
        legend.text = element_text( size = 18) )+
  ylab('')+xlab('')+labs(fill=str_wrap("Average Log2FC", 10))
p8


#---------------------------------------------------------------------------------------------------#
# Figure 4B: Enrichment scores of typical hallmarks within flows across four distinct patterns. 
# GSEA HALLMARK
gmt <- read.gmt("h.all.v2022.1.Hs.symbols.gmt")
geneList<- bulk.use[, 1] 
names(geneList) <- rownames(bulk.use)
geneList <- sort(geneList, decreasing = T) #从高到低排序
GSEA_results <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1) #GSEA分析
GSEA_results <- as.data.frame(GSEA_results)
GSEA_results$Flow <- colnames(bulk.use)[1]
for (i in 2:ncol(bulk.use)){
  geneList<- bulk.use[, i] 
  names(geneList) <- rownames(bulk.use)
  geneList <- sort(geneList, decreasing = T) #从高到低排序
  GSEA_tmp <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1) #GSEA分析
  GSEA_tmp <- as.data.frame( GSEA_tmp)
  GSEA_tmp$Flow <- colnames(bulk.use)[i]
  GSEA_results <- rbind( GSEA_results,  GSEA_tmp)
}
write.csv(GSEA_results, './Data/Figure4/Acquired_Resistance_GSEA_HALLMARK.csv')

# GSEA KEGG
gmt <- read.gmt("c2.cp.kegg.v2022.1.Hs.symbols.gmt")
geneList<- bulk.use[, 1] 
names(geneList) <- rownames(bulk.use)
geneList <- sort(geneList, decreasing = T) #从高到低排序
GSEA_results <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1) #GSEA分析
GSEA_results <- as.data.frame(GSEA_results)
GSEA_results$Flow <- colnames(bulk.use)[1]
for (i in 2:ncol(bulk.use)){
  geneList<- bulk.use[, i] 
  names(geneList) <- rownames(bulk.use)
  geneList <- sort(geneList, decreasing = T) #从高到低排序
  GSEA_tmp <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1) #GSEA分析
  GSEA_tmp <- as.data.frame( GSEA_tmp)
  GSEA_tmp$Flow <- colnames(bulk.use)[i]
  GSEA_results <- rbind( GSEA_results,  GSEA_tmp)
}
write.csv(GSEA_results, './Data/Figure4/Acquired_Resistance_GSEA_KEGG.csv')


# Score
gmt1 <- read.gmt("./Analysis/Data/h.all.v2022.1.Hs.symbols.gmt")
gmt2 <- read.gmt("./Analysis/Data//c2.cp.kegg.v2022.1.Hs.symbols.gmt")
gmt <- rbind(gmt1, gmt2)

GSEA_results.h <- read.csv('./Data/Figure4/Acquired_Resistance_GSEA_HALLMARK.csv', header = T, row.names = 1, check.names = F)
GSEA_results.k <- read.csv('./Data/Figure4/Acquired_Resistance_GSEA_KEGG.csv', header = T, row.names = 1, check.names = F)
GSEA_results <- rbind(GSEA_results.h, GSEA_results.k)

for (p in unique(gmt$term)){
  path <- paste0(p, '_NES')
  GSEA_results_part <- GSEA_results[which(GSEA_results$ID == p), ]
  cluster_bulk@meta.data[, path] <- GSEA_results_part$NES[match(colnames(cluster_bulk), GSEA_results_part$Flow)]
  cluster_bulk@meta.data[, path] <- GSEA_results_part$NES[match(colnames(cluster_bulk), GSEA_results_part$Flow)]
}

library(ggrepel)
for (p in unique(gmt$term)){
  path <- paste0(p, '_NESrawAdj')
  score <- cluster_bulk@meta.data[, paste0(p, '_NES')]
  score <-  ifelse( score >= 0, 1, -1)* (abs(score) - min(abs(score)))/(max(abs(score)) - min(abs(score)))
  cluster_bulk@meta.data[, path] <- score
}


raw_name <- c(
  # Quiescence
  'HALLMARK_HYPOXIA_NESrawAdj', 'HALLMARK_ANGIOGENESIS_NESrawAdj',
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_NESrawAdj', 'HALLMARK_P53_PATHWAY_NESrawAdj',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB_NESrawAdj',
  # Cell Cycle I
  'KEGG_CELL_CYCLE_NESrawAdj', 'HALLMARK_DNA_REPAIR_NESrawAdj', 'HALLMARK_MYC_TARGETS_V1_NESrawAdj',
  'HALLMARK_FATTY_ACID_METABOLISM_NESrawAdj', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY_NESrawAdj', 
  # EMT
  'HALLMARK_KRAS_SIGNALING_UP_NESrawAdj', 'HALLMARK_INTERFERON_ALPHA_RESPONSE_NESrawAdj', 
  'HALLMARK_COAGULATION_NESrawAdj', 
  # Cell Cycle II
  'HALLMARK_MTORC1_SIGNALING_NESrawAdj', 'HALLMARK_CHOLESTEROL_HOMEOSTASIS_NESrawAdj'
  
  
)
replace_name <- c( 'Hypoxia', 'Angiogenesis', 'EMT', 'P53 Pathway', 'TNFa Signaling via NFkB', 
                   'Cell Cycle','DNA Repair', 'MYC Targets V1', 
                   'Fatty Acid Metabolism', 'Reactive Oxygen Species Pathway', 
                   'KRAS Signaling Up','Interferon Alpha Response',  
                   'Coagulation', 'MTORC1 Signaling', 'Cholesterol Homeostasis'
)
path.names.df <- data.frame(Raw = raw_name, Replace =replace_name )
path.names.df <- path.names.df[c(6, 7, 8, 14, 15, 9, 3, 12, 11,  4,  1, 5), ]


path.plot.df <- cluster_bulk@meta.data[, c(path.names.df$Raw, 'Mode')]

colnames(path.plot.df) <- path.names.df$Replace[match(colnames(path.plot.df), path.names.df$Raw)]
colnames(path.plot.df)[13] <- 'Mode'
path.plot.df$Flow <- rownames(path.plot.df)
path.plot.df <- melt(path.plot.df, value.name = "Mode")
colnames(path.plot.df) <- c('Mode', 'Flow', 'Phenotype', 'Value')

mat <- dcast(path.plot.df, Mode+Flow ~ Phenotype, value.var = 'Value')

options(repr.plot.width=4, repr.plot.height=4)
set.seed(2023)
options(repr.plot.width=13, repr.plot.height=9)

# plot Figure 4B
p <- ggplot(path.plot.df, aes(x = Mode, y = Value, fill = Mode)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.25 ) +
  geom_boxplot(width = 0.65, outlier.size = 0)+
  geom_jitter(color="black", size=1, alpha=0.9, height = 0, width = 0.2) +
  
  scale_fill_manual(values = clu.colors[c(2,4,3,1)])+
  xlab("")+ylab('Enrichment Score')+ labs(fill='Pattern')+
  theme_bw()+
  theme(
    text = element_text('sans'),
    plot.title = element_text(size = 20), 
    axis.text.y = element_text(size = 18), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20, angle = 90), 
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size=20), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(0, "lines"),
    legend.position = 'right')+
  facet_wrap(~Phenotype, nrow = 3, dir="v", labeller = label_wrap_gen(width=20))
p


g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))

fills <- rep( c(colorRampPalette(colors = c(clu.colors[2], "white"))(10)[5],
                colorRampPalette(colors = c(clu.colors[4], "white"))(10)[5],
                colorRampPalette(colors = c(clu.colors[3], "white"))(10)[5],
                colorRampPalette(colors = c(clu.colors[1], "white"))(10)[5]), 3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

options(repr.plot.width=14, repr.plot.height=9)
grid.draw(g)


#---------------------------------------------------------------------------------------------------#
# Figure 4E: UMAP plot showing log2FC vectors from three data types: bulk samples from the PUMCH cohort (Bulk), 
# scRNA-seq samples, and flows from four patterns. 

library(ComplexHeatmap)
library(circlize)

# process PUMCH Cohort
sam.anno <- read.table("./Data/PUMCH_Cohort/design.txt", header = T, check.names = F)
rownames(sam.anno) <- sam.anno$sample1

surv.info <- read.table("./Data/PUMCH_Cohort/survival-info.txt", header = T, sep = "\t", check.names = F)
rownames(surv.info) <- surv.info$PDX_ID


sam.anno$surv_time <- surv.info[sam.anno$Model_NO, ]$Survival_time_month
sam.anno$drug <- "other"
gemc.pat <- c("PAX-SH-115", "PAX-BJ-014", "PAX-SH-054", "PAX-SH-129", "PAX-SH-110", 
              "PAX-BJ-181", "PAX-BJ-153", "PAX-BJ-172", "PAX-BJ-168", "PAX-SH-095")
sam.anno$drug[sam.anno$Model_NO %in% gemc.pat] <- "Gemcitabine"


pre.data <- read.table("./Data/PUMCH_Cohort/gene_expression.preTreat.GRCH37.txt", 
                       header = T, check.names = F)

post.data <- read.table("./Data/PUMCH_Cohort/gene_expression.postTreat.GRCH37.txt",
                        header = T, check.names = F)

all.equal(gsub("PRE", "", colnames(pre.data)), gsub("POS", "", colnames(post.data)))

post.expr <- post.data[ , grep("TPM$", colnames(post.data))]
pre.expr <- pre.data[ , grep("TPM$", colnames(pre.data))]

colnames(post.expr) <- substr(colnames(post.expr) , 5, 10)
colnames(pre.expr) <- substr(colnames(pre.expr) , 5, 10)
colnames(post.expr) <-  gsub('-', '', colnames(post.expr))
colnames(pre.expr) <-  gsub('-', '', colnames(pre.expr))

patient.id <- intersect(rownames(sam.anno), colnames(pre.expr))
post.expr <- post.expr[, patient.id]
pre.expr <- pre.expr[, patient.id]
anno <- sam.anno[patient.id, ]
genes.delta <- log2((post.expr+1) / (pre.expr+1))
rownames(pre.expr) <- pre.data$gene_id
rownames(post.expr) <- pre.data$gene_id
rownames(genes.delta) <- pre.data$gene_id

# process scRNA-seq to bulk
use_features <- intersect(pre.data$gene_name, rownames(cluster_bulk))
vec1.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P1-Ctrl')], 1, sum)
vec1.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P1-Gemc')], 1, sum)
vec1 <- log2(vec1.2 + 1) - log2(vec1.1 + 1)

vec2.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P2-Ctrl')], 1, sum)
vec2.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P2-Gemc')], 1, sum)
vec2 <- log2(vec2.2 + 1) - log2(vec2.1 + 1)

vec3.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P3-Ctrl')], 1, sum)
vec3.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P3-Gemc')], 1, sum)
vec3 <- log2(vec3.2 + 1) - log2(vec3.1 + 1)

vec4.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P4-Ctrl')], 1, sum)
vec4.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P4-Gemc')], 1, sum)
vec4 <- log2(vec4.2 + 1) - log2(vec4.1 + 1)

vec5.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P5-Ctrl')], 1, sum)
vec5.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P5-Gemc')], 1, sum)
vec5 <- log2(vec5.2 + 1) - log2(vec5.1 + 1)

vec6.1 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P6-Ctrl')], 1, sum)
vec6.2 <- apply(GetAssayData(smp, slot = 'counts')[use_features, which(smp$orig.ident == 'P6-Gemc')], 1, sum)
vec6 <- log2(vec6.2 + 1) - log2(vec6.1 + 1)
expr3 <- cbind(vec1, vec2, vec3, vec4, vec5, vec6)

# create Seurat object
use_features <- intersect(pre.data$gene_name, rownames(cluster_bulk))
gene_id <- pre.data$gene_id[match(use_features, pre.data$gene_name)]
gene.delta.use <- genes.delta[gene_id, ]
rownames(gene.delta.use) <- pre.data$gene_name[match(rownames(gene.delta.use) , pre.data$gene_id)]
expr2 <- as.matrix(gene.delta.use)

expr1 <- as.matrix(GetAssayData(cluster_bulk, slot = 'data')[use_features, ])
colnames(expr1) <- paste0(cluster_bulk$Mode, ' ',colnames(expr1) )

dim(expr1)
dim(expr2)
dim(expr3)

colnames(expr1) <- paste0('SC_', colnames(expr1) )
colnames(expr2) <- paste0('BU_', colnames(expr2) )
colnames(expr3) <- paste0('P', 1:6)
colnames(expr3) <- paste0('CO_', colnames(expr3) )

expr <- cbind(expr1, expr2, expr3)
expr[is.na(expr)] <- 0


delta.obj <- CreateSeuratObject(counts = expr , project = "cluster", min.cells = 1)
dim(delta.obj)
delta.obj@assays$RNA@data <- delta.obj@assays$RNA@counts
delta.obj <- FindVariableFeatures(delta.obj, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(delta.obj)
delta.obj <- ScaleData(delta.obj, features = all.genes, vars.to.regress = "orig.ident")

delta.obj <- RunPCA(delta.obj, features = VariableFeatures(object = delta.obj), npcs = 50 )
delta.obj <- FindNeighbors(delta.obj, dims = 1:30, k.param =8)
delta.obj <- FindClusters(delta.obj, resolution = 0.5)
delta.obj <- RunUMAP(delta.obj, dims = 1:30)
delta.obj$Sample <- colnames(delta.obj)
delta.obj$Group <- substr( colnames(delta.obj),1 , 2)


delta.obj$Mode <- substr(colnames(delta.obj), 4, 17)
delta.obj$Mode[which(delta.obj$orig.ident == 'BU')] <- 'Bulk'
delta.obj$Mode[which(delta.obj$orig.ident == 'CO')] <- 'scRNA-seq'
delta.obj$Mode[which(delta.obj$Mode %in% c('Cell Cycle I P' ))] <- 'Cell Cycle I'
delta.obj$Mode[which(delta.obj$Mode %in% c('Quiescence P1:', 'Quiescence P3:', 'Quiescence P4:', 
                                           'Quiescence P5:', 'Quiescence P6:'))] <- 'Quiescence'
delta.obj$Mode[which(delta.obj$Mode %in% c('Cell Cycle II '))] <- 'Cell Cycle II'
delta.obj$Mode[grep('^EMT ', delta.obj$Mode)] <- 'EMT'
delta.obj$Mode <- factor(delta.obj$Mode, levels = c('Bulk', 'scRNA-seq', 'Cell Cycle I', 'Cell Cycle II', 'EMT', 'Quiescence'))
table(delta.obj$Mode)


# plot Figure 4E
delta.obj$UMAP1 <- Embeddings(object = delta.obj[["umap"]])[,1]
delta.obj$UMAP2 <- Embeddings(object = delta.obj[["umap"]])[,2]

options(repr.plot.width=5.5, repr.plot.height=7)

p <- delta.obj@meta.data  %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Mode)) + 
  geom_point( size =5.5) + 
  scale_color_manual(values  =  c(colors_use[6], 'black', clu.colors[c(2, 4, 3, 1)]))+
  #geom_text_repel(aes(label= Flow ), color = 'black')+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill='transparent', color='white'), 
        #plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"), 
        plot.title = element_text(hjust = 0.5,size=16), 
        legend.title =  element_text(colour = "black",size=20),
        axis.text = element_text(colour = "black",size=14), 
        axis.title = element_text(colour = "black",size=14), 
        axis.line = element_line(colour = "black",size = 0.5), 
        axis.ticks = element_line(colour= "black",size=0.5),
        legend.key = element_rect(fill='transparent', color='transparent'), 
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(colour = "black", face="plain", size=20))+
  labs(color='Type')+
  guides(color=guide_legend(nrow=3,byrow=TRUE, title.position = "left"))
print(p)



#---------------------------------------------------------------------------------------------------#
# Figure 4F-G: Simpson indices for bulk samples from the PUMCH cohort. 
delta.obj$Sample <- as.character(delta.obj$orig.ident)
delta.obj$Sample[which(delta.obj$Sample == 'BU')] <- substr(colnames(delta.obj)[which(delta.obj$Sample == 'BU')], 4, 8)
delta.obj$Sample[which(delta.obj$Sample == 'CO')] <- substr(colnames(delta.obj)[which(delta.obj$Sample == 'CO')], 4, 5)
table(delta.obj$Sample)

delta.obj$TGI <- 0
delta.obj$TGI <- sam.anno$`TGI(%)`[match(delta.obj$Sample, rownames(sam.anno))]/100
delta.obj$TGI[which(delta.obj$Sample %in% paste0('P', 1:6))] <- c(0.7391, 0.6349, 0.3963, 0.5143, 0.7765, 0.8037) 


umap.df <- Embeddings(object = delta.obj[["umap"]])
dist.mat <- dist(as.matrix(umap.df), method = "euclidean")
dist.mat <- as.matrix(dist.mat)
dim(dist.mat)
dist.mat <- dist.mat[colnames(delta.obj)[which(delta.obj$orig.ident != 'SC')], colnames(delta.obj)[which(delta.obj$orig.ident == 'SC')]]
outliers <- c('SC_Quiescence P5: C6->G0',  'SC_EMT P2: C3->G4', 'SC_Cell Cycle I P6: C1->G3')
dist.mat <- dist.mat[, setdiff(colnames(dist.mat), outliers)]

max(dist.mat)
dim(dist.mat)

for (ii in 1:nrow(dist.mat) ){
  k_dist <- dist.mat[ii, order(dist.mat[ii, ], decreasing = F)[3]]
  dist.mat[ii, ] <- exp( -((dist.mat[ii, ]^2)/(k_dist))^1)
}



simpson.df <- data.frame(Sample = rownames(dist.mat), Simpson = rep(0, nrow(dist.mat)), 
                         TGI = delta.obj$TGI[match( rownames(dist.mat), colnames(delta.obj) )] , 
                         CellCycleI = rep(0, nrow(dist.mat)),
                         CellCycleII = rep(0, nrow(dist.mat)),
                         P2specific = rep(0, nrow(dist.mat)),
                         Quiescence = rep(0, nrow(dist.mat)))
nrow(simpson.df)
simpson <- c()
num <- 3

for (ii in 1:nrow(dist.mat) ){
  point <- rownames(dist.mat)[ii]
  # Cell Cycle I
  cols <- setdiff(colnames(delta.obj)[which(delta.obj$Mode == 'Cell Cycle I')], outliers)
  neighbor.cc1 <- cols[order(dist.mat[ii, cols], decreasing = T)][1:num]
  # Cell Cycle II
  cols <- colnames(delta.obj)[which(delta.obj$Mode == 'Cell Cycle II')]        
  neighbor.cc2 <- cols[order(dist.mat[ii, cols], decreasing = T)][1:num]
  # EMT
  cols <- setdiff(colnames(delta.obj)[which(delta.obj$Mode == 'EMT')], outliers)
  neighbor.emt <- cols[order(dist.mat[ii, cols], decreasing = T)][1:num]
  # Quiescence
  cols <- setdiff(colnames(delta.obj)[which(delta.obj$Mode == 'Quiescence')], outliers)
  neighbor.quie <- cols[order(dist.mat[ii, cols], decreasing = T)][1:num]
  
  prob <- c(sum(dist.mat[ii, neighbor.cc1]), 
            sum(dist.mat[ii, neighbor.cc2]), 
            sum(dist.mat[ii, neighbor.emt]), 
            sum(dist.mat[ii, neighbor.quie]))
  names(prob) <- c('Cell Cycle I', 'Cell Cycle II', 'EMT', 'Quiescence')
  
  prob <- prob/sum(prob)
  simpson.df$CellCycleI[ii] <- prob[1]
  simpson.df$CellCycleII[ii]  <- prob[2]
  simpson.df$P2specific[ii]  <- prob[3]
  simpson.df$Quiescence[ii]  <- prob[4]
  prob <- 1 - sum(prob^2)
  simpson <- c(simpson, prob)
  simpson.df$Simpson[ii] <- prob
}
print(cor(simpson.df$TGI, simpson.df$Simpson, method = 'pearson'))
print(cor(simpson.df$TGI, simpson.df$Simpson, method = 'spearman'))
delta.obj$Simpson <- simpson.df$Simpson[match(colnames(delta.obj), simpson.df$Sample)]


# Figure 4F
delta.obj$Simpson <- simpson.df$Simpson[match(colnames(delta.obj), simpson.df$Sample)]

options(repr.plot.width=4.5, repr.plot.height=4)
plot <- delta.obj@meta.data  %>%
  ggplot(aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = Simpson), size = 4) + 
  scale_color_continuous(na.value="grey80")+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill='transparent', color='white'), 
        plot.title = element_text(hjust = 0.5,size=20), 
        legend.title =  element_text(colour = "black",size=12),
        axis.text = element_text(colour = "black",size=14), 
        axis.title = element_text(colour = "black",size=14), 
        axis.line = element_line(colour = "black",size = 0.5), 
        axis.ticks = element_line(colour= "black",size=0.5),
        legend.text = element_text(colour = "black", face="plain", size=14))+
  labs('Score')+ labs(colour = "") 

print(plot)



# density plot for simpson index
g <- simpson.df %>%
  ggplot( aes(x=Simpson)) +
  geom_density(alpha=0.6, linewidth =1.5, colour = 'grey20' )+
  ylab("Density") +
  xlab("Simpson Index")
print(g)

# threshold for simpson index
dens <- layer_data(g, 1)
dens <- dens[order(dens$x),]
rle <- rle(diff(as.vector(dens$y)) < 0)
starts <- cumsum(rle$lengths) - rle$lengths + 1
maxima_id <- starts[!rle$values]
maxima1 <- dens[maxima_id[2],]
print(maxima1)

# Supplementary Figure 7B
g2 <- g + geom_vline(data = maxima1,
                     aes(xintercept = x), colour = 'red')+
  theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16), 
    strip.text = element_text(size=16),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1,
      lineend=1)
  )
options(repr.plot.width=4, repr.plot.height=3)
print(g2)


# Figure 4G
simpson.df$Simpson_Group <- 'High'
simpson.df$Simpson_Group[which(simpson.df$Simpson < maxima1$x)] <- 'Low'
simpson.df$Simpson_Group <- factor(simpson.df$Simpson_Group, levels =  c('High', 'Low'), ordered = T)

p <- ggplot(simpson.df, aes(x = Simpson_Group, y = TGI, fill = Simpson_Group)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.25 ) +
  geom_boxplot(width = 0.6, outlier.size = 0)+
  geom_jitter(color="black", size=1, alpha=0.9, height = 0, width = 0.2) +
  
  scale_fill_manual(values = c("#DE77AE", "#7FBC41"))+
  xlab("")+
  ylab('TGI %')+ labs(fill='Modes')+
  theme_bw()+
  theme(
    text = element_text('sans'),
    plot.title = element_text(size = 18), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 0), 
    axis.title.y = element_text(size = 20, angle = 90), 
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    strip.text = element_text(size=16), 
    axis.ticks.x=element_blank(),
    panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(0, "lines"),
    legend.position = 'none')+
  stat_compare_means(method = 'wilcox',label = "p", label.x.npc = 0.17, 
                     label.y=c(1.8), size =7.5)+
  ylim(c(-0.2, 2))

options(repr.plot.width=4, repr.plot.height=3)
print(p)



#---------------------------------------------------------------------------------------------------#
# Figure 4H, Survival analysis
# Calculate scores of gene signatures for bulk RNA-seq samples 
AddModuleScore_bulk <- function (object, features, nbin = 24, ctrl = 100, 
                                 k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = 1)
  
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = rownames(x = object))
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ", 
              paste(missing.features, collapse = ", "), 
              ifelse(delta.obj = FALSE, yes = ", attempting to find updated synonyms", 
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


TCGA_survival.df <- read.csv( 'TCGA_PAAD_Survival_info.csv', header = T, row.names = 1, check.names = F)
TCGA_pheno.df <- read.csv( 'TCGA_PAAD_Phenotype_info.csv', header = T, row.names = 1, check.names = F)
TCGA_expr.mat <- read.csv( 'TCGA_PAAD_FPKM.csv', header = T, row.names = 1, check.names = F)
TCGA_pheno.df <- TCGA_pheno.df[which(TCGA_pheno.df$primary_diagnosis.diagnoses == 'Infiltrating duct carcinoma, NOS' &
                                       TCGA_pheno.df$sample_type.samples == 'Primary Tumor'), ]

TCGA_survival.df <- TCGA_survival.df[TCGA_pheno.df$submitter_id.samples, ]
TCGA_expr.mat <- TCGA_expr.mat[, TCGA_pheno.df$submitter_id.samples]
colnames(TCGA_survival.df) <- c('sample', 'OS_STATUS', 'Patient', 'OS_DAYS')
TCGA_survival.df$OS_MONTHS <- round(TCGA_survival.df$OS_DAYS / 30, 4)

# Extract important genes for each pattern
genes.dot.df$Mode0_mean_rank <- rank(-abs(genes.dot.df$Mode0_Mean), ties.method = 'min')
genes.dot.df$Mode0_percent_rank <- rank(-abs(genes.dot.df$Mode0_Percent), ties.method = 'min')
genes.dot.df$Mode0_rank <- genes.dot.df$Mode0_mean_rank  + genes.dot.df$Mode0_percent_rank 

genes.dot.df$Mode1_mean_rank <- rank(-abs(genes.dot.df$Mode1_Mean), ties.method = 'min')
genes.dot.df$Mode1_percent_rank <- rank(-abs(genes.dot.df$Mode1_Percent), ties.method = 'min')
genes.dot.df$Mode1_rank <- genes.dot.df$Mode1_mean_rank + genes.dot.df$Mode1_percent_rank

genes.dot.df$Mode2_mean_rank <- rank(-abs(genes.dot.df$Mode2_Mean), ties.method = 'min')
genes.dot.df$Mode2_percent_rank <- rank(-abs(genes.dot.df$Mode2_Percent), ties.method = 'min')
genes.dot.df$Mode2_rank <- genes.dot.df$Mode2_mean_rank + genes.dot.df$Mode2_percent_rank

genes.dot.df$Mode3_mean_rank <- rank(-abs(genes.dot.df$Mode3_Mean), ties.method = 'min')
genes.dot.df$Mode3_percent_rank <- rank(-abs(genes.dot.df$Mode3_Percent), ties.method = 'average')
genes.dot.df$Mode3_rank <- genes.dot.df$Mode3_mean_rank + genes.dot.df$Mode3_percent_rank


# Quiescence
genes.dot.df.up <- genes.dot.df[which(genes.dot.df$Mode0_Mean > 0 ), ]
dim(genes.dot.df.up)
genes_use <-  genes.dot.df.up$Gene[order(genes.dot.df.up$Mode0_rank, decreasing = F)][1:200]
genes_use[grep('^RP', genes_use)]
genes_use <- setdiff(genes_use, c('RP11-138A9.2', 'RP11-363E7.4', 'RP1-313I6.12', 'RP11-540A21.2', 'RP11-138A9.1', 'RP6-99M1.2'))

# Quiescence
scores <- AddModuleScore_bulk(TCGA_expr.mat, features= list(genes_use))
length(genes_use)
scores <- as.vector(scores)
names(scores) <- colnames(TCGA_expr.mat)

patients <- intersect(rownames(TCGA_survival.df),colnames(TCGA_expr.mat))  
tcga.df <- cbind(TCGA_survival.df[patients, c('OS_STATUS','OS_MONTHS')], scores[patients])
tcga.df$OS_STATUS <- as.character(tcga.df$OS_STATUS) 
tcga.df$OS_STATUS <- as.numeric(substr(tcga.df$OS_STATUS, 1, 1))
colnames(tcga.df)[3] <- 'Value'  

tcga.df$group <- "Median"
tcga.df$group[tcga.df$Value > quantile(tcga.df$Value, 0.5)] <- "High"
tcga.df$group[tcga.df$Value <= quantile(tcga.df$Value, 0.5)] <- "Low"
surv_object <- Surv(time = tcga.df$OS_MONTHS, event = tcga.df$OS_STATUS)
fit <- survfit(surv_object ~ group, data =  tcga.df)


options(repr.plot.width=3.1, repr.plot.height=4.2)

p.surv.tcga <- ggsurvplot(fit, pval = TRUE,
                          legend.labs=c("High" ,  "Low"),
                          title = paste0(''),
                          ggtheme = theme_bw(),
                          legend.title = "Group", 
                          conf.int = F,
                          risk.table = TRUE,
                          tables.height = 0.25,
                          tables.theme = theme_cleantable())
p.surv.tcga$plot <- p.surv.tcga$plot + 
  xlab('Months')+
  theme(
    legend.position=c(0.8, 0.78),
    plot.title = element_text(size = 14, vjust = 0.5), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12), 
    strip.text = element_text(size=14),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(#坐标轴刻度线的设置
      colour="red",
      size=.5,
      linetype=1))

print(p.surv.tcga) 

# The code for plotting Kaplan–Meier survival curves for the other three patterns is the same as for the quiescence pattern.
# Cell Cycle I
genes.dot.df.up <- genes.dot.df[which(genes.dot.df$Mode1_Mean > 0 ), ]
print(genes_use[grep('^RP', genes_use)])
genes_use <-  genes.dot.df.up$Gene[order(genes.dot.df.up$Mode1_rank, decreasing = F)][1:200]

# EMT
genes.dot.df.up <- genes.dot.df[which(genes.dot.df$Mode2_Mean > 0  ), ]
genes_use <-  genes.dot.df.up$Gene[order(genes.dot.df.up$Mode2_rank, decreasing = F)][1:200]
print(genes_use[grep('^RP', genes_use)])
genes_use <- setdiff(genes_use,  c("RP5-857K21.11", "RP11-174G6.5",  "RP11-362F19.1"))

# Cell Cycle II
genes.dot.df.up <- genes.dot.df[which(genes.dot.df$Mode3_Mean > 0 &  genes.dot.df$Mode3_Percent > 0.6 ), ]
genes_use <-  genes.dot.df.up$Gene[order(genes.dot.df.up$Mode3_rank, decreasing = F)][1:200]
print(genes_use[grep('^RP', genes_use)])
genes_use <- setdiff(genes_use, c('RP11-466H18.1'))
