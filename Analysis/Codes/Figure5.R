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


setwd('./Analysis/')
colors_use <- c("#4E7EAE", '#F24E4E', '#7DCE7C', '#7EB8CA', '#B38FBC', '#FDC7B2', '#C64E6B')
flow.meta <- read.csv('./Data/Mateiral/Flow_information.csv', header = TRUE, check.names = F)
smp <- readRDS('./Data/Combine_PDAC_Seurat_Object.RDS')


#------------------------------------------------------------------------------------------#
# Figure 5A, Genes ranked by average log2FC between gemcitabine-treated and control conditions across quiescence pattern flows

# Average expression of control and treated clusters for PDAC flows
markers <- read.csv( './Data/Figure4/Flow_Acquired_Resistance_DEGs_LogFC0_Pct10.csv')
genes_use <- unique(markers$Gene)
length(genes_use)

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
# write.csv(bulk.use, 'Acquired_Resistance_markers_bulk_logFC_matrix.csv')


# Average expression of control and treated clusters for PDAC flows
pre.mat <- matrix(0, nrow = nrow(smp), ncol = length(unique(markers$Flow)))
pos.mat <- matrix(0, nrow = nrow(smp), ncol = length(unique(markers$Flow)))

rownames(pre.mat) <- rownames(smp)
colnames(pre.mat) <- unique(markers$Flow)
rownames(pos.mat) <- rownames(smp)
colnames(pos.mat) <- unique(markers$Flow)
for (clus in colnames(pre.mat)){
  print(clus)
  compare.1 <- paste0(substr(clus, 1, 2),  '_', substr(clus, 5, 6))
  compare.2 <- paste0(substr(clus, 1, 2),  '_', substr(clus, 9, 10))
  bulk.1 <- smp@assays$RNA@data[, which(smp$Clusters %in% compare.1)]
  bulk.2 <- smp@assays$RNA@data[, which(smp$Clusters %in% compare.2)]
  pre.mat[rownames(bulk.1),clus] <- rowMeans(bulk.1)
  pos.mat[rownames(bulk.2),clus] <- rowMeans(bulk.2)
}
dim(pre.mat)
dim(pos.mat)


# Genes' average log2FC in quiescence pattern
cluster_bulk <- readRDS('./Data/Figure4/Acquired_Resistance_Cluster_bulk.RDS')
bulk.use <- read.csv(  './Data/Acquired_Resistance_markers_bulk_logFC_matrix.csv', row.names = 1, check.names = F)

markers_up <- markers[which(markers$Flow %in% flow.meta$Flow[flow.meta$Drug_Response_Type == 'Acquired Resistance'] & markers$p_val_adj < 0.05 & markers$avg_log2FC > 0), ]
markers_down <- markers[which(markers$Flow %in% flow.meta$Flow[flow.meta$Drug_Response_Type == 'Acquired Resistance'] & markers$p_val_adj < 0.05 & markers$avg_log2FC < 0), ]

flows.tmp <- colnames(cluster_bulk)[which(cluster_bulk$Mode == 'Quiescence' )]
bulk.flow <- bulk.use[, flows.tmp]
length(flows.tmp)

markers.flow.up <- markers_up[which(markers_up$p_val_adj < 0.05 & markers_up$avg_log2FC > 0 & markers_up$Flow %in% flows.tmp), ]
markers.flow.down <- markers_down[which(markers_down$p_val_adj < 0.05 & markers_down$avg_log2FC < 0 & markers_down$Flow %in% flows.tmp), ]

mat_up <- dcast(markers.flow.up[,c('Gene', 'Flow', 'avg_log2FC')], Gene~Flow)
rownames(mat_up) <- mat_up$Gene
mat_up<-mat_up[, -1]
mat_up[is.na(mat_up)] <- 0

mat_down <- dcast(markers.flow.down[,c('Gene', 'Flow', 'avg_log2FC')], Gene~Flow)
rownames(mat_down) <- mat_down$Gene
mat_down<-mat_down[, -1]
mat_down[is.na(mat_down)] <- 0

times_up <- table(markers.flow.up$Gene)
times_down <- table(markers.flow.down$Gene)

mean_delta <- apply(bulk.flow, 1, mean)
mean_percent_up <-  apply(bulk.flow,  1, function(x) sum(x>0))  / length(flows.tmp)
mean_percent_down <-  apply(bulk.flow,  1, function(x) sum(x<0)) / length(flows.tmp)


# paired Wilcoxon test for identifying consistently up- or down-regulated genes 
options(warn=-1)
genes_score$Quiescence_wilcox <- 1

mode.pre.mat <- pre.mat[, flows.tmp]
mode.pos.mat <- pos.mat[, flows.tmp]

for (gene in genes.use){
  df.tmp1<- mode.pre.mat[gene, ]
  df.tmp2 <- mode.pos.mat[gene, ]
  pval <- wilcox.test(as.vector(unlist(df.tmp1)), as.vector(unlist(df.tmp2)),  paired = T)        
  genes_score[gene, 'Quiescence_wilcox'] <- pval$p.value
}
genes_score$Quiescence_wilcox.adj <- p.adjust(genes_score$Quiescence_wilcox, method = 'BH')
genes_score <- genes_score[order(genes_score$Quiescence_avg_log2FC, decreasing = T), ]
write.csv(genes_score, paste0( './Data/Quiescence_genes_rank.csv'))



# plot Figure5A
genes.df <- read.csv('./Data/Quiescence_genes_rank.csv', 
                     row.names = 1)
genes.df$Type <- 'None'
genes.df$Type[which(genes.df$Quiescence_wilcox.adj < 0.05 & genes.df$Quiescence_avg_log2FC > 0 )] <- 'Up'
genes.df$Type[which(genes.df$Quiescence_wilcox.adj < 0.05 & genes.df$Quiescence_avg_log2FC < 0)] <- 'Down'
table(genes.df$Type)
genes.df$Rank <- rank(genes.df$Quiescence_avg_log2FC)
write.csv(genes.df, paste0( './Data/Quiescence_genes_rank.csv'))


genes.df.sub <- genes.df[which(genes.df$Quiescence_wilcox.adj < 0.05), ]
genes.df.sub <- genes.df.sub[genes.df.sub$Gene %in% c('GDF15', 'PHLDA2', 'LMO7', 'PIM1', 'FTL', 'ZBTB43', 
                                                      'CENPF', 'MKI67', 'TUBA1B', 'ASPM', 'TOP2A', 'CDKN3'), ]
genes.df.sub <- genes.df.sub[order(genes.df.sub$Quiescence_avg_log2FC),]


library(ggplot2)
library(ggrepel)

options(repr.plot.width=8.5, repr.plot.height=5.5)
genes.df$Type <- factor(genes.df$Type, levels = c('Up', 'None', 'Down'), ordered = T)

p1 <- ggplot(genes.df, aes(x=Rank, y=Quiescence_avg_log2FC, color=Type)) +
  geom_point(size = 3.5, alpha= 0.8) + # Show dots
  scale_color_manual(labels = c('Up regulated (98)', "Not changed","Down regulated (141)" ), values  = c("#CD2626", "grey", '#1874CD'))+
  xlab('Rank') +ylab('Average Log2FC')+labs(colour="Type")+
  theme_bw()+
  theme(
    legend.key.height=unit(1.8, "line"),
    legend.position = 'right',
    legend.key.size = unit(0.7, "cm"),
    plot.title = element_text(size = 18), 
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20, angle = 90, vjust= 1.5), 
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 0), 
    strip.text = element_text(size=16),
    axis.ticks.x=element_line(colour="black", linewidth = 1),
    axis.ticks.y=element_line(colour="black", linewidth = 1),
    axis.ticks=element_line(colour="red", size=.5,linewidth=1))+
  geom_text_repel( data = genes.df.sub[7:11, ],
    aes(x=Rank, y=Quiescence_avg_log2FC), 
    size = 5.5, color = 'black',fontface = "italic", hjust = "left",
    nudge_y = genes.df.sub$Quiescence_avg_log2FC[7:11]-0.5,
    nudge_x = -3500, direction = "y", label= genes.df.sub$Gene[7:11], max.overlaps = 10)+
  geom_text_repel(data = genes.df.sub[1:6, ], aes(x=Rank, y=Quiescence_avg_log2FC), 
    size = 5.5, color = 'black',fontface = "italic", hjust = "left",
    nudge_y = genes.df.sub$Quiescence_avg_log2FC[1:6]+0.75, nudge_x = 1500,
    direction = "y", label= genes.df.sub$Gene[1:6], max.overlaps = 10)+
  geom_text_repel( data = genes.df.sub[1, ], aes(x=Rank, y=Quiescence_avg_log2FC), 
    size = 5.5, color = 'black',fontface = "italic",
    label= genes.df.sub$Gene[1], max.overlaps = 10, nudge_y = -0.02)+
  geom_text_repel(
    data = genes.df.sub[12, ],
    aes(x=Rank, y=Quiescence_avg_log2FC), 
    size = 5.5, color = 'black',fontface = "italic",
    label= genes.df.sub$Gene[12], max.overlaps = 10, 
    nudge_y = 0.07)+
  xlim(c(0, 14000)) + ylim(c(-1.1, 0.815))

print(p1)



#------------------------------------------------------------------------------------------#
# Figure 5B
cluster_bulk <- readRDS("./Data/Figure4/Acquired_Resistance_Cluster_bulk.RDS")
scbulk_pre <- read.csv('./Data/Figure4/Flow_pre_mean_expression.csv', row.names = 1, check.names = F)
scbulk_pos <- read.csv('./Data/Figure4/Flow_pos_mean_expression.csv',  row.names = 1, check.names = F)

flow.meta$Modes <- flow.meta$Drug_Response_Type

# scRNA-seq modes
library(ggpubr)
library(ggtext)

gene <- 'GDF15'

df.tmp1<- scbulk_pre[gene, ]
df.tmp2 <- scbulk_pos[gene, ]
df.tmp <- rbind(df.tmp1, df.tmp2)
rownames(df.tmp) <- c('Ctrl', 'Gemc')
df.tmp <- as.data.frame(t(as.matrix(df.tmp)))

df.tmp$Group <- flow.meta$Mode[match(rownames(df.tmp), flow.meta$Flow)]
df.tmp$Group <- paste0(rownames(df.tmp), '_', df.tmp$Group)
df.tmp <- melt(df.tmp)

df.tmp$Flow <- substr(df.tmp$Group, 1, 10)
df.tmp$Mode <- substr(df.tmp$Group, 12, 45)
table(df.tmp$Mode)

df.tmp <- df.tmp[which(! df.tmp$Mode %in% c('Absolutely Intrinsic Resistance', 'Relatively Intrinsic Resistance', 'Sensitivity')), ]
df.tmp$Modes <- cluster_bulk$Mode[match(df.tmp$Flow, cluster_bulk$Sample)]
colnames(df.tmp)[c(2, 3)] <- c('Condition', 'Value')

p1.sc <- ggpaired(df.tmp, x = "Condition", y = "Value", id = 'Group',
                  fill = "Condition", line.color = "black", line.size = 0.2, lwd=3, width = 0.6,
                  palette = "jco", facet.by = "Modes", short.panel.labs = T, nrow = 1)+
  theme_bw()+
  ylab('GDF15 expression')+xlab('')+
  stat_compare_means( size = 6,paired = TRUE, method = "wilcox.test", label.y = rep(2.9, 4), label = "p",  label.x.npc = 0.1)+
  theme(axis.text.x = element_text( size=20),
        axis.text.y = element_markdown(size =18),
        axis.title.y = element_text( size=20, vjust = 2), 
        legend.position="none", 
        plot.title = element_text(size=20, hjust = 0.5),
        panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(size = 20) )+
  scale_fill_manual(values =c('#4E7EAE', '#F24E4E'))+
  ylim(c(0.1, 3.05))

options(repr.plot.width=8, repr.plot.height=4.5)
print(p1.sc)


g <- ggplot_gtable(ggplot_build(p1.sc))
stripr <- which(grepl('strip-t', g$layout$name))
stripr
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
grid.draw(g)




#------------------------------------------------------------------------------------------#
# Figure 5C
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
rownames(pre.expr) <- pre.data$gene_id
rownames(post.expr) <- pre.data$gene_id

# bulk PDX 
options(repr.plot.width=4, repr.plot.height=4)
gene <- 'GDF15'
gene.id <- pre.data$gene_id[match(gene, pre.data$gene_name)]
print(paste0(gene, '   ', gene.id))
df.tmp1<- pre.expr[gene.id, ]
df.tmp2 <- post.expr[gene.id, ]

df.tmp <- rbind(df.tmp1, df.tmp2)
rownames(df.tmp) <- c('Ctrl', 'Gemc')
df.tmp <- as.data.frame(t(as.matrix(df.tmp)))
df.tmp$Group <- rownames(df.tmp)
df.tmp <- melt(df.tmp)
colnames(df.tmp)[c(2, 3)] <- c('Condition', 'Value')
df.tmp$Value <- log2(df.tmp$Value + 1)
df.tmp$Modes <- 'Bulk PDX'
p1.bulk <- p1.sc <- ggpaired(df.tmp, x = "Condition", y = "Value", id = 'Group',
                             fill = "Condition", line.color = "black", line.size = 0.2, lwd=3, width = 0.6,
                             palette = "jco", short.panel.labs = T, nrow = 1, facet.by = "Modes")+
  theme_bw()+
  ylab('GDF15 expression')+xlab('')+
  stat_compare_means( size = 6,paired = TRUE, method = "wilcox.test", label.y = 10, label = "p",  label.x.npc = 0.15)+
  theme(axis.text.x = element_text( size=20),
        axis.text.y = element_markdown(size =18),
        axis.title.y = element_text( size=20, vjust = 2), 
        legend.position="none", 
        plot.title = element_text(size=20, hjust = 0.5),
        strip.text.x = element_text(size = 20) )+
  scale_fill_manual(values =c('#4E7EAE', '#F24E4E'))+
  ylim(c(3.9, 10.5))

dim( df.tmp)
options(repr.plot.width=3.2, repr.plot.height=5)
grid.draw(p1.bulk)




#------------------------------------------------------------------------------------------#
# Figure 5D
obj <- readRDS('/home/glab/lxq/PDAC/public_data/GSE186960/Seurat_for_Plot.RDS')
obj$gse <- 'GSE186960'
options(repr.plot.width=4.5, repr.plot.height=6)

my_comparisons = list(c1 = c('Ctrl', 'Gemc'), c2 = c('Ctrl', 'GemcResis'), c3 = c('Gemc', 'GemcResis'))

table(obj$orig.ident)
p2 <- VlnPlot(obj, 'GDF15', pt.size = 0.01)+
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "solid",color="red", size=1)+
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox',p.adjust.method="BH",label.x=1, label.y = c(4.5, 5.4, 5), size = 5)+
  theme_bw()+
  ylab('GDF15 expression')+xlab('')+
  theme(
    #axis.text.x = element_text( size=18, angle = 30, vjust = 1, hjust = 1),
    axis.text.x = element_text( size=20, angle = 0),
    axis.text.y = element_text(size =20),
    axis.title.y = element_text( size=20, vjust = 2), 
    legend.position="none", 
    plot.title = element_text(size=0, hjust = 0.5),
    strip.text.x = element_text(size = 22) )+
  scale_fill_manual(values =c('#4E7EAE', '#7DCE7C', '#F24E4E'))+
  ylim(c(-0.01, 5.8))+                                
  scale_x_discrete(labels = str_wrap(c('Ctrl', 'Gemc 24h', 'Gemc resistant'), width = 5))+
  facet_grid(~obj$gse)
print(p2)




#------------------------------------------------------------------------------------------#
# Figure 5E, # NCI-60
library(readr)
expr_raw <- read.table('./NCI-60/GSE116444_series_matrix.txt.gz',fill = T,
                       comment.char = "!",header = T) 

sample <- read_tsv('./NCI-60/sample.tsv')
ss <- str_split_fixed(sample$Title, "_", n=4)
sample$Cell_line <- ss[,1]
sample$Drug <- ss[,2]
sample$Concentration <- ss[,3]
sample$Treat_time <- ss[,4]

reftop <- read_tsv('./NCI-60/GSE116444.top.table.tsv')

dataset <- 'GSE116444'
gene <- 'GDF15'
id <- reftop$ID[which(reftop$Gene.symbol == gene)]
expr <- expr_raw[which(expr_raw$ID_REF %in% id),]
ss <- apply(expr[2:ncol(expr)], 1, sum)


# subplot 1
data <-  data.frame(Cell_line = sample$Cell_line,
                    Time = sample$Treat_time,
                    Concentration = sample$Concentration,
                    Treat = paste0(sample$Concentration, '-', sample$Treat_time),
                    Gene = as.vector(unlist(expr[2,2:ncol(expr)])))
data$Treat <- factor(data$Treat , levels = c('0nM-2h', '0nM-6h', '0nM-24h', '200nM-2h',   '200nM-6h',  '200nM-24h', 
                                             '2000nM-2h',  '2000nM-6h',  '2000nM-24h'), ordered = T)
data$Time <- factor(data$Time , levels = c('2h', '6h', '24h'), ordered = T)
data$Concentration <- factor(data$Concentration , levels = c('0nM', '200nM', '2000nM'), ordered = T)
data$Sample <- paste0(data$Cell_line, '_', data$Time)

my_comparisons = list(c1 = c('0nM', '200nM'), c2 = c('200nM', '2000nM'),  c3 = c('0nM', '2000nM'))

complete <- names(table(data$Sample))[table(data$Sample) == 3]
aa <- data[which(data$Sample %in% complete),]

data$gse <- data$Cell_line
data$gse <- 'NCI-60: Concentration'
table(data$gse)

options(repr.plot.width=9, repr.plot.height=6.3)

p1 <- ggplot(data = data, aes(x=Concentration, y=Gene)) +
  facet_grid(~Time)+
  geom_boxplot(aes(fill = Concentration),  outlier.size = 0, width = 0.6,  lwd=0.8)+
  geom_point( size=0.5)+
  scale_fill_manual(values = c('#4E7EAE', '#FDAE61','#F24E4E'))+
  ylab(paste0(gene, ' expression')) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text( size=20, angle = 20, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size =20),
        axis.title.x = element_text( size=20), 
        axis.title.y = element_text( size=20), 
        plot.title = element_text(size=20, hjust = 0.5),
        strip.text.x = element_text(size = 22), 
        legend.position = 'none', 
        panel.spacing = unit(0, "lines")
  )+
  geom_line(data = aa, mapping  =aes(x=Concentration, y=Gene, group = Sample), color = '#8F8F8F',  alpha =0.4)+
  stat_compare_means(mapping = aes(x=Concentration, y=Gene, group = Concentration), 
                     data = aa, method = "wilcox.test", comparisons = my_comparisons,  paired=T, label.y = c(13.4, 14.2, 15.1), label = "p", size =6)+
  xlab('')+
  ylim(c(min(data$Gene - 0.5), 16))
print(p1)



# subplot 2
data <-  data.frame(Cell_line = sample$Cell_line,
                    Time = sample$Treat_time,
                    Concentration = sample$Concentration,
                    Treat = paste0(sample$Concentration, '-', sample$Treat_time),
                    Gene = as.vector(unlist(expr[2,2:ncol(expr)])))
data$Treat <- factor(data$Treat , levels = c('0nM-2h', '0nM-6h', '0nM-24h', '200nM-2h',   '200nM-6h',  '200nM-24h', 
                                             '2000nM-2h',  '2000nM-6h',  '2000nM-24h'), ordered = T)
data$Time <- factor(data$Time , levels = c('2h', '6h', '24h'), ordered = T)
data$Concentration <- factor(data$Concentration , levels = c('0nM', '200nM', '2000nM'), ordered = T)
data$Sample <- paste0(data$Cell_line, '_', data$Time)

data$gse <- data$Cell_line
data$gse <- 'NCI-60: Treated time'
table(data$gse)

my_comparision <- list(c("2h", "6h"), c("6h", "24h"), c("2h", "24h"))
data$Sample <- paste0(data$Cell_line, '_', data$Concentration)

complete <- names(table(data$Sample))[table(data$Sample) == 3]
aa <- data[which(data$Sample %in% complete),]

options(repr.plot.width=9, repr.plot.height=6)
p1 <- ggplot(data = data, aes(x=Time, y=Gene)) +
  facet_grid(~Concentration)+
  geom_boxplot(aes(fill = Time),  outlier.size = 0, width = 0.6,  lwd=0.8)+
  geom_point(size=0.5)+
  scale_fill_manual(values = c('#4E7EAE', '#FDAE61','#F24E4E'))+
  ylab(paste0(gene, ' expression')) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text( size=20,hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(size =20),
        axis.title.x = element_text( size=20), 
        axis.title.y = element_text( size=20), 
        plot.title = element_text(size=20, hjust = 0.5),
        strip.text.x = element_text(size = 22), 
        legend.position = 'none', 
        panel.spacing = unit(0, "lines")
  )+
  geom_line(data = aa, mapping  =aes(x=Time, y=Gene, group = Sample), color = '#8F8F8F',  alpha =0.4)+
  stat_compare_means(mapping = aes(x=Time, y=Gene, group = Time), 
                     data = aa, comparisons = my_comparision,  paired=T, 
                     method = "wilcox.test", label.y = c(13.4, 14.2, 15.1), label = "p", size =6)+
  xlab('')+
  ylim(c(min(data$Gene - 0.5), 16))

print(p1)


