library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(pheatmap)
library(ggpubr)


# Figure 6B
expr.data <- read.table('./Data/Figure6/expr_counts_final.txt',
                        header = T, sep = "\t", stringsAsFactors = F, check.names = F)
selGenes <- read.table('./Data/Figure6/Gene_annotation.txt',
                       header = T, sep = "\t", stringsAsFactors = F, check.names = F)

genes.use <- selGenes$gene_name[which(selGenes$gene_biotype == 'protein_coding')]
genes.use <- intersect(genes.use, rownames(expr.data))
expr.data <- expr.data[genes.use, ]

meta.data <- data.frame(Sample = colnames(expr.data), nCount = rep(1, 6), nFeature = rep(1, 6), percent.mito = rep(1, 6))
meta.data$nCount <- apply(expr.data, 2, sum)
meta.data$nFeature <- apply(expr.data > 0 , 2, sum)

Group <- c(rep('Knockdown', 3), rep('Control', 3))
meta.data$group <- Group

# DEseq2
dds <- DESeqDataSetFromMatrix( countData = expr.data, colData = meta.data, design = ~ group) 
dds <- dds[rowSums(counts(dds)) > 1, ] 
dds <- DESeq(dds) 
res <- results(dds)
sum(res$pvalue < 0.05, na.rm=TRUE)
summary(res)
#write.table(res, paste0('KDvsC_DEgenes_all.txt'), sep = "\t", quote = F, row.names = T)

resOrdered <- res[order(res$pvalue), ]
resOrdered$gene <- rownames(resOrdered)
resOrdered <- resOrdered[, c("gene", colnames(resOrdered)[1:(ncol(resOrdered)-1)])]
resOrdered <- subset(resOrdered, pvalue < 0.05 & !is.na(res$pvalue))


# Volcano plot

# The function to scale the coordinate axes
library(ggrepel)
squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  # return the transformation
  return(trans_new("squished", trans, inv))
}


res.selected <-  as.data.frame(subset(res,!is.na(res$padj)))
dim(res.selected)

res.selected$Direction <- 'NS'
res.selected$Direction[which(res.selected$log2FoldChange > 0.25 & res.selected$padj < 0.05)] <- 'Up'
res.selected$Direction[which(res.selected$log2FoldChange < -0.25 & res.selected$padj < 0.05)] <- 'Down'
res.selected$Name <- rownames(res.selected)
res.selected$Direction <- factor(res.selected$Direction, levels = c('Down', 'NS', 'Up'), ordered = T)
res.selected$padj.log <- abs(log10(res.selected$padj))

res.selected2 <- res.selected[which(res.selected$Direction %in% c('Up', 'Down')), ]
min(abs(res.selected2$log2FoldChange))

library(scales)

options(repr.plot.width=9, repr.plot.height=6)
pp <- ggplot(res.selected, aes(x = log2FoldChange, y = padj.log, colour=Direction)) +
  geom_point(alpha=0.4, size=3) +
  scale_y_continuous(trans = squish_trans(50, 75, 12),
                     breaks = c(0, 10, 20, 30, 40, 50,  75, 80))+
  scale_color_manual(values=c("#546de5", "grey","#ff4757"))+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +
  labs(x="log2(FoldChange, Knockdown/Ctrl)", y="-log10(p.adj)", color = "")+
  scale_color_manual(labels = c("Down regulated (297)", "Not changed",  'Up regulated (708)'), values  = c('#1874CD',  "grey", "#CD2626"))+
  geom_text_repel(data = res.selected2, 
                  aes(x = log2FoldChange, y = padj.log, 
                      label = Name), color = 'black',
                  size = 4,box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.2, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE, 
                  max.overlaps = 20)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  theme(axis.title.y = element_text(size = 14, face = "plain"),
        axis.title.x = element_text(size = 14, face = "plain"),
        axis.text.x = element_text(size = 14, face = "plain"),
        axis.text.y = element_text(size = 14, face = "plain",))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right") +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.key.height = unit(0.8, "cm"), 
        legend.key.width = unit(0.3, "cm"))


print(pp)



#-----------------------------------------------------------------------------#
# Figure 6C
gmt <- read.gmt("./Data/h.all.v2022.1.Hs.symbols.gmt")
# up
gene_df <- res.selected[which(res.selected$Direction == 'Up'), ]
genes <- gene_df$Name
egmt_h1 <- enricher(genes, TERM2GENE = gmt, 
                    pvalueCutoff =0.05, qvalueCutoff = 1)
hall.up <- as.data.frame(egmt_h1)
hall.up$Direction <- 'Up'
# down
gene_df <- res.selected[which(res.selected$Direction == 'Down'), ]
genes <- gene_df$Name
egmt_h2 <- enricher(genes, TERM2GENE = gmt, 
                    pvalueCutoff =0.05, qvalueCutoff = 1)
hall.down <- as.data.frame(egmt_h2)
hall.down$Direction <- 'Down'


plot.df <- rbind(hall.up, hall.down)
plot.df$ID  <- c('Mitotic Spindle', 'PI3K/AKT/mTOR signaling', 'TNFa signaling via NFkB', 'Interferon gamma response', 
                 'Interferon alpha response', 'P53 pathway')

plot.df$log10padj <- -log10(plot.df$p.adjust)
plot.df$ID <- factor(plot.df$ID , levels = rev(plot.df$ID[c(3, 4, 5, 6, 1, 2)]), ordered = T)
plot.df$Direction <- factor(plot.df$Direction , levels = c('Down', 'Up'), ordered = T)

p <- ggplot(data=plot.df, aes( x=ID,y=log10padj,fill = Direction))+
  coord_flip()+xlab('hallmark ID')+
  
  geom_bar(stat="identity")+theme_minimal()+scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme(axis.title.y = element_text(size=20, face="plain", colour = 'black'),
        axis.text.x = element_text(size=20, face="plain", colour = 'black'),
        axis.text.y = element_text(size=20, face="plain", colour = 'black'),  
        axis.title.x = element_text(size=20, face="plain", colour = 'black'))+
  theme(legend.text=element_text(size=20, face="plain", colour = 'black'), 
        legend.title=element_text(size=20, face="plain", colour = 'black'))+
  ylab('|log10(p.adj)|') + xlab('')

options(repr.plot.width=8.2, repr.plot.height=3.2)
print(p)



#-----------------------------------------------------------------------------#
# Figure 6F

sam.anno <- read.table("./Data/PUMCH_Cohort/design.txt", header = T, check.names = F)
rownames(sam.anno) <- sam.anno$sample1

surv.info <- read.table("./Data/PUMCH_Cohort/survival-info.txt", header = T, sep = "\t", check.names = F)
rownames(surv.info) <- surv.info$PDX_ID

sam.anno$surv_time <- surv.info[sam.anno$Model_NO, ]$Survival_time_month

pre.data <- read.table("./Data/PUMCH_Cohort/gene_expression.preTreat.GRCH37.txt", 
                       header = T, check.names = F)

post.data <- read.table("./Data/PUMCH_Cohort/gene_expression.postTreat.GRCH37.txt",
                        header = T, check.names = F)

post.expr <- post.data[ , grep("TPM$", colnames(post.data))]
pre.expr <- pre.data[ , grep("TPM$", colnames(pre.data))]

colnames(post.expr) <- substr(colnames(post.expr) , 5, 10)
colnames(pre.expr) <- substr(colnames(pre.expr) , 5, 10)
colnames(post.expr) <-  gsub('-', '', colnames(post.expr))
colnames(pre.expr) <-  gsub('-', '', colnames(pre.expr))

patient.id <- intersect(rownames(sam.anno)[which( (! is.na(sam.anno$surv_time)) & sam.anno$surv_time > 2)], colnames(pre.expr))
post.expr <- post.expr[, patient.id]
pre.expr <- pre.expr[, patient.id]
anno <- sam.anno[patient.id, ]

genes.delta <- log2((post.expr+1) / (pre.expr+1))

rownames(pre.expr) <- pre.data$gene_id
rownames(post.expr) <- pre.data$gene_id
rownames(genes.delta) <- pre.data$gene_id



# Survival analysis for GDF15 logFC
options(repr.plot.width=3.5, repr.plot.height=4.3)

gene <- 'GDF15'
gene.id <- pre.data$gene_id[match(gene, pre.data$gene_name)]
gene.id
# replace genes.delta with  pre.expr to obtain the relationship between pre-treatment expression of GDF15 and survival periods
surv.df <-data.frame(Patient = colnames(genes.delta), Log2FC = as.vector(unlist(genes.delta[gene.id, ])))
surv.df$surv_time <- sam.anno[colnames(genes.delta), 'surv_time']
surv.df$response <- sam.anno[colnames(genes.delta), ]$condition
surv.df$Drug <- sam.anno[colnames(genes.delta), ]$drug

surv.df$group <- "Median"
surv.df$group[surv.df$Log2FC > quantile(surv.df$Log2FC, 0.75)] <- "High"
surv.df$group[surv.df$Log2FC < quantile(surv.df$Log2FC, 0.25)] <- "Low"
surv.df$group <- factor(surv.df$group, levels = c('High', 'Median', 'Low'), ordered = T)
surv.df$surv_status <- 2

surv_object <<- Surv(time = surv.df$surv_time, event = surv.df$surv_status)
fit <- survfit(surv_object ~ group, data =  surv.df)
survdiff(formula = Surv(surv_time, surv_status) ~ group, data = surv.df)


p.surv <- ggsurvplot(fit, pval = F,
                     legend.labs=c("High" ,"Median",  "Low"),
                     palette = c("#DE77AE","#FFB90F", "#7FBC41"),
                     ggtheme = theme_bw(),
                     legend.title = "log2FC", 
                     conf.int = F,
                     risk.table = TRUE,
                     tables.height = 0.25,
                     tables.theme = theme_cleantable())
p.surv$plot <- p.surv$plot + 
  xlab('Months')+
  theme(
    legend.position=c(0.52, 0.78),
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

options(repr.plot.width=3.6, repr.plot.height=4.2)
print(p.surv)


# obtain p values
surv.df2 <- surv.df[which(surv.df$group != 'High'), ]
surv_object <<- Surv(time = surv.df2$surv_time, event = surv.df2$surv_status)
fit <- survfit(surv_object ~ group, data =  surv.df2)
survdiff(formula = Surv(surv_time, surv_status) ~ group, data = surv.df2)
p.surv <- ggsurvplot(fit, pval = T,  ggtheme = theme_bw(), conf.int = F, risk.table = TRUE,  tables.height = 0.25, tables.theme = theme_cleantable())
p.surv


surv.df2 <- surv.df[which(surv.df$group != 'Low'), ]
surv_object <<- Surv(time = surv.df2$surv_time, event = surv.df2$surv_status)
fit <- survfit(surv_object ~ group, data =  surv.df2)
survdiff(formula = Surv(surv_time, surv_status) ~ group, data = surv.df2)
p.surv <- ggsurvplot(fit, pval = T,  ggtheme = theme_bw(), conf.int = F, risk.table = TRUE,  tables.height = 0.25, tables.theme = theme_cleantable())
p.surv

