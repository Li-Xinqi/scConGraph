library(Matrix)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(rgl)
library(car)
library(ggsci)


setwd('F:/PDAC_writing/github/upload/')
flow.meta <- read.csv("F:/PDAC_writing/github/upload/Flow_information.csv", header = T, check.names = F)
colnames(flow.meta) <- c("Flow", "Sample", "Source", 
                         "Target", "Probability_of_Flow_in_Control_Cluster",  "Probability_of_Flow_in_Treated_Cluster", 
                         "Percent_of_Flow_in_Control_Sample", "Percent_of_Flow_in_Treated_Sample", 
                         "STC", "RAC_score", "AAC_score", "Drug_Response_Type" )
table(flow.meta$Drug_Response_Type)

flow.meta$color <- 'black'
flow.meta$color[which(flow.meta$Drug_Response_Type == "Acquired Resistance")] <- '#E377C2'
flow.meta$color[which(flow.meta$Drug_Response_Type == "Absolutely Intrinsic Resistance")] <- '#FFC300'
flow.meta$color[which(flow.meta$Drug_Response_Type == "Relatively Intrinsic Resistance")] <- '#EE2C2C'
flow.meta$color[which(flow.meta$Drug_Response_Type == "Sensitivity")] <- '#1ABC9C'

flow.meta$Drug_Response_Type <- factor(flow.meta$Drug_Response_Type, levels = c('Acquired Resistance', 
                                                          'Absolutely Intrinsic Resistance', 
                                                          'Relatively Intrinsic Resistance', 
                                                          'Sensitivity'), ordered = T)



colors <- pal_lancet("lanonc", alpha = 0.8)(9)
flow.meta$Patient_colors <- flow.meta$Sample
flow.meta$Patient_colors[which(flow.meta$Sample == 'P1')] <- colors[1]
flow.meta$Patient_colors[which(flow.meta$Sample == 'P2')] <- colors[2]
flow.meta$Patient_colors[which(flow.meta$Sample == 'P3')] <- colors[3]
flow.meta$Patient_colors[which(flow.meta$Sample == 'P4')] <- colors[4]
flow.meta$Patient_colors[which(flow.meta$Sample == 'P5')] <- colors[5]
flow.meta$Patient_colors[which(flow.meta$Sample == 'P6')] <- colors[6]

###----------------------------------------------------------------------------------###
# Figure 2D
# 3D scatter
library(scatterplot3d)
STC <- flow.meta$STC
AAC_score <- flow.meta$AAC_score
RAC_score <- flow.meta$RAC_score


pdf('3d_scatter_plot.pdf',onefile=TRUE,width=4.8,height=4.5)

scatterplot3d( x= STC, y= log2(RAC_score+1), z = log2(AAC_score+1),
               type='p',color = flow.meta$color,col.axis= 'grey', 
               highlight.3d=F,angle=210,grid=T,box=T,scale.y=1,
               cex.symbols=2,pch=20,col.grid='grey')
dev.off()



# legend
pdf('3d_scatter_plot_legend.pdf',onefile=TRUE,width=5,height=4)
plot.new()
legend(x = "center", title.adj = 0.09,
       bg = rgb(1, 1, 1, alpha = 1),
       cex = 1,title = 'Drug response', 
       c('Acquired Resistance', 
         'Absolutely Intrinsic Resistance', 
         'Relatively Intrinsic Resistance', 
         'Relative Sensitivity'),
       #fill=c( '#A852CC', '#EE2C2C', "#FF7F24","#FFC300"),box.col="grey")
       fill=c("#E377C2", '#EE2C2C', "#FFC300", "#1ABC9C"),box.col="grey")
dev.off()



###----------------------------------------------------------------------------------###
# Figure 2E

pdf('3d_scatter_plot_sample.pdf',onefile=TRUE,width=4.8,height=4.5)

scatterplot3d( x= STC, y= log2(RAC_score+1), z = log2(AAC_score+1),
               type='p',color = flow.meta$Patient_colors,col.axis= 'grey', 
               highlight.3d=F,angle=210,grid=T,box=T,scale.y=1,
               cex.symbols=2,pch=20,col.grid='grey', xlab = '', ylab = '', zlab = '')
dev.off()



pdf('3d_scatter_plot_sample_legend.pdf',onefile=TRUE,width=3,height=4)
plot.new()
legend(x = "center",
       bg = rgb(1, 1, 1, alpha = 1),
       #inset = c(-0.35, 0),
       paste0('P', 1:6),
       fill=colors[1:6],box.col="grey")
dev.off()


###---------------------------------------------------------------------------------------###
# Figure 2F
drtype <- c('Acquired Resistance', 
            'Absolutely Intrinsic Resistance', 
            'Relatively Intrinsic Resistance', 
            'Sensitivity')
df <- data.frame(Sample = c(), Drug_Response_Type = c(), Percent = c())
for (pa in paste0('P', 1:6)){
  for (dr in drtype){
    sub <- flow.meta[which(flow.meta$Sample == pa & flow.meta$Drug_Response_Type == dr), ]
    if (nrow(sub) > 0){
      ctrl.perct <- sum(sub$Percent_of_Flow_in_Control_Sample)
      gemc.perct <- sum(sub$Percent_of_Flow_in_Treated_Sample)
      df.sub <- data.frame(Sample = paste0(pa, c('_Ctrl', '_Gemc')), 
                           Drug_Response_Type = rep(dr, 2), Percent = c(ctrl.perct, gemc.perct))
      df <- rbind(df, df.sub)
    }else{
      df.sub <- data.frame(Sample = paste0(pa, c('_Ctrl', '_Gemc')), 
                           Drug_Response_Type = rep(dr, 2), Percent = rep(0, 2))
      df <- rbind(df, df.sub)
    }
    
  }
}

df$Drug_Response_Type <- factor(df$Drug_Response_Type, levels = c('Acquired Resistance', 
                                            'Absolutely Intrinsic Resistance', 
                                            'Relatively Intrinsic Resistance', 
                                            'Sensitivity'), ordered = T)
df$Condition <- substr(df$Sample, 4, 10)
df$Patient <- substr(df$Sample, 1, 2)


p1 <- ggplot(df, aes(x = Condition, y = Percent*100, fill = Drug_Response_Type)) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "", y = "Percent", fill = "Drug response") +
  theme_bw() +
  scale_fill_manual(values=c("#E377C2", '#EE2C2C', "#FFC300", "#1ABC9C"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=12), 
        legend.title = element_text(size=12))+
  facet_wrap(~Patient,ncol = 6)
p1

pdf('Barplot_Drug_response_proportion.pdf',onefile=TRUE,width=6,height=4)
print(p1)
dev.off()
