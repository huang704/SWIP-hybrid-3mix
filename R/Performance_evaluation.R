library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pROC)
library(RColorBrewer)
library(MSstats)

methods <- c("moderatedCMS", "moderatedCS", "moderatedCM", "moderatedC", "limmaCS", "limmaCMS")
res <- list()
count = 0

for(k in seq_along(methods)){
  count <- count + 1
  if(k <= 4){
    load(paste0('TMT-Swip-PD-by-', methods[k], ".testing-EB.rda"))
    data.res <-data.res$ComparisonResult
  }
  else{
    load(paste0('TMT-Swip-PD-by-', methods[k], ".testing.rda"))
  }
  data.res <- data.res %>% filter(!is.na(adj.pvalue))
  
  # MSstats::groupComparisonPlots(data=data.res,
  #                               type="VolcanoPlot",
  #                               logBase.pvalue=10,
  #                               width = 5,
  #                               height = 5,
  #                               address=paste0(methods[k], "-"))
  
  data.res <- data.res[, c("Protein", "Label", "log2FC", "SE", "DF", "pvalue",  "adj.pvalue")]
  
  data.res <- as.data.frame(data.res)
  data.res$Label <- as.character(data.res$Label)
  data.res$Protein <- as.character(data.res$Protein)
  data.res$Method <- methods[k]
  res[[count]] <- data.res
}


res1 <- data.table::rbindlist(res)
res1$log2FC <- as.numeric(as.character(res1$log2FC))
res1$adj.pvalue <- as.numeric(as.character(res1$adj.pvalue))
res1 <- within(res1, Method <- factor(Method, levels = methods))

res1$pred <- ifelse(res1$adj.pvalue <= 0.05, 1, 0)

num_testable_prots <- res1 %>% group_by(Method, Label) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Method, n)
num_testable_prots

perf <- res1 %>% filter(adj.pvalue <= 0.05) %>% 
  group_by(Method, Label) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Label, n)

########################################################################
# Distribution of SE
blues_colors <- brewer.pal(n = 9, name = "Blues")
purple_colors <- brewer.pal(n = 9, name = "PuRd")
colors <- c(purple_colors[9], blues_colors[6], blues_colors[8], purple_colors[7], purple_colors[4])

sub_res1 <- res1 %>% filter(Label %in% c("Mutant-Control", "Mutant.F4-Control.F4", "Mutant-ControlF10-F4") & Method != "moderatedCMS")

sub_res1[sub_res1$Label == "Mutant-Control", "Label"] <- "MUT-WT"
sub_res1[sub_res1$Label == "Mutant.F4-Control.F4", "Label"] <- "MUT_BF1-WT_BF1"
sub_res1[sub_res1$Label == "Mutant-ControlF10-F4", "Label"] <- "(MUT_BF7-MUT_BF1)-(WT_BF7-WT_BF1)"

sub_res1$Label <- factor(sub_res1$Label, levels = c("MUT-WT", "MUT_BF1-WT_BF1", "(MUT_BF7-MUT_BF1)-(WT_BF7-WT_BF1)"))

pdf("Compare_different_methods_se.pdf",width=8,height=5)
g4 <- sub_res1 %>% 
  ggplot(aes(x=Label, y= SE, color= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_color_manual(values=colors) + 
  ylim(0,0.3)+
  labs(x=(""), y = ("Standard error"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=30, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'right')

g4
dev.off()

pdf("Compare_different_methods_df.pdf",width=8,height=5)
g4 <- sub_res1 %>% 
  ggplot(aes(x=Label, y= DF, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_fill_manual(values=colors) +
  labs(x=(""), y = ("Degree of freedoms"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=30, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'right')

g4
dev.off()