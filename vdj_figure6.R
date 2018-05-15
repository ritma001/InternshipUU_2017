plotAffinity_epi <- function(df, title.label){
  library(ggplot2)
  library(ggpubr)
  # remove repetitive CDR3 across individual
  df <- df[!duplicated(df[, "CDR3"]), ]
  
  p <- ggplot(df, aes(x = dominant, y = affinity, fill = dominant)) +
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white", color = "lightgray") +
    geom_violin(trim = F, alpha = 0.4, adjust = 2) +
    labs(x = "Immune response", y = "Binding energy (kB)", title = title.label) +
    scale_fill_manual(name = "", values = c("lightblue", "darkorange2")) +
    custom_theme +
    theme(legend.position = "none",
          axis.text = element_text(face = "bold", size = 15), 
          plot.title = element_text(face = "bold", size = 20)) +
    stat_compare_means(aes(label = ..p.signif..), hide.ns = TRUE, #comparisons = list(c("ID", "SD")),
                       method = "wilcox.test", paired = FALSE,
                       label.x = 1.4, label.y.npc = 'top', size = 10)
  return(p)
}

# group the plot of GIL, DEN1 and DEN34
library(ggarrange)

gilAffinity <- plotAffinity_epi(a2GIL, "IAV") + theme(legend.position = "top")
den1Affinity <- plotAffinity_epi(a2DEN1, "DEN1") + theme(axis.title = element_blank())
den34Affinity <- plotAffinity_epi(a2DEN34, "DEN3/4") + theme(axis.title = element_blank())


ggarrange(gilAffinity, ggarrange(den1Affinity, den34Affinity, nrow = 2, labels = c("B", "C")), ncol = 2, 
          labels = "A")
