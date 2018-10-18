################################################################################
# figure 3 #
# CDR3 hydropathy #
################################################################################
# fig 3A#
## add properties
addProperties.noepi <- function(vdj){
  library("alakazam")
  library("Peptides")
  
  # add the column of boman value
  vdj$boman <- sapply(vdj$CDR3, boman)
  
  # add 10 columns contataining 10 Kidera factors (KF)
  kf <- sapply(vdj$CDR3, kideraFactors)
  kf.df <- as.data.frame(t(matrix(unlist(kf), nrow = length(kf[[1]]))))
  colnames(kf.df) <- paste0("KF", 1:10)
  vdj.kf <- cbind(vdj, kf.df)
  
  # add 8 properties calculated from amino-acid contaning CDR3 seq
  vdj.prop <- aminoAcidProperties(vdj.kf, property = c("gravy", "bulk", "aliphatic",
                                                       "polarity", "charge", "basic", "acidic", "aromatic"), "CDR3")
  return(vdj.prop)
}

colnames(ncd8b.nr)[8] <- "CDR3"
colnames(ncd8a.nr)[8] <- "CDR3"
ncd8b.nr <- addProperties.noepi(ncd8b.nr)
ncd8a.nr <- addProperties.noepi(ncd8a.nr)

# merge non-naive and naive data of alpha and beta chains
mergeNaiveChain <- function(nna, nnb, na, nb, parm){
  col <- c("CDR3", parm)
  mergeNaive <- do.call("rbind", list(nna[,col], nnb[, col], na[, col], nb[, col]))
  mergeNaive$dataset <- c(rep("Non-naive", (dim(nna)[1] + dim(nnb)[1])), 
                          rep("Naive", (dim(na)[1] + dim(nb)[1]))) 
  mergeNaive$chain <- c(rep("alpha", dim(nna)[1]), rep("beta", dim(nnb)[1]),
                        rep("alpha", dim(na)[1]), rep("beta", dim(nb)[1]))
  return(mergeNaive)
}

# generate merge dataset
cd8ab.n <- mergeNaiveChain(cd8a.nr, cd8b.nr, ncd8a.nr, ncd8b.nr, "CDR3_AA_GRAVY")

################################################################################
# violin plot naive vs non-naive of CDR3 alpha and beta chains
F3A <- ggplot(cd8ab.n, aes(x = dataset, y = CDR3_AA_GRAVY, fill = dataset)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "lightgray") +
  geom_violin(trim = F, alpha = 0.4, adjust = 2) +
  labs(x = " ", y = "CDR3 hydropathy") +
  scale_fill_manual(values = c("#6E016B", "lightgray")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 10) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  facet_wrap(~chain, scales = "free", labeller = label_parsed) +
  theme(legend.position = "None",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm"), size = 25, face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"))
################################################################################
F3A <- F3A + 
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# fig 3B #
# compare the CDR3 hydropthy between ID vs SD for top 6 HLA
################################################################################

mhc <- names(sort(summary(cd8b.nr$MHC.A.new), decreasing = TRUE))[1:6]
t.mhc <- droplevels(cd8b.nr[cd8b.nr$MHC.A.new %in% mhc, ])


t.mhc$MHC.locus <- c(ifelse(t.mhc$MHC.A.new %in% c("HLA-A*01", "HLA-A*02", "HLA-A*03", "HLA-A*11"), 
                            "A", "B"))

t.cd8 <- cd8b0.01nr
t.cd8$MHC.locus <- c(ifelse(cd8b.nr$MHC.A.new %in% c("HLA-A*01", "HLA-A*02", "HLA-A*03", "HLA-A*11"), 
                            "A", "B"))

# violin plot
violinHydropathy <- function(df){
  ggplot(df, aes(x = dominant, y = CDR3_AA_GRAVY)) + # x= dominant
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
    geom_violin(aes(fill = dominant), trim = F, alpha = 0.4, adjust = 2) + # MHC.locus to dominant
    labs(x = " ", y = "CDR3 hydropathy") + #HLA locus
    scale_fill_manual(values = c("lightblue", "darkorange2")) +
    scale_y_continuous(limits = c(-2.5, 2.5)) +
    custom_theme +
    stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                       method = "wilcox.test", size = 8) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8) +
    facet_wrap(~MHC.A.new, scales = "free") + #~MHC.A.new for the MHC, ~Epitope 
    theme(legend.position = "None",
          strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
          strip.background = element_rect(colour = "white", fill = "white"))
}
################################################################################
# plot
F3B <- violinHydropathy(t.mhc) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# fig 3C #
# compare ID vs SD in top 3 epitopes
################################################################################

a2epitope <- names(sort(summary(droplevels(pull(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02", 
                                                           "Epitope"]))), decreasing = TRUE))[1:6]
# plot
F3C <- ggplot(droplevels(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02" & 
                               cd8b.nr$Epitope %in% a2epitope[1], ]), 
       aes(x = dominant, y = CDR3_AA_GRAVY)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = dominant), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = " ", y = "CDR3 hydropathy") +
  scale_fill_manual(values = c("lightblue", "darkorange2")) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  facet_wrap(~Epitope, scales = "free") +
  theme(legend.position = "None",
        axis.text.x = element_text(angle = 0),
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white"))

################################################################################
F3C <- F3C +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))
################################################################################

## combine the plots ##
library(grid)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

# Arrange the plots
print(F3A, vp = define_region(row = 1, col = 1:2))   # Span over two columns
print(F3B, vp = define_region(row = 2, col = 1:3))
print(F3C, vp = define_region(row = 1, col = 3))




























