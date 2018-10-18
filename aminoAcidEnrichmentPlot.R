################################################################################
# figure 5 #
# pLES #
################################################################################
## clustering ##
################################################################################
clusterpLES <- function(df){
  lapply(c("gplots", "RColorBrewer"), library, character.only =TRUE)
  # get the matrix
  m <- selectMSA(df, c(10:15), "pLES")
  # label CDR3 position
  colnames(m) <- c(as.character(104:111), "112.1", "111.1",as.character(112:118))
  # cluster #! exlclude low diversity column to prevent dist = NA so, the column 
  # dendogrm is constructed
  heatmap.2(m[, -c(1,2,17)], scale = "none", dendrogram = "both", na.color = "lightgray", 
            col= colorRampPalette(c("steelblue", "white", "red"))(n = 1000), 
            sepcolor = "gray", sepwidth=c(0.01,0.01),colsep=1:ncol(m), 
            rowsep=1:nrow(m), density.info = "none", trace = "none", key = T, 
            cexCol=1.2, cexRow=1.2, srtCol = 0, key.title = "", key.xlab = "pLES", 
            key.par=list(cex=1.0), adjCol=c(0.5,0))
}

# add hclustfun= function(x) hclust(x, method="average") and change scale = "column" 
# for different presentation

clusterpLES_aa <- function(df){
  lapply(c("gplots", "RColorBrewer"), library, character.only =TRUE)
  # get the matrix
  m <- selectMSA(df, c(10:15), "pLES")
  # label CDR3 position
  colnames(m) <- c(as.character(104:111), "112.1", "111.1",as.character(112:118))
  # cluster #! exlclude low diversity column to prevent dist = NA so, the column 
  # dendogrm is constructed
  heatmap.2(m, Colv = NA, scale = "none", dendrogram = "row", na.color = "lightgray", 
            col= colorRampPalette(c("steelblue", "white", "red"))(n = 1000), 
            sepcolor = "gray", sepwidth=c(0.01,0.01),colsep=1:ncol(m), 
            rowsep=1:nrow(m), density.info = "none", trace = "none", key = T, 
            cexCol=1.2, cexRow=1.2, srtCol = 0, key.title = "", key.xlab = "pLES", 
            key.par=list(cex=1.0), adjCol=c(0.5,0))
}
################################################################################
F5A <- clusterpLES_aa(a2GIL)
F5B <- clusterpLES_aa(a2GLC)
F5C <- clusterpLES_aa(a2NLV)

library(ggseqlogo)
library(ggplot2)

F5D <- ggseqlogo(selectMSA(a2GIL, c(10:15), "pLES"),  
                 method='custom', seq_type = "aa") + 
  scale_y_continuous(limits = c(0, NA)) +
  invisible(scale_x_discrete(limits = c(as.character(104:111), "112.1", "111.1",
                                        as.character(112:118)))) +
  labs(title = " ") +
  custom_theme +
  theme(#axis.text.x = element_text(angle = 45),
    legend.position = "top",
    axis.ticks.x = element_line(color = "black"))

F5E <- ggseqlogo(selectMSA(a2GLC, c(10:15), "pLES"),  
                 method='custom', seq_type = "aa") + 
  scale_y_continuous(limits = c(0, NA)) +
  invisible(scale_x_discrete(limits = c(as.character(104:111), "112.1", "111.1",
                                        as.character(112:118)))) +
  labs(title = " ") +
  custom_theme +
  theme(#axis.text.x = element_text(angle = 45),
    legend.position = "none",
    axis.ticks.x = element_line(color = "black"))

F5F <- ggseqlogo(selectMSA(a2NLV, c(10:15), "pLES"),  
                 method='custom', seq_type = "aa") + 
  scale_y_continuous(limits = c(0, NA)) +
  invisible(scale_x_discrete(limits = c(as.character(104:111), 
                                        "112.1", "111.1",
                                        as.character(112:118)))) +
  labs(title = " ") +
  custom_theme +
  theme(#axis.text.x = element_text(angle = 45),
    legend.position = "none",
    axis.ticks.x = element_line(color = "black"))


################################################################################
library(grid)

# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))

# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
# Arrange the plots
print(F5A, vp = define_region(row = 1, col = 1))   # Span over two columns
print(F5B, vp = define_region(row = 2, col = 1))
print(F5C, vp = define_region(row = 3, col = 1))
print(F5D, vp = define_region(row = 1, col = 2))
print(F5E, vp = define_region(row = 2, col = 2))
print(F5F, vp = define_region(row = 3, col = 2))





F5D <- F5D + 
  ggtitle("E") +
  theme(plot.title = element_text(hjust = -0.02, vjust = 0.2, size = 20))

F5E <- F5E + 
  ggtitle("F") +
  theme(plot.title = element_text(hjust = -0.02, vjust = 0.2, size = 20))

F5F <- F5F + 
  ggtitle("G") +
  theme(plot.title = element_text(hjust = -0.02, vjust = 0.2, size = 20))


ggarrange(F5D, F5E, F5F, nrow = 3)


# arrang Fig1A and 1B
rl <- lapply(list("F5A_gil.png", "F5B_glc.png", "F5C_nlv.png"), png::readPNG)
gl <- lapply(rl, grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl, ncol = 2)

################################################################################
# non VJ top 3 epitopes

t <- as.data.frame(do.call(cbind, sapply(c("GILGFVFTL", "NLVPMVATV", "GLCTLVAML"), function(epi) {
  cdr3.nr <- cd8b.nr[!duplicated(cd8b.nr[, "CDR3"]),]
  count <- countAmino(pull(unique(cdr3.nr[cdr3.nr$Epitope == epi & cdr3.nr$dominant == "ID", "nonVJ"])))$count
  list(count*100/sum(count))
})))

# alternative multiple ks test (Anderson-Darling) for the AA distributions fron the 3 epitopes
t.test <- t(t)
colnames(t.test) <- amino.acid

if(!require(kSamples)){
  install.packages(kSamples, dependencies = TRUE)
}

ad.test(t.test[1,], t.test[2,], t.test[3,])

################################################################################
# reshape data for plotting
t$amino.acid <- amino.acid
t <- reshape::melt(t, idvar = "amino.acid", value.name = "frequency")

 +# plot the amino acid count 
ggplot(t, aes(x = amino.acid, y = value,  fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Amino acid", y = "Relative frequency (%)") + 
  scale_fill_manual(name = "Epitope\n", labels = c("GIL", "GLC", "NLV"),
                    values = c("lightgray", "skyblue", "steelblue4")) +
  scale_y_continuous(breaks = seq(0, max(t$value), 5)) +
  custom_theme +
  theme(legend.position = c(0.8, 0.7))
################################################################################


