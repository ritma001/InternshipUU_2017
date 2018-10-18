# fig 4A#
F4A <- plotAmino(cd8b.nr, "nonVJ") + 
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill=alpha('white', 0.5)),
        strip.background = element_rect(color = NA, fill = "white")) +
  facet_wrap(~hd, scales = "free") 
################################################################################
F4A <- F4A +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))
################################################################################
# count Amino acid and test the difference
t <- countAB(cd8b.nr[cd8b.nr$dominant == "ID",], cd8b.nr[cd8b.nr$dominant == "SD",])
ks.test(t[1, c("D", "E", "K", "N", "Q", "R")], t[2, c("D", "E", "K", "N", "Q", "R")])
ks.test(t[1, c("A", "C", "F", "I", "L", "M", "V", "W")], t[2, c("A", "C", "F", "I", "L", "M", "V", "W")])
ks.test(t[1, c("G", "H", "P", "S", "T", "Y")], t[2, c("G", "H", "P", "S", "T", "Y")])

################################################################################

logEnrichment <- function(t){
  # log comparisons
  tperc <- apply(t, 2, FUN = function(x) x*100/rowSums(t))
  enrichment <- log(tperc["df1",]/tperc["df2",])
  tperc <- as.data.frame(t(rbind(tperc, enrichment)))
  tperc$aa <- amino.acid
  ##return(tperc)
  # plot
  ggplot(tperc, aes(x = reorder(aa, +enrichment), y = enrichment)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Log enrichment ratio") +
    custom_theme
}

################################################################################
F4B <- logEnrichment(countAB(cd8b.nr, cd8a.nr)) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################

# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))

# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
# Arrange the plots
print(F4A, vp = define_region(row = 1, col = 1:2))   # Span over two columns
print(F4B, vp = define_region(row = 2, col = 1))



