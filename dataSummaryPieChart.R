################################################################################
# figure 1 #
# VDJ db summary: Absolute count of CDR3 per virus species and MHC
################################################################################
# count and combine rare HLA (frequency < 1.0%) as "Others"
countMHC <- function(df){
  invisible(lapply(c("tidyr", "dplyr"), library, character.only = TRUE))
  # summarise data
  sumDat <- df %>% 
    group_by(MHC.A.new) %>% 
    dplyr::summarise(count = n()) %>% 
    mutate(frequency = round(count*100/sum(count), 1))
  
  # replace MHC.A.new as others if the frequency < 1%
  sumDat$MHC.A.new <- as.character(sumDat$MHC.A.new)
  sumDat[which(sumDat$frequency < 1, arr.ind = TRUE), "MHC.A.new"] <- "Others"
  return(sumDat)
}


barPlot.hla <- function(df){
  # import necessary packages
  invisible(lapply(c("tidyr", "dplyr", "reshape", "ggplot2", "ggpubr", "RColorBrewer"), 
                   library, character.only = TRUE))
  
  # indentidy non-rare MHC
  MHC <- dplyr::pull(countMHC(df), "MHC.A.new")
  
  # find relative frequency of ID and SD in each MHC subset
  t <- df %>% 
    group_by(Epitope.species, MHC.A.new) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(total = sum(count))
  
  t$MHC.A.new <- as.character(t$MHC.A.new)
  t[!t$MHC.A.new %in% MHC, "MHC.A.new"] <- "Others"
  
  levels(t$Epitope.species) <- c(levels(t$Epitope.species), "YFV", "IAV")
  t$Epitope.species[t$Epitope.species == 'YellowFeverVirus'] <- 'YFV'
  t$Epitope.species[t$Epitope.species == 'InfluenzaA'] <- 'IAV'
  
  # sort factor level by total count of the MHC to control the colors
  mhcCount <- aggregate(t["count"], t["MHC.A.new"], sum)
  reorderLevel <- mhcCount[with(mhcCount, order(count, decreasing = T)), "MHC.A.new"]
  t$MHC.A.new <- factor(t$MHC.A.new, levels = reorderLevel)
  
  p <- ggplot(t, aes(x= reorder(Epitope.species, -total), y =  count)) +
    geom_bar(stat = "identity", aes(fill = MHC.A.new)) + 
    geom_text(aes(Epitope.species, total + 70, label= total), size = 7) +
    labs(x = "Virus species", y = "Absolute count of CDR3 sequences") +
    scale_fill_brewer(name = "MHC alleles", palette = "Set3") +
    custom_theme + 
    theme(axis.text.x = element_text(angle = 45),
          legend.position = "right")
  
  return(p)
}
################################################################################
fig1 <- barPlot.hla(cd8b.nr)

################################################################################
# figure 2 #
# proportion of MHC alleles in VDJ db #
################################################################################
pie3D.hla <- function(df){
  # import necessary packages
  invisible(lapply(c("tidyr", "dplyr", "reshape", "plotrix", "RColorBrewer"), 
                   library, character.only = TRUE))
  
  # format the data
  t <- df %>% 
    group_by(MHC.A.new) %>% 
    dplyr::summarise(count = n()) %>% 
    mutate(frequency = round(count*100/sum(count), 1))
  
  # replace MHC.A.new as others if the frequency < 1%
  t$MHC.A.new <- as.character(t$MHC.A.new)
  t[which(t$frequency < 1, arr.ind = TRUE), "MHC.A.new"] <- "Others"
  t <- aggregate(t[,2:3], t[,1], sum)
  
  # sort frequency column in a decreasing order
  t <- t[order(t$frequency, decreasing = TRUE),]
  
  # plot pie 3D chart
  pie <- pie3D(t$frequency, #main = "Proportion of the different MHC molecules", 
               explode= 0.15, theta = 0.7, radius = 1.6, start = 0, 
               mar = c(2.5,2.5,2.5,2.5), border = "gray", font = 2,
               col = brewer.pal(12, "Set3")[1:dim(t)[1]])
  
  pie.pos <- pie # a vector of label position from pie
  pie.pos[8] <- 6.0
  pie.pos[9] <- 6.20
  pie.pos[10] <- 6.40
  pie3D.labels(pie.pos, labels = paste(t$MHC.A.new, "\n", paste0(t$frequency,"%")), 
               labelcex = 1.1, labelrad = 1.8, minsep = 0.3) #1.1 1.5fig1
  return(pie)
}
################################################################################
pie <- pie3D.hla(cd8b.nr)

################################################################################
# arrang Fig1A and 1B
rl <- lapply(list("1A_vdj_epitope_count.png", "2A_vdj_mhc_count.png"), png::readPNG)
gl <- lapply(rl, grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl)

################################################################################
Fig1A <- png::readPNG("1A_vdj_epitope_count.png")
Fig2A <- png::readPNG("2A_vdj_mhc_count.png")
F1A <- gridExtra::arrangeGrob(grid::rasterGrob(Fig1A), top = textGrob("A", x = unit(0.3, "npc")
                                    , y   = unit(0, "npc"), #just=c("left","top"),
                                    gp=gpar(col="black", fontsize=18)))


F1B <- gridExtra::arrangeGrob(grid::rasterGrob(Fig2A), top = textGrob("B", x = unit(0.3, "npc")
                                                    , y   = unit(0, "npc"), #just=c("left","top"),
                                                    gp=gpar(col="black", fontsize=18)))

gridExtra::grid.arrange(F1A, F1B, ncol = 1)
################################################################################



















