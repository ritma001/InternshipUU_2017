################################################################################
# figure 2 #
# CDR3 length analysis #
################################################################################
# fig 2A #
# compare alpha and beta length
################################################################################
summariseDatI <- function(df){
  # generate summary of CDR3 length of CD8+ T cells of vdj, naive and public CDR3
  length.df <- summariseDat(df, NA, "CDR3.length") # change Vseq to nonVJ length or else
  colnames(length.df)[1] <- "dataset"
  length.df[is.na(length.df)] <- deparse(substitute(df))
  return(length.df)
}


plotLengthAB <- function(sum.dat, fill.var, fill.name, legend.label, psignif){
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  
  p <- ggplot(data = sum.dat, aes(x= as.factor(CDR3.length), y = frequency)) + 
    geom_bar(aes_string(fill = fill.var), stat = "identity", position = "dodge") + 
    labs(x = "CDR3 length (AA)", y = "Relative frequency (%)") +
    scale_fill_manual(name = fill.name, labels = legend.label,
                      values = c("navy", "orange")) + 
    scale_y_continuous(limits = c(0, max(sum.dat$CDR3)), expand = c(0,0)) +
    theme(legend.position = c(0.8,0.8)) +
    custom_theme +
    annotate("text", x = 8, y = max(sum.dat$frequency), 
           label = psignif, size = 10)
  
  return(p)
}

################################################################################
# test the distribution using Komogorov-Smirnov test then add p-signif to the plot
ks.test(cd8a.nr$CDR3.length, cd8b.nr$CDR3.length)

# plot and add p-signif from ks
F2A <- plotLengthAB(rbind(summariseDatI(cd8a.nr),summariseDatI(cd8b.nr)), 
                    "dataset", "TCR chain", c(expression(paste("   ", alpha)), 
                                              expression(paste("   ", beta))), "****") +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# fig 2 B #
# Violin plot of CDR3 length ID vs SD
################################################################################
violinLength <- function(df){
  library(ggplot2)
  library(ggpubr)
  p <- ggplot(df, aes(x = dominant, y = CDR3.length, fill = dominant)) + 
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white", color = "gray") +
    geom_violin(trim = F, alpha = 0.4, adjust = 2) +
    
    labs(x = "Immune response", y = "CDR3 length (AA)") +
    scale_fill_manual(values = c("lightblue", "darkorange1")) +
    custom_theme +
    theme(legend.position = "None") +
    stat_compare_means(aes(x = dominant, y = CDR3.length, label = ..p.signif..),
                       data = df, label.x = 1.5, method = "wilcox.test", size = 10,
                       hide.ns = TRUE) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8, hjust = 0.2)
  return(p)
}

################################################################################
F2B <- violinLength(cd8b.nr) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))


################################################################################
# fig 2C #
# distribution of non-VJ length ID vs SD
################################################################################
plotLength.nonvj <- function(vdj, fill.var, fill.name){
  invisible(lapply(c("tidyr", "dplyr", "ggplot2"), library, character.only = TRUE))
  
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either ID or SD
  t <- vdj %>%
    group_by_(fill.var, ~nonVJ.length) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(rel.freq = n*100/sum(n)) %>%
    filter(n > 0)
  
  df <- as.data.frame(t)
  names(df)[1] <- "variable" 
  
  df <- complete(df, variable, nonVJ.length, fill = list(n = 0, rel.freq = 0))
  
  # suppress the warning message from the color extrapolation
  options(warn = -1)

  # plot the data 
  p <- ggplot(data = df, aes(x= as.factor(nonVJ.length), y = rel.freq, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(#title = expression(paste("Distribution of non-VJ length of CD8+T-cell ", beta, " chain")),
         x = "Non-VJ length (AA)", y = "Relative frequency (%)") +
    scale_fill_manual(name = fill.name, values = c("lightblue", "darkorange1")) +
    scale_y_continuous(breaks = seq(0, max(df$rel.freq), 5),
                       limits = c(0, max(df$rel.freq)), expand = c(0,0)) +
    custom_theme +
    theme(legend.position = c(0.8,0.6))
  # reset the warning message to the default value
  options(warn = 0)
  return(p)
}
################################################################################
library(tidyr)
library(dplyr)
# run ks test
ks.test(pull(cd8b.nr[cd8b.nr$dominant == "ID", "nonVJ.length"]), 
        pull(cd8b.nr[cd8b.nr$dominant == "SD", "nonVJ.length"]))
# plot 
F2C <- plotLength.nonvj(cd8b.nr, "dominant", "Immune response") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# fig 3D #
# compare the CDR3 length between non-naive (ID and SD combined) with the naive #
################################################################################

plotLength.summary <- function(sumDat, fill.var, fill.name, legend.label, psignif){
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  
  p <- ggplot(data = sumDat, aes(x= as.factor(CDR3.length), y = frequency)) + 
    geom_bar(aes_string(fill = fill.var), stat = "identity", position = "dodge") + 
    labs(#title = expression(paste("Distribution of CDR3 length of CD8+T-cell ", beta, " chain")),
      x = "CDR3 length (aa)", y = "Relative frequency (%)") +
    scale_fill_manual(name = fill.name, labels = legend.label,
                      values = c("#6E016B", "lightgray")) + # green "#E7B800", yellow "#E7B800"
    scale_y_continuous(limits = c(0, max(sumDat$CDR3)+2), expand = c(0,0)) +
    theme(legend.position = c(0.8,0.8)) +
    custom_theme +
    annotate("text", x = 8, y = max(sumDat$frequency), label = psignif, size = 10)
  return(p)
}


violinLength.summary <- function(sumDat){
  library(ggplot2)
  library(ggpubr)
  p <- ggplot(sumDat, aes(x = dataset, y = CDR3.length, fill = dataset)) + 
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white", color = "gray") +
    geom_violin(trim = F, alpha = 0.4, adjust = 2) +
    labs(x = "Cell population", y = "CDR3 length (AA)") +
    scale_fill_manual(values = c("#6E016B", "#B3DE69")) +
    custom_theme +
    theme(legend.position = "None") +
    stat_compare_means(aes(x = dataset, y = CDR3.length, label = ..p.signif..),
                       data = sumDat, label.x = 1.5, method = "wilcox.test", size = 10,
                       hide.ns = TRUE) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8, hjust = 0.2)
  return(p)
}

################################################################################
# summarise data #
t <-rbind(summariseDatI(cd8b.nr), summariseDatI(ncd8b.nr))
t$dataset <- ifelse(t$dataset == "cd8b.nr", "Non-naive", "Naive")

# frequency distribution #
ks.test(cd8b.nr$CDR3.length, ncd8b.nr$CDR3.length, alternative= "l")
  ## if significance meaning that the cumulative freq of non-naive > naive cd8 is 
  ## to the left 
F2D <- plotLength.summary(t ,"dataset", "Cell population", c("Non-naive", "Naive"), "**") +
  ggtitle("D") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))


# violin plot #
violinLength.summary(t)

# check distribution
length <- wideFormat(t, "dataset", "CDR3.length")

# nonVJ (AA)
ncd8b.ins <- ncd8b.nr[ncd8b.nr$nonVJ.length >= 3,]

# using mapply to pass getInsCodon to increse speed
InsCodecd8b <- mapply(getInsCodon, ncd8b.ins$NTseq, ncd8b.ins$Ins)
################################################################################

# extract only the added nucleotide sequences codon -> the codon combining between
# part of V/J with added nucleotides are excluded
getInsCodon <- function(seq,ins){
  if (!require(stringr)){ 
    install.packages(stringr) 
  } 
  # get start and end index of added nucleotide
  ind <- str_locate(as.character(seq), as.character(ins))
  if(is.na(ind[1])){
    return(NA)
  } else {
    # extract the begining of codon dereived from added nucleotides
    getStartCodon <- function(ind){
      if(ind > 3 & ind%%3 == 1){
        return(ind)
      } else {
        while(ind%%3 != 1){
          ind = ind + 1
        }
        return(ind)
      }
    }
    # extract codon sequence of added nucleotides
    subseq <- substr(seq, getStartCodon(ind[1]), ind[2])
    return(subseq)
  }
}


write.table(InsCodecd8b, file = "naive_pt.txt", 
            append = FALSE, quote = FALSE, sep = "\n",eol = "\n", 
            row.names = FALSE, col.names = FALSE)
# import python output and add as the new column
addedAmino <- read.table("D:/InternshipUU_2017/VDJdb/R_VDJdb/naive_pt_out.txt")
################################################################################
# compare the length of non-VJ
ks.test(cd8b.nr$nonVJ.length, sapply(as.character(addedAmino[,1]), nchar, USE.NAMES = F))
################################################################################
# fig 2 E #
# CDR3 length comparison between different epitope length #
################################################################################
plotLength.epitope <- function(df){
  lapply(c("ggplot2", "RColorBrewer", "ggpubr"), library, character.only = TRUE)
  ##fitLinear <- lm(scale(CDR3.length)~scale(Epitope.length)-1, data = df) # -1 to drop intercept 
  r <- round(cor(df$CDR3.length, df$Epitope.length, method = "spearman"),2)
  
  p <- ggplot(df, aes(x = Epitope.length, y = as.numeric(CDR3.length))) +
    stat_boxplot(aes(group = Epitope.length), geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(aes(group = Epitope.length), width = 0.1, fill = "white", color = "gray") +
    geom_violin(aes(group = Epitope.length), fill = "white", alpha = 0.5, 
                 color = "lightgray", trim = F, alpha = 0.4, adjust = 2, lwd = 1.0) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dotted", size = 1.2) +
    scale_x_continuous(breaks = c(8:11)) +
    labs(x = "Epitope length (AA)", y = "CDR3 length (AA)") +
    stat_summary(aes(group = Epitope.length, label=round(..y..,2)), fun.y=mean, geom="text", 
                 color = "black", size = 8) +
    custom_theme 
  
  p <- p + 
    annotate("text", x = 8:10, y = 23, label = "****", size = 8) +
    annotate("text", x = 11, y = 24, label = paste("r = ", r), size = 6, fontface = "italic")
  
  return(p)
}
################################################################################
t <- cd8b.nr[!duplicated(cd8b.nr[, "CDR3"]) & cd8b.nr$MHC.A.new != "HLA-A*02", ]

F2E <- plotLength.epitope(t) +
  ggtitle("E") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))
################################################################################

# combine all Figure 2 plots (2A to 2E)
library(ggpubr)

ggarrange(F2A, ggarrange(F2B, F2C, nrow = 2, labels = c("B", "C")), ncol = 2, labels = "A")

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
print(F2A, vp = define_region(row = 1, col = 1:2))   # Span over two columns
print(F2B, vp = define_region(row = 2, col = 1))
print(F2C, vp = define_region(row = 2, col = 2))
print(F2D, vp = define_region(row = 3, col = 1))
print(F2E, vp = define_region(row = 3, col = 2))


















































