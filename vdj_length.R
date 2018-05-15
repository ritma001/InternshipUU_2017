################################################################################
#
# Author: Wannisa Ritmahan 
# Aim: CDR3 length analysis on the TCR derived from the VDJdb
# 
################################################################################
## custom theme ##
################################################################################
custom_theme <- theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "white"),
                      legend.position = c(0.8,0.6),
                      legend.key = element_rect(fill = "white"),
                      legend.background = element_rect(fill=alpha("white", 0.5)),
                      legend.text = element_text(size=17),
                      legend.title = element_text(size=18,face="bold"),
                      axis.text = element_text(size=17, face = "bold"),
                      axis.title = element_text(size=19,face="bold"),
                      axis.line = element_line(colour = "black"),
                      strip.text = element_text(size=18, face="bold"),
                      plot.title = element_text(size=16, face="bold", hjust = 0.5))

custom_theme2 <- theme(panel.background = element_rect(fill = "white"),
                      legend.position = c(0.2,0.7),
                      legend.key = element_rect(fill = "white"),
                      legend.background = element_rect(fill=alpha("white", 0.5)),
                      legend.text = element_text(size=17),
                      legend.title = element_text(size=18,face="bold"),
                      axis.text = element_text(size=17, face = "bold"),
                      axis.title = element_text(size=19,face="bold"),
                      axis.line = element_line(colour = "black"),
                      strip.text = element_text(size=18, face="bold"),
                      plot.title = element_text(size=16, face="bold", hjust = 0.5))


################################################################################
## extract p-value from permutation test ##
################################################################################
extractPvalue <- function(df){
  library(coin)
  # extract p-value from permutation test
  permutation <- independence_test(as.numeric(CDR3.length) ~ factor(MHC.A.new), 
                                   data = df)
  invisible(strOptions(permutation))
  pv.perm <- permutation@distribution@pvalue(permutation@statistic@teststatistic)
  return(pv.perm)
}

################################################################################
## plot CDR3 length distribution ##
################################################################################
plotLength <- function(vdj, fill.var, fill.name){
  invisible(lapply(c("tidyr", "dplyr", "ggplot2", "RColorBrewer"), library,  
         character.only = TRUE))
  
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
  
  # select the color regarding the number of variables
  color.value <- length(levels(factor(unlist((cd8b[, fill.var])))))
  if(color.value < 12 & color.value > 2){color <- brewer.pal(n = color.value, name = "Set3")
  } else if (color.value == 2){color <- c("dimgray", "darkorange")
    } else {color <- colorRampPalette(brewer.pal(color.value, "Set3"))(color.value)}
  
  # plot the data 
  p <- ggplot(data = df, aes(x= as.factor(nonVJ.length), y = rel.freq, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(title = expression(paste("Distribution of non-VJ length of CD8+T-cell ", beta, " chain")),
         x = "Non-VJ length (aa)", y = "Relative frequency (%)") +
    scale_fill_manual(name = fill.name, values = color) +
    scale_y_continuous(breaks = seq(0, max(df$rel.freq), 5),
                       limits = c(0, max(df$rel.freq)), expand = c(0,0)) +
    custom_theme +
    theme(legend.position = c(0.8,0.6))
  return(p)
  # reset the warning message to the default value
  options(warn = 0)
}

# call function #
plotLength(cd8b.nr, "dominant", "Immune response")
plotLength(cd8b.nr[cd8b.nr$Epitope.species == "HCV", ], "MHC.A.new", "HLA-A class")
plotLength(cd8b.nr[!(as.numeric(cd8b.nr$MHC.A.new) %in% which(table(cd8b.nr$MHC.A.new) < 46)),],
           "MHC.A.new", "")
plotLength(cd8b.nr[cd8b.nr$Epitope.species == "HCV",] %>% filter(frequency > 1), "Epitope.species", "Epitope species")
################################################################################
# boxplot of CDR3 length
################################################################################
boxplotLength <- function(vdj){
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  #vdj <- summariseDat(vdj, "dominant", "CDR3.length") %>% filter(count > 1)
  p <- ggplot(data = vdj, aes(x = factor(dominant), y = CDR3.length, fill = dominant)) + 
    geom_boxplot() +
    labs(title = "The CDR3 length of CD8+T-cells",
         x = "Immune response", y = "CDR3 length") +
    scale_fill_manual(name = "Immnue response", values = c("dimgray", "darkorange")) +
    theme(legend.position = "none") +
    custom_theme +
    stat_compare_means(size = 5, label.y = max(vdj$CDR3.length) + 2, 
                       label.x.npc = "center", method = 'wilcox.test', paired = F) + 
                       #method.args = list(alternative = "l", var.equal = F)) +
    stat_compare_means(aes(label = ..p.format..), comparison = list(c("ID", "SD")),
                       label.y = max(vdj$CDR3.length)) +
    stat_summary(fun.y=mean, geom="point",color = "darkgray") +
    stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", 
                 color = "black", size=5)
  return(p)
}

# call function #
t <- summariseDat(rbind(cd8a.nr,cd8b.nr), "Gene", "CDR3.length") %>% filter(count > 1)
boxplotLength(t)

names(cd8b.nr)[names(cd8b.nr) == "Gene"] <- "chain"


################################################################################
# compare alpha and beta chains
boxplotLength.CDR3 <- function(df){
  invisible(lapply(c("ggplot2", "ggpubr"), library, character.only = TRUE))
  p <- ggplot(df, aes(x = factor(Gene), y = as.numeric(CDR3.length), color = Gene)) +
    geom_jitter(size = 3, width = 0.3) +
    geom_boxplot(alpha = 0.5, width = 0.7, color = "darkgray") +
    scale_color_manual(values = c("navy", "orange2")) + 
    scale_x_discrete(labels=c("TRA" = expression(alpha), "TRB" = expression(beta))) + 
    scale_y_continuous(breaks = seq(min(df$CDR3.length), max(df$CDR3.length), 2)) +
    labs(x = "CD8+TCR chain", y = paste0("\n", "CDR3 length (aa)")) +
    custom_theme +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 25)) +
    stat_compare_means(aes(x = factor(Gene), y = as.numeric(CDR3.length),
                           label = ..p.signif..), label.x = 1.4, 
                       method = "wilcox.test", paired = F, size = 10) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8) +
    facet_grid(~Epitope.length)
  return(p)
}

################################################################################
# boxplot CDR3 length vs Epitope length
################################################################################

boxplotLength.epitope <- function(df){
  # import necessary packages
  lapply(c("ggplot2", "RColorBrewer", "ggpubr"), library, character.only = TRUE)
  
  # transform charcter to factor vector
  df$Epitope.length <- factor(df$Epitope.length)
  comparisons <- list(c("8", "9"), c("8", "10"), c("8", "11"), 
                       c("9", "10"), c("10", "11")) 
  # plot
  p <- ggplot(df, aes(x = Epitope.length, y = as.numeric(CDR3.length))) + 
    geom_jitter(width = 0.3, size = 3)+ #aes(color = MHC.A.new), , color = "orange2") + 
    geom_boxplot(alpha = 0.7,  width = 0.7, color = "darkgray") +
    geom_violin(trim = F , color = "gray", alpha = 0.7, width = 0.8, 
                scale = "width", size = 2) +
    labs(x = "Epitope length (aa)", y = "CDR3 length (aa)") + 
    custom_theme +
    theme(legend.position = "bottom") + 
    scale_color_manual(values = colorRampPalette(brewer.pal(color.value, "Set3"))(16)) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 7) +
    stat_compare_means(aes(label = ..p.signif..), label.y = max(df$CDR3.length), 
                       label.x = 2, size = 7) 
  
  return(p)
}

# call the function
violinPlotLength.hla(cd8b.nr)


################################################################################
# violin plot of CDR3 length
################################################################################

violinLength <- function(df){
  library(ggplot2)
  library(ggpubr)
  #df.p <- summariseDat(df, "dominant", "CDR3.length") %>% filter(count > 1)
 
  p <- ggplot(df, aes(x = dominant, y = CDR3.length, fill = dominant)) + 
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white") +
    geom_violin(trim = F, alpha = 0.4, adjust = 2) +
    
    labs(x = "Immune response", y = "CDR3 length (AA)") +
    scale_fill_manual(values = c("dimgray", "darkorange")) +
    custom_theme +
    theme(legend.position = "None") +
    stat_compare_means(aes(x = dominant, y = CDR3.length, label = ..p.signif..),
                       data = df, label.x = 1.5, method = "wilcox.test", size = 10,
                       hide.ns = TRUE) + #, 
                              #method.args = list(alternative = "l", var.equal = F)) +
    #facet_wrap(~MHC.A.new)
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8, hjust = 0.2)
  return(p)
}

################################################################################
# statistic test length
################################################################################
mean.diff <- function(df, parm){
  library(dplyr)
  sd.mean <- mean(df[df$dominant == "SD", parm])
  id.mean <- mean(df[df$dominant == "ID", parm])
  mean.diff <- sd.mean - id.mean
  return(mean.diff)
}

################################################################################
permutation <- function(df, n, iter, mean.diff){
  library(dplyr)
  # initialise the vector to keep mean difference in each iteration
  rand.dist <- c()
  
  for (i in seq(iter)){
    # set pseudo random number generator for reproducibility
    set.seed(1+i)
    rand1 <- sample(df, n, F)
    rand2 <- sample(df, n, F) #pull(sample_n(df, n, F))

    # calculate mean difference
    rand.dist <- c(rand.dist, mean(rand1)-mean(rand2))
  }
  # generate histogram of permutation 
  hist(rand.dist)
  # print permutation p value
  print(sum(abs(rand.dist) > abs(mean.diff))/iter)
}

cd8b.id <- cd8b[cd8b$dominant == "ID", ]
cd8b.id <- as.character(unique(cd8b.id$CDR3))
cd8b.sd <- cd8b[cd8b$dominant == "SD", ]
cd8b.sd <- as.character(unique(cd8b.sd$CDR3))

# compare id and sd unqiue CDR3
cd8b.unique <- data.frame(unique.cdr3 = c(cd8b.id, cd8b.sd), 
                             CDR3.length = sapply(c(cd8b.id, cd8b.sd), nchar),
                             dominant = c(rep("ID", length(cd8b.id)), rep("SD",length(cd8b.sd))))

meandiff <- mean.diff(cd8b.unique, "length")
permutation(cd8b.unique$length, 1000, 10000, meandiff)


# test of length (after proof of nalmality by shapiro)
t.test(cd8b.nr$CDR3.length, sample_n(ncd8b.nr, 5000)$CDR3.length, alternative = 'l')

# test of difference 
wilcox.test(cd8a.nr$CDR3.length, cd8b.nr$CDR3.length)
wilcox.test(cd8b.nr[cd8b$dominant == "ID", "CDR3.length"], 
            cd8b.nr[cd8b$dominant == "SD","CDR3.length"], alternative = "l")

wilcox.test(sapply(cd8b.nr[cd8b.nr$dominant == "ID", "Vseq"], nchar), sapply(cd8b.nr[cd8b.nr$dominant == "SD", "Vseq"], nchar))

# test of normality
library(data.table)
shapiro.test(sample_n(ncd8b.nr, 5000)$CDR3.length) # max for shapiro is 5000

################################################################################
# test if 2 data are from the same distribution
# one-sample KS test by comparing the ID value with the discrete distribution of SD length
library(dgof) # use ks.test from this librayry for a discrete dist comparison
ks.test(cd8b.nr[cd8b.nr$dominant == "ID", "CDR3.length"], 
        ecdf(cd8b.nr[cd8b.nr$dominant == "SD", "CDR3.length"]))


# test if the treatments have an effect on the outcomes
cd8bnr.length <- cd8b.nr %>%
  group_by_("dominant", ~CDR3.length) %>%
  dplyr::summarise(n = n()) %>%
  filter(n > 1) %>%
  mutate(rel.freq = n/sum(n))

################################################################################
## calculate count/perc  table
################################################################################

getCountTable <- function(parm1.list, parm2.list, format = "perc"){
  library(dplyr)
  library(plyr)
  parm1 <- summary(factor(pull(parm1.list)))
  parm2 <- summary(factor(pull(parm2.list)))
  
  if(format != "perc"){
    comb.parm <- rbind.fill.matrix(t(parm1), t(parm2))
  } else {
    # convert count to percentage
    parm1.perc <- parm1*100/sum(parm1)
    parm2.perc <- parm2*100/sum(parm2)
    # merge 2 vectors to matrix
    comb.parm <- rbind.fill.matrix(t(parm1.perc), t(parm2.perc))
  }
  # add zero count/perc to the data with no observation
  comb.parm[is.na(comb.parm)] <- 0
  rownames(comb.parm) <- c("A", "B")
  return(comb.parm)
}

## chisuare test of the count table oibtained from getCountTable
chisqTest <- function(countData){
  t <- cbind(t(countData), t(apply(t(countData), 1, function(x) {
    if(sum(x) == 0){
      ch <- c(5, 5)
    } else {
      ch <- chisq.test(x)#, correct = TRUE)#, simulate.p.value = TRUE) # to increase power in p-value estimation of low count
      c(unname(ch$statistic), ch$p.value)}
    })))
  
  colInd <- dim(t)[2]
  colnames(t)[c(colInd-1, colInd)] <- c("X-square statistic", "p-value")
  # add column of p significance
  t <- cbind(t, ifelse(t[, "p-value"] > 0.05, "", 
               ifelse(t[, "p-value"] < 0.05 & t[, "p-value"] > 0.01, "*", 
                      ifelse(t[, "p-value"] < 0.01 & t[, "p-value"] > 0.001, "**", "***"))))
  
  colnames(t)[colInd + 1] <- "p-signif"
  return(t)
}


## call fucntion ##
t <- getCountTable(cd8a.nr[, "CDR3.length"], cd8b.nr[, "CDR3.length"], "perc")
chisqTest(t)

################################################################################
# calculate the different frequency between alpha and beta
diffFrequency <- function(df1, df2, length){
  df1.freq <- sum(df1$CDR3.length == length)/dim(df1)[1]
  df2.freq <- sum(df2$CDR3.length == length)/dim(df2)[1]
  diffFreq <- abs(df1.freq - df2.freq)
  return(diffFreq)
}

# permutation length frequency 
permutationLength <- function(df, n, iter, length, mean.diff){
  library(dplyr)
  # calculate the actual difference of the input data
  
  # initialise the vector to keep mean difference in each iteration
  rand.dist <- c()
  
  for (i in seq(iter)){
    # set pseudo random number generator for reproducibility
    set.seed(1+i)
    rand1 <- sum(sample(df, n, F) == length, na.rm = TRUE) / n
    rand2 <- sum(sample(df, n, F) == length, na.rm = TRUE) / n
    # count only the sequence than calculate the frequency 
    
    # calculate mean difference
    rand.dist <- c(rand.dist, rand1-rand2)
  }
  # generate histogram of permutation 
  ##hist(rand.dist)
  # print permutation p value (2-way test)
  print(sum(abs(rand.dist) > abs(mean.diff))/iter)
}

################################################################################
## Suumarise data ##
################################################################################

# get summary of the data based on specified parm1 such as immune response 
# (dominant, Gene-chain) and parm2 (CDR3.length, amino acid) etc.

summariseDat <- function(vdj, parm1, parm2){
  # import necessary library
  library(tidyr)
  library(dplyr)
  # generate a dataset containing only the non-repetitive CDR3
  #vdj.nr <- excludeRepeat(vdj)
  # get summary of the data grouped by argument parm
  vdj.sel <- vdj %>%
    group_by_(parm1, parm2) %>%
    dplyr::summarise(count = n()) %>% 
    mutate(frequency = count*100/sum(count))
  
  # transform column and the fill the missing data with zero
  length.df <- data.frame(vdj.sel)
  # exclude sequence in which count == 1
  #length.df <- length.df[length.df$count != 1, ]
  return(length.df)
}

cd8abnr.length <- summariseDat(rbind(cd8a.nr, cd8b.nr), "Gene", "CDR3.length")
sumlen <- summariseDat(cd8b.nr, "dominant", "CDR3.length")
################################################################################

summariseDatI <- function(df){
  # generate summary of CDR3 length of CD8+ T cells of vdj, naive and public CDR3
  length.df <- summariseDat(df, NA, "Vseq") # change Vseq to nonVJ length or else
  colnames(length.df)[1] <- "dataset"
  length.df[is.na(length.df)] <- deparse(substitute(df))
 return(length.df)
}

summariseDatII <- function(df, length){
  # generate summary of CDR3 length of CD8+ T cells of vdj, naive and public CDR3
  length.df <- summariseDat(df, NA, length)
  colnames(length.df)[1] <- "dataset"
  length.df[is.na(length.df)] <- deparse(substitute(df))
  return(length.df)
}
cd8bnr.mergen <- rbind(summariseDatI(cd8b.nr), summariseDatI(ncd8b.nr))

cd8bnr.mergep <- rbind(summariseDatI(cd8b.nr), summariseDatI(pcd8b.nr))

#, summariseDatI(ncd8b.nr)
##do.call(rbind, sapply(list(cd8b.nr,ncd8b.nr,pcd8b.nr), summariseDatI, simplify = F))
################################################################################

plotLength.summary <- function(sum.dat, fill.var, fill.name, legend.label, t, y){
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  
  p <- ggplot(data = sum.dat, aes(x= as.factor(CDR3.length), y = frequency)) + 
    geom_bar(aes_string(fill = fill.var), stat = "identity", position = "dodge") + 
    labs(#title = expression(paste("Distribution of CDR3 length of CD8+T-cell ", beta, " chain")),
         x = "CDR3 length (aa)", y = "Relative frequency (%)") +
    scale_fill_manual(name = fill.name, labels = legend.label,
                      values = c("navy", "orange")) + 
    scale_y_continuous(limits = c(0, max(sum.dat$CDR3)), expand = c(0,0)) +
    # alpha color = c("#88419D",  "#66C2A5"), purple-green = c("#6E016B", "#B3DE69")
    # alpha/beta ("navy", "orange2")
    theme(legend.position = c(0.8,0.8)) +
    custom_theme # +
    annotate("text", x = c(rownames(t)), y = y + 0.2, label = t[, "p-signif"], size = 7)
    
  return(p)
}




# call function #
plotLength.summary(rbind(summariseDatI(cd8b.nr),summariseDatI(pcd8b.nr)),
                         "dataset",  expression(paste("CD8+T-cell ", beta, " chain")),
                   c("Unique", "Public"))

plotLength.summary(rbind(summariseDatI(cd8a.nr),summariseDatI(ncd8a.nr)), 
                   "dataset", expression(paste("CD8+T-cell ", beta, " chain")), 
                   c("Non-naive", "Naive"))

plotLength.summary(rbind(summariseDatI(cd4b.nr),summariseDatI(ncd4b.nr)), 
                   "dataset", "Cell population", c("Non-naive", "Naive"))

plotLength.summary(wideFormat(summariseDatI(cd8bnr), "dominant", "CDR3.length"), "dominant", 
                   "Immune response", c("ID", "SD"))
plotLength.summary(summariseDat(cd8b[cd8b$MHC.A.new == "HLA-A*02",], 
                                "dominant", "CDR3.length"),
                   "dominant", "Immune response", c("ID", "SD"))

################################################################################
wideFormat <- function(sum.long, parm1, parm2){
  library(reshape2)
  if(!require(fifer))
    {installed.packages("fifer")}
  # converse the result of summariseLength to wide format
  df.wide  <- dcast(data = sum.long, value.var = "count", 
                         formula = sprintf("%s ~ %s", parm1, parm2))
  # remove the column containg NA due to no count at certain length
  df.wide <-  df.wide[, apply(df.wide, 2, function(col) !any(is.na(col)))]
  
  # first column to rownames
  df.wide <- data.frame(df.wide[,-1], row.names = df.wide[,1])
  # assign proper column names (i.e. length)
  ##colnames(df.wide) <- sapply(colnames(df.wide), 
                                  ##function(string) as.numeric(strsplit(string, "\\D+")[[1]][-1]))

  return(df.wide)
}

################################################################################
# statistic test length frequency
################################################################################
# Chi-square test
chisq.test(wideFormat(merge.df, "dataset", "CDR3.length"))
chisq.test(wideFormat(rbind(summariseDatI(cd8a.nr),summariseDatI(ncd8a.nr)),"dataset", "CDR3.length"))


# Miscellaneous
# test color from Rcolorbrewer (# '9' = can be change to cover the length(brewer.pal(n = 9, name = "BuPu")))
plot(1:9, 1:9, col = brewer.pal(n = 9, name = "BuPu"), pch=19, cex=3, xlab="", ylab="")
brewer.pal(n = 9, name = "BuPu") 

################################################################################
# boxplot regarding MHC molecule and epitope
################################################################################
boxplotLengthFacet <- function(df){
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  
  p <- ggplot(data = df, aes(x = dominant, y = CDR3.length, fill = dominant)) + 
    geom_boxplot() +
    labs(title = "The CDR3 length of CD8+T-cells",
         x = "Immune response", y = "CDR3 length") +
    scale_fill_manual(name = "Immnue response", values = c("dimgray", "darkorange")) +
    theme(legend.position = "none") +
    custom_theme +
    stat_compare_means(size = 5, label.y = max(df$CDR3.length) + 2, 
                       label.x.npc = "center") +
    stat_compare_means(aes(label = paste("p-value", ..p.format..)), comparison = list(c("ID", "SD")),
                       label.y = max(df$CDR3.length)) +
    stat_summary(fun.y=mean, geom="point",color = "darkgray") +
    stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", 
                 color = "black", size=5)
  p <- p + facet_grid(~MHC.A.new)
  return(p)
}

################################################################################
# random sampling plot
################################################################################
randSampling <- function(vdj.chain, nsize, parm){
  library(dplyr)
  library(tidyr)
  
  # select the 80 CDR3 sequences from each bin of hydropathy.epi.scale
  t <- vdj.chain %>%
    group_by_(parm) %>%
    sample_n(size = nsize, replace = F)#do(sample_n(., nsize))
  #result.t <- t[, c("Epitope", "Epitope.species", "CDR3.length", 
                    #"Epitope.length", "dominant")]
  return(t)
}  

kfoldSampling <- function(vdj.chain, nsize, iter, parm){
  # to iterate randSampling and keep the sampled data in the same dataframe
  df.list = list()
  for (i in seq(iter)){
    # set pseudo random nyumber generator for a sake of reproducibility
    set.seed(1+i)
    # call random sampling and keep the result in the datalist
    rand.df<- randSampling(vdj.chain, nsize, parm)
    rand.df$iter <- i
    df.list[[i]] <- rand.df 
  }
  # generate the dataframe after iterations 
  df.final <- do.call(rbind, df.list)
  df.final$iter <- as.factor(df.final$iter)
  return(df.final)
}


plotLengthRand <- function(vdj.chain, nsize, kfold, parm){
  library(tidyr)
  library(dplyr)
  library(plotrix)
  library(ggplot2)
  # to plot length vs epitope length in the random dataset #
  vdj.rand <- kfoldSampling(vdj.chain, nsize, kfold, parm)
  # summarise the data and calculate mean, sd and std.err
  t <- vdj.rand %>% group_by(iter, Epitope.length) %>%  # ,dominant, Epitope.length) %>% 
    dplyr::summarise(mean=mean(CDR3.length), sd=sd(CDR3.length), 
                     se = std.error(CDR3.length))
 
  bp <- ggplot(t, aes(x = factor(Epitope.length), y = as.numeric(mean))) + 
    geom_jitter(aes(colour = iter), size = 3, position = position_dodge(width = 0.25)) +
    # geom_line(aes(colour = iter, group = iter), position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymax = as.numeric(mean) + se, ymin = as.numeric(mean) - se, 
                      colour = iter), position = position_dodge(width = 0.25)) +
    #geom_boxplot(colour = "dimgray", alpha = 0.5, width = 0.4) + 
    scale_color_manual(values = rep("gray", kfold)) + 
    labs(x = "Epitope length (amino acids)",  y = "CDR3 length (amino acids)") +
    theme_bw() +
    custom_theme + theme(legend.position = "none") +
    stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", color = "black", size=5)
  
  return(bp)
}
m <- lm(c(12.82, 12.39, 12.40, 11.42)~c(8:11))
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)

plotLengthRand(cd8b.nr, 100,1000, "Epitope.length") + 
  geom_abline(intercept = a, slope = b, colour = "darkorange")


plhlab1 <- plotLengthRand(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02",], 20, 12, "Epitope.length")
plhlab2 <- plotLengthRand(vdj.1b.hlab2, 20, 12, "Epitope.length")
plhlab3 <- plotLengthRand(vdj.1b.hlab3, 20, 12, "Epitope.length")
plhlab4 <- plotLengthRand(vdj.1b.hlab4, 20, 12, "Epitope.length")




textlab <- paste("y = ",b,"x + ",a, sep="")
print(textlab)
annotate(geom = "text", x = 10, y = 16.5, label = linear(lm(Epitope.length ~mean, t)))

################################################################################
## Barplot of junctional region ##
################################################################################
plotLength.nonVJ <- function(vdj){
  invisible(lapply(c("tidyr", "dplyr", "ggplot2", "RColorBrewer", "broman"), 
                   library, character.only = TRUE))
  # change TRA and TRB
  vdj$Gene <- ifelse(vdj$Gene == "TRA", "alpha", "beta")
  
  #select non-reducdant CDR3
  cdr3.nr <- vdj[!duplicated(vdj[,c("CDR3", "dominant")]),]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either ID or SD
  vdj <- droplevels(cdr3.nr)
  t <- vdj %>%
    group_by(dominant, T.cell, Gene, nonVJ.length) %>%
    dplyr::summarise(n = n()) %>%
    mutate(rel.freq = n*100/sum(n)) %>%
    filter(n > 1)
  
  rel.t <- as.data.frame(t)
  
  rel.t <- complete(rel.t, dominant, T.cell, Gene, nonVJ.length, fill = list(n = 0, rel.freq = 0))
  
  p <- ggplot(data = rel.t, aes(x= as.factor(nonVJ.length))) + 
    geom_bar(aes(y = rel.freq, fill = dominant), stat = "identity", position = "dodge") + 
    labs(title = expression(paste("Distribution of junctional length of CD8+CDR3 ", 
                                  alpha , " and ", beta, " chain")),
                            x = "non-VJ length (aa)", y = "Relative frequency") +
    custom_theme + 
    theme(legend.position = "bottom", 
          strip.text.y = element_text(face = "bold", size = 25),
          strip.switch.pad.grid = unit(0, "cm")) +
    scale_fill_manual(name = "Immune response", values = c("dimgray", "darkorange")) +
    facet_grid(Gene~., labeller = label_parsed) 
    
  return(p)
}

plotLength.nonVJ(rbind(cd8a.nr,cd8b.nr))


################################################################################
## plot CDR3.length vs Epitope.length  reregarding Epitope.length ##
################################################################################

# plot random sampling data
plotLength.epitope <- function(df, parm, hla){
  lapply(c("tidyr", "plotrix", "ggplot2", "ggpmisc", "ggpubr"), 
         library, character.only = TRUE)
  
  sum.data <- df %>% group_by(Epitope.length) %>% 
    dplyr::summarise(mean=mean(CDR3.length), sd=sd(CDR3.length), 
                     se = as.numeric(std.error(CDR3.length)))
  
  # plot 
  p <- ggplot(df, aes(x = Epitope.length, y = as.numeric(CDR3.length))) + 
    geom_jitter(aes(fill = Epitope.length, group = Epitope.length),size = 5, 
                position = position_dodge(width = 0.25),
                color = "lightgray") +
    geom_errorbar(aes(ymax = as.numeric(mean) + se, ymin = as.numeric(mean) - se, 
                      fill = Epitope.length), position = position_dodge(width = 0.25),
                  color = "lightgray", lty = 3) +
    geom_smooth(data = df, aes(x = Epitope.length, y = CDR3.length), 
                method = "lm",se = F,formula = y~x, size = 1, lty = 2, color = "orange1") +
    #geom_text(x = 0.6, y = 0.13, label = lm_eqn(vdj.rand), size = 5, parse = TRUE) +
    #geom_text(x = 0.6, y = 0.09, label = paste0("p-value ",compare$p.format), size = 5) +
    
    labs(x = "Epitope length (AA)", y = "CDR3 length (AA") + 
    theme_bw() +
    custom_theme + 
    stat_summary(data = df, aes(x =  Epitope.length, y = CDR3.length, label=round(..y..,2)),
                 fun.y=mean, geom="text", color = "black", size=5)
  
  return(p)
}

plotLength.epitope(rbind(ncd8b.nr, ncd8a.nr), "chain","")
plotLength(rbind(cd8b.nr[,-49], cd8a.nr), "chain","")

################################################################################
#Miscellaneous
count.hla <- table(cd8b.nr$MHC.A.new)
select.hla <- count.hla[count.hla > 45]

# remove the rare CDR3 length (select the length that appear > 1% of the whole data)
cd8bnr.length <- cd8b.nr %>%
  group_by(CDR3.length) %>%
  dplyr::summarise(n = n()) %>%
  mutate(frequency = n*100/sum(n)) %>%
  filter(frequency > 1)

cd8b.nr.filter <- cd8b.nr[cd8b.nr$MHC.A.new %in% names(select.hla) & 
                            cd8b.nr$CDR3.length %in% cd8bnr.length$CDR3.length,]

# group HLA-B based on the similarity of anchor residues
## HLA-B alleles ##
hlab1 <- droplevels(cd8b.nr.filter[cd8b.nr.filter$MHC.A.new %in% c("HLA-B*51", 
              "HLA-B*53", "HLA-B*42", "HLA-B*35", "HLA-B*81", "HLA-B*07"),])
hlab2 <- droplevels(cd8b.nr.filter[cd8b.nr.filter$MHC.A.new %in% c("HLA-B*58", "HLA-B*57"),])

hlab3 <- droplevels(cd8b.nr.filter[cd8b.nr.filter$MHC.A.new == "HLA-B*08",])

hlab4 <- droplevels(cd8b.nr.filter[cd8b.nr.filter$MHC.A.new == "HLA-B*27",])

# select the CDR3 and MHC class
ggplot(cd8b.nr.filter, aes(MHC.A.new, CDR3.length)) + 
  geom_boxplot() + 
  stat_compare_means()

################################################################################
## Miscellaneous ##
################################################################################

# generate the datasets with non-repetitive sequences
## write table 
#exportTable(spc_table, args[4])

## delete all objects except the function 
# rm(list = setdiff(ls(), lsf.str()))

## delete object by pattern
#rm(list = ls(pattern = "^b"))


