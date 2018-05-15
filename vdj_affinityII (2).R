# 6-23 Feb 2018
# personalised analysis
# select the ID with > 40% then study the  individual response
#
################################################################################
# explore the length of Epitope and CDR3
################################################################################
## getEpiProp: get an overview of propeties (hydropathy, size, amino acids)
## of epitopes (at the same length) in each MHC molecule
getEpiProp <- function(df, MHC.class, prop){
  # extract the data
  lapply(c("tidyr", "dplyr"), library, character.only = TRUE)
  df.epi <- df %>% 
    select(MHC.A.new, Epitope.length, Epitope, Epitope.species) %>%
    distinct(.keep_all = TRUE) %>%
    filter(MHC.A.new == MHC.class)
  
  # identify aa at each position
  for(i in unique(df.epi$Epitope.length)){
    print(i)
    df.epi2 <- df.epi[df.epi$Epitope.length == i, ]
    epi.msa <- moveGap(sapply(as.character(df.epi2$Epitope), addGap))
    epi.frequency <- positionalCount(epi.msa)
    # plot heatmap
    print(HM.prop(epi.frequency, prop))
  }
}

################################################################################
## plot
plotParmHla <- function(df, parm){
  # select unique sequence in id and sd datasets
  unique.df <- rbind(df[df$CDR3 == unique(as.character(df[df$dominant == "ID", "CDR3"])),],
                 df[df$CDR3 == unique(as.character(df[df$dominant == "SD", "CDR3"])),])
  ggplot(unique.df, aes(x = dominant)) + geom_boxplot(aes_string(y = parm))
}

################################################################################
# test hydropathy in different HLA groups
################################################################################
# 1. permutation test using coin package
library(coin)
independence_test(CDR3_AA_GRAVY ~ factor(dominant), 
                  data = cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02", ])

#2. Mann-whitney U test (Wilcox rank sum test) using unpaired wilco.test 
wilcox.test(CDR3_AA_GRAVY ~ dominant, data = cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02", ])

################################################################################
plotHydropathy.HLA <- function(df, hla){
  # import necessary packages
  invisible(lapply(c("ggplot2", "ggpubr", "coin"), library, character.only = TRUE))
  
  # subset data base on HLA molecule
  df <- droplevels(df[df$MHC.A.new == hla, ])
  
  # extract p-value from permutation test
  permutation <- independence_test(CDR3_AA_GRAVY ~ factor(dominant), data = df)
  invisible(str(permutation))
  pv.perm <- permutation@distribution@pvalue(permutation@statistic@teststatistic)
  pv <- round(pv.perm, 4)
  
  # plot 
  p <- ggplot(df, aes(x = dominant, y = CDR3_AA_GRAVY)) + 
    geom_jitter(aes(color = dominant), size = 3.0, width = 0.3) +
    geom_boxplot(alpha = 0.4, color = "gray", width = 0.65) +
    annotate(geom = "text", x = 1.5, y = max(df$CDR3_AA_GRAVY),  vjust = 0.5, 
              label = paste0("p-value ", pv), size = 7) +
    scale_color_manual(values = c("dimgray", "darkorange")) +
    labs(x = "Immune response", y = "CDR3 hydropathy") + 
    custom_theme +
    stat_summary(aes(label=round(..y.., 2)), fun.y=mean, geom = "text", 
                 color = "black", size = 7)
  
  return(p)
}

################################################################################
# pie chart
################################################################################

# pie chart
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
  pie <- pie3D(t$frequency, main = "Proportion of the different MHC molecules", 
        explode= 0.15, theta = 0.5, radius = 1.6, start = 0, 
        mar = c(2.5,2.5,2.5,2.5), border = "gray", font = 2,
        col = brewer.pal(12, "Set3")[1:dim(t)[1]])
  return(pie)
}

# call pie 3D
pie <- pie3D.hla(cd8b.nr)

################################################################################

# bold label
par(font = 2)

# add label maunally
pie.pos <- pie # a vector of label position from pie
pie.pos[8] <- 6.0
pie.pos[9] <- 6.20
pie.pos[10] <- 6.40

pie3D.labels(pie.pos, labels = paste(t$MHC.A.new, "\n", paste0(t$frequency,"%")), 
             labelcex = 0.9, labelrad = 1.9, minsep = 0.3)


################################################################################
# bar chart of vdj db summary #
################################################################################

barPlot.hla <- function(df){
  # import necessary packages
  invisible(lapply(c("tidyr", "dplyr", "reshape", "ggplot2", "ggpubr", "RColorBrewer"), 
                   library, character.only = TRUE))
  
  # exclude the MHC with low frequency
  df.new <- droplevels(df[df$MHC.A.new %in% 
                  names(table(df$MHC.A.new))[table(df$MHC.A.new) > 50] , ])
  
  # find relative frequency of ID and SD in each MHC subset
  t <- df.new %>% 
    group_by(Epitope.species, MHC.A.new) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(total = sum(count))
  
  levels(t$Epitope.species) <- c(levels(t$Epitope.species), "YFV", "IAV")
  t$Epitope.species[t$Epitope.species == 'YellowFeverVirus'] <- 'YFV'
  t$Epitope.species[t$Epitope.species == 'InfluenzaA'] <- 'IAV'
  
  p <- ggplot(t, aes(x= reorder(Epitope.species, -total), y =  count, fill = MHC.A.new)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(Epitope.species, total + 30, label= total), size = 8) +
    labs(x = "Epitope species", y = "Absolute count") +
    scale_fill_brewer(name = "MHC alleles", palette = "Set3") +
    custom_theme + 
    theme(axis.text.x = element_text(angle = 45),
          legend.position = "right")
  
  return(p)
}

barPlot.hla(cd8b.nr)


# get data summary
summary(combine.naive.nn %>% filter(dataset == "naive") %>% select(CDR3_AA_GRAVY))
sd((cd8b.nr %>% filter(frequency > 1) %>% select(CDR3_AA_GRAVY))[,1])

################################################################################
# select the tope ten epitope
library(data.table)
epitope.count <- sort(table(droplevels(cd8b.nr[cd8b.nr$frequency > 1, "Epitope"])), 
                      decreasing = T)
  
# select the top 10 highest epitope
epitope.topten <- names(epitope.count[1:10])

# plot length to compare ID and SD
for(epi in epitope.topten){
  df <- cd8b.nr[cd8b.nr$Epitope == epi, ] %>% filter(frequency > 1)
  p <- boxplotLength(df) + theme(legend.position = "none") 
  # only significant in "GTSGSPIVNR"; epitope.topten[7]
  print(p)
}

# explore the data
table(cd8b.nr[cd8b.nr$frequency > 1, "CDR3.length"])
table(ncd8b.nr[ncd8b.nr$Naive.count > 1, "CDR3.length"])

sd(cd8b.nr[cd8b.nr$frequency > 1, "CDR3.length"]-2)
mean(cd8b.nr[cd8b.nr$frequency > 1, "CDR3.length"]-2)

# check the column of naive data 
cdr3.col <- which(colnames(ncd8b.nr) == "AAseq", arr.ind = T)
colnames(ncd8b.nr)[cdr3.col] <- "CDR3"



# plot length of all Epitope
plotLength(droplevels(cd8b.nr[cd8b.nr$Epitope %in% epitope.topten, ]), 
           "Epitope", "Virus species")

# Plot cd8 alpha and beta datasets
cd8ab <- rbind(cd8a.nr, cd8b.nr)
plotLength.summary(summariseDat(cd8ab, "Gene", "CDR3.length"),
                   "Gene",  expression(paste("CD8+T-cell ", alpha, "and", beta, " chain")),
                   c(expression(paste("  ", alpha, " chain")), expression(paste("  ", beta, " chain"))))


################################################################################
# affinity model II: acute pathogenesis
################################################################################
# get the summary of mhc molecules
lapply(c("tidyr", "dpyr"), library, character.only = TRUE)
cd8bnr.mhc <- cd8b.nr %>% 
  filter(frequency > 1) %>% 
  group_by(MHC.A.new, Epitope, Epitope.species) %>% 
  dplyr::summarise(count = n()) 

# Dengue virus binding energies
a2DEN1 <- droplevels(cd8b.nr[cd8b.nr$Epitope == "GTSGSPIVNR", ])
a2DEN1 <- a2DEN1[!duplicated(a2DEN1[, "CDR3"]),]

a2DEN34 <- droplevels(cd8b.nr[cd8b.nr$Epitope == "GTSGSPIINR", ])
a2DEN34 <- a2DEN1[!duplicated(a2DEN1[, "CDR3"]),]

# Yellow fever virus
a2YFV <- droplevels(cd8b.nr1[cd8b.nr1$Epitope == "LLWNGPMAV", ])

################################################################################
plotAffinity <- function(df, epitope, mj.mat, dataset.label){
  library(ggplot2)
  library(ggpubr)
  df <- df[!duplicated(df[,c("CDR3")]),]
  df$CDR3.se <- substr(as.character(df$CDR3), 5, df$CDR3.length-5) #CDR3 beta 109(6) - 112(9)
  df <- df[df$CDR3.se != "",]
  df$affinity <- sapply(df$CDR3.se, calMotifAffinity, epitope, mj.mat)
  
  # plot
  ggplot(df, aes(x = dominant, y = affinity)) +
    geom_jitter(size = 3, width = 0.2) +
    geom_boxplot(alpha = 0.3, width = 0.4, color = "darkgray") +
    labs(x = "Immune response", y = "Binding energy (kB)", title = dataset.label) +
    custom_theme +
    theme(legend.position = "none",
          axis.text = element_text(face = "bold", size = 13)) +
    stat_compare_means(aes(label = paste0("p-value = ", ..p.format..)), 
                       method = "wilcox.test", paired = F, 
                       label.x = 1.4, label.y.npc = 'top', size = 5)
}

################################################################################
# load package for immuno repertoire analysis
install.packages("devtools")
devtools::install_github("imminfo/tcr", build_vignettes = FALSE)
install.packages("tcR", dependencies = TRUE)

# affinity of DEN1 and DEV3/4
aff.gts1 <- plotAffinity.denv(a2DEN1, "GSPIV", mj.mat, "GTS1-A11")
aff.gts34 <- plotAffinity.denv(a2DEN34, "GSPII", mj.mat, "GTS3/4-A11")

# affinity of YFV
aff.yfv <- plotAffinity(a2YFV, "NGPMA", mj.mat, "LLW-A2")

################################################################################
# GILGFVFTL variation
################################################################################
# significant analysis
aff.gil <- plotAffinity(a2GIL %>% filter(frequency > 1), "GFVFT", mj.mat, "GIL-A2")

# 
plotAffinity(a2GIL %>% filter(frequency > 1), "GIL", mj.mat, "GIL-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "ILG", mj.mat, "ILG-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "LGF", mj.mat, "LGF-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "GFV", mj.mat, "GFV-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "FVF", mj.mat, "FVF-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "VFT", mj.mat, "VFT-A2")
plotAffinity(a2GIL %>% filter(frequency > 1), "FTL", mj.mat, "FTL-A2")


################################################################################
# compare affinity across epitope reagrdless of immunodominance
t <- droplevels(cd8b.nr[cd8b.nr$Epitope %in% 
                           c("LLWNGPMAV", "NLVPMVATV", "GLCTLVAML", "GILGFVFTL"), ])

t$CDR3.se <- substr(as.character(t$CDR3), 5, t$CDR3.length-5)
t$epi.se <- substr(as.character(t$Epitope), 3, 8)
t <- t[t$CDR3.se != "",]
  
t$aff <- sapply(t$CDR3.se, calAffinity, t$epi.se, mj.mat)
  
ggplot(t, aes(x = Epitope, y = aff)) + 
  geom_jitter(size = 3, width = 0.2, color = "darkgray") +
  geom_boxplot(alpha = 0.3, width = 0.4, color = "black") +
  labs(x = "", y = "Binding energy (kB)") +
  custom_theme + stat_compare_means(aes(label = paste("p-value = ",..p.format..)), 
                                    label.x = 2, label.y = -1 ,size = 5)


################################################################################
# add affinity value to the dataframe 
addAffinity <- function(df, epitope){
  df$CDR3.se <- substr(as.character(df$CDR3), 
                       df$v.end + 1, df$CDR3.length - 5)
  df <- df[df$CDR3.se != "",]
  df$affinity <- sapply(df$CDR3.se, calMotifAffinity, epitope, mj.mat)
  return(df)
}

a2GIL <- addAffinity(a2GIL, "GFVFT")
a2GLC <- addAffinity(a2GLC, "TLVAM")
a2NLV<- addAffinity(a2NLV, "PMVAT")
a2DEN1 <- addAffinity(a2DEN1, "GSPIV")
a2DEN34 <- addAffinity(a2DEN34, "GSPII")

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


t <- addAffinity(cd8a0.01nr[cd8a0.01nr$Epitope.species == "InfluenzaA",],  "GFVFT")
t <-  rbind(t[, c("Gene", "affinity")], a2gilall[, c("Gene", "affinity")])
ggplot(t, aes(Gene, affinity)) + geom_boxplot() + custom_theme
