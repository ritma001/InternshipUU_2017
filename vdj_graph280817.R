#!usr/bin/r
################################################################################
# Author: Wannisa Ritmahan
# Established Date: 28/08/2017 Last Modification: 11/10/2017
#
# Aim: To study common feature of TCR motifs and epitopes from VDJdb
# and represent the result graphically
#
################################################################################
## Amino acid length plot ##
################################################################################
#
# plot the amino acid length distribution and highlight the different filter 
# parmrameter i.e. MHC.class (parm)
plotLength1 <- function(vdj, chain){
  library(ggplot2)
  library(ggpubr)
  p <- ggplot(droplevels(vdj), aes(x = T.cell, y = CDR3.length, fill = T.cell, color = T.cell)) +  
    geom_violin(trim = FALSE, size = 3) + 
    scale_fill_brewer(palette = "Set3") +
    scale_color_brewer(palette = "Set3") +
    #geom_jitter(shape = 19, alpha = 0.1, size = 3, width = 0.15) + 
    #geom_boxplot(alpha = 0.3, width = 0.4) + 
    ylab('CDR3 length') + xlab("T-cell population") + 
    ggtitle(paste("Distribution of CDR3 length of", chain, sep = " " )) +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15,face="bold"),
          strip.text = element_text(size=16, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          legend.text = element_text(size=14),
          legend.position = "none") + 
    facet_grid(~Gene)
  p <- p + stat_compare_means(label.y = 23, label.y.npc = "center", size = 5) 
  p <- p + stat_compare_means(aes(label = ..p.signif..),comparison = list(c("CD4+", "CD8+")), label.y = 21)
  p <- p + stat_summary(fun.y=mean, geom="point",color = "darkgray") +
    stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", color = "black", size=5)
  
    #facet_grid(MHC.class~dominant)
  ##Add p-value
  # p + stat_compare_means(method = "anova" or "ttest" if desired)
  # p + scale_y_continuous(breaks=seq(7,25,step = 2), labels=seq(7,25,step = 2),limits=c(7,25))
  return(p)
}

## call function ##
plotLength1(vdj,  "TCR alpha and beta chains in CD4+ and CD8+ T-cells")

################################################################################
plotLength2 <- function(vdj, chain){
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(broman)
  # select non-reducdant CDR3
  cdr3.nr <- vdj[!duplicated(vdj[,c("CDR3", "dominant")]),]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either IDD or SDD
  vdj <- droplevels(cdr3.nr)
  t <- vdj %>%
    group_by(dominant, T.cell, CDR3.length) %>%
    dplyr::summarise(n = n()) %>%
    mutate(rel.freq = n/sum(n))
  
  rel.t <- as.data.frame(t)
  
  rel.t <- complete(rel.t, dominant, T.cell, CDR3.length, fill = list(n = 0, rel.freq = 0))
  
  p <- ggplot(data = rel.t, aes(x= as.factor(CDR3.length))) + 
    geom_bar(aes(y = rel.freq, fill = dominant), stat = "identity", position = "dodge") + 
    labs(title = paste("The diversity of CDR3 length of", chain, sep = " "), x = "CDR3 length", y = "Relative frequency") +
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=13,face="bold"),
          strip.text = element_text(size=14, face="bold"),
          plot.title = element_text(size=14, face="bold", hjust = 0.5)) +
    #scale_fill_manual(values = c("dimgray", "gray")) +
    facet_wrap(~T.cell, scales = "free_y")
  return(p)
}

# call the function
pl2.a <- plotLength2(rbind(vdj.1a,vdj.2a), "alpha chain")
pl2.b <- plotLength2(rbind(vdj.1b,vdj.2b), "beta chain")

# test significant diffrence of CDR3 between IDD and SDD
mw.1a <- wilcox.test(vdj.1a[vdj.1a$dominant == "IDD","CDR3.length"], vdj.1a[vdj.1a$dominant == "SDD","CDR3.length"])
mw.1b <- wilcox.test(vdj.1b[vdj.1b$dominant == "IDD","CDR3.length"], vdj.1b[vdj.1b$dominant == "SDD","CDR3.length"])
mw.2a <- wilcox.test(vdj.2a[vdj.2a$dominant == "IDD","CDR3.length"], vdj.2a[vdj.2a$dominant == "SDD","CDR3.length"])
mw.2b <- wilcox.test(vdj.2b[vdj.2b$dominant == "IDD","CDR3.length"], vdj.2b[vdj.2b$dominant == "SDD","CDR3.length"])

library(ggpubr) ## to be done 
figure <- ggarrange(pl2.a, pl2.b, ncol = 1, nrow = 2) #, labels = LETTERS[1:6])
annotate_figure(figure, bottom = text_grob("CDR3 length", color = "black",
                                   hjust = 1, x = 0.5, face = "bold", size = 13),
                left = text_grob("Differential ratio of relative frequency between SDD and IDD response", color = "black",face = "bold", size = 13, rot = 90)
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)

################################################################################
dpVJ <- function(vdj.chain){
  library(ggplot2)
  library(ggpubr)
  library(magrittr)
  p <-  ggplot(vdj.chain, aes(x = dominant, y = CDR3.length)) +
    geom_jitter(aes(color = dominant), shape = 19, alpha = 0.2, size = 3, width = 0.15)  +
    labs(title = "Distribution of CDR3 length", x = "Immune response", y = "CDR3 length") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15,face="bold"),
          strip.text = element_text(size=16, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = "bottom") +
    facet_grid(Gene~T.cell) 
  ## add the chuck below to show the crossbar of mean ##
  p <- p + stat_summary(fun.y = "mean", geom = "crossbar", mapping=aes(ymin=..y.., ymax=..y..),
                   position=position_dodge(),show.legend = FALSE,color = "gray", width = 0.5)
  p <- p + stat_summary(fun.y=mean, geom="point") +
    stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", color = "white", size=5)
   
  return(p)
} 

################################################################################
# exploring the length of peptide binding in each T cell 
relFreq.epi <- function(t){
  rel.t <- t %>% 
  group_by(dominant,CDR3.length, Epitope.species) %>%
  dplyr::summarise(n = n()) %>%
  mutate(rel.freq = n/sum(n))

  #rel.t <- complete(rel.t, dominant, CDR3.length, Epitope.length, fill = list(n = 0, rel.freq = 0))
  return(rel.t)
}

plotrelFreq3 <- function(t, chain){
  library(tidyr)
  library(dplyr)
  p <- ggplot(data = t, aes(x= as.factor(CDR3.length))) + 
    geom_bar(aes(y = rel.freq, fill = Epitope.species), stat = "identity") + 
    labs(title = paste("The diversity of  CDR3 length vs epitope species of ", chain, sep = " "), x = "CDR3 length", y = "Relative frequency") +
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=13,face="bold"),
          strip.text = element_text(size=14, face="bold"),
          plot.title = element_text(size=14, face="bold", hjust = 0.5),
          legend.position = "bottom") +
    facet_wrap(~dominant, scales = "free_x")
  return(p)
}
 
# call fucntion
plotrelFreq3(relFreq.epi(vdj.1a), "alpha chain CD8+ T-cells")
################################################################################
## Epitope plot ##
################################################################################
#
# plot epitope species regarding the CDR3 length
epitope <- sort(summary(droplevels(vdj2$Epitope)),decreasing = TRUE)
#
# print output
write.table(paste(names(epitope),epitope, sep = "\t"), file = "epitope frequency.txt", 
            append = FALSE, quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
#
# plot diversity species for each epitope species 
p <- ggplot(data = vdj2.HLA2, aes(x=factor(Epitope.species), fill=Epitope)) + geom_bar()+
  ylab('Relative Frequency') + xlab("CDR3 length") + 
  ggtitle("The distribution of Epitopes for species") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
p + facet_wrap(~dominant, scales = "free_y" ,drop = TRUE)+ 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(13))
#
# pie chart of Epitope species presenting in dataset
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
epi <- as.data.frame(summary(droplevels(vdj$Epitope.species)))
pie.epi <- epi[epi != 0, ,FALSE]
pie.labels <- c(paste(row.names(epi),epi[,1],sep='\n'))
pie(epi[,1], col = brewer.pal(8,"Greys"), labels = pie.labels, main = 'Epitope species of CDR3 motifs included in this analysis')
#
################################################################################
## analysis VJ combination ##
################################################################################
#
# Plot relative frequency of V or J loci distribution of selected (alpha or beta) Chain
plotVJ <- function(vdj, gene, chain, main){
  # data processing
  library(dplyr)
  library(tidyr)
  # select either alpha or beta chain
  vdj <- vdj[vdj$Gene == chain,]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either IDD or SDD
  vdj <- droplevels(vdj)
  vdj.t <- vdj %>%
    group_by_(~dominant, gene) %>%
    dplyr::summarise(n = n()) %>%
    mutate(rel.freq = n/sum(n))
  
  # complete the table zero observation of certain V/J segments
  if(gene == "J"){
    vdj.t <- complete(vdj.t, dominant, J, fill = list(n = 0, rel.freq = 0))
  } else {
    vdj.t <- complete(vdj.t, dominant, V, fill = list(n = 0, rel.freq = 0))
  }
  vdj.t <- as.data.frame(vdj.t)
  vdj.t[order(vdj.t[, 4]),] 
  # plot
  library(ggplot2)
  p <- ggplot(vdj.t, aes_string(x = gene)) + 
    geom_bar(aes(y = rel.freq, fill = dominant), stat = "identity", position = "dodge") + 
    labs(title = main, x = "Gene loci", y = "Relative frequency") + 
    scale_fill_manual(values = c(gray(0.4),"gray")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")
  return(p)
}  

# call the function
plotVJ(vdj1, "V", "TRB", "Distribution of V loci in TCR beta chain of CD8+ T-cells")
plotVJ(vdj1, "V", "TRA", "Distribution of V loci in TCR alpha chain of CD8+ T-cells")
plotVJ(vdj1, "J", "TRA", "Distribution of J loci in TCR alpha chain of CD8+ T-cells")
plotVJ(vdj1, "J", "TRB", "Distribution of J loci in TCR beta chain of CD8+ T-cells")

plotVJ(vdj2, "V", "TRB", "Distribution of V loci in beta chain CDR3 of MHC-II")
plotVJ(vdj2, "V", "TRA", "Distribution of V loci in alpha chain CDR3 of MHC-II")
plotVJ(vdj2, "J", "TRA", "Distribution of J loci in alpha chain CDR3 of MHC-II")
plotVJ(vdj2, "J", "TRB", "Distribution of J loci in beta chain CDR3 of MHC-II")

################################################################################
plotVJ2 <- function(vdj.chain, gene, title.name){
  # data processing
  library(dplyr)
  library(tidyr)
  library(reshape)
  library(ggplot2)
  
  # select non-reducdant CDR3
  cdr3.nr <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant")]),]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either IDD or SDD
  vdj <- droplevels(cdr3.nr)
  vdj.t <- vdj %>%
    group_by_(~dominant, gene) %>%
    dplyr::summarise(n = n()) %>%
    mutate(rel.freq = n/sum(n))

  vdj.t <- as.data.frame(vdj.t)
  colnames(vdj.t) <- c("response", "gene", "count", "rel.freq")
  
  # complete the table zero observation of certain V/J segments
  ##vdj.tc <- complete(vdj.t, response, gene, fill = list(count = 0, rel.freq = 0))
    
  #compare the difference between IDD and SDD response
  t <- cast(vdj.t, gene~response, value = "rel.freq")
  t$diff <- log2(t$IDD/t$SDD)
  t <- as.data.frame(t)
  t <- t[complete.cases(t), ]
  print(t)
 
  # select only the the difference morethan 1.0
  t <- t[abs(t$diff) > 1 ,]

  #plot*gray(0.4) is darken than gray*
  p <- ggplot(droplevels(t), aes(x = reorder(gene,diff), y = diff)) +
    geom_bar(stat = "identity", width = 0.5) + 
    scale_fill_manual(values = c(gray(0.4), "gray")) + 
    labs(x = NULL, y = "log2(IDD/SDD)") + ggtitle(title.name) + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text = element_text(size=15, face="bold"),
          plot.title = element_text(size=15, face="bold", hjust = 0.5),
          legend.position = "none")
  return(p)
}

v.1b <- plotVJ2(vdj.1b, "V.new", "Distribution of V loci of CD8+ T-cell beta chain") 
v.1a <- plotVJ2(vdj.1a, "V.new", "Distribution of V loci of CD8+ T-cell alpha chain") 
v.2b <- plotVJ2(vdj.2b, "V.new", "Distribution of V loci of CD4+ T-cell beta chain") 
v.2a <- plotVJ2(vdj.2a, "V.new", "Distribution of V loci of CD4+ T-cell alpha chain")
  
j.1b <- plotVJ2(vdj.1b, "J", "Distribution of J loci of CD8+ T-cell beta chain") 
j.2b <- plotVJ2(vdj.2b, "J", "Distribution of J loci of CD4+ T-cell beta chain")
j.1a <- plotVJ2(vdj.1a, "J.new", "Distribution of J loci of CD8+ T-cell alpha chain")
j.2a <- plotVJ2(vdj.2a, "J.new", "Distribution of J loci of CD4+ T-cell alpha chain") 

## arrange the plots 
library(ggpubr)
figure.va <- ggarrange(v.1a, v.2a,j.1a, j.2a, ncol = 2, nrow = 2) 
figure.vb <- ggarrange(v.1b, v.2b,j.1b, j.2b, ncol = 2, nrow = 2) #, labels = LETTERS[1:6])
annotate_figure(figure, text_grob("Comparison of V and J loci distribution between IDD and SDD responses", color = "black", face = "bold", size = 14),
                bottom = text_grob("Gene subgroup", color = "black",
                                   hjust = 1, x = 0.5, face = "bold", size = 13),
                left = text_grob("Differential ratio of relative frequency between SDD and IDD response", color = "black",face = "bold", size = 13, rot = 90)
                #fig.lab = "Figure 1", fig.lab.face = "bold"
                  )

## inspect the V_J bias in related to CDR3 length
vjtable <- function(vdj.chain){
  cdr3.nr <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant")]),]
  t <- cdr3.nr %>% 
    group_by(V.new, dominant, CDR3.length) %>% 
    dplyr::summarise(n = n()) %>% 
    mutate(rel.freq = n/sum(n))
  return(t)
}

t.1b <- vjtable(vdj.1b)
t.2b <- vjtable(vdj.2b)

reduceRed <- function(vdj.chain, region){
  t <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant")]),]
  species <- summary(droplevels(t[t$V.new == region, "Epitope.species"]))
  return(species)
}

table <- function(vdj.chain, gene){
  cdr3.nr <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant")]),]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either IDD or SDD
  vdj <- droplevels(cdr3.nr)
  vdj.t <- vdj %>%
    group_by_(~dominant, gene) %>%
    dplyr::summarise(n = n()) %>%
    mutate(rel.freq = n/sum(n))
  vdj.t <- as.data.frame(vdj.t)
  colnames(vdj.t) <- c("response", "gene", "count", "rel.freq")
  # complete the table zero observation of certain V/J segments
 # vdj.tc <- complete(vdj.t, response, gene, fill = list(count = 0, rel.freq = 0))

  #compare the difference between IDD and SDD response
  t <- cast(vdj.t, gene~response, value = "rel.freq")
  t$diff <- log2(t$IDD) - log2(t$SDD)
  t <- t[!is.infinite(rowSums(t[,-1])),]
  return(t)
}

table2 <- function(vdj.chain){
  t <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant")]),]
  result <- t %>% 
    group_by(Epitope.species, dominant, CDR3.length) %>% 
    summarise(n = n())
  return(result)
}


################################################################################
# exclude the segments with relative frequency less than 1% (0.01)

plotV <- function(t, chain, dominant){
  library(ggplot2)
  library(RColorBrewer)
  rel.freq <- prop.table(summary(droplevels(t$V)))
  sel.V <- rel.freq[as.vector(rel.freq) > 0.01] # 0.1 for alpha and 0.2 for beta
  t <- t[t$V %in% names(sel.V),]
  
  #show the highest freq at the top by rearranging the order of V seqments
  t$V <- factor(t$V, levels = names(sort(summary(droplevels(t$V)),decreasing = TRUE)))
  #extrapolate the Rbrewer color palette
  colourCount = length(levels(t$V))
  getPalette = colorRampPalette(brewer.pal(12, "Set3"))
  
  p <- ggplot(data = t, aes(x=factor(CDR3.length), fill=V, order = -as.numeric(V))) + geom_bar(stat = "identity")+
    ylab('Frequency') + xlab("CDR3 length") + ggtitle(paste0(paste0("The distribution of V seqments for each CDR3 length in", chain), dominant))+
    theme(plot.title = element_text(hjust = 1))
  p + facet_wrap(~dominant, scales = "free_y" ,drop = TRUE)+ scale_fill_manual(values = getPalette(colourCount))
}

# call function
plotV(vdj2.a," alpha chain of ", "IDD and SDD")
plotV(vdj2.b," beta chain of ", "IDD and SDD")

plotV2 <- function(t,chain, dominant){
  #show the highest freq at the top by rearranging the order of V seqments
  t$V <- factor(t$V, levels = names(sort(summary(droplevels(t$V)),decreasing = TRUE)))
  #extrapolate the Rbrewer color palette
  p <- ggplot(data = t, aes(x=factor(V), fill = dominant)) + geom_bar(stat = "identity", position = "dodge")+
    ylab('Relative Frequency') + xlab("V segment") + ggtitle(paste0(paste0("The distribution of V seqments in ", chain), dominant)) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  p + facet_wrap(~dominant, scales = "free_y" ,ncol = 1)
  p
}

plotV2(vdj2.a," alpha chain of ", "IDD and SDD") 
plotV2(vdj2.b," beta chain of ", "IDD and SDD")

plotV2(HLA2a.idd, " alpha chain of ", "IDD")


################################################################################
plotJ <- function(t,chain,dominant){
  library(ggplot2)
  library(RColorBrewer)

  rel.freq <- prop.table(summary(droplevels(t$J)))
  sel.J <- rel.freq[as.vector(rel.freq) > 0.01]
  t <- t[t$J %in% names(sel.J),]
  
  #show the highest freq at the top by rearranging the order of V seqments
  t$J <- factor(t$J, levels = names(sort(summary(droplevels(t$J)),decreasing = TRUE)))
  
  #extrapolate the Rbrewer color palette
  colourCount = length(levels(t$J))
  getPalette = colorRampPalette(brewer.pal(12, "Set3"))
  
  # calculate frequency for each group 
  
  # plot
  p <- ggplot(data = t, aes(x=factor(CDR3.length), fill=J, order = -as.numeric(J))) + geom_bar()+
    ylab('Relative Frequency') + xlab("CDR3 length") + ggtitle(paste0(paste0("The distribution of J seqments for each CDR3 length in", chain), dominant))+
    theme(plot.title = element_text(hjust = 0.5))
  p + facet_wrap(~dominant, scales = "free_y" ,drop = TRUE)+ scale_fill_manual(values = getPalette(colourCount))
}

#call function 

plotJ(vdj2.a," alpha chain of ", "IDD and SDD")
plotJ(vdj2.b," beta chain of ", "IDD and SDD")

#if no facet_wrap and p+ scale_fill_manual(...)

plotJ(sdd.b, " beta chain of ","SDD")
plotJ(sdd.a, " alpha chain of ","SDD")
plotJ(idd.a, " alpha chain of ","IDD")
plotJ(idd.b, " beta chain of ","IDD")

################################################################################
## amino acid distribution ##
################################################################################
#
# count amino acids (1 mer) of cdr3 sequencing: input =  a vector of CDR3 sequencing 
countAmino <- function(cdr){
  t <- data.frame("amino_acid" = c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E"),
                  "properties" = c(rep("hydrophobic",8),rep("hydrophilic",8),rep("charged",4)), "count" = rep(0, 20))
  #chemical properties from http://www.proteinstructures.com/Structure/Structure/amino-acids.html
  for(a in cdr){
    count <- table(strsplit(toupper(as.character(a)), ''))
    for(e in names(count)){
      t[t$amino_acid == e,3] <- t[t$amino_acid == e,3] + as.vector(count[e])
    }
  }
  #calculate relative frequency and add to the last column
  rel.freq <- t$count/sum(t$count)
  t$frequency <- rel.freq
  return(t)
}
#
## call function
cdr.sdd.a <- countAmino(vdj2.a[vdj2.a$dominant == "SDD", "CDR3"])
cdr.idd.a <- countAmino(vdj2.a[vdj2.a$dominant == "IDD", "CDR3"])
cdr.sdd.b <- countAmino(vdj2.b[vdj2.b$dominant == "SDD", "CDR3"])
cdr.idd.b <- countAmino(vdj2.b[vdj2.b$dominant == "IDD", "CDR3"])

################################################################################
## V-J distribution ##
################################################################################

findOverlap <- function(vdj.chain, vj.segment){
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  vdj.idd <- vdj.chain[vdj.chain$dominant == "IDD", vj.segment]
  vdj.sdd <- vdj.chain[vdj.chain$dominant == "SDD", vj.segment]
  int <- intersect(vdj.idd, vdj.sdd)
  unique.idd <- levels(droplevels(vdj.idd[!vdj.idd %in% int]))
  unique.sdd <- levels(droplevels(vdj.sdd[!vdj.sdd %in% int]))
  # make the dataframe by joining (un)equal length of vectors: int, unique.idd & .sdd
  pre.df <- list(int,unique.idd, unique.sdd)
  # check max nr. of rows
  max.row <- max(sapply(pre.df, length))
  #insert the needed `NA`s to each dataframe
  new.pre.df <- lapply(pre.df, function(x) {as.data.frame(x)[1:max.row,]})
  # construct a datafraframe 
  df <- as.data.frame(matrix(unlist(new.pre.df), nrow=max.row))
  colnames(df) <- c("common.segment", "unique.idd", "unique.sdd")
  return(df)
}

# make the venn diagram
library('VennDiagram')
grid.newpage()
vd <- draw.pairwise.venn(57+5,7+5,5, category = c("SDD", "IDD"), lty = rep("blank", 2), fill = c("light blue", "gray"), 
                         alpha = rep(0.5, 2), cat.pos = c(180, 180), euler.d = T, cex = 2, cat.cex = 2,
                         sep.dist = 0.03, scaled = F)
grid.text("Non-redundance of 12-amio acid long CDR3 of CD8+ T-cell alpha chain", y= 0.9, gp=gpar(cex=2))

################################################################################
## Part V: Calulate the relative frequency based on the idd and sdd ##
################################################################################
## input: the dataset contains both IDD and SDD, 
## column with rel. frequency for each dominant type
#
relFreq <- function(table, parm){
  rel.freq <- c()
  # calculate the total number of IDD and SDD
  sum.IDD <- nrow(table[table$dominant == "IDD", ])
  sum.SDD <- nrow(table[table$dominant == "SDD", ])
  
  # ensure the the collumn is a factor 
  table[, parm] <- factor(table[, parm])
  
  # row-wise iteration and keep the frequency in the rel.freq vector
  for(i in 1:nrow(table)){
    if(table[i,"dominant"] == "IDD"){
      level.freq <- length(table[table$dominant == "IDD" & table[,parm] == table[i,parm], parm])/sum.IDD
    } else {
      level.freq <- length(table[table$dominant == "SDD" & table[,parm] == table[i,parm], parm])/sum.SDD
    }
    rel.freq <- c(rel.freq, level.freq)
  }
  table$rel.freq <- rel.freq
  return(table)
}
# plot the 2plot top-bottom of IDD and SDD with black color
library(ggplot2)
plotrelFreq <- function(table, var, main){
  #table[,var] <- factor(table[,var], levels = names(sort(summary(droplevels(table[,var])),decreasing = TRUE)))
  p <- ggplot(data = table, aes_string(x=var)) + 
    geom_bar(aes(y=rel.freq, fill = dominant),stat = "identity", position = "dodge") + 
    labs(title = main, y = "Relative frequency") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p + facet_wrap(~, nrow =2, scales = "free_y")
  #return(p)
}
# call function
plotrelFreq(relFreq(vdj2.a, "V"), "V" , "V distribution of alpha chain")
plotrelFreq(relFreq(vdj2.b, "V"), "V" , "V distribution of beta chain")

# plot 2 bar next to each other
library(ggplot2)
library(tidyr)
library(reshape2)

plotrelFreq2 <- function(table, var, main){
  table[,var] <- factor(table[,var], levels = names(sort(summary(droplevels(table[,var])),decreasing = TRUE)))
  p <- ggplot(data = table, aes_string(x=var)) + 
    geom_bar(aes(y=rel.freq, fill = Epitope),stat = "identity", position = "dodge") + 
    labs(title = main, y = "Relative frequency") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p + facet_wrap(~dominant)
}
# call function
plotrelFreq2(relFreq(vdj.2a, "V"), "V" , "V distribution of alpha chain")
plotrelFreq2(relFreq(vdj.2b, "V"), "V" , "V distribution of beta chain")
