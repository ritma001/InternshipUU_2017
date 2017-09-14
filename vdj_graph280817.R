#!usr/bin/r
################################################################################
# Author: Wannisa Ritmahan
# Established Date: 28/08/2017 Last Modification: 12/09/2017
#
# Aim: To study common feature of TCR motifs and epitopes from VDJdb
# and represent the result graphically (require library ggplot2, reshape2?, scales?)
#
################################################################################
## Amino acid length plot ##
################################################################################
#
# plot the amino acid length distribution and highlight the different filter 
# parmrameter i.e. MHC.class (parm)

plotLength <- function(t){
  p <- ggplot(data = t, aes(x=factor(CDR3.length))) + 
    geom_bar(aes(y = ..count../sum(..count..), fill = dominant)) + 
    labs(title = "CDR3 length Distribution", x = "CDR3 length", y = "Relative frequency") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    facet_wrap(~dominant, scales = "free_y")
  return(p)
}
#
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
epi <- as.data.frame(summary(droplevels(vdj$Epitope.species)))
pie.epi <- epi[epi != 0, ,FALSE]
pie.labels <- c(paste(row.names(epi),epi[,1],sep='\n'))
pie(epi[,1],col = rainbow(length(epi[,1]),0.3), labels = pie.labels, main = 'Epitope species 
    of Human\' s CDR3 motifs from VDJdb database' )
#
################################################################################
## analysis VJ combination ##
################################################################################
#
# Plot relative frequency of V or J loci distribution of selected (alpha or beta) Chain
plotVJ <- function(vdj, gene, chain, main){
  # data processing
  library(plyr)
  library(dplyr)
  library(tidyr)
  # select either alpha or beta chain
  vdj <- vdj[vdj$Gene == chain,]
  # calculate relative frequency based on V or J loci normalised with the total 
  # loci in either IDD or SDD
  vdj <- droplevels(vdj)
  vdj.t <- vdj %>%
    group_by_(~dominant, gene) %>% # mind the order domianant and V
    summarise(n = n()) %>%
    mutate(rel.freq = n/sum(n))
  
  # complete the table zero observation of certain V/J segments
  if(gene == "J"){
    vdj.t <- complete(vdj.t, J, fill = list(n = 0, rel.freq = 0))
  } else {
    vdj.t <- complete(vdj.t, V, fill = list(n = 0, rel.freq = 0))
  }
  vdj.t <- as.data.frame(vdj.t)
  vdj.t[order(vdj.t[, 4]),] 
  
  # plot
  library(ggplot2)
  p <- ggplot(vdj.t, aes_string(x = gene)) + 
    geom_bar(aes(y = rel.freq, fill = dominant), stat = "identity", position = "dodge") + 
    labs(title = main, x = "Gene loci", y = "Relative frequency") + 
    scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")
  return(p)
}  

# call the function
plotVJ(vdj2, "V", "TRA", "Distribution of V loci in alpha chain CDR3")

################################################################################

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
  
  p <- ggplot(data = t, aes(x=factor(CDR3.length), fill=V, order = -as.numeric(V))) + geom_bar()+
    ylab('Frequency') + xlab("CDR3 length") + ggtitle(paste0(paste0("The distribution of V seqments for each CDR3 length in", chain), dominant))+
    theme(plot.title = element_text(hjust = 0.5))
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

