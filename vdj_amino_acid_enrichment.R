################################################################################
#
# Author: Wannisa Ritmahan
# Aim: inspect an influence of amino acid arrangement and composition of CDR3 on
# the immunodominance
# 
################################################################################
## An overview of amino acid distribution ##
################################################################################
countNucleotide <- function(ins.base){
  # convert factor to character #
  if(is.factor(ins.base) == TRUE){
    ins.list <- as.character(levels(ins.base))[ins.base]
  } else {
    ins.list <- ins.base
  }
  # count amino acid in each cdr3 #
  nucleotides <- c("A","T","C","G")
  
  # count the occurance of each amino acid #
  count <- sapply(nucleotides, function(x) x<-sum(x==unlist(strsplit(ins.list,""))))
  count.df <- data.frame("nucleotides" = c(names(count)), "count" = c(as.numeric(count)))
  count.df$frequency <- count.df$count*100/(sum(count.df$count))
  return(count.df)
}

################################################################################
countAmino <- function(cdr3.list){
  # convert factor to character #
  if(is.factor(cdr3.list) == TRUE){
    cdr3.list <- as.character(levels(cdr3.list))[cdr3.list]
  }
  # count amino acid in each cdr3 #
  amino.acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
  
  # count the occurance of each amino acid #
  count <- sapply(amino.acid, function(x) x<-sum(x==unlist(strsplit(cdr3.list,""))))
  count.df <- data.frame("amino.acid" = c(names(count)), "count" = c(as.numeric(count)))
  return(count.df)
}

################################################################################
plotAmino <- function(vdj, seq){
  # get count table of ID and SD from combining countAmino
  invisible(lapply(c("dplyr", "tidyr", "ggplot2", "ggsignif"), library ,
            character.only = TRUE))
  aa.id <- countAmino(as.character(vdj[vdj$dominant == "ID", seq]))
  aa.sd <- countAmino(as.character(vdj[vdj$dominant == "SD", seq]))
  aa.count <- rbind(aa.id, aa.sd)
  aa.count$dominant <- c(rep("ID", dim(aa.id)[1]), rep("SD",dim(aa.sd)[1]))

  # calculate the frequency
  aa.freq <- aa.count %>% 
    group_by(dominant) %>%
    dplyr::mutate(freq = count*100/sum(count))

  # group the amino acid based on hydropathy
  aa <- c("I", "V", "L", "F", "C","M", "A", "W", 
          "G", "T", "S", "Y", "P", "H", 
          "N", "D", "Q", "E", "K", "R")
  hd <- c(rep("Hydrophobic", 8), rep("Neutral",6), rep("Hydrophilic",6))
  hd.index <- data.frame(t(rbind(aa, hd)))
  

  # addd column to aa.freq dataframe
  aa.freq$hd <- hd.index$hd[match(aa.freq$amino.acid, hd.index$aa)]
  
  # plot
  p <- ggplot(aa.freq, aes(x = amino.acid, y = freq, fill = dominant)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(name = "Immnue response", values = c("lightblue", "darkorange1")) +
    labs(x = "Amino acid", y = "Relative frequency (%)", fill = NULL) + 
    scale_y_continuous(breaks = seq(0, max(aa.freq$freq) + 5, 5), 
                       limits = c(0, max(aa.freq$freq) + 5), expand = c(0,0)) +
    custom_theme +
    theme(legend.position = "bottom") +
    facet_grid(~hd, scales = "free_x", space = "free_x")
  
  return(p)
}

################################################################################
# plot nonVJ, V, J #
plotAmino(cd8b.nr, "nonVJ") + 
  labs(title = "Relative frequency of added amino acids at the junctional region") + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        strip.background = element_rect(color = "black"))

cdr3V <- plotAmino(vdj.1b, "Vseq")
cdr3J <- plotAmino(vdj.1b, "Jseq")

plot <- ggarrange(cdr3V, cdr3J, nrow = 2, ncol = 1, 
                            labels = c("CDR3", "non VJ")) 

# calculate the enrichment score for each amino acid #
require(reshape)
aa.enrichment <- dcast(aa.freq, amino.acid ~dominant)
aa.enrichment$ratio <- aa.enrichment[,2]/aa.enrichment[,3]
aa.enrichment$log.ratio <- log(aa.enrichment$ratio)

################################################################################
## An overview of amino acid distribution in naive TCR (peter dataset) ##
################################################################################
# extract only the added nucleotide sequences codon
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

# select only the sequence with nonVJ.length >= 3
ncd8b.ins <- ncd8b.nr[ncd8b.nr$nonVJ.length >= 3,]

# using mapply to pass getInsCodon to increse speed
InsCode <- mapply(getInsCodon, ncd8b.ins$NTseq, ncd8b.ins$Ins)

# export this sequence to be used as input for Python code
write.table(InsCode, file = "naive_pt.txt", 
            append = FALSE, quote = FALSE, sep = "\n",eol = "\n", 
            row.names = FALSE, col.names = FALSE)
# import python output and add as the new column
addedAmino <- read.table("D:/InternshipUU_2017/VDJdb/R_VDJdb/naive_pt_out.txt")
ncd8b.ins$nonVJ <- addedAmino
ncd8b.nr$nonVJ.length <- sapply(ncd8b.nr$nonVJ)

################################################################################
# plot naive added Amino acid #
################################################################################
plotAminoI <- function(peptide.list, data.label, color){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggsignif)
  aa.count <- countAmino(peptide.list)

  # calculate the frequency
  aa.freq <- aa.count %>% 
    mutate(freq = count*100/sum(count))
  print(aa.freq)
  
  # plot the amino acid count 
  p <- ggplot(aa.freq, aes(x = amino.acid, y = freq)) + 
    geom_bar(stat = "identity", position = "dodge", fill = color) +
    labs(x = "Amino acid", y = "Frequency") + 
    annotate("text", x = 18, y = max(aa.freq$freq) - mean(aa.freq$freq)+3, 
             label = data.label, size = 5) + 
    scale_y_continuous(expand = c(0,0)) +
    #title = "Amino acid abundance related to the immune responses") +
    custom_theme 
  return(p)
}
# call function #
plotAminoI(addedAmino[,1], " dataset", "blue")

################################################################################
# plot naive and effector cell 
################################################################################
plotAminoII <- function(naive, vdj){
  invisible(lapply(c("dplyr", "tidyr", "ggplot2", "ggsignif"), library,
            character.only = TRUE))
  # count the occurance of each amino acids
  aa.naive <- countAmino(naive)
  aa.vdj <- countAmino(vdj)
  aa.count <- rbind(aa.naive, aa.vdj)
  
  # add unique name for each datasets
  aa.count$type <- c(rep("alpha", dim(aa.naive)[1]), rep("beta",dim(aa.vdj)[1]))
  
  # calculate count frequency
  aa.freq <- aa.count %>% 
    group_by(type) %>%
    mutate(frequency = count*100/sum(count))
  
  #(chisqDat(aa.freq, "type", "amino.acid"))
  
  # plot the amino acid count 
  p <- ggplot(aa.freq, aes(x = amino.acid, y = frequency,  fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Amino acid", y = "Relative frequency (%)") + 
    scale_fill_manual(name = "T-cell population", labels = c("naive", "non-naive"), 
                      values = c("#6E016B", "#B3DE69")) +
    #title = "Amino acid abundance related to the immune responses") +
    custom_theme +
    theme(legend.position = "bottom")
  return(aafreq)
}

plotAminoII(addedAmino[,1], cd8b.nr$nonVJ) + 
  annotate(geom = "text", x = 18, y = 0.2, label = "non VJ region", size = 5)
plotAminoII(naive.pt8b.ins$AAseq, vdj.1b$CDR3) + 
  annotate(geom = "text", x = 18, y = 0.2, label = "CDR3 sequences", size = 5)


################################################################################
## An overview of amino acid distribution ##
################################################################################
# Calculate enrichment score #
calFrequency <- function(peptide.list){
  # convert factor to character #
  if(is.factor(peptide.list) == TRUE){
    peptide.list <- as.character(levels(peptide.list))[peptide.list]
  } 
  # calculate enrichment score for each position #
  amino.acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y",
                  "C","W","R","K","D","E")
  # split aa in a peptide string to each position (column)
  count.aa <- t(sapply(peptide.list, function(s){unlist(strsplit(s,""))})) # aa x position
  # count aa for each position of peptide #
  count.all <- apply(count.aa, 2, function(aa) sapply(amino.acid, 
                                                      function(x) x<-sum(x==unlist(aa,""))))
  # calculate the percentage and round to 2 decimal point #
  count.perc <- round(prop.table(count.all, margin = 2)*100,1)
  return(count.perc)
}

# calulate log enrichment score #
calLogScore <- function(vdj.hla,length, sd.cutoff, response){
  # select the length of peptide #
  vdj.hla <- vdj.hla[vdj.hla[,'CDR3.length'] == length, ]
  
  # discard the SD with low frequency
  vdj.sel <- vdj.hla[vdj.hla[,'frequency'] >= sd.cutoff, ]
  print(summary(vdj.hla[vdj.hla[,'frequency'] >= sd.cutoff,"dominant"]))
  
  # find perc peptide in ID and SD 
  perc <- calFrequency(vdj.sel[vdj.sel$dominant == response, "CDR3"])
  # select the SD with frequency more than sd.cutoff
  ##sd.perc <- calFrequency(vdj.sel[vdj.sel$dominant == "SD", type])
  ##log.ratio <- log(id.perc/sd.perc)

  return(perc)
}

################################################################################
# heatmap #
################################################################################
calRatio <-  function(vdj.hla, length, sd.cutoff){
  id.freq <- calLogScore(vdj.hla, length, sd.cutoff, "ID")
  sd.freq <- calLogScore(vdj.hla, length, sd.cutoff, "SD")
  log.ratio <- round(log(id.freq/sd.freq),2)
  # replace inf and -inf 
  log.ratio[log.ratio == -Inf] <- -5  # extreme minus value when aa only in SD
  log.ratio[log.ratio == Inf] <- 5    # extreme minus value when aa only in ID
  log.ratio[is.nan(log.ratio)] <- NA
  return(log.ratio)
}

################################################################################
# add a column of prop; hydropathy/size 
identAminoHydropathy <- function(aa.vector){
  if(aa.vector %in% c("A", "F", "C", "I", "L", "M","V","W")){
     hd <- "Hydrophobic"
    } else if(aa.vector %in% c("G", "H", "P", "S","T","Y")){
      hd <- "Neutral"
    } else if(aa.vector %in% c("D", "E","K", "N", "Q", "R")){
      hd <- "Hydrophilic"
    } else {
      hd <- "-"
    }
  return(hd)
}

# sapply(aa.vector, identHydropathy)

identAminoSize <- function(aa.vector){
  if(aa.vector %in% c("F", "Y", "W")){
    hd <- "Very large"
  } else if(aa.vector %in% c("I", "K", "L", "M", "R")){
    hd <- "Large"
  } else if(aa.vector %in% c("E", "H", "Q", "V")){
    hd <- "Medium"
  } else if(aa.vector %in% c("C", "D", "N", "P", "T")){
    hd <- "Small"
  } else if(aa.vector %in% c("A", "G", "S")){
    hd <- "Very small"
  } else {
    hd <- "-"
  }
  return(hd)
}

HM.prop <- function(ratio, prop){
  # import required packages
  invisible(lapply(c("data.table", "ggplot2", "RColorBrewer", "reshape"), 
         library, character.only = TRUE))

  # transform matrix to dataframe format
  df.ratio <- as.data.frame(ratio, row.names = c(gsub("-*","", row.names(ratio))))
  # set amino acid (rowname) as the first column
  df.ratio2 <- setDT(df.ratio, keep.rownames = T)
  colnames(df.ratio2) <- c("CDR3", colnames(df.ratio)[-1])
  #sapply(seq(1,dim(df.ratio2)[2]-1), as.character)
  # reshape data from wide to long format
  df.ratio3 <- melt(df.ratio2, id.vars = "CDR3", 
                    variable.name = "Position", value.name = "Amino.acid")
  df.ratio3$Amino.acid <- unlist(df.ratio3$Amino.acid)
  # define color 
  if(prop == "hydropathy"){
    df.ratio3$Amino.acid <- factor(sapply(df.ratio3$Amino.acid, identAminoHydropathy))
    sfm <- scale_fill_manual(name = "Hydropathy", values = c("-" =  "lightgray", "Hydrophilic" = "#BEBADA", 
                                 "Hydrophobic" = "#66C2A5", "Neutral" = "#FFFFB3"))
  } else if(prop == "size") {
    df.ratio3$Amino.acid <- factor(sapply(df.ratio3$Amino.acid, identAminoSize))
    sfm <- scale_fill_manual(name = "Size", values = c("-" = "lightgray", "Very large" = "gold",
                                  "Large" = "orange", "Medium" = "#B3DE69", "Small" = "#810F7C",
                                  "Very small" = "#4D004B"))
  } else if(prop == "amino.acids") {
    sfm <- scale_fill_manual(values = c(sort(colorRampPalette(brewer.pal(11, "Spectral"))(11)), 
                                        brewer.pal(9, "Set3"), "lightgray"))
  }
  # create heatmap using geom_tile()
  p <- ggplot(df.ratio3, aes(Position, CDR3, fill = Amino.acid)) + 
    geom_tile(aes(fill = Amino.acid), color = "white", size = 0.01) + 
    #scale_fill_manual(name = "Properties", values = color.code) + 
    scale_x_discrete(labels = sapply(1:dim(df.ratio)[2], as.character)) +
    sfm +
    labs(title = "Amino acid distribution", y = NULL) +
    custom_theme +
    theme(axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x=element_blank())
  return(p)
}

t <- positionalCount(a2GILhighID.msa)
enrich.sd <- rownames(t)[which(t[,8] == "R")]
# remove gap in CDR3
sd.cdr3 <- sapply(enrich.sd, function(x) gsub('-*','',x))
enrichsd.df <- cd8b.nr[cd8b.nr$CDR3 %in% sd.cdr3, ]
HM.prop(t, "size")

################################################################################
permutationTest <- function(df, iter, parm){
  library(dplyr)
  parm.id <- pull(df[df$dominant == "ID", parm])
  mean.id <- mean(parm.id)
  
  # generate SD subset
  parm.sd <- pull(df[df$dominant == "SD", parm])
  # define the  number of random sampling (eqaul to ID data)
  n <- length(parm.id)
  # initialise a vector for keeing the mean difference value
  rand.dist <- c()
  for (i in seq(iter)){
    # set pseudo random nyumber generator for a sake of reproducibility
    set.seed(1+i)
    # call random sampling and keep the result in the datalist
    rand.parm <- sample(parm.sd, n, F)
    
    # calculate mean difference
    mean.sd <- mean(rand.parm)
    rand.dist <- c(rand.dist, mean.id-mean.sd)
  }

  p.value <- sum(rand.dist > (mean(parm.id) - mean(parm.sd)))/iter

  return(p.value)
}

################################################################################
HM <- function(log.ratio){
  require(data.table)
  require(ggplot2)
  require(RColorBrewer)
  require(reshape)
  # transform matrix to dataframe format
  df.ratio <- data.frame(log.ratio)
  # set amino acid (rowname) as the first column
  df.ratio2 <- setDT(df.ratio, keep.rownames = TRUE)[]
  colnames(df.ratio2) <- c("Amino acids", colnames(log.ratio))
  #sapply(seq(1,dim(df.ratio2)[2]-1), as.character)
  # reshape data from wide to long format
  df.ratio3 <- melt(df.ratio2, id.vars = "Amino acids", 
                    variable.name = "Position", value.name = "Log enrichment Score")
  
  # reoreder amino acid factor based on hydropathy
  df.ratio3$`Amino acids` <- factor(df.ratio3$`Amino acids`, 
         levels = c("A", "F", "C", "I", "L", "M","V","W", 
                    "G", "H", "P", "S","T","Y", 
                    "D", "E","K", "N", "Q", "R"))
  
  #return(df.ratio3)}
  # create heatmap using geom_tile()
  p <- ggplot(df.ratio3, aes(`Amino acids`, Position)) + 
    geom_tile(aes(fill = `Log enrichment Score`), colour = "white") + 
    scale_fill_gradient2(name = "pLES\n", low = "steelblue",mid = "white", high = "red", 
                         na.value = "snow3" ,limits=c(-5,5)) + 
    labs(title = "Log enrichment score per position (pLES) of CDR3 peptide") +
    custom_theme
  
  return(p)
}

HM(log(positionalFrequency(cd8bhlaa1ID.msa)/positionalFrequency(cd8bhlaa1SD.msa)))
HM(log(positionalFrequency(cd8bhlaa2ID.msa, T)/positionalFrequency(cd8bhlaa2SD.msa, T)))

HM(log(positionalFrequency(a2NLVID.msa, T)/positionalFrequency(a2NLVSD.msa, T))) + 
  labs(title = "NLV")
HM(log(positionalFrequency(a2GILID.msa, T)/positionalFrequency(a2GILSD.msa, T))) + 
  labs(title = "GIL")
HM(log(positionalFrequency(a2GLCID.msa, T)/positionalFrequency(a2GLCSD.msa, T))) + 
  labs(title = "GLC")

################################################################################
# plot heatmap #
plotHM <- function(vdj.hla, length, cutoff){
  vdj.nr <- vdj.hla[!duplicated(vdj.hla[,c("CDR3", "dominant")]),]
  vdj.les <- calRatio(vdj.nr,length, cutoff)
  return(HM(vdj.les))
}

# call for heatmap (finally!)
plotHM(vdj.1b[!duplicated(vdj.1b[,c("CDR3", "dominant")]),],
       length = 13, 0.5)

################################################################################
# calculate the sum of log enrichment score
################################################################################
# log enrichment freq of each amino acid regardless position in CDR3 
# (unlike calLogScore)
calLogScore.I <- function(vdj.hla){
  # remove non-reduncdant CDR3
  vdj.nr <- vdj.hla[!duplicated(vdj.hla[,c("CDR3", "dominant")]),]
  # sum occurance of each amino acid in sd and id response
  aa.id <- countAmino(vdj.nr[vdj.nr$dominant == "ID", "CDR3"])
  aa.id$perc <- round(aa.id$count*100/sum(aa.id$count),1)
  print(aa.id)
  aa.sd <- countAmino(vdj.nr[vdj.nr$dominant == "SD", "CDR3"])
  aa.sd$perc <- round(aa.sd$count*100/sum(aa.sd$count),1)
  print(aa.sd)
  log.score <- data.frame("Amino acid" = aa.id[,1], 
                          "Enrichment" = round(aa.id$perc/aa.sd$perc,2),
                          "Log enrichment" = round(log(aa.id$perc/aa.sd$perc),2))
  return(log.score)
}

################################################################################
# fit a continuous distribution to infer the importance of each position in CDR3
# since the max CDR3 length = 23 and position 6-8 are higly detected in contact #
fitChisq <- function(df.value, vjust.value, max){
  require(ggplot2)
  df <- data.frame(x = c(0:max), y = dchisq(seq(0,max), df = df.value))
  p <- ggplot(df , aes(x = x, y = y, label = y)) +
    stat_function(fun = dchisq, args = list(df = df.value)) +
    geom_point(color = "red", size = 3) +
    geom_text(aes(label = round(y,3), vjust = vjust.value)) + 
    #geom_segment(data=df, aes(xend=-Inf, yend=y)) +
    geom_segment(data=df, aes(xend=x, yend=-Inf),
                 lty = "dotted", color = "red", size = 1.0) +
    #geom_vline(aes(xintercept = 0:23), lty = "dotted", color = "red", size = 0.5) +
    scale_x_continuous(breaks = seq(0,max), labels = seq(0,max)+1) +
    labs(x = "Position of amino acids in CDR3 peptide", y = "Density",
         title = "Chi-square distribution (df = ?)") +
    custom_theme +
    theme(legend.position = "None")
    
  return(p)
}

# call for fitting #
fitChisq(8, 0.01, 22)


################################################################################
# calculate frequency per position to use as weight #
################################################################################
# saperate each amino acid based on the position in CDR3 sequences
positionalCount <- function(cdr3){
  count <- list()
  max.length <- max(sapply(as.character(cdr3), nchar))
  for(i in seq(max.length)){
    count[[i]] <- as.list(sapply(cdr3, substr, i, i, simplify = FALSE))
  }
  ##count.mat <- t(do.call(rbind, count)) # give a matrix of list instead  of character
  count.mat <- matrix(unlist(count), nrow = length(cdr3))
  rownames(count.mat) <- cdr3
  return(count.mat)
}

positionalFrequency <- function(cdr3, filter.lowoccupancy){
  count.mat <- positionalCount(cdr3)
  
  # count aa for each position of peptide #
  amino.acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S",
                  "T","Y","C","W","R","K","D","E")
  
  count.all <- apply(count.mat, 2, function(aa) sapply(amino.acid, 
                                                    function(x) x<-sum(x==unlist(aa,""))))
  # assgin the colume names (== position of amino acids)
  colnames(count.all) <- as.character(seq(dim(count.all)[2]))
  
  if(filter.lowoccupancy == T){
    # remove the column with occupancy < 0.50
    count.all <- count.all[,colSums(count.all)/length(cdr3) > 0.50]
  }
  ## print occupancy in each column-position
  ##print(apply(count.all, 2, sum)/length(cdr3))

  # calculate the percentage and round to 2 decimal point #
  count.tb <- prop.table(count.all, margin = 2) #*100 after =2 ) to get percentage
  return(count.tb)
}

pLES <- function(vdj, f){
  pf.id <- positionalFrequency(vdj[vdj$dominant == "ID", "CDR3"], 0)
  pf.sd <- positionalFrequency(vdj[vdj$dominant == "SD", "CDR3"], f)
  
  # pLES
  log.ratio <- round(log(pf.id/pf.sd),2)
  
  # replaace -inf (log(id == 0)/(sd != 0)) and inf (log(id != 0)/(sd == 0))
  log.ratio[log.ratio == -Inf] <- -5  # extreme minus value when aa only in SD
  log.ratio[log.ratio == Inf] <- 5    
  
  return(log.ratio)
}

t <- pLES(vdj.1b)
t <- rbind(t, apply(t,2,sum, na.rm = T))
# 
################################################################################

positionalBar <- function(freq, prop){
  require(data.table)
  require(ggplot2)
  require(RColorBrewer)
  require(reshape)
  # transform matrix to dataframe format
  df <- data.frame(freq)
  # set amino acid (rowname) as the first column
  df2 <- setDT(df, keep.rownames = TRUE)[]
  colnames(df2) <- c("Amino acids", colnames(freq))
  # reshape data from wide to long format
  df3 <- melt(df2, id.vars = "Amino acids", 
                    variable.name = "Position", value.name = "Frequency")
  if(prop == "hydropathy"){
    df3$`Amino acids` <- factor(df3$`Amino acids`, levels = c("A", "F", "C", "I", "L",
                                                              "M","V","W", "G", "H", "P",
                                                              "S","T","Y", "D", "E","K",
                                                              "N", "Q", "R"))
    color.code <- c(rep("#9f4867", 8), rep("#5fb951", 6), rep("#a25acd",6))
    name.code <- c("Hydrophobic", "Neutral", "Hydrophilic")
  } else if(prop == "size") {
    df3$`Amino acids` <- factor(df3$`Amino acids`, levels = c("F", "Y", "W", "I", "K",
                                                              "L","M","R", "E", "H", "Q",
                                                              "V","C","D", "N", "P","T",
                                                              "A", "G", "S"))
    color.code <- c(rep("gold", 3), rep("orange", 5), rep("#B3DE69",4),
                    rep("#810F7C", 5), rep("#4D004B", 3))
    name.code <- c("Very large", "Large", "Medium", "Small", "Very small")
  } else if(prop == "amino.acids") {
    color.code <- sort(sort(c(colorRampPalette(brewer.pal(11, "Spectral"))(11), brewer.pal(9, "Set3"))))
  }
  # extrapolate the color
  ##getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  # create stacked bar plot
  p <- ggplot(df3, aes(x = Position, y = Frequency, fill = `Amino acids`)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = color.code) + #, label = name.code) + 
    labs(title = "Amino acid frequency in each IMGT position of CDR3 peptides") +
    custom_theme
  
  return(p)
}

positionalBar(positionalFrequency(cd8b.nr$CDR3))
################################################################################

sumLogEnrichment <- function(cdr3, log.score){
  # specify weigth 
  weight.df <- data.frame(x = c(0:38)+1, y = dchisq(seq(0,38), 18))
  # calculate sum of log score * weight of CDR3 peptide
  score <- sum(sapply(sapply(cdr3, function(x) {unlist(strsplit(x,""))}), function(a) 
    {log.score[log.score$Amino.acid == a, "Log.enrichment"]}) * weight.df[1:nchar(cdr3),"y"])
  
  return(score)
}

################################################################################
# cross-validation (CV) fn taking a list of train index at a time and return the
# log enrichment score (les) of the test set
cv <- function(train.index, vdj.nr){
  # devide data to train and test set
  train.set <- vdj.nr[-train.index,]
  test.set <- vdj.nr[train.index,]
  
  # calculate score taken into account the log enrichment ratio and weight matrix
  log.score <- calLogScore.I(train.set)
  sum.score <- sapply(as.character(test.set$CDR3), sumLogEnrichment, log.score, simplify = F)
  test.set$les.score <- as.numeric(sum.score)
  
  return(test.set)
}

################################################################################

modelTest.chi <- function(vdj.hla, kfold){
  library(caret)
  library(klaR)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  
  # remove repetitive CDR3
  ##vdj.nr <- vdj.hla[!duplicated(vdj.hla[,c("CDR3", "dominant")]),]
  
  # convert from fcator to character
  vdj.hla[, c("CDR3", "dominant")] <- sapply(vdj.hla[, c("CDR3", "dominant")], as.character)
  
  # create 5-fold CV dataset
  folds <- createFolds(vdj.hla[, "dominant"], k = kfold, list = TRUE, returnTrain = FALSE)
  
  ##trainIndex <- createDataPartition(vdj.nr[,"dominant"], p = 0.80, list = F
  
  # train the model and show the les of testset each fold
  sapply(folds, function(f) {
    test.set <- cv(f, vdj.hla)
    p <- ggplot(test.set, aes(x = dominant, y = les.score)) + 
      geom_jitter(aes(color = dominant), size = 3, width = 0.20) +
      geom_boxplot(alpha = 0.4, width = 0.40) +
      stat_compare_means(label.x = 1.4, size = 5) +
      labs(x = "Immune response", y =  "Sum of log enrichment scores") +
      custom_theme
    print(p)
    })
}

################################################################################
sumpLES <- function(cdr3, pLES.mat){
  # calculate the pLES score for  each cdr3 using pLES.mat 
  cdr3.list <- unlist(strsplit(cdr3,''))
  sum.pLES <- sum(sapply(cdr3.list, function(a) (pLES.mat[a,match(a, cdr3.list)]), simplify = T))
  return(sum.pLES)
}

ples.matA2 <- log(positionalFrequency(cd8bhlaa2ID.msa)/positionalFrequency(cd8bhlaa2SD.msa))


modelTest.pLES <- function(vdj, kfold){
  library(caret)
  library(klaR)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  
  # convert from fcator to character
  vdj[, c("CDR3", "dominant")] <- sapply(vdj[, c("CDR3", "dominant")], as.character)
  
  # create 5-fold CV dataset
  folds <- createFolds(vdj[, "dominant"], k = kfold, list = TRUE, returnTrain = FALSE)
  
  # train the model and show the les of testset each fold
  sapply(folds, function(f){
    # calculate ples matrix of from the training data
    msa.id <- alignGap(vdj[-f,], "ID")
    msa.sd <- alignGap(vdj[-f,], "SD")
    pLES.mat <- log(positionalFrequency(msa.id, F)/positionalFrequency(msa.sd, F))
    pLES.mat[pLES.mat == -Inf] <- -5  # extreme minus value when aa only in SD
    pLES.mat[pLES.mat == Inf] <- 5 
    # add  row containing zero for gap "-"
    pLES.mat <- rbind(pLES.mat, rep(0, dim(pLES.mat)[2]))
    rownames(pLES.mat)[21] <- "-"

    # test with the test set
    test.set <- vdj[f,]
    test.msa <- sapply(as.character(test.set[, "CDR3"]),addGap)
    test.msa <- moveGap(test.msa)
    test.set$les.score <- sapply(as.character(test.msa), sumpLES, pLES.mat)
    
    p <- ggplot(test.set, aes(x = dominant, y = les.score)) + 
      geom_jitter(aes(color = dominant), size = 3, width = 0.20) +
      geom_boxplot(alpha = 0.4, width = 0.40) +
      #stat_compare_means(label.x = 1.4, size = 5) +
      labs(x = "Immune response", y =  "Sum of log enrichment scores") +
      custom_theme +
      theme(legend.position = "None")
    
    print(permutationTest(na.omit(test.set), 1000, "les.score"))
  })
}

################################################################################
# compare alpha and beta chains
assignHydropathy <- function(aa){
  if(aa %in% c("A", "F", "C", "I", "L", "M","V","W")){
    return("Hydrophobic")
  } else if(aa %in% c("G", "H", "P", "S", "T", "Y")){
    return("Neutral")
  } else 
    return("Hydrophilic")
}
compare_cd8ab$hydropathy <- sapply(compare_cd8ab$amino.acid, assignHydropathy)

ggplot(compare_cd8ab, aes(x = reorder(amino.acid, log) , y = log, fill = hydropathy)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(name = "Hydropathy", values = c("#66C2A5", "#FFFFB3", "#BEBADA")) +
  labs(x = "Amino acid", y = "Log enrichment ratio") +
  custom_theme +
  theme(legend.position = c(0.8,0.3),
        legend.background = element_rect(fill=alpha('white', 0.5)),
        legend.text = element_text(size=13))











