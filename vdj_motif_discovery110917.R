#!usr/bin/r
# Author: Wannisa Ritmahan
# Established Date: 11/09/2017 Last Modification: 14/09/2017
# Aim: To explore the distincive motifs in IDD and SDD datasets
#
# Method: After filtering data with missing frequency and excluding human 
# epitopes, the remained dataset contains 6272 samples. The duplicated CDR3 
# sequences are removed then, the non-redundant sequences are search for motifs,
# test for statistic significance and the frequency of motifs in the CDR3 sequences. 
# 
################################################################################
## count (2-4)-mer amino acid ##
# %ToDO (vecterising instead of looping CDR3 to increase speed)
################################################################################
# generate a list of non-reduncdant CDR3 sequences the vdj table should only 
# contain either the "IDD" or "SDD" dataset of specific chain (alpha or beta)
nrCDR3 <- function(vdj){
  cdr3 <- levels(droplevels(vdj$CDR3))
  return(cdr3)
}
################################################################################
# function to generate the matrix or array to recorde the n-mer motif
countTable <- function(n.mer){
  amino.acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
  count.t <- array(0L, dim = rep(20,n.mer))
  dimnames(count.t) <- rep(list(amino.acid),n.mer)
  return(count.t)
}
################################################################################
# get row index (in df-which is the table to keep frequency) of the n-mer motif 
searchMotif <- function(motif, nmer, df){
  i = 1
  inds <- c()
  while(i < nchar(motif)){
    mer <- substr(motif,i,(i+nmer-1))
    ind <- match(mer, df)
    if(!is.na(ind)){
      inds <- c(inds, ind)
    }
    i = i + 1
  }
  return(inds)
}
################################################################################
# count n-mer (2-4) and keep the number in the table generated from countTable()
countMer <- function(vdj, n.mer){
  cdr <- nrCDR3(vdj)
  # count n-mer motifs
  count.ar <- countTable(n.mer)
  # convert matrix to dataframe
  df <- adply(count.ar,seq(n.mer))
  df$motif <- factor(apply(df[,-ncol(df)],1,paste,collapse = ""))
  # get the row index of motif
  motif.ind <- lapply(as.character(cdr), searchMotif, nmer = n.mer, df = df$motif)
  # count while looping
  for(ind in motif.ind){
    ind.ul <- unlist(ind)
    df[ind.ul, "V1"] <- df[ind.ul, "V1"] + 1
  }
  row.remove <- which(df[,(n.mer+1)] == 0)
  result <- df[-row.remove,]
  result <- result[order(-xtfrm(result[,(n.mer+1)])),]
  prop <- result$V1/sum(result$V1)
  result$prop <- prop
  return(result)
}

##countMer gives accurate result (check : sum(vdj2.a[vdj2.a$dominant == "IDD", "CDR3.length"]-(nmer-1)))
# only with 4mer not 2 or 3 ????
# call function

################################################################################
## plot n-mer
################################################################################
# Plot the motif
nameMotif <- function(idd,sdd, i){
  # generate the list of top 20 motifs for both IDD and SDD
  motifs <- c(as.character(idd$motif[i:(i+19)]),as.character(sdd$motif[i:(i+19)]))
  combined.motifs <-  unique(motifs)
  return(combined.motifs)
}

plotMotif <- function(idd,sdd, i, main){
  library(reshape2)
  # generate the list of top 20 motifs for both IDD and SDD
  motifs <-sort(nameMotif(idd,sdd,i))
  idd <- idd[order(idd$motif),]
  sdd <- sdd[order(sdd$motif),]
  cdr4 <- data.frame(motif = motifs, idd = idd[idd$motif %in% motifs, "prop"], sdd = sdd[sdd$motif %in% motifs, "prop"])
  ## print(cdr4)
  cdr4.m <- melt(cdr4, id = "motif")
  #plot
  library(ggplot2)
  cdr4.p <- ggplot(cdr4.m, aes(x=motif,y=value, fill = variable)) + geom_bar(stat = "identity", position = "dodge" ,width = 0.5) + 
  ggtitle(main) + xlab("4-mer amino acid motif") + ylab("Relative Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = "top")
  cdr4.p <- cdr4.p# + facet_grid(~variable)
  
  return(cdr4.p)
}

# call function
  # alpha chain
cdr3a.idd.nr.4m <- countMer(vdj2.a[vdj2.a$dominant == "IDD", ], 4)
cdr3a.sdd.nr.4m <- countMer(vdj2.a[vdj2.a$dominant == "SDD", ], 4)

plotMotif(cdr3a.idd.nr.4m, cdr3a.sdd.nr.4m, 1, "Distribution of 4-mer motifs of non-reduncdant CD3 of alpha chain")

  # beta cahin
cdr3b.idd.nr.4m <- countMer(vdj2.b[vdj2.b$dominant == "IDD", ], 4)
cdr3b.sdd.nr.4m <- countMer(vdj2.b[vdj2.b$dominant == "SDD", ], 4)

plotMotif(cdr3b.idd.nr.4m, cdr3b.sdd.nr.4m, 4, "Distribution of 4-mer motifs of non-reduncdant CD3 of beta chain")

################################################################################
## export sequence to fasta file format
################################################################################
exportSeq <- function(sequences,filename){
  export <- paste(paste(">",1:length(sequences),sep="") ,sequences, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}
################################################################################
## extract motifs
################################################################################
## Alpha chain ##
# extarct the sequences containing certain motifs
a.sdd <- grep("GSQGNLIF", cdr3a.sdd.nr) # 103 from 580 sequences 
motif.a.sdd <- cdr3a.sdd.nr[a.sdd]
exportSeq(motif.a.sdd)

motif.a.sdd.df <- vdj2a.sdd[vdj2a.sdd$CDR3 %in% motif.a.sdd, ]


## Beta chain ##
# extract sequences containing certain motifs
b.idd1 <- grep("GNTIYF", cdr3b.idd.nr) 
motif.b.idd1 <- cdr3b.idd.nr[b.idd1]
motif.b.idd.df1 <- vdj2b.idd[vdj2b.idd$CDR3 %in% motif.b.idd1, ]

b.idd2 <- grep("CASSF", cdr3b.idd.nr) 
motif.b.idd2 <- cdr3b.idd.nr[b.idd2]
motif.b.idd.df2 <- vdj2b.idd[vdj2b.idd$CDR3 %in% motif.b.idd2, ]


