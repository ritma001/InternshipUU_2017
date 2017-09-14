#!usr/bin/r
# Author: Wannisa Ritmahan
# Established Date: 05/09/2017 Last Modification: 11/09/2017
# Aim: To infer statistically significance of the 4-mer motif by pair permutation
#
# Method: randomly dividing the data into datasets with the size equals to 544 and 
# 5728 samples (to resemble the size of IDD and SDD datasets, respectively)
################################################################################
# install packages
#
#install.packages(plyr, dependencies = TRUE)
# create 4-mer dataframe
genMotif <- function(){
  library(plyr, lib.loc = "/home/wannisa/R/x86_64-pc-linux-gnu-library/3.4/")
  # provide a list of amino acids
  amino_acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
  # generate 20*20*20*20 4D-array 
  ar <- array(rep(0,160000), dim = c(20,20,20,20), dimnames = list(amino_acid,amino_acid,amino_acid, amino_acid))
  # convert array to dataframe
  df <- adply(ar,c(1,2,3,4))
  # concatenate the 4-amino acid motifs
  four.mer <- apply(df[,-ncol(df)],1,paste,collapse = "")
  
  return(four.mer)
}

# pair permutation test
pairPerm <- function(table, chain , dominant, count.mat, four.mer){
  # calculate the total number of IDD, SDD of specified TCR chain
  sum.dominant <- nrow(table[table$dominant == dominant & table$Gene == chain, ])

  # selected the dataframe containing specified chain
  df <- table[table$Gene == chain,]
  
  print("Start looping")
  #  run count4mer2 function for each dataset
  iter = 1
  while(iter < 1001){
    # randomly divide the data to 2 set 
    sel.row <- sample(1:nrow(df), sum.dominant)
    df1 <- df[sel.row,"CDR3" ]
    # remove redundant sequences 
    cdr3 <- levels(droplevels(df1))
    
    # keep the frequency of each iteration
    for(a in cdr3){
      i = 1
      while(i < nchar(a)-3){
        mer <- substr(a,i,i+3)
        # get the corresponding motif index
        ind <- match(mer, four.mer)
        count.mat[ind, 1] <- count.mat[ind, 1] +1  
        # change[ind, 1] to [ind, iter] to keep count from each iteration in separate column
        i = i+1 
      }
    }
    print(paste("end of loop ", iter))
    iter = iter + 1
  }
  # delete the row containing no count and make the dataframe with the corresponding motif added 
  names <- which(rowSums(count.mat)!=0)
  df <- data.frame(motif = factor(four.mer[names]), count = count.mat[names,])
  rownames(df) <- names
  
  # write the output as csv file 
  return(df)
}

# call these outside the pp function
# generate matrix to keep data
count.mat <- matrix(0L, nrow = 160000, ncol = 1)
# generate all possilbe 4-mer motif
four.mer <- genMotif()

# call pairPerm() and export the table in text file
args = commandArgs(trailingOnly = TRUE) 
vdj2 <- read.csv(args, header = TRUE, row.names = 1)

ppa.idd <- pairPerm(vdj2, "TRA", "IDD", count.mat, four.mer)
write.csv(ppa.idd, file = "ppA_IDD.csv")

ppb.idd <- pairPerm(vdj2, "TRB", "IDD", count.mat, four.mer)
write.csv(ppb.idd, file = "ppB_IDD.csv")

ppa.sdd <- pairPerm(vdj2, "TRA", "SDD", count.mat, four.mer)
write.csv(ppa.sdd, file = "ppA_SDD.csv")

ppb.sdd <- pairPerm(vdj2, "TRB", "SDD", count.mat, four.mer)
write.csv(ppb.sdd, file = "ppB_SDD.csv")
