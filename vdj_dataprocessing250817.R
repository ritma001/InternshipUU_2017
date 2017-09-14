#!usr/bin/r 
# 
# Author: Wannisa Ritmahan 
# Established Date: 28/08/2017  Last modification: 03/09/2017
# Aim: R code for manipulating data from the VDJdb for Immunodominance studies
# 
# To be improved: 
# 1. check consistency with of word i.e. table, t, vdj etc.
# 2. incorperate the code with bash script to reduce repetitive work
#
################################################################################
## Part I: Data processing ##
################################################################################

## read vdj.csv or .txt file
readTable <- function(file){
  t <- read.table(file,header = TRUE, sep = '\t',quote = "\"",dec = ".")
  return(t)
}

## exclude the data containing missing frequency
freqFilter <- function(vdj){
  ind <- c()
  t <- readTable(vdj)
  for(i in 1:nrow(t)){
    logi <- pmatch('{frequency: ,',t[i,'Method'])
    if(is.na(logi)){
      ind = c(i,ind)
    }
  }
  ft <- t[ind,]
  return(ft)
}

## extract frequency
parseFreq <- function(table){
  freq = c()
  for(i in table[,14]){
    f_value <- gsub(".frequency: |,.*","",i)
    freq = c(freq,f_value)
  }
  return(freq)
}

## extract subject id
parseSubject <- function(table){
  id = c()
  for(i in table[,15]){
    subject_id <- gsub(".*subject.id: |,.*","",i)
    if(subject_id != ""){
      id = c(id,subject_id)
    } else {
      id = c(id,NA)
    }
  }
  return(id)
}

## check frequency format and convert into number
freqFormat <- function(str){ 
  s1  <- sapply(strsplit(str,'/|%'),as.numeric)
  if(length(s1) == 1){
    number = s1[1]
  } else {
    number = s1[1,1]/s1[2,1]
  }
  return(number)
}

## generate the CDR3 & epitope length vector then, add to the table
aminoLength <- function(table){
  lenCDR3 <- nchar(as.character(table$CDR3))
  table[,'CDR3.length'] <- lenCDR3
  
  lenEpitope <- nchar(as.character(table$Epitope))
  table[,'Epitope.length'] <- lenEpitope
  return(table)
}

################################################################################
## Part II: Identification of immunodominant- & subdominant epitopes ##
################################################################################

## Identify IDD and SDD for each subject_id
identIDD <- function(table){
  table[,'dominant'] <- NA
  # iterate for all samples in each subject to find the max frequency
  for(i in levels(table$subject.id)){
    # get the row indicies foreach subject ID
    nrow <- which(table$subject.id == i)
    IDD <- max(table$frequency[nrow])
    # get the row index containing max frequency for each subject
    indmax <- which(table$frequency[nrow] == IDD)
    # assign IDD and SDD in the 'dominant column'
    table[nrow[indmax],'dominant'] <- 'IDD'
    table[nrow[-indmax],'dominant'] <- 'SDD'
  }
  table[,'dominant'] <- as.factor(table[,'dominant'])
  return(table)
}

################################################################################
## Part III: Create several datasets with different criteria ##
################################################################################

## filter the epitope.species and MHC class (optional)
spcFilter <- function(table, species, MHC){
  t <- table
  if(missing(MHC)){
    MHCI <- t[t$Epitope.species == species & t$MHC.class == 'MHCI',]
    MHCII <- t[t$Epitope.species == species & t$MHC.class == 'MHCII',]
    st <- rbind(MHCI,MHCII)
  } else {
    st <- t[t$Epitope.species == species & t$MHC.class == MHC,]
  }
  return(st)
}

## randomly select the CDR3 squences for MSA (less than 1000 sequences are accepted)
randCDR <- function(vdj, species,n){
  #library('dplyr')
  t <- spcFilter(vdj, species, MHC)
  rand <- sample_n(t, as.numeric(n))
  return(rand)
}

## create the 'IDD' and 'SDD ' datasets with option to exclude human sample
parseISDD <- function(table, epitope){
  df <- table[table$dominant == epitope & table$Epitope.species != 'HomoSapiens',]
  return(df)
}

## create the different dataset using gene (TRA-alpha or TRB-beta chain) as criteria
parseAB <- function(table, chain){
  df <- table[table$Gene == chain,]
  return(df)
}

## combine the genetic regions of V, J and MHC.A after the * ex. HLA-A*02 & HLA-A*02:01 
## are grouped to HLA-A*02
combineRegion <- function(table, parm){
  # convert factor to character column
  table[, parm] <- as.character(table[, parm])
  result = c()
    for(i in table[,parm]){
      if(parm == "V" | parm == "J"){
        new.value <- gsub("[\\*|/|-].*$","",i) # for V and J loci
        result = c(result,new.value)
      } else {
        new.value <- gsub(":.*$","",i) # for MHC.A
        result = c(result,new.value)
      }
    }
  table[, parm] <- factor(result)
  return(table)
}

################################################################################
## Part IV: Write CDR3 with custom header as .fasta file ##
################################################################################

## write the table as in fasta format I
exportTable <- function(t,filename){
  export <- paste(paste(paste(paste(">",as.character(row.names(t)),sep=""),t$dominant, sep = ""),t$subject.id,sep="") ,t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

## write the table in fasta format II (clearer headder than format I)
exportTable2 <- function(t,filename){
  export <- paste(paste(paste(">",t$dominant, sep = ""),1:length(t$dominant),sep="") ,t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

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
    geom_bar(aes(y=rel.freq),stat = "identity", position = "dodge") + 
    labs(title = main, y = "Relative frequency") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p + facet_wrap(~dominant, ncol = 1, scales = "free_y")
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
    geom_bar(aes(y=rel.freq, fill = dominant),stat = "identity", position = "dodge") + 
    labs(title = main, y = "Relative frequency") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p
}
# call function
plotrelFreq2(relFreq(vdj2.a, "V"), "V" , "V distribution of alpha chain")
plotrelFreq2(relFreq(vdj2.b, "V"), "V" , "V distribution of beta chain")

########################s########################################################
## Execution of the code above ##
################################################################################

## pass arguments from the command lines
args = commandArgs(trailingOnly = TRUE) 
# args[1] == VDJ data, args[2] == species, args[3] == MHC class

## Filter rawdata by species and MHC
table <- freqFilter('vdj220817.txt') #'vdj220817.txt #args[1]'

## extract frequency and subject.id 
table$Method <- parseFreq(table)
table$Meta <- as.factor(parseSubject(table))
colnames(table)[14:15] <- c('frequency','subject.id')

## format frequency to numeric value
for(i in 1:nrow(table)){
  table[i,14] <- freqFormat(table[i,14])
}
table <- na.omit(table)

## identify IDD and SDD
table <- identIDD(table)

## check the length of CDR3 motif and epitopes
table <- aminoLength(table)

## select species and MHC class/ no human epitopes
#spc_table <- spcFilter(table, args[2], args[3])
vdj2 <- table[table$Epitope.species != "HomoSapiens",]

# group V, J and MHC region 
vdj2.j <- combineRegion(vdj2, "J")
vdj2.v <- combineRegion(vdj2.j, "V")
vdj2 <- combineRegion(vdj2.v, "MHC.A")

# generae 2 datasets based on TCR chain
vdj2.a <- parseAB(vdj2,"TRA")
vdj2.b <- parseAB(vdj2,"TRB")

vdj2a.idd <- vdj2.a[vdj2.a$dominant == "IDD",]
vdj2a.sdd <- vdj2.a[vdj2.a$dominant == "SDD",]
vdj2b.idd <- vdj2.b[vdj2.b$dominant == "IDD",]
vdj2b.sdd <- vdj2.b[vdj2.b$dominant == "SDD",]

# generate 2 data sets based on dominant types
vdj2.idd <- vdj2[vdj2$dominant == "IDD",]
vdj2.sdd <- vdj2[vdj2$dominant == "SDD",]

################################################################################
## get a vector of non-redundnat CDR3 sequences
################################################################################
# generate the datasets with non-repetitive sequences
cdr3a.idd.nr <- levels(droplevels(vdj2a.idd$CDR3))
cdr3a.sdd.nr <- levels(droplevels(vdj2a.sdd$CDR3))
cdr3b.idd.nr <- levels(droplevels(vdj2b.idd$CDR3))
cdr3b.sdd.nr <- levels(droplevels(vdj2b.sdd$CDR3))

## write table 
#exportTable(spc_table, args[4])





