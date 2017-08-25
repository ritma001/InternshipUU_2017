#!usr/bin/r 
# Author: Wannisa Ritmahan 
# Date: 23-25 August,2017
# Aim: to learn to write R code to manipulate data retrieved from the VDJdb
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

## Identify IDD and SDD for each subject_id
identIDD <- function(table){
  table[,'dominant'] <- NA
  # iterate for all samples in each subject to find the max frequency
  for(i in levels(table$subject.id)){
    #get the row indicies foreach subject ID
    nrow <- which(table$subject.id == i)
    IDD <- max(table$frequency[nrow])
    #get the row index containing max frequency for each subject
    indmax <- which(table$frequency[nrow] == IDD)
    #assign IDD and SDD in the 'dominant column'
    table[nrow[indmax],'dominant'] <- 'IDD'
    table[nrow[-indmax],'dominant'] <- 'SDD'
  }
  table[,'dominant'] <- as.factor(table[,'dominant'])
  return(table)
}

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

## randomly select the CDR3 squences for MSA
randCDR <- function(vdj, species,n){
  #library('dplyr')
  t <- spcFilter(vdj, species, MHC)
  rand <- sample_n(t, as.numeric(n))
  return(rand)
}

## write the table as in fasta format
exportTable <- function(t,filename){
  export <- paste(paste(paste(paste(">",as.character(row.names(t)),sep=""),t$dominant, sep = ""),t$subject.id,sep="") ,t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}


################################################################################
## pass arguments from the command lines
args = commandArgs(trailingOnly = TRUE) 
# args[1] == VDJ data, args[2] == species, args[3] == MHC class

## Filter rawdata by species and MHC
table <- freqFilter(args[1]) #'vdj220817.txt'

## extract frequency and subject.id 
table$Method <- parseFreq(table)
table$Meta <- factor(parseSubject(table))
colnames(table)[14:15] <- c('frequency','subject.id')

## format frequency to numeric value
for(i in 1:nrow(table)){
  table[i,14] <- freqFormat(table[i,14])
}
table <- na.omit(table)

## identify IDD and SDD
table <- identIDD(table)

## select species and MHC class
spc_table <- spcFilter(table, args[2], args[3])

## write table 
exportTable(spc_table, args[4])
















