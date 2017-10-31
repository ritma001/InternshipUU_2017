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
  for(i in table[,"Method"]){
    f_value <- gsub(".frequency: |,.*","",i)
    freq = c(freq,f_value)
  }
  return(freq)
}

## extract subject id, the replica.id and subject.cohort represented in "Meta" column in json format 
parseSubject.meta <- function(Meta){
    #subject <- sub(".*?subject.id: (.*?),.*","\\1",Meta)
    clone <- sub(".*?clone.id: (.*?),.*", "\\1", Meta)
    replica <- sub(".*?replica.id: (.*?),.*", "\\1", Meta)
    epitope <- sub(".*?epitope.id: (.*?),.*", "\\1", Meta)
    #samples <- sub(".*?samples.found: (.*?),.*", "\\1", Meta)
    #studies <- sub(".*?studies.found: (.*?),.*", "\\1", Meta)
    cohort <- sub(".*?subject.cohort: (.*?),.*", "\\1", Meta)
    
    id = paste(paste(paste(clone, replica, sep = "*"),epitope, sep = "*"), cohort, sep = "*")
    #if(subject_id != ""){
     # id = paste(paste(paste(subject, replica, sep = "*"), cohort, sep = "*"), clone, sep = "*")
    #} else {
     # id = NA
    #}
  return(id)
}

parseSubject.method <- function(Method){
  identification <- sub(".*?identification: (.*?),.*","\\1",Method)
  verification <- sub(".*?verification: (.*?)}", "\\1", Method)
  method.id = paste(identification, verification, sep = "*")
  return(method.id)
}

## extract the cell.subset to identify either the effector or memory T cell is present 
parseSubset <- function(Meta){
    cell.subset <- sub(".*cell.subset: (.*?),.*","\\1",Meta)
    if(!cell.subset == ""){
      subset = cell.subset
    } else {
      subset = "unavailable"
    }
  return(subset)
}


## extract the jStart and vEnd positon of given CDR3
parseVJ <- function(table){
  table$v.end <- as.numeric(sapply(table$CDR3fix, function(x) sub(".*?vEnd: (.*?),.*","\\1",x)))
  table$j.start <- as.numeric(sapply(table$CDR3fix, function(x) sub(".*?jStart: (.*?),.*","\\1",x)))
  return(table)
}


## parse non-VJ
parse.nonvj <- function(table){
  table$non.VJ <- substr(as.character(table[,"CDR3"]), table$v.end +1, table$j.start-1)
  return(table)
}


## check frequency format and convert into number
freqFormat <- function(str){ 
    s1  <- sapply(strsplit(str,'/|%|//'),as.numeric)
    if(length(s1) == 1){
      number = s1[1]
    } else {
      number = s1[1,1]*100/s1[2,1]
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
## add the denominator to the subject 
getDenom <- function(table){
  # check multiple measurement of the same sample
  row.name = grep(paste(c("/", "//"),collapse="|"), table$raw.frequency)
  table$subject.id <- as.character(table$subject.id)
  # add the denominator to subject id
  for(i in row.name){
    denominator = unlist(strsplit(table[i, "raw.frequency"],"/|//"))
    id = paste(as.character(table$subject.id[i]), denominator[2], sep = "*")
    table$subject.id[i] <- id 
  }
  table$subject.id <- factor(table$subject.id)
  return(table)
}

mergeFreq <- function(table){
  library(dplyr)
  # deal with  subjec.id in which total frequency > 100%
  # so, compile based on the CDR3 seq 
  checksum <- aggregate(table$frequency, by=list(subject.id=table$subject.id), sum)
  merge.subject <- as.character(checksum[checksum$x > 100, "subject.id"])
  drop.row = c()
  for(i in merge.subject){
    #print(i)
    sum.CDR3 <- aggregate(frequency ~ CDR3, table[table$subject.id == i,], sum)
    merge.row <- table[table$subject.id == i,]
    merge.row <- merge.row[!duplicated(merge.row$CDR3), ]
    merge.row$frequency <- (sum.CDR3$frequency*100)/sum(sum.CDR3$frequency)
    nrow <- which(table$subject.id == i)
    # remove the old rows from the dataframe then add the new one
    drop.row <- c(drop.row, nrow)
    table <- rbind(table, merge.row)
  }
  if(length(drop.row) == 0){
    result <- table
  } else {
    result <- table[-drop.row, ]
  }
  return(result)
}

## Identify IDD and SDD for each subject_id
identISD <- function(table){
  table[, "dominant"] <- 0L
  # iterate for all samples in each subject to find the max frequency
  for(i in levels(droplevels(table$subject.id))){
    # get the row indicies foreach subject ID
    id.row <- which(table$subject.id == i)
    # check the max frequency 
    max.freq <- max(table$frequency[id.row])
    # get the row index containing max frequency for each subject
    max.row <- which(table$frequency == max.freq)
    #select the row with max.freq defined of the certain object
    sel.row <- intersect(id.row, max.row)
    # assign IDD and SDD in the 'dominant column'
    table[sel.row,'dominant'] <- 'ID'
    table[id.row[!id.row %in% sel.row],'dominant'] <- 'SD'
  }
  table[,'dominant'] <- as.factor(table[,'dominant'])
  return(table)
}

# remove the subject id containing only one frequency value 
#(to reduce the bias since IDD is only defined in this subject)
remLowFreq <- function(table){
  t <- table %>% group_by(subject.id) %>% dplyr::summarise(n = n()) %>% filter(n == 1)
  table <- table[!table$subject.id %in% t$subject.id,]
  return(table)
}

# remove the subject that contain only IDD
remISD <- function(vdj.chain){
  library(dplyr)
  library(reshape2)
  library(data.table)
  library(tidyr)
  t <- vdj.chain %>% group_by(subject.id, dominant) %>% dplyr::summarise(n = n())
  t <- complete(t, subject.id, dominant, fill = list(n = 0))
  t2 <- dcast(t, subject.id~dominant, value.var = "n", fun.aggregate = sum, drop = FALSE)
  # remove the row containing either IDD zero
  t2.idd <- setDT(t2)[, .SD[!any(.SD[, -1, with = F] == 0)], by = subject.id]
  #Subset the subject after removing
  vdj <- vdj.chain[vdj.chain$subject.id %in% t2.idd$subject.id, ]
  return(vdj)
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
  name <- paste(parm, "new", sep = ".")
  table[, name] <- factor(result)
  return(table)
}

################################################################################
## Part IV: Write CDR3 with custom header as .fasta file ##
################################################################################

## write the table as in fasta format I
exportTable <- function(t,filename){
  export <- paste(paste(paste(paste(">",as.character(row.names(t)),sep=""),
                              t$dominant, sep = ""),t$subject.id,sep="") ,t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, 
              quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

## write the table in fasta format II (clearer headder than format I)
exportTable2 <- function(t,filename){
  export <- paste(paste(paste(">",t$dominant, sep = ""),1:length(t$dominant),
                        sep=""), t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, 
              quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

########################s########################################################
## Execution the code above ##
################################################################################

## pass arguments from the command lines
#*args = commandArgs(trailingOnly = TRUE) 
#*args[1] == VDJ data, args[2] == species, args[3] == MHC class

## Filter rawdata by species and MHC ##
raw.table <- freqFilter('D:/InternshipUU_2017/VDJdb/R_VDJdb/vdj241017.txt') #'vdj220817.txt #args[1]'

## reduce the resolution of V and J loci ##
#group V, J and MHC region 
table.j <- combineRegion(raw.table, "J")
table.v <- combineRegion(table.j, "V")
table <- combineRegion(table.v, "MHC.A")
rm(list = ls(pattern = "^table."))

## extracct VJ position ##
table <- parseVJ(table)

## check the length of CDR3 motif and epitopes ##
table <- aminoLength(table)

## extract cell subset from the VDJ db ##
# cell.subset
table$cell.subset <- as.factor(sapply(table$Meta, parseSubset))
# remove the Naive T-cell population
table <- table[!table$cell.subset == "CD45RA+CCR7+CD8+", ]
# add the column to indicate type of T-cell which indirectly infers from MHC class
table$T.cell <- as.factor(sapply(table$MHC.class, function(x) ifelse(x == "MHCI", "CD8+", "CD4+")))

## extract frequency ##
# extract frequency from the Meta column
table$raw.frequency <- parseFreq(table)
# format frequency to numeric value
frequency <- sapply(table$raw.frequency, freqFormat)
# check the format of frequency if it is in x/x, x% otherwise multiple 100 to the fraction
rowname <- unique(grep(paste(c("/", "%", "//"),collapse="|"), table$Method))
frequency[-rowname] <- frequency[-rowname]*100
# add the frequency column in the table
table$frequency <- frequency

## remove the row (subject) with the frequency = 100% (so, only dominance is defined) ##
table <- table[table$frequency != 100, ]

## select only CDR3 non-specifc to the human epitopes ##
table <- table[table$Epitope.species != "HomoSapiens",]

## add the cloumns of non-VJ sequences and their lengths ##
table$nonVJ <- substring(table$CDR3, table$v.end+1, table$j.start-1)
table$nonVJ.length <- sapply(table$nonVJ, nchar)

## remove some study ##
# remove PMID:15849183 due to different experimental setup; modified epitope #
table <- table[which(table$Reference != "PMID:15849183"),]
# remove PMID:12555663 due to (T-cell clone frequency from 4 donors not each subject)
table <- table[which(table$Reference != "PMID:12555663"),]
# reomove repetitive submission PMID:19017975 
table <- table[!rownames(table) %in% c(sapply(seq(7302,7310),as.character)), ]

################################################################################
## separate data based on the MHC class ##
vdj2 <- table[table$MHC.class == "MHCII",]
vdj1 <- table[table$MHC.class == "MHCI",]

# generae 4 datasets based on TCR chain  and MHC. calss then, identify the IDD and SDD
vdj.2a <- parseAB(vdj2,"TRA")
vdj.2b <- parseAB(vdj2,"TRB")
vdj.1a <- parseAB(vdj1,"TRA")
vdj.1b <- parseAB(vdj1,"TRB")

################################################################################
## add subject.id ##
# add column of raw.subject.id containing only the subject.id from the database
getID <- function(table){
  raw.subject.id <- sapply(table$Meta, function(x) sub(".*?subject.id: (.*?),.*","\\1", x))
  replica.id <- sapply(table$Meta, function(x) sub(".*?replica.id: (.*?),.*", "\\1",x))
  other.id <- paste(table$Epitope, table$Reference, sep = "*")
  
  #combine all the information and use as the subject id 
  subject.id <- paste(paste(raw.subject.id, replica.id, sep = "*"), other.id,sep = "*")
  return(subject.id)
}

## call function ##
addID <- function(vdj.chain){
  vdj.chain$subject.id <- getID(vdj.chain)
  vdj.chain <- getDenom(vdj.chain)
  return(vdj.chain)
}

vdj.2a <- addID(vdj.2a)
vdj.2b <- addID(vdj.2b)
vdj.1a <- addID(vdj.1a)
vdj.1b <- addID(vdj.1b)

################################################################################
## Immunodominance/subdominance assignment ##
################################################################################

assignDominance <- function(vdj.chain){
  # compile the total frequency more than 100% 
  vdjm <- mergeFreq(vdj.chain)
  # remove the smaple with only IDD defined
  vdjmr <- remLowFreq(vdjm)
  # assign the type of immune response(ID/SD)
  result <- remISD(identISD(vdjmr))
  return(result)
}

## call function ##
vdj.2a <- assignDominance(vdj.2a)
vdj.2b <- assignDominance(vdj.2b)
vdj.1a <- assignDominance(vdj.1a)
vdj.1b <- assignDominance(vdj.1b)

## final check at the total T-cell clone freq for each subject ##
freq.2a <- aggregate(vdj.2a$frequency, by = list(subject.id = vdj.2a$subject.id), sum)
freq.2b <- aggregate(vdj.2b$frequency, by = list(subject.id = vdj.2b$subject.id), sum)
freq.1a <- aggregate(vdj.1a$frequency, by = list(subject.id = vdj.1a$subject.id), sum)
freq.1b <- aggregate(vdj.1b$frequency, by = list(subject.id = vdj.1b$subject.id), sum)

################################################################################
## get a vector of non-redundnat CDR3 sequences
################################################################################
# generate the datasets with non-repetitive sequences
## write table 
#exportTable(spc_table, args[4])

## delete all objects except the function 
# rm(list = setdiff(ls(), lsf.str()))

