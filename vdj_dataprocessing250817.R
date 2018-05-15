# Author: Wannisa Ritmahan 
# Created on: 28/08/2017
# Aim: R code for manipulating data from the VDJdb for Immunodominance studies
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
    clone <- sub(".*?clone.id: (.*?),.*", "\\1", Meta)
    replica <- sub(".*?replica.id: (.*?),.*", "\\1", Meta)
    epitope <- sub(".*?epitope.id: (.*?),.*", "\\1", Meta)
    cohort <- sub(".*?subject.cohort: (.*?),.*", "\\1", Meta)
    
    id = paste(paste(paste(clone, replica, sep = "*"),epitope, sep = "*"), cohort, sep = "*")
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
  table <- table[table$j.start != -1, ]
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
  lenCDR3 <- nchar(as.character(table$CDR3))-2
  table[,'CDR3.length'] <- lenCDR3
  
  lenEpitope <- nchar(as.character(table$Epitope))
  table[,'Epitope.length'] <- lenEpitope
  return(table)
}

################################################################################
## Part II: Identification of immunodominant- & subdominant CDR3 ##
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
  if(length(merge.subject == 0)){
    return(table)
  } else {
    for(i in merge.subject){
      sum.CDR3 <- aggregate(frequency ~ CDR3, table[table$subject.id == i,], sum)
      # select only the rows containing the subject.id == i
      merge.row <- table[table$subject.id == i,]
      # remove repetitive CDR3
      merge.row <- merge.row[!duplicated(merge.row$CDR3), ]
      
      # combine freq and recalculate the total frequency
      sum.CDR3$corr.frequency <- (sum.CDR3$frequency*100)/sum(sum.CDR3$frequency)
      # correct the freq in merge.row  by replacing with correponding corr.freq value
      merge.row[match(sum.CDR3$CDR3, merge.row$CDR3), "frequency"] <- sum.CDR3$corr.frequency
      
      # remove the old row of subject i before adding since the merge row contains the same rowname 
      rem.row <- rownames(table[which(table$subject.id == i),])
      # remove the old rows from the dataframe then add the new one
      table <- table[!rownames(table) %in% rem.row, ]
      # add merge row
      table <- rbind(table, merge.row)
    }
    return(table)
  }
}

## Identify ID and SD for each subject.id
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
    # assign ID and SD in the 'dominant column'
    table[sel.row,'dominant'] <- 'ID'
    table[id.row[!id.row %in% sel.row],'dominant'] <- 'SD'
  }
  return(table)
}

# remove the subject id containing only one frequency value
remLowFreq <- function(table){
  t <- table %>% group_by(subject.id) %>% dplyr::summarise(n = n()) %>% filter(n == 1)
  table <- table[!table$subject.id %in% t$subject.id,]
  return(table)
}

# remove the subject that contain only ID
remISD <- function(vdj.chain){
  library(dplyr)
  library(reshape2)
  library(data.table)
  library(tidyr)
  t <- vdj.chain %>% group_by(subject.id, dominant) %>% dplyr::summarise(n = n())
  t <- complete(t, subject.id, dominant, fill = list(n = 0))
  t2 <- dcast(t, subject.id~dominant, value.var = "n", fun.aggregate = sum, drop = FALSE)
  # remove the row containing either ID zero
  t2.id <- setDT(t2)[, .SD[!any(.SD[, -1, with = F] == 0)], by = subject.id]
  #Subset the subject after removing
  vdj <- vdj.chain[vdj.chain$subject.id %in% t2.id$subject.id, ]
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

## create the 'ID' and 'SD ' datasets with option to exclude human sample
parseISD <- function(table, epitope){
  df <- table[table$dominant == epitope & table$Epitope.species != 'HomoSapiens',]
  return(df)
}

## create the different dataset using gene (TRA-alpha or TRB-beta chain) as criteria
parseAB <- function(table, chain){
  df <- table[table$Gene == chain,]
  return(df)
}


## remove the suffix of of V, J and MHC.A to enable comparison 
## at all studies at higher resolution
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
## Part IV: Write CDR3 with custom header as .fasta file if required ##
################################################################################

## write the table as in fasta format I
exportFasta<- function(t,filename){
  export <- paste(paste(paste(paste(">",as.character(row.names(t)),sep=""),
                              t$dominant, sep = ""),t$subject.id,sep="") ,t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, 
              quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

## write the table in fasta format II (clearer headder than format I)
exportFasta2 <- function(t,filename){
  export <- paste(paste(paste(">",t$dominant, sep = ""),1:length(t$dominant),
                        sep=""), t$CDR3, sep="\n")
  write.table(export, file = paste(filename,"fasta",sep="."), append = FALSE, 
              quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

########################s########################################################
## Part V: combined functions previously created to generate a flow of data 
## processing and for an ease for code execution
################################################################################
# extrat prefix of V,J (TR[A,B][V,J]xx) and MHC (HL(A,B)-*xx)
reduceSolution <- function(raw.table){
  table.j <- combineRegion(raw.table, "J")
  table.v <- combineRegion(table.j, "V")
  table <- combineRegion(table.v, "MHC.A")
  
  return(table)
}

extractInfo <- function(table){
  ## extract VJ position ##
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
  
  ## add the cloumns of non-VJ sequences and their lengths ##
  table$nonVJ <- substring(table$CDR3, table$v.end+1, table$j.start)
  table$nonVJ.length <- sapply(table$nonVJ, nchar)
  
  table$Vseq <- substring(table$CDR3, 1, table$v.end)
  table$Jseq <- substring(table$CDR3, table$j.start, table$CDR3.length)

  return(table)
}

removeData <- function(table){
  ## select only CDR3 non-specifc to the human epitopes ##
  table <- table[table$Epitope.species != "HomoSapiens",]
  
  ## remove some study ##
  # remove PMID:15849183 due to different experimental setup; modified epitope #
  table <- table[which(table$Reference != "PMID:15849183"),]
  # remove PMID:28146579 due to the memory T-cell study
  table <- table[which(table$Reference != "PMID:28146579"),]
  # remove PMID:12555663 due to (T-cell clone frequency from 4 donors not each subject)
  table <- table[which(table$Reference != "PMID:12555663"),]
  # reomove repetitive submit PMID:19017975 
  ##table <- table[!rownames(table) %in% c(sapply(seq(808,816),as.character)), ]
  
  return(table)
}

calculateFrequency <- function(table){
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
  
  return(table)
}

parseTCR <- function(table, mhc, chain){
  ## separate data based on the MHC class and TCR chain ##
  table.mhc <- table[table$MHC.class == mhc,]
  table.chain <- table.mhc[table.mhc$Gene == chain,]
  
  return(table.chain)
}

## recalculate frequency for selecting rare clonotype < 0.01% ##
recalculateFreq <- function(df){
  invisible(lapply(c("tidyr", "dplyr"), library, character.only = TRUE))
  # sort dataframe by subject.id column 
  df.sort <- df[order(df$subject.id),]
  # recalculate the frequency and add this column (recalFreq) to the df.sort
  df <- df.sort %>% 
    group_by(subject.id) %>% 
    mutate(recalFreq = frequency*100/sum(frequency))
  return(df)
}

## add subject.id ##
################################################################################
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

## add the frequency of multiple clones in each individual
remMultiClone <- function(df){
  library(dplyr)
  # check if more then one CDR3 in individual then add the frequency of that CDR3
  df.new <- df %>%
    group_by(subject.id, CDR3, V.new, J.new) %>%
    mutate(frequency=paste(sum(frequency),collapse='')) %>% ungroup()
  
  # collapse calumn in individual 
  result <- df.new[!duplicated(df.new[,c("subject.id", "CDR3", "V.new", "J.new")]),]
  result$frequency <- as.numeric(result$frequency)
  return(result)
}

## Immunodominance/subdominance assignment ##
################################################################################
assignDominance <- function(vdj.chain){
  # compile the total frequency more than 100% 
  vdjm <- mergeFreq(vdj.chain)
  # remove the subjec that contains only 1 sequence
  vdj.mr <- remLowFreq(vdjm)
  # assign the type of immune response(ID/SD)
  vdj.ISD <- identISD(vdj.mr)
  # remove the subject containing only ID
  result <- remISD(vdj.ISD)
  # add corrected frequency
  result <- recalculateFreq(result)
  return(result)
}

## assign ID/SD
addDominance <- function(table){
  # calculate freq #
  table.freq <- calculateFrequency(table)
  # add subject.id and denom
  table.id <- addID(table.freq)
  # add frequency of the same CDR3 seq in each individual
  table.rem <- remMultiClone(table.id)
  # assign the immune response type
  table.ID <- assignDominance(table.rem)
  
  return(table.ID)
}

## create a dataset containing only the unique CDR3 sequences #
excludeRepeat <- function(vdj.chain){
  cdr3.nr <- vdj.chain[!duplicated(vdj.chain[,c("CDR3", "dominant", "Epitope")]),]
  return(cdr3.nr)
}

# add peptide properties calculated from Peptides and alakazam packages
################################################################################
#* Kidera factor *#
# calculate 10 Kidera factors, a physico-chemical properties related to structure
callLibrary <- function(lib){
  if(!require(noqoute(lib)))
  {installed.packages(lib)}
  library(lib)
}

addProperties <- function(vdj){
  library("alakazam")
  library("Peptides")
  
  # add 10 columns contataining 10 Kidera factors (KF)
  kf <- sapply(vdj$CDR3, kideraFactors)
  kf.df <- as.data.frame(t(matrix(unlist(kf), nrow = length(kf[[1]]))))
  colnames(kf.df) <- paste0("KF", 1:10)
  vdj.kf <- cbind(vdj, kf.df)
  
  # add 8 properties calculated from amino-acid contaning CDR3 seq
  vdj.prop <- aminoAcidProperties(vdj.kf, property = c("gravy", "bulk", "aliphatic",
                       "polarity", "charge", "basic", "acidic", "aromatic"),"CDR3")
  # add grand hydropathy (gravy) of epitope
  vdj.prop.epi <- aminoAcidProperties(vdj.prop, property = c("gravy"), "Epitope")
  return(vdj.prop.epi)
}

## final check at the total T-cell clone freq for each subject ##
#checkFreq <- aggregate(cd8a$frequency, by = list(subject.id = cd8a$subject.id), sum)

################################################################################
## process raw data of VDJ db ##
################################################################################

## import tab-delimited table downloaded from VDJ db web server
raw.table <- freqFilter('D:/InternshipUU_2017/VDJdb/R_VDJdb/vdj190318.tsv')

# extract info and remove some studies
cleanData <- function(raw.vdj){
  library(tidyr)
  library(dplyr)
  vdj <- reduceSolution(raw.vdj) %>%
    extractInfo %>%
    removeData
  return(vdj)
}

processData <- function(vdj, mhc, chain){
  # divide data based on MHC (MHCI or MHCII) and TCR chain (TRA or TRB)
  vdj.chain <- parseTCR(vdj, mhc,chain) %>%
    addDominance %>%                 # defined ID/SD
    filter(recalFreq > 0.01) %>%     # exclude rare CDR3 with low frequency
    excludeRepeat %>%                # exclude repetitive sequences
    as.data.frame %>%                # convert grouped_df to data.frame
    addProperties                    # add CDR3 pro
  return(vdj.chain)
}

################################################################################
## process raw data
vdj <- cleanData(raw.table)

## generate datasets based on TCR chain  and MHC
cd8a.nr <- processData(vdj, "MHCI","TRA")
cd8b.nr <- processData(vdj, "MHCI","TRB")

################################################################################
# select the top three epitopes
a2GIL <- cd8b.nr[cd8b.nr$Epitope== "GLCTLVAML", ] 
a2GIL <- cd8b.nr[cd8b.nr$Epitope== "GILGFVFTL", ]  
a2NLV <- cd8b.nr[cd8b.nr$Epitope== "NLVPMVATV", ] 


################################################################################
## select public TCR ##
################################################################################
selectPublicTCR <- function(table, cutoff){
  require(tidyr)
  require(dplyr)
  t <- table %>% 
    group_by(CDR3) %>% 
    dplyr::summarise(n = n()) %>%
    filter(n < cutoff) 
  table <- table[!table$CDR3 %in% t$CDR3,]
  return(table)
}

# select sequence that presents > 1
pcd8b <- selectPublicTCR(cd8b, 2)
# change InfluenzaA to Influenza-A
pcd8b[pcd8b$Epitope.species == "InfluenzaA", "Epitope.species"] <- "Influenza-A"
# extract unique CDR3 based on CDR3 seq, dominant and epitope
pcd8b.nr <- excludeRepeat(pcd8b)
pcd8b.nr <- addProperties(pcd8b.nr)

################################################################################
## import naive T-cell data ##
################################################################################
selectCDR3 <- function(d, cell, chain, id){
  library(dplyr)
  library(tidyr)
  # select the CDR3 with only Naive count and zero Non-naive count (NN count) to
  # avoid technical error from cell contamination
  naive <- d[which(d[, 1] != 0),]
  naive <- naive[which(naive[, 2] == 0),]
  # add the CD3 length, nonV-J length, column
  naive <- naive %>% 
    mutate(CDR3.length = sapply(as.character(naive$AAseq), nchar)-2,
           cell.subset = rep(cell, dim(naive)[1]),
           chain = rep(chain, dim(naive)[1]),
           subject.id = rep(id, dim(naive)[1]),
           nonVJ.length = sapply(as.character(naive$Ins), nchar))
  return(naive)
}

# export CDR3 to be used as an input for python script to extract CDR3 AA
exportCDR3.naive <- function(table, filename, column){
  # the arg column can be either "CDR3" or "Epitope"
  cdr3.table <- paste(paste0(">", as.character(row.names(table))), 
                      table[,column], sep = "\n")
  write.table(cdr3.table , file = paste0(filename, ".txt"), append = FALSE, 
              quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}

################################################################################
# read csv file; 1 file/donor of CDR3 beta chain

pt_8b1 <- read.csv("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/classified8b1.csv")
pt_8b2 <- read.csv("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/classified8b2.csv")

naive.pt8b.1 <- selectCDR3(pt_8b1, "CD8+", "Beta", "1")
naive.pt8b.2 <- selectCDR3(pt_8b2, "CD8+", "Beta", "2")

# combined data from 2 donors and select the CDR3 clone with absolute count > 1
naive.pt8b <- rbind(naive.pt8b.1 , naive.pt8b.2)
ncd8b.nr <- naive.pt8b[!duplicated(naive.pt8b[, "AAseq"]) & naive.pt8b$Naive.count > 1,]

# add nonVJ AA
exportCDR3.naive(naive.pt8b, "peter_naivecd8b", "AAseq")
naive.hydropathy.cdr3 <- read.table("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/peter_naivecd8b_hd.txt")
naive.pt8b$hydropathy.cdr3 <- naive.hydropathy.cdr3[,1]

################################################################################
# read csv file; 1 file/donor of CDR3 alpah chain

pt_8a1 <- read.csv("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/classified8a1.csv")
pt_8a2 <- read.csv("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/classified8a2.csv")

naive.pt8a.1 <- selectCDR3(pt_8a1, "CD8+", "Alpha", "1")
naive.pt8a.2 <- selectCDR3(pt_8a2, "CD8+", "Alpha", "2")

naive.pt8a <- rbind(naive.pt8a.1 , naive.pt8a.2)
ncd8a.nr <- naive.pt8a[!duplicated(naive.pt8a[, "AAseq"]) & naive.pt8a$Naive.count > 1,]

# add nonVJ AA
exportCDR3.naive(naive.pt8a, "peter_naivecd8a", "AAseq")
naive.hydropathy.cdr3 <- read.table("D:/InternshipUU_2017/Naive_Tcell/Peter_NT/peter_naivecd8a_hd.txt")
naive.pt8a$hydropathy.cdr3 <- naive.hydropathy.cdr3[,1]

################################################################################
# generate the CDR-H3 dataset containing CDR3 from the heavy chain of B-cell
# (downloaded from abYsis-human heavy chain)
################################################################################

readCDRH3 <- function(file){
  cdrh3 <- read.table(file, col.names = "CDRH3")
  cdrh3$insertion <- apply(cdrh3, 1, function(seq) substr(seq, 7, nchar(seq)-2))
  return(cdrh3)
}

abysis1 <- readCDRH3("ABYSIS_CDRH3.txt")
abysis401 <- readCDRH3("ABYSIS401_CDRH3.txt")
abysis801 <- readCDRH3("ABYSIS801_CDRH3.txt")
abysis1201 <- readCDRH3("ABYSIS1201_CDRH3.txt")

bcr.cdrh3 <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                    list(abysis1, abysis401,abysis801, abysis1201))



