################################################################################
#
# Author: Wannisa Ritmahan
# Aim: MSA alignment based on IMGT position
# 
################################################################################
# map IMGT position to CDR3 sequences
################################################################################
# initialise amino acids and gap
#amino.acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T",
 #               "Y","C","W","R","K","D","E", "-")

defineGap <- function(cdr3.length){
  #cat("Only valid for the CDR3 length from 7 to 14 amino-acid long\n")
  # check if the length is even or odd number
  if(cdr3.length %% 2 == 0){
    i = 14
    exc = c(8)
    for(i in c(14, 12, 10, 8)){
      inc = (14-i)/2
      #print(inc)
      exc = c(exc, 8+inc, 8-inc)
      if(i == cdr3.length){
        result <- sort(unique(exc))
        return(result)
      } else {
        i = i+2
      }
    }
  } else if (cdr3.length %% 2 == 1){
    i = 13
    exc = c(8, 9)
    for(i in c(13, 11, 9, 7)){
      inc = (13-i)/2
      #print(inc)
      exc = c(exc, 8+inc, 8-inc, 9+inc, 9-inc)
      if(i == cdr3.length){
        result <- sort(unique(exc))
        return(result)
      } else {
        i = i+2
      }
    }
  }
}

countAminoImgt <-function(cdr3){
  aa <- c()
  # check the CDR3 length
  cdr3.length <- nchar(cdr3)
  # specified the postion to be filled with zero according to the length 
  zero.ind <- defineGap(cdr3.length)
  # initialise a list ot count amino acid at certain position
  for(i in 1:15){
    if(i < suppressWarnings(min(zero.ind))){
      aa <- c(aa, as.list(sapply(cdr3, substr, i, i)))
    } else if(i %in% zero.ind){
      aa <- c(aa, NA)
    } else {
      inc <- length(zero.ind)
      aa <- c(aa, as.list(sapply(cdr3, substr, i-inc, i-inc)))
    }
  }
  result <- as.list(aa)
  return(result)
}

locateImgtPos1 <- function(cdr3.list){
  count <- sapply(cdr3.list, countAminoImgt, simplify = FALSE)
  
  count.mat <- do.call(rbind, count)
  #count aa for each position of peptide #
  count.all <- apply(count.mat, 2, 
                     function(aa) sapply(amino.acid, function(x) x<-sum(x==unlist(aa,""), na.rm = T)))
  # calculate the percentage and round to 2 decimal point #
  count.perc <- round(prop.table(count.all, margin = 2)*100,1)
  # replace NA with zero
  count.perc[is.na(count.perc)] <- 0
  # assgin the colume names (== position of amino acids)
  colnames(count.perc) <- as.character(seq(dim(count.perc)[2]))
  return(count.perc)
}

################################################################################
# call function
t <- locateImgtPos1(as.character(cd8b.nr[cd8b.nr$CDR3.length < 15, "CDR3"]))

################################################################################
## Fill gap of CDR3 sequence to produce MSA following imgt position ##
################################################################################

addGap <- function(cdr3){
  if(nchar(cdr3) > 15){
    # # gap position of CDR3 length 16-33
    gp <- unique(c(seq(0,8), seq(0,-8)))
    gp <- gp[order(abs(gp))] + 17

    # get the index to be added gap ("-")
    l <- nchar(cdr3)
    gap <- 33 - l
    gap.ind <- sort(gp[1:gap])
    
    # add gap(s) to CDR3 of length 16-33
    gap.cdr3 <- paste(c(substr(cdr3,1,min(gap.ind)-1), rep("-",length(gap.ind)), 
                        substr(cdr3, min(gap.ind), 33)), collapse = "")
  } else if(nchar(cdr3) == 15){
    # add gap(s) to CDR3 of length == 15
    gap.cdr3 <- paste(c(substr(cdr3,1,8), rep("-",18), 
                        substr(cdr3,9,15)), collapse = "")
    # add gap(s) to CDR3 of length < 15
  } else {
    gap.ind <- defineGap(nchar(cdr3))
    gap.cdr3 <- paste(c(substr(cdr3,1,min(gap.ind)-1), rep("-",length(gap.ind)+18), 
                        substr(cdr3, min(gap.ind), nchar(cdr3))), collapse = "")
  }
  # gap.cdr3 will have the length 33 regardless initial cdr3 length
  return(gap.cdr3)
}

# delete middle gap
moveGap <- function(gap.cdr3){
  # assign each seq in different cell of matrix row
  m <- matrix(unlist(sapply(gap.cdr3, strsplit, "")), ncol = 33, byrow = T)
  # check column containing only "-" (gap)
  gap.col <- colSums(m == "-")
  names(gap.col) <- 1:33
  rem.col <- as.numeric(names(gap.col[gap.col == dim(m)[1]]))
  # remove the column cotaining only "-"
  m <- m[, -rem.col]
  # concatenate each row of the matrix
  result <- apply(m, 1, paste, collapse = "")
  return(result)
}
#c(5:8,13:15)
################################################################################
# export seqwuence in fasta format #
################################################################################
# export gap.cdr3 to be used as an input for sequence logo
exportFastaGap <- function(t, seq, name){
  export <- paste(paste(paste(">",seq,sep=""),
                                    t$subject.id,sep=""),seq, sep="\n")
  
  write.table(export, paste0(name,".fasta"), append = FALSE, quote = FALSE, 
              sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
}
################################################################################
alignGap <- function(d, response){
  d.subset <- sapply(d[d$dominant == response, "CDR3"], as.character)
  gap.cdr3 <- sapply(d.subset,addGap)
  
  # remove centrally located gap
  gap.cdr3.mv <- moveGap(gap.cdr3)
  return(gap.cdr3.mv)
}

alignGap1 <- function(d){
  d.subset <- as.character(d)
  gap.cdr3 <- sapply(d.subset,addGap)
  
  # remove centrally located gap
  gap.cdr3.mv <- moveGap(gap.cdr3)
  return(gap.cdr3.mv)
}

################################################################################
# export gap.cdr3 of ID and SD in fasta format
writeFastaGap <- function(d, response, filename){
  exp.seq <- alignGap(d, response)
  # call exportFastaGap
  exportFastaGap(d[d$dominant == response, ], exp.seq, paste0(filename, response))
}
################################################################################
# call function
# *** first check the summary(factor(cd8b.nr$CDR3.length)) then, decide min and max 
# length from 9 and 20- 20 is the max length that exist in both ID and SD responses
cd8b.sellen <- droplevels(cd8b.nr[cd8b.nr$CDR3.length > 9 & cd8b.nr$CDR3.length < 20,])
cd8b.hlaa2 <- droplevels(cd8b.sellen[cd8b.sellen$MHC.A.new == "HLA-A*02",])
cd8b.hlaa2.id <- cd8b.hlaa2[cd8b.hlaa2$dominant == "ID", ]
cd8b.hlaa2.sd <- cd8b.hlaa2[cd8b.hlaa2$dominant == "SD", ]

cd8b.hlaa1 <- droplevels(cd8b.sellen[cd8b.sellen$MHC.A.new == "HLA-A*01",])
cd8b.hlaa1.id <- cd8b.hlaa1[cd8b.hlaa1$dominant == "ID", ]
cd8b.hlaa1.sd <- cd8b.hlaa1[cd8b.hlaa1$dominant == "SD", ]

# generate aligned cdr3 based on imgt position with gap
cd8bhlaa2ID.msa <- alignGap(cd8b0.01nr, "ID")
cd8bhlaa2SD.msa <- alignGap(cd8b0.01nr, "SD")

cd8bhlaa1ID.msa <- alignGap(cd8b.hlaa1, "ID")
cd8bhlaa1SD.msa <- alignGap(cd8b.hlaa1, "SD")


cd8bhlaa11ID.msa <- alignGap(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*11",], "ID")
cd8bhlaa11SD.msa <- alignGap(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*11",], "SD")


################################################################################
# export CDR3 of hlaa2

writeFastaGap(cd8b.hlaa2, "ID", "cdr3hlaa2_vdj04jan18_5.13")
writeFastaGap(cd8b.hlaa2, "SD", "cdr3hlaa2_vdj04jan18_5.13")

write.table(c(cd8bhlaa2ID.msa, cd8bhlaa2SD.msa), file = "cd8ba2_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)

# export epitope seq
write.table(levels(cd8b.msa.hlaa2$Epitope), file = "epihlaa2_vdj04jan18.fasta", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)


write.table(cd8bhlaa2ID.msa, file = "cdr3a2ID_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(cd8bhlaa2SD.msa, file = "cdr3a2SD_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)

# export CDR3 of hlaa1

writeFastaGap(cd8b.hlaa1, "ID", "cdr3hlaa1_vdj04jan18")
writeFastaGap(cd8b.hlaa1, "SD", "cdr3hlaa1_vdj04jan18")

write.table(cd8bhlaa1ID.msa, file = "cdr3a1ID_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(cd8bhlaa1SD.msa, file = "cdr3a1SD_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)

################################################################################
# epitope subset
a2GLCID.msa <- alignGap(a2GLC, "ID")
a2GLCSD.msa <- alignGap(a2GLC, "SD")

a2GILID.msa <- alignGap(a2GIL, "ID")
a2GILSD.msa <- alignGap(a2GIL, "SD")

a2NLVID.msa <- alignGap(a2NLV, "ID")
a2NLVSD.msa <- alignGap(a2NLV, "SD")

# export sequence alignment
writeFastaGap(a2GLC, "ID", "a2GLCID_vdj04jan18")
writeFastaGap(a2GLC, "SD", "a2GLCSD_vdj04jan18")

write.table(a2GLCID.msa, file = "a2GLCID.msa_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(a2GLCSD.msa, file = "a2GLCSD.msa_vdj04jan18.txt", append = FALSE, 
            quote = FALSE, sep = "\n",eol = "\n", row.names = FALSE, col.names = FALSE)


################################################################################

cd8bhlaa2ID.msa <- alignGap(cd8b0.01nr[cd8b0.01nr$CDR3.length %in% 10:15, ], "ID")
cd8bhlaa2SD.msa <- alignGap(cd8b0.01nr[cd8b0.01nr$CDR3.length %in% 10:15, ], "SD")

(positionalFrequency(cd8bhlaa2ID.msa, F)/positionalFrequency(cd8bhlaa2SD.msa, F))

