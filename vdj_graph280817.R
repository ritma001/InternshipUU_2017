#!usr/bin/r
# Author: Wannisa Ritmahan
# Established Date: 28/08/2017 Last Modification: 03/09/2017
# Aim: To study common feature of TCR motif and epitopes from VDJdb

################################################################################
#layout(matrix(1:2, nrow=1))
library(ggplot2)
library(plyr)

## calculate prop  of certain CDR3 length
vdj.new <-ddply(vdj,.(dominant),summarise,
                prop = prop.table(table(dominant)),
                Frequency = as.integer(table(CDR3.length)),
                CDR3.length = names(table(CDR3.length))) # fix the code here 


ggplot(vdj.new, aes(MHC.class, Frequency,fill=dominant)) +
  geom_bar(stat = "identity", position = "dodge") + labs(title = 'MHC.class', xlab="MHC.class")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(vdj.new, aes(MHC.class,  Frequency,fill=dominant)) +
  geom_bar(stat = "identity", position = "dodge") + labs(title = 'MHC.class', xlab="MHC.class")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  facet_grid(~dominant)


################################################################################
# pie chart of Epitope species presenting in dataset

epi <- as.data.frame(summary(droplevels(vdj$Epitope.species)))
pie.epi <- epi[epi != 0, ,FALSE]
pie.labels <- c(paste(row.names(epi),epi[,1],sep='\n'))
pie(epi[,1],col = rainbow(length(epi[,1]),0.3), labels = pie.labels, main = 'Epitope species 
    of Human\' s CDR3 motifs from VDJdb database' )

################################################################################
# analysis VJ combination 
 
# barplot of VJ for IDD and SDD for alpha and beta chain
plotAB <- function(t, chain, main){
  df <- sort(summary(droplevels(t[, chain]))/dim(t)[1], decreasing = TRUE)
  bp <- barplot(df[df>0.01],xaxt = "n", ylab = "Frequency", main = main)
      #,col = rainbow(length(droplevels(t[, chain])))[t$V])
  text(bp, par("usr")[3], labels = names(df[df>0.01]), srt = 45, 
       adj = c(1.1,1.1),xpd= TRUE, cex = 0.8)
}

# call the function
plotAB(idd.alpha, "V", "V segment of alpha chain in IDD data")
plotAB(sdd.alpha, "V", "V segment of alpha chain in SDD data")
plotAB(idd.beta, "V", "V segment of beta chain in IDD data")
plotAB(sdd.beta, "V", "V segment of beta chain in SDD data")


plotAB(idd.alpha, "J", "J segment of alpha chain in IDD data")
plotAB(sdd.alpha, "J", "J segment of alpha chain in SDD data")
plotAB(idd.beta, "J", "J segment of beta chain in IDD data")
plotAB(sdd.beta, "J", "J segment of beta chain in SDD data")

################################################################################

# sort the factor by level frequency
sortFactor <- function(t, chain){
  t <- sort(summary(droplevels(t[, chain])),decreasing = TRUE)
  return(t)
}

################################################################################
## amino acid distribution: the input is the vector extract from the VDJ table ##
## input: cdr <- vdj$CDR3

## count amino acids (1 mer) of cdr3 sequencing: input =  a vector of CDR3 sequencing 
countAmino <- function(cdr){
  t <- data.frame("amino_acid" = c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E"),
                  "properties" = c(rep("hydrophobic",8),rep("hydrophilic",8),rep("charged",4)), "count" = rep(0, nrow(t)))
  #chemical properties from http://www.proteinstructures.com/Structure/Structure/amino-acids.html
  for(a in cdr){
    count <- table(strsplit(toupper(as.character(a)), ''))
    for(e in names(count)){
      t[t$amino_acid == e,3] <- t[t$amino_acid == e,3] + as.vector(count[e])
    }
  }
  #calculate relative frequency and add to the last column
  relative_freq <- t$count/sum(t$count)
  t$frequency <- relative_freq
  return(t)
}

## call function
cdr.sdd.a <- countAmino(sdd.alpha$CDR3)
cdr.idd.a <- countAmino(idd.alpha$CDR3)
cdr.idd.b <- countAmino(idd.beta$CDR3)
cdr.sdd.b <- countAmino(sdd.beta$CDR3)

################################################################################
## count 2-mer amino acid 

count2mer <- function(cdr,dominant){
  # initialise the 20*20 matrix to store frequency
  amino_acid<- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
  m <- matrix(rep(0,400), nrow = 20, ncol=20, byrow = TRUE)
  dimnames(m) <- list(amino_acid,"amino_acid" = amino_acid)
  # count the 2 mer of amino acid
  for(a in cdr){
    i = 1
    while(i < nchar(a)-1){
      mer <- substr(a,i,i+1)
      split.mer <- unlist(strsplit(mer,split = ""))
      m[split.mer[1],split.mer[2]] <- m[split.mer[1],split.mer[2]] + 1
      i = i+1 # i = i+2 if the non-ovelapping sequence is analysed...which one to be chosen?
    }
  }
  # convert matrix to dataframe
  df <- as.data.frame(as.table(m))
  # remove the row containing zero
  row.remove <- which(df[,ncol(df)] == 0)
  result <- df[-row.remove,]
  # order the value based on the frequency values
  result <- result[order(-xtfrm(result[,ncol(result)])),]
  # combine mer keep the 2-mer and frequency column only
  result$two.mer <- factor(apply(result[,-ncol(result)],1,paste,collapse = ""))
  #result[,1:2] <- NULL if the first two column is designed to be deleted after being combined
  # calculate relative frequency
  prop <- result$Freq/nrow(result)
  result$prop <- prop
  return(result)
}

## call function
cdr2.idd.a <- count2mer(idd.alpha$CDR3)
cdr2.sdd.a <- count2mer(sdd.alpha$CDR3)
cdr2.idd.b <- count2mer(idd.beta$CDR3)
cdr2.sdd.b <- count2mer(sdd.beta$CDR3)

################################################################################
## count 3-mer amino-acid

count3mer <- function(cdr){
  library(plyr)
 amino_acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
 ar <- array(rep(0,8000), dim = c(20,20,20), dimnames = list(amino_acid,amino_acid,amino_acid))
 for(a in cdr){
   i = 1
   while(i < nchar(a)-2){
     mer <- substr(a,i,i+2)
     split.mer <- unlist(strsplit(mer,split = ""))
     ar[split.mer[1],split.mer[2],split.mer[3]] <- ar[split.mer[1],split.mer[2],split.mer[3]] + 1
     i = i+1 # i = i+2 if the non-ovelapping sequence is analysed
   }
 }
 
 # convert array to dataframe
 df <- adply(ar,c(1,2,3))
 row.remove <- which(df[,ncol(df)] == 0)
 result <- df[-row.remove,]
 result <- result[order(-xtfrm(result[,ncol(result)])),]
 result$three.mer <- factor(apply(result[,-ncol(result)],1,paste,collapse = ""))
 prop <- result$V1/nrow(result)
 result$prop <- prop
 return(result)
}

# call function
cdr3.idd.a <- count3mer(idd.alpha$CDR3)
cdr3.sdd.a <- count3mer(sdd.alpha$CDR3)
cdr3.idd.b <- count3mer(idd.beta$CDR3)
cdr3.sdd.b <- count3mer(sdd.beta$CDR3)

# ar3 <- count3mer(sdd.beta$CDR3)
# (ar3[which(!ar3==0)]) ; check non-zero values
# df3 <- ; 
# df3 <- df3[df3[, ncol] != 0] ; remove the 3-mer that are not represented

################################################################################
## count 4-mer 

count4mer <- function(cdr){
  library(plyr)
  amino_acid <- c("A","I","L","M","F","V","P","G","Q","N","H","S","T","Y","C","W","R","K","D","E")
  ar <- array(rep(0,160000), dim = c(20,20,20,20), dimnames = list(amino_acid,amino_acid,amino_acid, amino_acid))
  for(a in cdr){
    i = 1
    while(i < nchar(a)-3){
      mer <- substr(a,i,i+3)
      split.mer <- unlist(strsplit(mer,split = ""))
      ar[split.mer[1],split.mer[2],split.mer[3],split.mer[4]] <- ar[split.mer[1],split.mer[2],split.mer[3],split.mer[4]] + 1
      i = i+1 # i = i+2 if the non-ovelapping sequence is analysed...which one to be chosen?
    }
  }
  # convert array to dataframe
  df <- adply(ar,c(1,2,3,4))
  row.remove <- which(df[,ncol(df)] == 0)
  result <- df[-row.remove,]
  result <- result[order(-xtfrm(result[,ncol(result)])),]
  result$four.mer <- factor(apply(result[,-ncol(result)],1,paste,collapse = ""))
  prop <- result$V1/nrow(result)
  result$prop <- prop
  return(result)
}

# call function
cdr4.idd.a <- count4mer(idd.alpha$CDR3)
cdr4.sdd.a <- count4mer(sdd.alpha$CDR3)
cdr4.idd.b <- count4mer(idd.beta$CDR3)
cdr4.sdd.b <- count4mer(sdd.beta$CDR3)

################################################################################
## plot n-mer
plot1mer <- function(idd,sdd){
  cdr1 <- rbind(idd,sdd)
  cdr1$dominant <- factor(c(rep("IDD",20),rep("SDD",20)), levels = c("IDD","SDD"))
  library(ggplot2)
  cdr1.p <- ggplot(cdr1, aes(x=amino_acid,y=frequency,fill=properties, order=-as.numeric(properties))) + geom_bar(stat = "identity")
  cdr1.p + facet_wrap(~dominant) + ggtitle("Amino acid distribution in IDD and SDD alpha chain") + xlab("1-mer amino acid") + ylab("Relative Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  }

plot3mer <- function(idd,sdd){
  cdr3 <- rbind(idd[2:21,],sdd[2:21,])
  cdr3$dominant <- factor(c(rep("IDD",20),rep("SDD",20)), levels = c("IDD","SDD"))
  #plot
  library(ggplot2)
  cdr3.p <- ggplot(cdr3, aes(x=three.mer,y=prop, fill = dominant)) + geom_bar(stat = "identity", position = "dodge" ,width = 0.5)
  cdr3.p <-  cdr3.p + ggtitle("CDR3 Alpha Chain") + xlab("3-mer amino acid motif") + ylab("Relative Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  return(cdr3.p)
}

plot4mer <- function(idd,sdd){
  cdr4 <- rbind(idd[2:21,],sdd[2:21,])
  cdr4$dominant <- factor(c(rep("IDD",20),rep("SDD",20)), levels = c("IDD","SDD"))
  #plot
  library(ggplot2)
  cdr4.p <- ggplot(cdr4, aes(x=four.mer,y=prop, fill = dominant)) + geom_bar(stat = "identity", position = "dodge" ,width = 0.5)
  cdr4.p <-  cdr4.p + ggtitle("CDR3 Beta Chain") + xlab("4-mer amino acid motif") + ylab("Relative Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  return(cdr4.p)
}

################################################################################
# create dataset

vdj.a <- vdj[vdj$Gene == "TRA",]
vdj.b <- vdj[vdj$Gene == "TRB",]

idd.a <- vdj.a[vdj.a$dominant=='IDD',]
idd.b <- vdj.b[vdj.b$dominant=='IDD',]
sdd.a <- vdj.a[vdj.a$dominant=='SDD',]
sdd.b <- vdj.b[vdj.b$dominant=='SDD',]
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
    ylab('Relative Frequency') + xlab("CDR3 length") + ggtitle(paste0(paste0("The distribution of V seqments for each CDR3 length in", chain), dominant))+
    theme(plot.title = element_text(hjust = 0.5))
  p + facet_wrap(~dominant, scales = "free_y" ,drop = TRUE)+ scale_fill_manual(values = getPalette(colourCount))
}

# call function
plotV(vdj.a," alpha chain of ", "IDD and SDD")
plotV(vdj.b," beta chain of ", "IDD and SDD")

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

plotJ(vdj.a," alpha chain of ", "IDD and SDD")
plotJ(vdj.b," beta chain of ", "IDD and SDD")


#if no facet_wrap and p+ scale_fill_manual(...)

plotJ(sdd.b, " beta chain of ","SDD")
plotJ(sdd.a, " alpha chain of ","SDD")
plotJ(idd.a, " alpha chain of ","IDD")
plotJ(idd.b, " beta chain of ","IDD")

################################################################################


