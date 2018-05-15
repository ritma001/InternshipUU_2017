# 6-23 Feb 2018
# personalised analysis
# select the ID with > 40% then study the  individual response
#
################################################################################
##cd8b.exc <- cd8b[cd8b$Reference != "PMID:28423320",]

# divide the data (cd8b.exc = cd8b excluded PMID:28423320)
a2GLC <- cd8b[cd8b$Epitope== "GLCTLVAML", ] 
a2GIL <- cd8b[cd8b$Epitope== "GILGFVFTL", ]  
a2NLV <- cd8b[cd8b$Epitope== "NLVPMVATV", ] 

# a function that returns the subject.id with high frequency of ID (log(ID/SD) >0.5)
filterLowFreq <- function(df, cutoff.freq){
  library(tidyr)
  library(dplyr)
  library(reshape)
  sum.df <- df %>% 
    group_by(subject.id, dominant) %>% 
    dplyr::summarise(count = max(frequency))
  
  sel.subject <- cast(sum.df, subject.id~dominant, value = "count") %>% 
    filter(log(ID/SD) > cutoff.freq) 
  return(as.character(sel.subject$subject.id))
}

################################################################################
## public clone ID and SD ##
################################################################################
identDominantClone <- function(df, n){
  library(tidyr)
  library(dplyr)
  library(reshape)
  sum.df <- df %>% 
    group_by(CDR3, nonVJ, V.new, J.new, v.end, CDR3.length, dominant) %>% 
    dplyr::summarise(count = n()) %>%
    filter(count > n)
  result <- cast(sum.df, CDR3+nonVJ+V.new+J.new+v.end+CDR3.length~dominant,
                 value = "count", fun.aggregate = sum)
  return(result)
}

GLC.sum <- identDominantClone(a2GLC.high1.0, 0)
GIL.sum <- identDominantClone(a2GIL.high1.0, 0)
NLV.sum <- identDominantClone(a2NLV.high1.0, 0)

################################################################################
# generate datasubset (high ID frequency and exclusion of low SD freq < 1%)
library(dplyr)
a2GLC.high1.0 <- a2GLC[a2GLC$subject.id %in% filterLowFreq(a2GLC, 1.0), ] %>% 
  filter(frequency > 1)

a2GIL.high1.0 <- a2GIL[a2GIL$subject.id %in% filterLowFreq(a2GIL, 1.0), ] %>% 
  filter(frequency > 1)

a2NLV.high1.0 <- a2NLV[a2NLV$subject.id %in% filterLowFreq(a2NLV, 1.0), ] %>% 
  filter(frequency > 1)

# alignment
a2GLChighID.msa <- alignGap(a2GLC.high1.0, "ID")
a2GLChighSD.msa <- alignGap(a2GLC.high1.0, "SD")

a2GILhighID.msa <- alignGap(a2GIL.high1.0, "ID")
a2GILhighSD.msa <- alignGap(a2GIL.high1.0, "SD")

a2NLVhighID.msa <- alignGap(a2NLV.high1.0, "ID")
a2NLVhighSD.msa <- alignGap(a2NLV.high1.0, "SD")

# check the CDR3 that present as ID and SD across individual
remRepISD <- function(sum.df, df){
  # check the CDR3 the are ID and SD
  rem.sd <- as.character(sum.df[sum.df$ID != 0 & sum.df$SD != 0, "CDR3"])
  rem.df <- rbind(df[df$dominant == "SD" & !df$CDR3 %in% rem.sd, ], df[df$dominant == "ID", ])
  return(rem.df)
}
# call function
droplevels(rempRepISD(GLC.sum, a2GLC.high1.0))
################################################################################

HM(log(positionalFrequency(a2NLVhighID.msa, T)/positionalFrequency(a2NLVhighSD.msa, T))) + 
  labs(title = "NLV")
HM(log(positionalFrequency(a2GILhighID.msa, T)/positionalFrequency(a2GILhighSD.msa, T))) + 
  labs(title = "GIL")
HM(log(positionalFrequency(a2GLChighID.msa, T)/positionalFrequency(a2GLChighSD.msa, T))) + 
  labs(title = "GLC")

################################################################################


# Length comparison
plotLength(a2GLC[a2GLC$subject.id %in% filterLowFreq(a2GLC, 1.0), ], 
           "dominant", "GLC 1.0") # change filter(n > 1) to (n > 0)

# V, J analysis
ggplot(a2GLC[a2GLC$subject.id %in% filterLowFreq(a2GLC, 1.0), ] %>% 
         filter(frequency > 1.0), aes(x = V.new, fill = subject.id)) + 
  geom_bar(position = "dodge") + theme(axis.text.x = element_text(angle = 45))

################################################################################
## Part V:individual frequency plot ID and SD ##
################################################################################

indSummary <- function(vdj, sd.cutoff){
  library(tidyr)
  library(dplyr)
  library(plotrix)
  # calculate the mean SD of each individual
  sum.t <- vdj %>%
    filter(frequency >= sd.cutoff) %>%
    group_by(subject.id, dominant) %>% 
    dplyr::summarise(harmonic.mean = length(CDR3_AA_GRAVY)/sum(1/(CDR3_AA_GRAVY)),
                     se = as.numeric(std.error(CDR3_AA_GRAVY)))
    sum.t[is.na(sum.t)] <- 0
  #spread(dominant, mean)'
  return(sum.t)
}

plotIndFreq <- function(sum.t){
  require(tidyr)
  require(dplyr)
  library(ggplot2)
  library(ggpubr)

  p <- ggplot(sum.t, aes(x = dominant, y = harmonic.mean)) +
    geom_jitter(aes(group = subject.id), size = 3.5,
               position = position_dodge(width = 0.25)) +
    geom_line(aes(group = subject.id), size = 1.0,
              position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymax = as.numeric(harmonic.mean) + se, 
                      ymin = as.numeric(harmonic.mean) - se, color = subject.id), 
                  position = position_dodge(width = 0.25)) +
    geom_boxplot(alpha = 0.3, width = 0.3, color = "darkgray") +
    custom_theme +
    theme(legend.position = "none",
          axis.text = element_text(face = "bold", size = 13)) +
    stat_compare_means(aes(x = dominant, y= harmonic.mean), 
                       method = "wilcox.test", label.x = 1.4,
                       size = 5)
  #stat_compare_means(sum.t, aes(x = dominant, y = mean), paired = TRUE, method = "wilcox.test")
  return(p)
} 

################################################################################
# call function #
plotIndFreq(indSummary(a2NLV.high1.0, 1)) + 
  labs(x = "Immune responses", y = "Harmonic mean of CDR3 length (aa)", 
       title = "CDR3 length comparison within subject in the NLV dataset")
plotIndFreq(indSummary(a2GLC[a2GLC$subject.id %in% filterLowFreq(a2GLC, 1.0), ], 1))


################################################################################
# check the prevalence of V-J combination in each data set
################################################################################
mapColorV <- function(df, brewerPal, n){
  library(RColorBrewer)
  vGenes <- c(names(sort(summary(droplevels(df$V.new1)), decreasing = T)), "Others")
  cols <- c(colorRampPalette(brewer.pal(n, name = brewerPal))(length(vGenes)-1), "beige")
  vColors <- cbind(vGenes, cols)
  return(vColors)
}

chordPlot <- function(df, brewerPal, df.color, n){
  lapply(c("tidyr", "dplyr", "circlize", "RColorBrewer"), library, character.only = TRUE)
  
  t2 <- data.frame(table((df[, "V.new1"]))*100/dim(df)[1], stringsAsFactors = F)
  t2$indicator <- ifelse(t2$Freq < 3, "Others", levels(t2$Var1))
  colnames(t2)[1] <- "V.new1"

  t1 <- df %>% 
    group_by(V.new1, J.new1) %>% 
    dplyr::summarise(n = n()) %>% 
    plyr::mutate(f = n*100/sum(n))
  
  t12 <- left_join(t1, t2[, -2], by = "V.new1")
  t <- reshape::cast(t12, indicator~J.new1, value = "n", fun.aggregate = sum)
  t.new <- data.frame(t[,-1], row.names = t[,1])
  t.mat <- as.matrix(t.new)
  # sort column 
  t.mat <- t.mat[names(sort(rowSums(t.mat), decreasing = TRUE)),]
  
  # set color for V segments
  fullColors <- as.data.frame(mapColorV(df.color, brewerPal, n))
  vColors <- sapply(rownames(t.mat), function(v) as.character(fullColors[fullColors$vGenes == v, 2]))
  
  
  #12, "Set3" for gil
  # 11, "Paired" for nlv
  # and 11, "Spectral" for glc
  # plot chord diagram
  cd <- chordDiagram(t.mat, annotationTrack = c("grid"), 
               annotationTrackHeight = uh(6, "mm"),
               row.col = vColors,
               grid.col = c(vColors, setNames(rep("lightgray", 2), c("J1", "J2"))))
  
  for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
                facing = "bending", niceFacing = TRUE, col = "black", 
                cex = 0.8, font = 2)}
  return(cd)
  circos.clear()
  # add % of V genes
  
}


# tranform df 
t <- reshape2::dcast(df, V.new1~J.new1, value.var = "dominant", 
                     fun.aggregate = base::length) 
# export top V
sort(reshape::cast(cdglc_id, V.new1~J.new1, value = "f", fun.aggregate = sum) %>% rowSums)



# call function
chordPlot(a2GIL[a2GIL$dominant == "ID",], "Set3", a2GIL, 12)
chordPlot(a2GIL[a2GIL$dominant == "SD",], "Set3", a2GIL, 12)
chordPlot(a2GLC[a2GLC$dominant == "ID",], "Spectral", a2GLC, 11)
chordPlot(a2GLC[a2GLC$dominant == "SD",], "Spectral", a2GLC, 11)
chordPlot(a2NLV[a2NLV$dominant == "ID",], "Paired", a2NLV, 11)
chordPlot(a2NLV[a2NLV$dominant == "SD",], "Paired", a2NLV, 11)


# combine png files of chord diagrams

rl.gil <- lapply(list("4A_gil_ID.png", "4A_gil_SD.png"), png::readPNG)
rl.glc <- lapply(list("4A_glc_ID.png", "4A_glc_SD.png"), png::readPNG)
rl.nlv <- lapply(list("4A_nlv_ID.png", "4A_nlv_SD.png"), png::readPNG)

gl <- lapply(list(rl.gil,rl.glc, rl.nlv) grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl)




################################################################################
# merge frequency of multiple count of a single clone
################################################################################

# code while developing remMutlClone
x <- droplevels(cd8b[cd8b$subject.id == "human_subject0011**GILGFVFTL*PMID:28636592*25",])
x.new <- x %>% group_by(CDR3) %>% mutate(frequency = sum(frequency)) %>% ungroup()
x.new2 <- data.table(x, key="CDR3") #x.new2 <- aggregate(.~CDR3, x, FUN = paste, collapse=",")
x.new3 <- aggregate(x[,-3], list(x[,3]), function(x) paste0(sum(x$frequency)))
x.new4 <- ddply(x,"CDR3",numcolwise(sum))

DT <- as.data.table(x)
# which columns are numeric 
numeric_cols <- which(sapply(DT, is.numeric))
DT[, lapply(.SD, sum), by = x, .SDcols = numeric_cols]

################################################################################
# calculate properties of CDR3 using library alakaszam and Peptides
################################################################################
callLibrary <- function(lib){
  if(!require(noqoute(lib)))
  {installed.packages(lib)}
  library(lib)
}

callLibrary("alakazam")
callLibrary("Peptides")

a2GIL <- addProperties(a2GIL) # 10 kidera factors and 1 broman
a2GIL <- aminoAcidProperties(a2GIL, seq = "CDR3") # 9 chem prop inclding length

a2GLC <- addProperties(a2GLC)
a2GLC <- aminoAcidProperties(a2GLC, seq = "CDR3") 

# run random forest; no difference between ID and SD but possible 
# to differentiate between GLC and GIL Epitope
t <- droplevels(rbind(a2GIL[, c(10, 18:19, 34:54)], a2GLC[, c(10, 18:19, 34:54)]))
vdj.rf[t,5]

################################################################################
# a2 high1.0 data 
################################################################################
#a2GILhigh1.0
plotAmino(a2GIL.high1.0[a2GIL.high1.0$dominant == "ID",], "nonVJ")

plotAmino(a2GIL.high1.0[a2GIL.high1.0$V == "TRBV19*01",], "nonVJ") + labs(title="GILhigh V19")
plotAmino(a2GIL.high1.0[a2GIL.high1.0$V != "TRBV19*01",], "nonVJ") + labs(title="GILhigh nonV19")

mean(dplyr::pull(a2GIL.high1.0[a2GIL.high1.0$V.new == "TRBV19",], nonVJ.length))

#a2GLChigh1.0
cat(a2GLChighID.msa, sep = "\n")
chordPlot(a2GLC.high1.0[a2GLC.high1.0$dominant == "ID",])



ggplot(cd8b1.0[cd8b1.0$Epitope == "NLVPMVATV",], aes(x = CDR3, y = CDR3_AA_GRAVY, colour = dominant)) + geom_point() + theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) + facet_grid(~CDR3.length)

# 
vdj.rf(a2GIL.high1.0[18:19, 23, 34])

################################################################################
# calculate affinity binding between CDR3 and Epitope using M-J matrix
################################################################################
library(BioPhysConnectoR) # containing code for calling M-J matrix

mj.mat<-as.matrix(read.table(system.file("mj1.txt", package = "BioPhysConnectoR")))
# mj.mat represent the positive value of original MJ matrix 
mj.aa <- c("C","M", "F", "I", "L", "V", "W", "Y", "A", "G", "T", "S", "N", "Q", "D", "E", "H", "R", "K", "P")
mj.ind <- cbind(mj.aa, 1:20)

# add extra row and column containing for gap("-")
mj.mat <- rbind(mj.mat, rep(0,20))
mj.mat <- cbind(mj.mat, rep(0,21))

#assignt rownames and colnames for mj matrix
rownames(mj.mat) <- c(mj.aa, "-")
colnames(mj.mat) <- c(mj.aa, "-")

calAffinity <- function(cdr3, epitope, mj.mat){
  cdr3.split <- unlist(strsplit(cdr3, ""))
  epi.split <- unlist(strsplit(epitope,""))
  # initilise a vector to store the energy
  affinity <- c()
  for(i in 1:nchar(cdr3)){
    affinity <- c(affinity, mj.mat[cdr3.split[i], epi.split[i]])
  }
  # sum the interaction energy and add mimus sign since the original value is negative
  affinity.eng <- -sum(affinity)
  return(affinity.eng) 
}

################################################################################
# GIL-flu #
# select the CDR3 of predominant V-J combination
chordPlot(a2GIL.high1.0)
chordPlot(a2GIL.high1.0[a2GIL.high1.0$dominant == "ID", ])

a2GIL.V19 <- a2GIL.high1.0[a2GIL.high1.0$V.new == "TRBV19",] #& a2GIL.high1.0$J.new == "TRBJ2", ]
gil58 <- "FVFT"

gilid.cdr3 <- unique(substr(as.character(a2GIL.V19[a2GIL.V19$dominant == "ID", "CDR3"]), 5, 8))
gilsd.cdr3 <- unique(substr(as.character(a2GIL.V19[a2GIL.V19$dominant == "SD", "CDR3"]), 5, 8))

# plot
plot(sapply(gilid.cdr3, calAffinity, gil58, mj.mat), pch = 19)
point(sapply(gilsd.cdr3, calAffinity, gil58, mj.mat), pch = 19, col = "red")

# calculate intra-individual affinity comparison between ID and SD responses
indAffinity <- function(vdj, epitope, mjmat){
  library(tidyr)
  library(dplyr)
  library(plotrix)
  # add 2 columns; contact.region & affinity
  contact.region <- sapply(as.character(vdj$CDR3), substr, 
                           vdj$v.end + 1, vdj$CDR3.length -5) 
  vdj$affinity <- sapply(contact.region, calAffinity, epitope, mjmat)
  # calculate the mean SD of each individual
  sum.t <- vdj %>%
    group_by(subject.id, dominant) %>% 
    dplyr::summarise(harmonic.mean = length(affinity)/sum(1/(affinity)),
                     se = as.numeric(std.error(affinity))) 
  sum.t[is.na(sum.t)] <- 0
  #spread(dominant, mean)
  return(sum.t)
}

# call function
t <- indAffinity(a2GIL.V19, gil58, mj.mat)#5, 8)
plotIndFreq(t) + 
  labs(x = "Immune response", y = "Harmonic mean of affinity", title = "GIL_V19")

################################################################################
# GLC-EBV epitope #
chordPlot(a2GLC.high1.0)

a2GLC.V20 <- a2GLC.high1.0[a2GLC.high1.0$V.new == "TRBV20",] 
a2GLC.V29 <- a2GLC.high1.0[a2GLC.high1.0$V.new == "TRBV29",] 
a2GLC.V7 <- a2GLC.high1.0[a2GLC.high1.0$V.new == "TRBV7",] 
glc58 <- "LVAM"

t <- indAffinity(a2GLC.high1.0, glc58, mj.mat, 7, 9)
plotIndFreq(t) + 
  labs(x = "Immune response", y = "Harmonic mean of affinity", title = "GIC_V20")

################################################################################
# generate the unique CDR3 from ID and SD of each data subsets (GIL, GLC, NLV)
gilid.cdr3 <- unique(as.character(a2GIL.high1.0[a2GIL.high1.0$dominant == "ID", "CDR3"]))
gilsd.cdr3 <- unique(as.character(a2GIL.high1.0[a2GIL.high1.0$dominant == "SD", "CDR3"]))
# remove the ID sequece in SD
gilsd.cdr3 <- gilsd.cdr3[!gilsd.cdr3 %in% gilid.cdr3]


glcid.cdr3 <- unique(as.character(a2GLC.high1.0[a2GLC.high1.0$dominant == "ID", "CDR3"]))
glcsd.cdr3 <- unique(as.character(a2GLC.high1.0[a2GLC.high1.0$dominant == "SD", "CDR3"]))
glcsd.cdr3 <- glcsd.cdr3[!glcsd.cdr3 %in% glcid.cdr3]

nlvid.cdr3 <- unique(as.character(a2NLV.high1.0[a2NLV.high1.0$dominant == "ID", "CDR3"]))
nlvsd.cdr3 <- unique(as.character(a2NLV.high1.0[a2NLV.high1.0$dominant == "SD", "CDR3"]))
nlvsd.cdr3 <- nlvsd.cdr3[!nlvsd.cdr3 %in% nlvid.cdr3]

# export the unique sequence with the similar name

write.table(glcsd.cdr3, file = "D:/TCR2Epitope-master/input/glcsd.txt",
            append = FALSE, quote = FALSE, sep = "\n",eol = "\n", 
            row.names = FALSE, col.names = FALSE)

################################################################################ 
# import 1x600 vector of cdr3 from Ewald's pepline
gilid.vt <- read.table("D:/TCR2Epitope-master/output/gilid_vt.txt")
gilsd.vt <- read.table("D:/TCR2Epitope-master/output/gilsd_vt.txt")

################################################################################ 
# affinity calculation mod"el
################################################################################ 

# select the peptide residues
gil.se <- substr("GILGFVFTL", 5, 7)

addGapMotif <- function(motifs){
  # add gap at the begining
  motif.gb <- gsub("^.", "-", motifs)
  # add gap at the end 
  motif.ge <- gsub(".$", "-", motifs)
  # add gap at the middle (this will automatically change motifs)
  substring(motifs,2,2) <- rep('-', length(motifs))
  # add all the motifs with 1 gap
  motifs.gap <- c(motif.gb, motif.ge, motifs)
  return(motifs.gap)
}

reverseString <- function(string){
  rev.s <- paste(rev(substring(string, 1:nchar(string), 1:nchar(string))), collapse = "")
  return(rev.s)
}

generateMotif <- function(cdr3){
  if(nchar(cdr3) == 1){
    motifs.gap <- c(paste0("--", cdr3), paste0("-", cdr3, "-"), paste0(cdr3, "--"))
    } else if(nchar(cdr3) == 2){
    motifs.gap <- c(paste0("-", cdr3), paste0(cdr3, "-"), gsub("(.)(.)", "\\1-\\2", cdr3))
    } else if(nchar(cdr3) > 2){
    motifs <- c()
    i = 1
    while(i < (nchar(cdr3)-3+1+1)){
      motif <- substr(cdr3, i, i+2)
      motifs <- c(motifs, motif)
      i = i + 1
    }
    motifs.gap <- c(motifs, addGapMotif(motifs))
    # add the reverse motif 
    ##rev.motifs <- unname(sapply(motifs.gap, reverseString))
    ##motifs.gap <- c(motifs.gap, rev.motifs)
  }
  return(motifs.gap)
}

checkLength <- function(m){
  m.split <- unlist(strsplit(m, ""))
  count <- length(m.split[m.split %in% LETTERS])
  return(count)
}

calMotifAffinity <- function(cdr3, epitope, mj.mat){
  motifs2 <-  sapply(generateMotif(epitope), 
                   function(e) sapply(generateMotif(cdr3), 
                                      function(m) calAffinity(m, e, mj.mat)/min(sapply(c(m, e), checkLength))))
  
  # print the pair motifs of CDR3 - peptide
  ##ind <- which(motifs2 == min(motifs2), arr.ind = T)
  ##print(paste(colnames(motifs2)[ind[2]], rownames(motifs2)[ind[1]], ","))
  result <- min(motifs2)
  
  return(result)
}

# call funtion
sapply(cdr3.se, calMotifAffinity, gil.se)


# plot 
plotAffinity.ind <- function(df, epitope, mj.mat){
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(plotrix)
  df$CDR3.se <- substr(as.character(df$CDR3), 
                      df$v.end + 1, df$CDR3.length - 5)
  # remove the row containing empty CDR3.se
  df <- df[df$CDR3.se != "",]
  df$affinity <- sapply(df$CDR3.se, calMotifAffinity, epitope, mj.mat)
  # summarise the data
  sum.df <- df %>%
    group_by(subject.id, dominant) %>% 
    dplyr::summarise(harmonic.mean = length(affinity)/sum(1/(affinity)),
                     se = as.numeric(std.error(affinity)))
  # plot
  ggplot(sum.df, aes(x = dominant, y = harmonic.mean)) +
    geom_jitter(aes(group = subject.id), size = 3.5,
                position = position_dodge(width = 0.25)) +
    geom_line(aes(group = subject.id), size = 1.0,
              position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymax = as.numeric(harmonic.mean) + se, 
                      ymin = as.numeric(harmonic.mean) - se, color = subject.id), 
                  position = position_dodge(width = 0.25)) +
    geom_boxplot(alpha = 0.3, width = 0.3, color = "darkgray") +
    custom_theme +
    theme(legend.position = "none",
          axis.text = element_text(face = "bold", size = 13)) +
    stat_compare_means(aes(x = dominant, y= harmonic.mean), 
                       method = "wilcox.test", label.x = 1.4,
                       size = 5)
}

plotAffinity <- function(df, epitope, mj.mat, dataset.label){
  library(ggplot2)
  library(ggpubr)
  df <- df[!duplicated(df[,c("CDR3")]),]
  df$CDR3.se <- substr(as.character(df$CDR3), 
                       df$v.end + 1, df$CDR3.length - 5)
  df <- df[df$CDR3.se != "",]
  df$affinity <- sapply(df$CDR3.se, calMotifAffinity, epitope, mj.mat)
  
  # plot
  ggplot(df, aes(x = dominant, y = affinity)) +
    geom_jitter(size = 3, width = 0.2) +
    geom_boxplot(alpha = 0.3, width = 0.4, color = "darkgray") +
    labs(x = "Immune response", y = "Binding energy (kB)", title = dataset.label) +
    custom_theme +
    theme(legend.position = "none",
          axis.text = element_text(face = "bold", size = 13)) +
    stat_compare_means(aes(label = paste0("p-value = ", ..p.format..)), 
                       method = "wilcox.test", paired = F, 
                       label.x = 1.4, label.y.npc = 'top', size = 5)
}
  
# call the plot
plotAffinity(a2GIL.high1.0, "GFVFT", mj.mat, "GIL-A2")
aff.nlv <- plotAffinity(a2NLV %>% filter(frequency > 1), "PMVAT", mj.mat, "NLV-A2")
aff.glc <- plotAffinity(a2GLC %>% filter(frequency > 1), "TLVAM", mj.mat, "GLC-A2")
aff.gil <- plotAffinity(a2GIL %>% filter(frequency > 1), "GFVFT", mj.mat, "GIL-A2")

library(ggpubr) ## to be done 
figure <- ggarrange(aff.gil, aff.glc, aff.nlv, ncol = 3, nrow = 1)
annotate_figure(figure, bottom = text_grob("Immune response", color = "black",
                                           hjust = 1, x = 0.55, face = "bold", size = 15),
                left = text_grob("Binding energy (kB)", color = "black",face = "bold", size = 15, rot = 90)
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)
  
  
# test set of gil CDR3 from PMID:28423320
b27HIV <- droplevels(cd8b1.0[cd8b1.0$Epitope == "KRWIILGLNK", ])
b57HIV <- droplevels(cd8b1.0[cd8b1.0$Epitope == "KAFSPEVIPMF", ])

plotAffinity(b27HIV, "LGN", mj.mat, "KRW-B27")

################################################################################ 
# contact epitope with CDR3 +/- gap
motifs <- c(sapply(generateMotif(cdr3), calAffinity, epitope, mj.mat)/3, 
            sapply(addGapMotif(generateMotif(cdr3)), calAffinity, epitope, mj.mat)/2)
m1 <- c(epitope, names(motifs)[which.min(motifs)], min(motifs))

# gapped epitope with contact CDR3
motifs2 <-  sapply(c(addGapMotif("FVF")), 
                   function(e) sapply(generateMotif(cdr3), 
                                      function(m) calAffinity(m, e, mj.mat)/2))
ind <- which(motifs2 == min(motifs2), arr.ind = T)
m2 <- c(colnames(motifs2)[ind[2]], rownames(motifs2)[ind[1]], min(motifs2))
print(m1)
print(m2)

################################################################################ 









































