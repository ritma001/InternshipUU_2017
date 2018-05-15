#
# Author: Wannisa Ritmahan
# Aim: supplementary figure 

################################################################################
# common funciton
################################################################################
invisible(lapply(c("ggplot2", "ggpubr", "tidyr", "dplyr", "reshape", "plotrix",
                   "grid"), library, character.only = TRUE))

# theme to combine with custom theme for nice figure
alt_theme <- theme(legend.position = "None",
         strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
         strip.background = element_rect(colour = "white", fill = "white"),
         plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

################################################################################
# figure 2: CDR3 length analysis
################################################################################
## import paired TCRv data downcloaded from vdj db
pairTCR <- readTable("D:/InternshipUU_2017/VDJdb/R_VDJdb/vdj220318_pair.tsv")

## select only CD8+ (MHC.class == "MHCI) and no human epitope included
pairTCR <- droplevels(pairTCR[pairTCR$MHC.class == "MHCI" & 
                                pairTCR$Epitope.species != "HomoSapiens", ])
pairTCR <- extractInfo(pairTCR)

## add common characteristics
pairTCR <- addProperties(pairTCR)

## extract high resolution of V-J and MHC molecules
pairTCR <- reduceSolution(pairTCR)


# remove missingpair from complex.id
t <- pairTCR %>% 
  select(complex.id, Epitope, Epitope.length, MHC.A.new, Gene, CDR3.length) %>%
  cast(complex.id + Epitope + Epitope.length + MHC.A.new ~ Gene, value = "CDR3.length")

removeComplex <-  t[is.na(t[,5:6]), 1]
pairDat <- droplevels(pairTCR[!pairTCR$complex.id %in% removeComplex, ])
t <- droplevels(t[!t$complex.id %in% removeComplex, ])
# compare the length between alpha and beta (-1 = CDR3A < CDR3B length and 1 otherwise)
t$comparison <- apply(t[, 5:6], 1, 
                      FUN = function(r) ifelse(r[1] == r[2], "0", ifelse(r[1] < r[2], "-1", "1")))

pairDat$comparison <- unlist(sapply(pairDat$complex.id, 
                                     function(x) t[t$complex.id == x, "comparison"]))
################################################################################
SF2A <- ggplot(pairDat, aes(x = Gene, y = CDR3.length, color = comparison)) +
  geom_jitter(aes(group = complex.id), size = 3.5,
              position = position_dodge(width = 0.25)) +
  geom_line(aes(group = complex.id), size = 0.5,
            position = position_dodge(width = 0.25)) +
  labs(x = "TCR chain", y = "CDR3 length (AA)") +
  scale_color_manual(values = c("red3", "lightgray", "yellow")) +
  scale_x_discrete(labels=c("TRA" = expression(alpha), "TRB" = expression(beta))) +
  custom_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 25)) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################

SF2B <- ggplot(rbind(cd8a.nr,cd8b.nr), aes(x = factor(Gene), y = nonVJ.length)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = Gene), trim = F, alpha = 0.4, adjust = 3) +
  labs(x = "TCR chain", y = "non-VJ length (AA)") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta))) +
  scale_fill_manual(values = c("navy", "orange2")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                     color = "black", size = 8) +
  theme(legend.position = "None", axis.text.x = element_text(size = 25)) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
SF2C <- ggplot(rbind(cd8a.nr,cd8b.nr), aes(x = factor(Gene), y = nchar(Vseq))) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = Gene), trim = F, alpha = 0.4, adjust = 3) +
  labs(x = "TCR chain", y = "V length (AA)") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta))) +
  scale_fill_manual(values = c("navy", "orange2")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  theme(legend.position = "None", axis.text.x = element_text(size = 25)) +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
SF2D <- ggplot(rbind(cd8a.nr,cd8b.nr), aes(x = factor(Gene), y = nchar(Jseq))) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = Gene), trim = F, alpha = 0.4, adjust = 3) +
  labs(x = "TCR chain", y = "J length (AA)") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta))) +
  scale_fill_manual(values = c("navy", "orange2")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  alt_theme +
  ggtitle("D") 


## combine plots ##
# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))

# Arrange the plots
print(SF2A, vp = define_region(row = 1, col = 1))   # Span over two columns
print(SF2B, vp = define_region(row = 2, col = 1))
print(SF2C, vp = define_region(row = 1, col = 2))
print(SF2D, vp = define_region(row = 2, col = 2))

################################################################################
# figure 3: CDR3 length vs epitope length
################################################################################
mhc <- names(sort(summary(cd8b.nr$MHC.A.new), decreasing = TRUE))[1:6]
t.mhc <- droplevels(cd8b.nr[cd8b.nr$MHC.A.new %in% mhc, ])

SF3A <- ggplot(t.mhc, aes(x = dominant, y = CDR3.length)) +
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = dominant), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = " ", y = "CDR3 length (AA)") +
  scale_fill_manual(values = c("lightblue", "darkorange2")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  facet_wrap(~MHC.A.new, scales = "free") + 
  theme(legend.position = "None",
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################

a2epitope <- names(sort(summary(droplevels(pull(cd8b.nr[cd8b.nr$MHC.A.new == "HLA-A*02", 
                                            "Epitope"]))), decreasing = TRUE))[1:6]

SF3B <- ggplot(droplevels(cd8b.nr[cd8b.nr$Epitope %in% a2epitope, ]), 
               aes(x = dominant, y = CDR3.length)) + 
    stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
    geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
    geom_violin(aes(fill = dominant), trim = F, alpha = 0.4, adjust = 2) +
    labs(x = " ", y = "CDR3 length (AA)") +
    scale_fill_manual(values = c("lightblue", "darkorange2")) +
    custom_theme +
    stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                       method = "wilcox.test", size = 8) +
    stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
                 color = "black", size = 8) +
    facet_wrap(~Epitope, scales = "free") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 0),
          strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
          strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))


################################################################################

t <- cd8b.nr[cd8b.nr$Epitope.species %in% c("EBV", "HIV-1") & 
               cd8b.nr$MHC.A.new != "HLA-A*02",]
cor <- c()
for(i in c("EBV", "HIV-1")){
  df <- t[t$Epitope.species == i,]
  r <- round(cor(df$CDR3.length, df$Epitope.length, method = "spearman"),2)
  cor <- c(cor, r)
}

# pair comparison 
pair_comparison <- list(c("8","9"), c("8","10"), c("8", "11"), c("9", "10"), 
                        c("9", "11"), c("10", "11"))

# position of cor

SF3C <- ggplot(t, aes(x = factor(Epitope.length), y = as.numeric(CDR3.length))) +
  stat_boxplot(aes(group = Epitope.length), geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(aes(group = Epitope.length), width = 0.1, fill = "white", color = "gray") +
  geom_violin(aes(group = Epitope.length), fill = "white", 
              color = "lightgray", trim = F, alpha = 0.4, adjust = 2) +
  labs(x = "Epitope length (AA)", y = "CDR3 length (AA)") +
  stat_compare_means(comparisons = pair_comparison, label = "p.signif", 
                     cex = 1.5, method = 'wilcox.test')+
  stat_summary(aes(group = Epitope.length, label=round(..y..,2)), fun.y=mean, geom="text", 
               color = "black", size = 8) +
  facet_wrap(~Epitope.species, scales = "free") +
  custom_theme +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  geom_text(aes(x, y, label=lab),
            data=data.frame(x=4, y=20, lab=c(paste0("r = ",cor[1]), paste0("r = ", cor[2])),
            Epitope.species=c("EBV", "HIV-1")), vjust=1, col = "red", size = 6, fontface = "italic") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))


################################################################################
# summarise data then, plot the sum of paired A/B CDR3 length vs epitope length
t <- pairDat %>% 
  filter(MHC.A.new != "HLA-A*02", !Epitope.length %in% c(12, 13)) %>%
  select(complex.id, Epitope, Epitope.length, MHC.A.new, Gene, CDR3.length) %>%
  cast(complex.id + Epitope + Epitope.length + MHC.A.new ~ Gene, value = "CDR3.length")

t <- t %>%
  mutate(sum.length = apply(t[,5:6], 1, sum))


SF3D <- ggplot(t, aes(x = factor(Epitope.length), y = as.numeric(sum.length))) +
  stat_boxplot(aes(group = Epitope.length), geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(aes(group = Epitope.length), width = 0.1, fill = "white", color = "gray") +
  geom_violin(aes(group = Epitope.length), fill = "white", 
              color = "lightgray", trim = F, alpha = 0.4, adjust = 2) +
  labs(x = "Epitope length (AA)", y = "CDR3 length (AA)") +
  stat_compare_means(label = "p.signif", label.x.npc = 0.4, label.y.npc = 1.0, size = 8)+
  stat_summary(aes(group = Epitope.length, label=round(..y..,2)), fun.y=mean, geom="text", 
               color = "black", size = 8) +
  custom_theme +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("D") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))


################################################################################

SF3E <- ggplot(pairDat[!pairDat$Epitope.length %in% c(12, 13),], 
               aes(x = Gene, y = as.numeric(CDR3.length))) +
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "gray") +
  geom_violin(aes(fill = Gene), color = "lightgray", trim = F, alpha = 0.4, adjust = 2) +
  labs(x = "", y = "CDR3 length (AA)") +
  scale_fill_manual(values = c("navy", "orange2")) +
  scale_x_discrete(labels=c("TRA" = expression(alpha), "TRB" = expression(beta))) +
  stat_compare_means(label = "p.signif", label.x.npc = 0.4, label.y.npc = 1.0, size = 8)+
  stat_summary(aes(group = Epitope.length, label=round(..y..,2)), fun.y=mean, geom="text", 
               color = "black", size = 8) +
  facet_wrap(~Epitope.length, scales = "free") +
  custom_theme +
  ggtitle("E") +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20),
        legend.position = "none") 


################################################################################
## combine plots ##
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(SF3D, vp = define_region(row = 1, col = 1))
print(SF3E, vp = define_region(row = 1, col = 2))

################################################################################
# figure 4: CDR3 hydropathy
################################################################################
# SF4A: hydropathy of top five epitopes from pair A/B CDR3 

pairEpitope <- names(sort(summary(droplevels(pull(pairDat["Epitope"]))), 
                          decreasing = TRUE))[1:5]

SF4A <-  ggplot(pairDat[pairDat$Epitope %in% pairEpitope,]
                , aes(x = Gene, y = CDR3_AA_GRAVY)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = Gene), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = " ", y = "CDR3 hydropathy", title = "") +
  scale_fill_manual(values = c("navy", "orange2")) +
  scale_x_discrete(labels=c("TRA" = expression(alpha), "TRB" = expression(beta))) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  facet_wrap(~Epitope, scales = "free") + 
  theme(legend.position = "None",
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white"),
        plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20)) +
  ggtitle("A") 

################################################################################
# SF3B: hydropathy of top five epitopes from pair A/B CDR3 

combineHLAB <- droplevels(cd8b.nr[cd8b.nr$MHC.A.new %in% 
      c("HLA-B*81", "HLA-B*51", "HLA-B*53", "HLA-B*42", "HLA-B*35", "HLA-B*07"), ])
  
SF4B <- ggplot(combineHLAB, aes(x = dominant, y = CDR3_AA_GRAVY)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = dominant), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = " ", y = "CDR3 hydropathy") +
  scale_fill_manual(values = c("lightblue", "darkorange2")) +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  theme(legend.position = "None",
        axis.text.x = element_text(angle = 0),
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# SF3C: CDR3 hydropathy between MHC molecules encoded from A locus 
# (HLA-A*01, -A*02 and -A*11) and  the B locus (HLA-B*07, -B*27 and -B*42)

t.mhc <- cd8b.nr[cd8b.nr$MHC.A.new %in% c("HLA-A*01", "HLA-A*02", "HLA-A*11",
                                          "HLA-B*07", "HLA-B*27", "HLA-B*42"),]
t.mhc$MHC.locus <- c(ifelse(t.mhc$MHC.A.new %in% c("HLA-A*01", "HLA-A*02", "HLA-A*11"), 
                            "A", "B"))

SF4C <- ggplot(t.mhc, aes(x = MHC.locus, y = CDR3_AA_GRAVY)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = MHC.locus), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = "HLA locus", y = "CDR3 hydropathy") +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  theme(legend.position = "None",
        axis.text.x = element_text(angle = 0),
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
# SF3D: Epitope hydropathy between A and B loci

SF4D <- ggplot(t.mhc, aes(x = MHC.locus, y = Epitope_AA_GRAVY)) + 
  stat_boxplot(geom = "errorbar", width = 0.05, size = 1) +
  geom_boxplot(width = 0.1, fill = "white", color = "darkgray") +
  geom_violin(aes(fill = MHC.locus), trim = F, alpha = 0.4, adjust = 2) +
  labs(x = "HLA locus", y = "Epitope hydropathy") +
  custom_theme +
  stat_compare_means(aes(label = ..p.signif..),label.x = 1.5, 
                     method = "wilcox.test", size = 8) +
  stat_summary(aes(label=round(..y.., 2)), fun.y = mean, geom = "text", 
               color = "black", size = 8) +
  theme(legend.position = "None",
        axis.text.x = element_text(angle = 0),
        strip.text.x = element_text(margin = margin(0.28,0,0.28,0, "cm")),
        strip.background = element_rect(colour = "white", fill = "white")) +
  ggtitle("D") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 0.2, size = 20))

################################################################################
## combine plots ##
grid.newpage()
# Create layout
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3)))
# Arrange the plots
print(SF4B, vp = define_region(row = 1, col = 1))
print(SF4C, vp = define_region(row = 1, col = 2))
print(SF4D, vp = define_region(row = 1, col = 3))


################################################################################
# figure 5: Relative freq of AA in the top 6 MHC #
################################################################################
# AA composition of nonVJ from the in the top six MHC ##

plotAA.mhc <- function(df, seq){
  # get count table of ID and SD from combining countAmino
  invisible(lapply(c("dplyr", "tidyr", "ggplot2", "ggsignif", "ggpubr"), library ,
                   character.only = TRUE))
  
  mhc <- names(summary(droplevels(df$MHC.A.new)))
  # iteratively add the amino acid count
  
  aa.count <- countAmino(as.character(df[df$MHC.A.new == mhc[1], seq]))
  aa.count$MHC <- rep(mhc[1], 20)
  for(i in mhc[-1]){
    count <- countAmino(as.character(df[df$MHC.A.new == i, seq]))
    count$MHC <- rep(i, 20)
    aa.count <- rbind(aa.count, count)
  }
  
  # calculate the frequency
  aa.freq <- aa.count %>% 
    group_by(MHC) %>%
    dplyr::mutate(freq = count*100/sum(count))
  
  
  # group the amino acid based on hydropathy
  aa <- c("I", "V", "L", "F", "C","M", "A", "W", 
          "G", "T", "S", "Y", "P", "H", 
          "N", "D", "Q", "E", "K", "R")
  hd <- c(rep("Hydrophobic", 8), rep("Neutral",6), rep("Hydrophilic",6))
  hd.index <- data.frame(t(rbind(aa, hd)))
  
  
  # addd column to aa.freq dataframe
  aa.freq$hd <- hd.index$hd[match(aa.freq$amino.acid, hd.index$aa)]
  # add locus 
  aa.freq$MHC.locus <- c(ifelse(aa.freq$MHC %in% c("HLA-A*01", "HLA-A*02", "HLA-A*03", "HLA-A*11"), 
                                "A", "B"))
  
  for(i in amino.acid){
    ks <- kruskal.test(freq ~ factor(MHC), data = dropleveles(aa.freq[aa.freq$amino.acid == i, ]))
    cat(paste(i, round(ks$p.value,5),"\n",  sep = "\t"))
  }
  
  # plot
  p <- ggplot(aa.freq, aes(x = amino.acid, y = freq, fill = MHC)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(name = "HLA allels", palette = "Set3") +
    labs(x = "Amino acid", y = "Relative frequency (%)", fill = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
    custom_theme +
    theme(legend.position = "bottom", 
          legend.background = element_rect(fill=alpha('white', 0.5)),
          strip.background = element_rect(color = NA, fill = "white")) +
    facet_wrap(MHC.locus~hd, scales = "free", shrink = TRUE) 
  
  return(p)
}

SF5A <- plotAA.mhc(t.mhc, "nonVJ")


################################################################################
# AA comparison between naive and non-naive of alpha and beta chains

SF6A <- plotpLES(selectMSA(cd8b.nr, c(10:15), "pLES")) +
  scale_y_discrete(labels =  c(as.character(104:111), "112.1", "111.1",
                               as.character(112:118)))

################################################################################
# figure 7: RF classifiers from 24 features
################################################################################
SF7A <- vdjRF(a2GIL[, c(18:19, 23, 28, 34, 36:53, 57)], 5, T)
SF7B <- vdjRF(a2GLC[, c(18:19, 23, 28, 34, 36:53, 57)], 5, T)
SF7C <- vdjRF(a2NLV[, c(18:19, 23, 28, 34, 36:53, 57)], 5, T)

################################################################################
# pca on the the dataset
vdj.pca <- function(vdj, color.factor){
  # import library #
  library(ggfortify)
  library(pca3d)
  
  # convert type factor to character and interger to numeric
  vdj <- convertType(is.integer, as.numeric, vdj)
  vdj <- convertType(is.factor, as.character, vdj)
  
  # pca #
  num.factors <- sapply(vdj, is.numeric)
  pca <- prcomp(na.omit(vdj[,num.factors]), center = T, scale. = T)
  
  p <- autoplot(pca, 
                data = na.omit(vdj), 
                colour = color.factor,
                size = 5) +
    scale_color_manual(values = c("lightblue", "darkorange")) +
    custom_theme +
    theme(legend.position = "bottom")

  return(p)
}

################################################################################

SF7D <- vdj.pca(a2GIL[, c(34, 36:53)], "dominant")
SF7E <- vdj.pca(a2GLC[, c(34, 36:53)], "dominant")
SF7F <- vdj.pca(a2NLV[, c(34, 36:53)], "dominant")
















