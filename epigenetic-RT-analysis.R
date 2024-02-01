#Analysis replication Timing naive vs formative
library(preprocessCore)
library(DNAcopy)
library(viridis) #needed on cluster
library(readr) #needed on cluster
library(dplyr)
library(ggplot2)
library(umap)
library(RColorBrewer)
library(ComplexHeatmap)
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10) #for GC content
library(chromoMap)

#Set working dir
setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/")
#Load delta RT data
dRT <- read.table("/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-20kb-bins.bed", header=T, sep ="\t")
dRT.gr <- makeGRangesFromDataFrame(dRT, keep.extra.columns = TRUE)

#Load and cleanup epigenetics data at H3K4me1 peaks from Boileau et al 2023 Genome Biol (www.github.com/)
supermatrix.raw <- read_delim("/blellochlab/data1/deniz/go-data/analysis/repli-seq-mll34/withv65/bg/R/epigenetic-rna-data-ryan/supermatrix.k4m1.txt","\t", escape_double = FALSE, trim_ws = TRUE)
supermatrix <- data.frame(supermatrix.raw)
#supermatrix <- supermatrix[abs(supermatrix$pk_fold_nf_wt_k4m1) > 1,] #Filter for dynamic H3K4me1 peaks (based on Boileau et al)
colnames(supermatrix)[1:3] <- c("Chr","Start","End")
supermatrix <- supermatrix[,-c(12:14,29:32,45:56,61:64,67:69,72:74,77:84,89)]
supermatrix$rna_fold_nf_dcd <- supermatrix$FdCD - supermatrix$NdCD
supermatrix$pk_fold_n_dcdwt_k4m1 <- supermatrix$N_dCD_k4m1 - supermatrix$N_WT_k4m1
supermatrix$pk_fold_n_kowt_k4m1 <- supermatrix$N_CKO_k4m1 - supermatrix$N_WT_k4m1
supermatrix.gr <- makeGRangesFromDataFrame(supermatrix, keep.extra.columns = TRUE)
colnames(supermatrix)
columnsToAverage <- c(13:20,23:55) #Define columns to be averaged per region

#Calculate median signal of each dataset per 50kb RT bin
ovlp <- findOverlaps(dRT.gr,supermatrix.gr,ignore.strand=T)
averagedValues  <- apply(supermatrix[subjectHits(ovlp),columnsToAverage], 2, function(x) tapply(x, factor(queryHits(ovlp)), median,na.rm=T))
dRT[,colnames(averagedValues)] <- NA
dRT[unique(queryHits(ovlp)),43:ncol(dRT)] <- averagedValues

#Cleanup columns names
colnames(dRT)[43:ncol(dRT)] <- c(
  'RNA.AVGnaive.WT', 'RNA.AVGformative.WT', 'RNA.AVGnaive.KO', 
  'RNA.AVGformative.KO', 'RNA.AVGnaive.dKO', 'RNA.AVGformative.dKO', 'RNA.AVGnaive.dCD', 'RNA.AVGformative.dCD', 
  'dRNA.WT','dRNA.dKO', 'dRNA.naive.dKOvsWT', 'H3K4me1.AVGnaive.WT', 'H3K4me1.AVGformative.WT', 
  'H3K4me1.AVGnaive.KO','H3K4me1.AVGformative.KO', 'H3K4me1.AVGnaive.dKO', 'H3K4me1.AVGformative.dKO', 'H3K4me1.AVGnaive.dCD', 
  'H3K4me1.AVGformative.dCD', 'H3K27ac.AVGnaive.WT', 'H3K27ac.AVGformative.WT', 'H3K27ac.AVGnaive.dKO', 'H3K27ac.AVGformative.dKO', 
  'ATAC.AVGnaive.WT', 'ATAC.AVGformativeWT', 'ATAC.AVGnaive.dKO', 'ATAC.AVGformative.dKO', 'dH3K4me1.WT', 
  'dH3K27ac.WT', 'dH3K4me1.dKO', 'dH3K27ac.dKO', 'dH3K4me1.naive.dKOvsWT', 'dH3K27ac.naive.dKOvsWT', 
  'dH3K4me1.KO','dH3K4me1.dCD', 'dATAC.WT', 'dATAC.dKO', 'dATAC.naive.dKOvsWT', 
  'dRNA.dCD', 'dH3K4me1.naive.dCDvsWT', 'dH3K4me1.naive.KOvsWT'
)

# 'chr', 'start', 'stop', 'RT.WT.formative.n1', 'RT.WT.formative.n2', 
# 'RT.WT.formative.n3', 'RT.dCD.formative.n1', 'RT.CD.formative.n2', 'RT.dCD.formative.n3', 'RT.dCD.naive.n1', 
# 'RT.dCD.naive.n2', 'RT.dCD.naive.n3', 'RT.dKO.formative.n1', 'RT.dKO.formative.n2', 'RT.dKO.formative.n3', 
# 'RT.dKO.naive.n1', 'RT.dKO.naive.n2', 'RT.dKO.naive.n3', 'RT.KO.formative.n1', 'RT.KO.formative.n2', 
# 'RT.KO.formative.n3', 'RT.KO.naive.n1', 'RT.KO.naive.n2', 'RT.KO.naive.n3', 'RT.WT.naive.n1', 
# 'RT.WT.naive.n2', 'RT.WT.naive.n3', 'RT.AVGformative.WT', 'RT.AVGnaive.WT', 'RT.AVGformative.dCD', 
# 'RT.AVGnaive.dCD','RT.AVGformative.dKO', 'RT.AVGnaive.dKO', 'RT.AVGformative.KO', 'RT.AVGnaive.KO', 
# 'dRT.WT', 'dRT.dCD', 'dRT.dKO', 'dRT.KO', 'dRT.naive.dKOvsKO', 
# 'dRT.naive.dCDvsKO', 'dRT.naive.dCDvsWT',

#Add relative RT for naive
dRT$dRT.naive.dKOvsWT <- dRT$AVGnaive.dKO - dRT$AVGnaive.WT
dRT$dRT.naive.KOvsWT <- dRT$AVGnaive.KO - dRT$AVGnaive.WT

#write.table(dRT, "deltaRT-50kbBins-epi-rna-binned.bed", row.names = F, quote = F, sep ="\t")
# dRT <- read.table("/blellochlab/data1/deniz/analysis/mll-rt-paper/deltaRT-50kbBins-epi-rna-binned.bed", header=T, sep ="\t")
# dRT.gr <- makeGRangesFromDataFrame(dRT, keep.extra.columns = TRUE)

#Heatmaps
##Naive heatmap CKO, DKO, dCD
#Create input and filter for bins that loose H3K4me1 more than s in dKO samples
#Use paste(shQuote(colnames(dRT)), collapse=", ") to quickly extract names to plot
colnames(dRT)
input <- dRT[,c('dRNA.WT','dRNA.dKO','dRNA.dCD','dRT.KO','dRT.dKO','dRT.dCD')]
input <- input[complete.cases(input),]
s=1
input <- dRT[dRT$dH3K4me1.naive.dKOvsWT-dRT$dH3K4me1.naive.KOvsWT < -s ,c(
                                               'chr',
                                               'dH3K4me1.naive.KOvsWT',  'dH3K4me1.naive.dKOvsWT', 'dH3K4me1.naive.dCDvsWT',
                                               'dRT.naive.KOvsWT','dRT.naive.dKOvsWT','dRT.naive.dCDvsWT'
                                               )]
input <- input[complete.cases(input),]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
annotation <- HeatmapAnnotation(df=data.frame("Feature"=c(rep("H3K4me1",3),rep("RT",3))),
                                annotation_height = unit(4,"mm"))

df <- input[region,-1]

library(circlize)
col_fun = colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
dim(df)
Heatmap(df, use_raster = T, top_annotation = annotation, cluster_columns = F, cluster_rows = T, km=5,
        #column_labels = c("Naive","Formative"),
        name = "log2(fc vs WT)", #title of legend
        show_row_dend = F,
        col= col_fun,
        column_labels = c(rep(c("Ctr","dKO","dCD"),2)),
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

#Group data by MLL3/4-dependent RT domains vs independent
s=1
input <- dRT[dRT$dH3K4me1.naive.dKOvsWT-dRT$dH3K4me1.naive.KOvsWT & dRT$dH3K4me1.naive.dCDvsWT-dRT$dH3K4me1.naive.KOvsWT < -s,c(
  'chr',
  'dH3K4me1.naive.KOvsWT',  'dH3K4me1.naive.dKOvsWT', 'dH3K4me1.naive.dCDvsWT',
  'dRT.naive.KOvsWT','dRT.naive.dKOvsWT','dRT.naive.dCDvsWT'
)]
input <- input[complete.cases(input),]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
annotation <- HeatmapAnnotation(df=data.frame("Feature"=c(rep("H3K4me1",3),rep("RT",3))),
                                annotation_height = unit(4,"mm"))

df <- input[region,-1]

groups <- rep("indep.",nrow(df))
groups[df$dRT.naive.dKOvsWT-df$dRT.naive.KOvsWT < -1 & df$dRT.naive.dCDvsWT-df$dRT.naive.KOvsWT < -1] <- "act.-dep."
groups[df$dRT.naive.dKOvsWT-df$dRT.naive.KOvsWT < -1 & df$dRT.naive.dCDvsWT-df$dRT.naive.KOvsWT >= -1] <- "prot-dep."

library(circlize)
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
dim(df)[1]/dim(dRT)[1] #10% of RT bins harbor MLL3/4-dependent K4me1
Heatmap(df, use_raster = T, top_annotation = annotation, cluster_columns = F, cluster_rows = T,
        row_split=groups,
        name = "log2(fc vs WT)", #title of legend
        show_row_dend = F,
        col= col_fun,
        column_labels = c(rep(c("Ctr","dKO","dCD"),2)),
        row_names_gp = gpar(fontsize = 1), # Text size for row names
        row_gap = unit(5, "mm"), border=T
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

##Naive to formative heatmap CKO, DKO, dCD
#Create input and filter for bins that loose H3K4me1 more than s in dKO samples
#Use paste(shQuote(colnames(dRT)), collapse=", ") to quickly extract names to plot
s=1
input <- dRT[abs(dRT$dH3K4me1.KO) > s,c(
  'chr',
  'dH3K4me1.KO', 'dH3K4me1.dKO', 'dH3K4me1.dCD',
  'dRT.KO', 'dRT.dKO', 'dRT.dCD'
)]
input <- input[complete.cases(input),]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]

annotation <- HeatmapAnnotation(df=data.frame("Feature"=c(rep("H3K4me1",3),rep("RT",3))),
                                annotation_height = unit(4,"mm"))

df <- input[region,-1]

library(circlize)
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

Heatmap(df, use_raster = T, top_annotation = annotation, cluster_columns = F,km=10,
        #column_labels = c("Naive","Formative"),
        name = "log2(fc vs WT)", #title of legend
        column_labels = c(rep(c("Ctr","dKO","dCD"),2)),
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))


#Group data by MLL3/4-dependent RT domains vs independent
s=1
input <- dRT[dRT$dH3K4me1.WT > s & dRT$dH3K4me1.KO > s & dRT$dH3K4me1.dKO < s & dRT$dH3K4me1.dCD < s,c(
  'chr',
  'dH3K4me1.KO', 'dH3K4me1.dKO', 'dH3K4me1.dCD',
  'dRT.KO', 'dRT.dKO', 'dRT.dCD'
)]
input <- input[complete.cases(input),]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
annotation <- HeatmapAnnotation(df=data.frame("Feature"=c(rep("H3K4me1",3),rep("RT",3))),
                                annotation_height = unit(4,"mm"))

df <- input[region,-1]

s=0.5
groups <- rep("indep.",nrow(df))
groups[df$dRT.KO > s & df$dRT.dKO < s & df$dRT.dCD < s] <- "act.-dep."
groups[df$dRT.KO > s & df$dRT.dKO < s & df$dRT.dCD > s] <- "prot-dep."

library(circlize)
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
dim(df)[1]/dim(dRT)[1] #3.5% of RT bins harbor MLL3/4-dependent K4me1
Heatmap(df, use_raster = T, top_annotation = annotation, cluster_columns = F, cluster_rows = T,
        row_split=groups,
        name = "log2(fc vs WT)", #title of legend
        show_row_dend = F,
        col= col_fun,
        column_labels = c(rep(c("Ctr","dKO","dCD"),2)),
        row_names_gp = gpar(fontsize = 1), # Text size for row names
        row_gap = unit(5, "mm"), border=T
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

# #Scale each assay
# df_scaled <- df %>% mutate(row = row_number()) %>% 
#   tidyr::pivot_longer(1:9,names_to="sample")
# df_scaled <- data.frame(df_scaled)
# df_scaled$assay <- ifelse(grep("H3K4",df_scaled$sample),"H3K4me1",ifelse(grep("RT",df_scaled$sample),"RT","RNA"))
# df_scaled <- df_scaled %>% group_by(assay) %>% mutate(zscore=scale(value))
# df_scaled <- data.frame(df_scaled[,c(1,2,5)])
# df_scaled <- df_scaled %>% tidyr::pivot_wider(names_from=sample,values_from=zscore,names_sep="_") %>% dplyr::select(-row)
# df_scaled <- data.frame(df_scaled)