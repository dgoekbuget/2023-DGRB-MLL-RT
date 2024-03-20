#This script is used to perform loess smoothening, segmentation, statistics, and other RT analyses done in manuscript
#The input is mm10 mapped counts turned into log2 ratio by dividing Early/Late counts as detailed in Marchal et al 2018 Nat Prot

#Load packages
library(readr)
library(preprocessCore)
library(DNAcopy)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(rafalib)
library(GenomicRanges)
library(RColorBrewer)
library(MASS)
library(ComplexHeatmap)
library(magick)
library(rafalib)
options(scipen = 0)

#Define file paths and other parameters
plotdir <- "/blellochlab/data1/deniz/analysis/mll-rt-paper/exports" #where to export plots
setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/mapping/bg-50kb-windows")

##Colors
geno_colors <- c("black","#999999","#CC3333","#996699")

#Import and prepare RT data
merge <- read_delim("merge_RT.txt", delim = "\t", 
                 escape_double = FALSE, trim_ws = TRUE)

merge_values<-as.matrix(merge[,4:ncol(merge)])
ad<-stack(merge[,4:ncol(merge)])$values
norm_data<-normalize.quantiles.use.target(merge_values,ad)
merge_norm<-data.frame(merge[,1:3],norm_data)
colnames(merge_norm)<-colnames(merge)

dat <- merge_norm
dat <- dat[complete.cases(dat),]
ind <- c(grep("chrM",dat$chr),grep("chrY",dat$chr))
dat <- dat[-ind,]
colnames(dat)[4:ncol(dat)] <- gsub(".bg","",colnames(dat)[4:ncol(dat)])
dat$start <- as.integer(dat$start)
dat$stop <- as.integer(dat$stop)

#Perform Loess smoothening and export bedgraph files for plotting (after bigwig conversion)
# # Needs to be done once only
# # Generate qnorm bedgraphs per sample
# for(i in 4:ncol(merge_norm)){write.table(merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)],
#                                          gsub(".bg" , ".qnorm.bg", colnames(merge_norm)[i]),
#                                           sep= "\t" ,row.names=FALSE, quote=FALSE, col.names = FALSE)}
# # Perform Loess smoothing
# ind <- c(grep("chrM",merge_norm$chr),grep("chrY",merge_norm$chr))
# chrs=unique(merge_norm$chr[-ind])
# AllLoess=list()
# for(i in 1:(ncol(merge_norm)-3)){
#   AllLoess[[i]]=data.frame();
#   cat("Current dataset:", colnames(merge_norm)[i+3], "\n");
#   for(Chr in chrs){
#     RTb=subset(merge_norm, merge_norm$chr==Chr);
#     lspan=300000/(max(RTb$start)-min(RTb$start));
#     cat("Current chrom:" , Chr, "\n");
#     RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
#     RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, merge_norm[which( merge_norm$chr==Chr & merge_norm$start %in% RTla$x),3],RTla$fitted);
#     colnames(RTl)=c("chr" , "start" , "end" ,colnames(RTb)[i+3]);
#     if(length(AllLoess[[i]])!=0){AllLoess[[i]]=rbind(AllLoess[[i]],RTl)};if(length(AllLoess[[i]])==0){AllLoess[[i]] = RTl}
#   }
# }
# #Write loess smoothened bedgraphs
# for(i in 1:length(AllLoess)){write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),],
#                                          gsub(".bg" , ".Loess.bg" , colnames(AllLoess[[i]]))[4], sep= "\t",
#                                          row.names=FALSE, quote=FALSE, col.names = FALSE)}

#Load Loess smoothened data
merge_loess <- read_delim("merge_Loess_norm_RT.txt", delim = "\t", 
                    escape_double = FALSE, trim_ws = TRUE)
dat_loess <- data.frame(merge_loess)
dat_loess <- dat_loess[complete.cases(dat_loess),]
colnames(dat_loess)[4:ncol(dat_loess)] <- gsub(".bg","",colnames(dat_loess)[4:ncol(dat_loess)])
dat_loess$start <- as.integer(dat_loess$start)
dat_loess$stop <- as.integer(dat_loess$stop)

#Reorder sample columns
dat <- dat_loess[,c(1:3,25:27,22:24,19:21,16:18,13:15,10:12,7:9,4:6)]
dat_loess <- dat_loess[,c(1:3,25:27,22:24,19:21,16:18,13:15,10:12,7:9,4:6)]

#Average replicates and generate the various RT comparison
dat$AVGnaive.WT <- (dat[,4] + dat[,5] + dat[,6])/3
dat$AVGformative.WT <- (dat[,7] + dat[,8] + dat[,9])/3
dat$AVGnaive.KO <- (dat[,10] + dat[,11] + dat[,12])/3
dat$AVGformative.KO <- (dat[,13] + dat[,14] + dat[,15])/3
dat$AVGnaive.dKO <- (dat[,16] + dat[,17] + dat[,18])/3
dat$AVGformative.dKO <- (dat[,19] + dat[,20] + dat[,21])/3
dat$AVGnaive.dCD <- (dat[,22] + dat[,23] + dat[,24])/3
dat$AVGformative.dCD <- (dat[,25] + dat[,26] + dat[,27])/3

dat_loess$AVGnaive.WT <- (dat_loess[,4] + dat_loess[,5] + dat_loess[,6])/3
dat_loess$AVGformative.WT <- (dat_loess[,7] + dat_loess[,8] + dat_loess[,9])/3
dat_loess$AVGnaive.KO <- (dat_loess[,10] + dat_loess[,11] + dat_loess[,12])/3
dat_loess$AVGformative.KO <- (dat_loess[,13] + dat_loess[,14] + dat_loess[,15])/3
dat_loess$AVGnaive.dKO <- (dat_loess[,16] + dat_loess[,17] + dat_loess[,18])/3
dat_loess$AVGformative.dKO <- (dat_loess[,19] + dat_loess[,20] + dat_loess[,21])/3
dat_loess$AVGnaive.dCD <- (dat_loess[,22] + dat_loess[,23] + dat_loess[,24])/3
dat_loess$AVGformative.dCD <- (dat_loess[,25] + dat_loess[,26] + dat_loess[,27])/3

dat$dRT.WT <- dat$AVGformative.WT - dat$AVGnaive.WT
dat$dRT.dCD <- dat$AVGformative.dCD - dat$AVGnaive.dCD
dat$dRT.dKO <- dat$AVGformative.dKO - dat$AVGnaive.dKO
dat$dRT.KO <- dat$AVGformative.KO - dat$AVGnaive.KO

dat$dRT.naive.dKOvsKO <- dat$AVGnaive.dKO - dat$AVGnaive.KO
dat$dRT.naive.dCDvsKO <- dat$AVGnaive.dCD - dat$AVGnaive.KO
dat$dRT.naive.dKOvsWT <- dat$AVGnaive.dKO - dat$AVGnaive.WT
dat$dRT.naive.dCDvsWT <- dat$AVGnaive.dCD - dat$AVGnaive.WT

dat_loess$dRT.WT <- dat_loess$AVGformative.WT - dat_loess$AVGnaive.WT
dat_loess$dRT.dCD <- dat_loess$AVGformative.dCD - dat_loess$AVGnaive.dCD
dat_loess$dRT.dKO <- dat_loess$AVGformative.dKO - dat_loess$AVGnaive.dKO
dat_loess$dRT.KO <- dat_loess$AVGformative.KO - dat_loess$AVGnaive.KO

dat_loess$dRT.naive.dKOvsKO <- dat_loess$AVGnaive.dKO - dat_loess$AVGnaive.KO
dat_loess$dRT.naive.dCDvsKO <- dat_loess$AVGnaive.dCD - dat_loess$AVGnaive.KO
dat_loess$dRT.naive.dKOvsWT <- dat_loess$AVGnaive.dKO - dat_loess$AVGnaive.WT
dat_loess$dRT.naive.dCDvsWT <- dat_loess$AVGnaive.dCD - dat_loess$AVGnaive.WT

#write.table(dat_loess,"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-50kb-bins-loess.bed", row.names = F, quote = F, sep ="\t")

#Circular binary segmentation
library(DNAcopy)
set.seed(178)

##Empiric determination of optimized segmentation paramters
input <- dat_loess[1:1000,]

mypar(4,5)
dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
for (i in seq(1,20,by=2)) {
  for(j in c(1e-15,1e-200)){
  seg.cna  = segment(dat.cna, nperm = 1000, alpha = j, undo.splits = "sdundo", 
                     undo.SD = i, verbose = 0)
  plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
       pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
       )
  legend("topright",legend=paste0("alpha=",j," undo.SD=",i),cex=0.8,box.lty=0,bg="transparent")
  }
}

mypar(4,5)
dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "delta.WT")
for (i in seq(1,20,by=2)) {
  for(j in c(1e-15,1e-200)){
    seg.cna  = segment(dat.cna, nperm = 1000, alpha = j, undo.splits = "sdundo", 
                       undo.SD = i, verbose = 0)
    plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
         pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
    )
    legend("topright",legend=paste0("alpha=",j," undo.SD=",i),cex=0.8,box.lty=0,bg="transparent")
  }
}
#Plot selected parameters (SD=5 and alpha 1e-15 for 50kb bins)
mypar(1,2)
dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
seg.cna  = segment(dat.cna, nperm = 1000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 5, verbose = 2)
plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
     pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
)
legend("topright",legend=paste0("alpha=",1e-15," undo.SD=",5),cex=0.8,box.lty=0,bg="transparent")

dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "delta.WT")
seg.cna  = segment(dat.cna, nperm = 1000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 5, verbose = 0)
plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
     pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
)
legend("topright",legend=paste0("alpha=",1e-15," undo.SD=",5),cex=0.8,box.lty=0,bg="transparent")


##Segmentation of naive and delta(naive-to-form)
input <- dat_loess

dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
seg.cna.naive  = segment(dat.cna, nperm = 10000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 5, verbose = 2)

dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "dRT.WT")
seg.cna.dRT  = segment(dat.cna, nperm = 10000, alpha = 1e-15, undo.splits = "sdundo", 
                         undo.SD = 5, verbose = 2)

##Descriptive statistics on segments
##Naive size
mypar(1,2)
hist(log10(seg.cna.naive$output$loc.end-seg.cna.naive$output$loc.start),breaks=20,
     main=paste0("Naive Segments (N=",nrow(seg.cna.naive$output),")"),
     xlab="Log10(bp)",xlim=c(5,7.5))
abline(v=median(log10(seg.cna.naive$output$loc.end-seg.cna.naive$output$loc.start)),lty="dashed",col="red")
##dRT
hist(log10(seg.cna.dRT$output$loc.end-seg.cna.dRT$output$loc.start),breaks=20,
     main=paste0("Diff. Segments (N=",nrow(seg.cna.dRT$output),")"),xlab="Log10(bp)",xlim=c(5,7.5))
abline(v=median(log10(seg.cna.dRT$output$loc.end-seg.cna.dRT$output$loc.start)),lty="dashed",col="red")


##PCA WT naive + formative
library(rafalib)
set.seed(178)
rt <- dat_loess[,grepl("WT.naive.n|WT.form.n",colnames(dat_loess))]
rt <- rt[complete.cases(rt),]
y <- rt - rowMeans(rt)
s <- svd(y)
label <- factor(colnames(y))

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100

##Plot variance explained
plot(pvar,ylab="Percent variability explained",cex=2,bg="gray",pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

##PC1 vs PC2
mycols <- rep(c('#8491b4ff','#00a087ff'),each=3)
set.seed(111)
# stripchart(PC1, method = "jitter",pch=21,bg=mycols,cex=2,vertical=T,
#            ylab=paste0("PC1 (",round(pvar[1],1),"% variance explained)"))

plot(PC1,PC2,xlab = paste0("PC1 (",round(pvar[1]),"% var)"), ylab = paste0("PC2 (",round(pvar[2]),"% var)"),
     bg = ifelse(grepl("naive",label),"#8491B4","#00a087ff"),cex = 2,pch=21)

##Heatmap RT WT naive + formative
set.seed(178)
input <- dat_loess[,grep("WT.naive.n|WT.form.n|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
km <-  kmeans(df, 10)
split <- paste0("C\n", km$cluster)
Heatmap(df, use_raster = T, cluster_columns = T,cluster_rows = F,split=split,
        column_labels = rep(c("Naive","Formative"),each=3),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        border=T, row_gap=unit(2,"mm"),
        row_names_gp = gpar(fontsize = 0) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

#Volcano plot RT dRT
set.seed(178)
pval.dRT <- apply(
  cbind(dat_loess$WT.form.n1.RT.Loess-dat_loess$WT.naive.n1.RT.Loess,
        dat_loess$WT.form.n2.RT.Loess-dat_loess$WT.naive.n2.RT.Loess,
        dat_loess$WT.form.n3.RT.Loess-dat_loess$WT.naive.n3.RT.Loess),1,function(row){
         ttest <- t.test(row,mu=0)
         pval <- ttest$p.value
         return(pval)
        }
)
pval.dRT.corr <- p.adjust(pval.dRT,method="fdr")
# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.WT,-log10(pval.dRT.corr), n = 250)
# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
mypar()
image(density_est, col = mycols, xlab = "dRT WT", ylab = "-log10(adj p-val)",xlim=c(-3.5,3.5),ylim=c(0,2),useRaster=T)
abline(h=-log10(0.05),lty=2)
legend("top",legend=paste0(round(100*(sum(pval.dRT.corr < 0.05)/nrow(dat_loess)),1),"% of genome"),cex=0.8,box.lty=0,bg="transparent")

#Sankey plots delta-RT (naive to form), WT
ind <- ind.p.wt <- pval.dRT.corr < 0.05  #define p-value cutoff
dat_sk <- dat[,c("AVGnaive.WT","AVGformative.WT","dRT.WT")]
dat_sk$naive.state <- ifelse(dat_sk$AVGnaive.WT > 0,"early","late")
dat_sk$transition.state <- ifelse(dat_sk$dRT.WT > 0 & ind,"earlier",
                                  ifelse(dat_sk$dRT.WT < 0 & ind,"later","unchanged"
                                  ))
table(dat_sk$transition.state)
dat_sk <- dat_sk[,4:5] %>% make_long(naive.state,transition.state)


dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#Extended Data Fig.3 - PCA WT, CKO, DKO - naive/form
set.seed(178)
rt <- dat_loess[,grepl("naive.n|form.n",colnames(dat_loess))]
rt <- rt[,!grepl("dCD",colnames(rt))]
y <- rt - rowMeans(rt)
s <- svd(y)
label <- factor(colnames(y))

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100
mycols <- brewer.pal(12,"Paired")
plot(pvar,ylab="Percent variability explained",cex=2,bg="gray",pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

plot(PC1,PC2,xlab = paste0("PC1 (",round(pvar[1]),"% var)"), ylab = paste0("PC2 (",round(pvar[2]),"% var)"),
     bg = ifelse(grepl("WT",label),"black",ifelse(grepl("dKO",label),"#CC3333",ifelse(grepl("dCD",label),"#996699","#999999"))),
     pch=ifelse(grepl("naive",label),21,23),cex = 2)

#Extended Data Fig.3 - PCA WT, CKO, DCD - naive
set.seed(178)
rt <- dat_loess[,grepl("naive.n|form.n",colnames(dat_loess))]
rt <- rt[,!grepl("dKO",colnames(rt))]
y <- rt - rowMeans(rt)
s <- svd(y)
label <- factor(colnames(y))

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100
mycols <- brewer.pal(12,"Paired")
plot(pvar,ylab="Percent variability explained",cex=2,bg="gray",pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

plot(PC1,PC2,xlab = paste0("PC1 (",round(pvar[1]),"% var)"), ylab = paste0("PC2 (",round(pvar[2]),"% var)"),
     bg = ifelse(grepl("WT",label),"black",ifelse(grepl("dKO",label),"#CC3333",ifelse(grepl("dCD",label),"#996699","#999999"))),
     pch=ifelse(grepl("naive",label),21,23),cex = 2)

##Heatmap RT WT,CKO,DKO,dCD naive
set.seed(178)
head(dat_loess)
input <- dat_loess[,grep("AVGnaive.WT|AVGnaive.KO|AVGnaive.dKO|AVGnaive.dCD|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
km <-  kmeans(df, 5)
split <- paste0("C\n", km$cluster)
Heatmap(df, use_raster = T, cluster_columns = F,cluster_rows = F,split=split,
        column_labels = c("WT","3KO","dKO","dCD"),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        border=T, row_gap=unit(2,"mm"),
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

##Heatmap RT DKO,dCD vs KO naive
set.seed(178)
input <- dat_loess[,grep("AVGnaive.WT|AVGnaive.KO|AVGnaive.dKO|AVGnaive.dCD|chr",colnames(dat_loess))]
input <- data.frame(chr=input$chr,
                    dKOvsWT=input$AVGnaive.dKO-input$AVGnaive.KO,
                    dCDvsWT=input$AVGnaive.dCD-input$AVGnaive.KO)
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue", "white", "orange"))
km <-  kmeans(df, 5)
split <- paste0("C\n", km$cluster)
Heatmap(df, use_raster = T, raster_by_magick = T, cluster_columns = F, cluster_rows = F, split=split,
        column_labels = c("dKO","dCD"),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        border=T, row_gap=unit(2,"mm"),
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

#Control RT of changing sites in naive
boxplot(dat_loess$AVGnaive.KO[abs(dat_loess$dRT.naive.dKOvsKO) > 1 & abs(dat_loess$dRT.naive.dCDvsKO) > 1],
        ylab="Naive control RT")

##Fold change PCA
set.seed(178)
rt <- dat_loess[,grepl("dRT.",colnames(dat_loess))]
rt <- rt[complete.cases(rt),]
y <- rt - rowMeans(rt)
y <- y[,c(1,4,3,2)] #reorder colnames by genotype
s <- svd(y)

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100
mycols <- brewer.pal(4,"Paired")
plot(pvar,ylab="Percent variability explained",cex=2,bg="gray",pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

mypar(1,1)
plot(PC1,PC2,xlab=paste0("PC1 (",round(pvar[1]),"% var)"), 
     ylab = paste0("PC2 (",round(pvar[2]),"% var)"),
     bg = geno_colors, pch=21,cex = 2)
legend("bottomleft",colnames(y),col=geno_colors,
       pch=16,box.lwd=1,cex = 1.5)

##Heatmap delta RT 3KO,dKO,dCD
set.seed(178)
input <- dat_loess[,grep("dRT.KO|dRT.dKO|dRT.dCD|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
df <- df[order(df$dRT.KO,decreasing = T),]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T,raster_by_magick = T,raster_resize_mat = T,cluster_rows = F,cluster_columns = T,
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        col = col_fun,
        border=T,
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

mypar()
col_fun = colorRamp2(c(-5, -2, 0, 2, 5), c("royalblue","dodgerblue", "white", "orange","orange3"))
col_fun(seq(0, 1, length = 20))
plot(NULL, xlim = c(-5, 5), ylim = c(0, 1))
x = seq(-5, 5, length = 20)
y = rep(0.5, 20)
points(x, y, pch = 16, col = col_fun(x), cex = 2)

#Correlation of RT changes between dCD and dKO in naive
mypar(1,1)
# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.naive.dKOvsKO,dat_loess$dRT.naive.dCDvsKO, n = 250)
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
image(density_est,
     useRaster=T,
     xlab="dKOvsKO",
     ylab = "dCDvsKO",
     col=mycols,
     xlim = c(-2.5,2.5),ylim=c(-2.5,2.5),
     main=paste0("Naive RT, rho=",round(cor(dat_loess$dRT.naive.dKOvsKO,dat_loess$dRT.naive.dCDvsKO,method="spearman"),3)))
abline(lm(dat_loess$dRT.naive.dKOvsKO~dat_loess$dRT.naive.dCDvsKO),col="dodgerblue",lty="dashed",lwd=3)

#Correlation of deltaRT changes dCD and dKO
x <- dat_loess$dRT.dKO-dat_loess$dRT.KO
y <- dat_loess$dRT.dCD-dat_loess$dRT.KO

# Calculate the 2D density using kde2d()
density_est <- kde2d(x, y, n = 250)

# Create a contour plot
image(density_est,
        useRaster=T,
        xlab="ddRT (dKO vs KO)",
        ylab = "ddRT (dCD vs KO)",
        col=mycols,
        xlim = c(-2.5,2.5),ylim=c(-2.5,2.5),
        main=paste0("ddRT, rho=",round(cor(x,y,method="spearman"),3)))
abline(lm(x~y),col="dodgerblue",lty="dashed",lwd=3)

#Histogram of RT distribution in naive
##compute densities
density1 <- density(dat_loess[,"AVGnaive.dKO"]-dat_loess[,"AVGnaive.KO"])
density2 <- density(dat_loess[,"AVGnaive.dCD"]-dat_loess[,"AVGnaive.KO"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-3,3),xlab="naive RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[3], lwd = 3)
lines(density2, col = geno_colors[4], lwd = 3)
legend("topright", legend = c("dKO", "dCD"), col = geno_colors[3:4], lty = 1, lwd = 2)


#Histogram of RT distribution in formative dKO vs ctr
##compute densities
density1 <- density(dat_loess[,"AVGformative.KO"])
density2 <- density(dat_loess[,"AVGformative.dKO"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-3,3),xlab="formative RT", main="Formative RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[2], lwd = 3)
lines(density2, col = geno_colors[3], lwd = 3)
legend("topright", legend = c("Ctr","dKO"), col = geno_colors[2:3], lty = 1, lwd = 2)

#Histogram of RT distribution in formative
##compute densities
density1 <- density(dat_loess[,"AVGnaive.KO"])
density2 <- density(dat_loess[,"AVGnaive.dKO"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-3,3),xlab="naive RT", main="Naive RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[2], lwd = 3)
lines(density2, col = geno_colors[3], lwd = 3)
legend("topright", legend = c("Ctr","dKO"), col = geno_colors[2:3], lty = 1, lwd = 2)

#Histogram of RT distribution in naive-to-form transition
##compute densities
density1 <- density(dat_loess[,"dRT.KO"])
density2 <- density(dat_loess[,"dRT.dKO"])
density3 <- density(dat_loess[,"dRT.dCD"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-3,3),ylim=c(0,1),xlab="dRT",title="")

# Plot density estimates as lines
lines(density1, col = geno_colors[2], lwd = 3)
lines(density2, col = geno_colors[3], lwd = 3)
lines(density3, col = geno_colors[4], lwd = 3)
legend("topright", legend = c("KO", "dKO", "dCD"), col = geno_colors[-1], lty = 1, lwd = 2)
ks.test(dat_loess[,"dRT.KO"], dat_loess[,"dRT.dKO"])
ks.test(dat_loess[,"dRT.KO"], dat_loess[,"dRT.dCD"])

#Volcano plot dRT
##dKO
set.seed(178)
pval.dRT <- apply(
  cbind(dat_loess$dKO.form.n1.RT.Loess-dat_loess$dKO.naive.n1.RT.Loess,
        dat_loess$dKO.form.n2.RT.Loess-dat_loess$dKO.naive.n2.RT.Loess,
        dat_loess$dKO.form.n3.RT.Loess-dat_loess$dKO.naive.n3.RT.Loess),1,function(row){
          ttest <- t.test(row,mu=0)
          pval <- ttest$p.value
          return(pval)
        }
)
pval.dRT.corr <- p.adjust(pval.dRT,method="fdr")

# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.dKO,-log10(pval.dRT.corr), n = 250)
# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
mypar()
image(density_est, col = mycols, xlab = "dRT dKO", ylab = "-log10(adj p-val)",xlim=c(-3.5,3.5),ylim=c(0,2),useRaster=T)
abline(h=-log10(0.05),lty=2)
legend("top",legend=paste0(round(100*(sum(pval.dRT.corr < 0.05)/nrow(dat_loess)),1),"% of genome"),cex=0.8,box.lty=0,bg="transparent")

#Sankey plots delta-RT (WT vs dKO)
ind <- pval.dRT.corr < 0.05  #define p-value cutoff
dat_sk <- dat[,c("dRT.WT","dRT.dKO")]
dat_sk$sign <- ind
dat_sk$wt.state <- ifelse(dat_sk$dRT.WT > 0 & ind.p.wt,"Earlier",ifelse(dat_sk$dRT.WT < 0 & ind.p.wt,"Later","Unchanged"))
dat_sk$transition.state <- ifelse(dat_sk$dRT.dKO > 0 & dat_sk$sign,"Earlier",
                                  ifelse(dat_sk$dRT.dKO < 0 & dat_sk$sign,"Later","Unchanged"
                                  ))
table(dat_sk$wt.state)
table(dat_sk$transition.state)
dat_sk <- dat_sk[,4:5] %>% make_long(wt.state,transition.state)
dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

##dCD
set.seed(178)
pval.dRT <- apply(
  cbind(dat_loess$dCD.form.n1.RT.Loess-dat_loess$dCD.naive.n1.RT.Loess,
        dat_loess$dCD.form.n2.RT.Loess-dat_loess$dCD.naive.n2.RT.Loess,
        dat_loess$dCD.form.n3.RT.Loess-dat_loess$dCD.naive.n3.RT.Loess),1,function(row){
          ttest <- t.test(row,mu=0)
          pval <- ttest$p.value
          return(pval)
        }
)
pval.dRT.corr <- p.adjust(pval.dRT,method="fdr")

# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.dCD,-log10(pval.dRT.corr), n = 250)
# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
mypar()
image(density_est, col = mycols, xlab = "dRT dCD", ylab = "-log10(adj p-val)",xlim=c(-3.5,3.5),ylim=c(0,2),useRaster=T)
abline(h=-log10(0.05),lty=2)
legend("top",legend=paste0(round(100*(sum(pval.dRT.corr < 0.05)/nrow(dat_loess)),1),"% of genome"),cex=0.8,box.lty=0,bg="transparent")

#Volcano plot dKO vs KO naive
set.seed(178)
pval.dRT <- apply(
  cbind(dat_loess$dKO.naive.n1.RT.Loess-dat_loess$KO.naive.n1.RT.Loess,
        dat_loess$dKO.naive.n2.RT.Loess-dat_loess$KO.naive.n2.RT.Loess,
        dat_loess$dKO.naive.n3.RT.Loess-dat_loess$KO.naive.n3.RT.Loess),1,function(row){
          ttest <- t.test(row,mu=0)
          pval <- ttest$p.value
          return(pval)
        }
)
pval.dRT.corr <- p.adjust(pval.dRT,method="fdr")
# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.naive.dKOvsKO,-log10(pval.dRT.corr), n = 250)
# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
mypar()
image(density_est, col = mycols, xlab = "naive RT (dKO vs Ctr)", ylab = "-log10(adj p-val)",xlim=c(-3.5,3.5),ylim=c(0,2),useRaster=T)
abline(h=-log10(0.1),lty=2)
abline(v=-1,lty=2)
abline(v=1,lty=2)
legend("top",legend=paste0(round(100*(sum(pval.dRT.corr < 0.1)/nrow(dat_loess)),1),"% of genome"),cex=0.8,box.lty=0,bg="transparent")

#Volcano plot dCD vs KO naive
set.seed(178)
pval.dRT <- apply(
  cbind(dat_loess$dCD.naive.n1.RT.Loess-dat_loess$KO.naive.n1.RT.Loess,
        dat_loess$dCD.naive.n2.RT.Loess-dat_loess$KO.naive.n2.RT.Loess,
        dat_loess$dCD.naive.n3.RT.Loess-dat_loess$KO.naive.n3.RT.Loess),1,function(row){
          ttest <- t.test(row,mu=0)
          pval <- ttest$p.value
          return(pval)
        }
)
pval.dRT.corr <- p.adjust(pval.dRT,method="fdr")
# Calculate the 2D density using kde2d()
density_est <- kde2d(dat_loess$dRT.naive.dCDvsKO,-log10(pval.dRT.corr), n = 250)
# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
mypar()
image(density_est, col = mycols, xlab = "naive RT (dCD vs Ctr)", ylab = "-log10(adj p-val)",xlim=c(-3.5,3.5),ylim=c(0,2),useRaster=T)
abline(h=-log10(0.1),lty=2)
abline(v=-1,lty=2)
abline(v=1,lty=2)
legend("top",legend=paste0(round(100*(sum(pval.dRT.corr < 0.1)/nrow(dat_loess)),1),"% of genome"),cex=0.8,box.lty=0,bg="transparent")


#Sankey plots delta-RT (WT vs dCD)
ind <- pval.dRT.corr < 0.05  #define p-value cutoff
dat_sk <- dat[,c("dRT.WT","dRT.dCD")]
dat_sk$sign <- ind
dat_sk$wt.state <- ifelse(dat_sk$dRT.WT > 0 & ind.p.wt,"Earlier",ifelse(dat_sk$dRT.WT < 0 & ind.p.wt,"Later","Unchanged"))
dat_sk$transition.state <- ifelse(dat_sk$dRT.dCD > 0 & dat_sk$sign,"Earlier",
                                  ifelse(dat_sk$dRT.dCD < 0 & dat_sk$sign,"Later","Unchanged"
                                  ))
table(dat_sk$wt.state)
table(dat_sk$transition.state)
dat_sk <- dat_sk[,4:5] %>% make_long(wt.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#Sankey plots RT (naive), DKO vs WT
s <- 0.5 #define log2 fold change cutoff
##Simple
dat_sk <- dat[,c("AVGnaive.WT","AVGnaive.dKO","dRT.naive.dKOvsWT")]
dat_sk$wt.state <- ifelse(dat_sk$AVGnaive.WT > 0,"early","late")
dat_sk$transition.state <- ifelse(dat_sk$dRT.naive.dKOvsWT > s,"Earlier",
                                  ifelse(dat_sk$dRT.naive.dKOvsWT < -s,"Later","Unchanged"
                                                              ))

dat_sk <- dat_sk[,4:5] %>% make_long(wt.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#Sankey plots RT (naive), dCD vs WT
s <- 0.5 #define log2 fold change cutoff
##Simple
dat_sk <- dat[,c("AVGnaive.WT","AVGnaive.dCD","dRT.naive.dCDvsWT")]
dat_sk$wt.state <- ifelse(dat_sk$AVGnaive.WT > 0,"early","late")
dat_sk$transition.state <- ifelse(dat_sk$dRT.naive.dCDvsWT > s,"Earlier",
                                  ifelse(dat_sk$dRT.naive.dCDvsWT < -s,"Later","Unchanged"
                                  ))

dat_sk <- dat_sk[,4:5] %>% make_long(wt.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#Sankey plots delta-RT WT vs DKO
s <- 0.5 #define log2 fold change cutoff
dat_sk <- dat[,c("AVGnaive.WT","AVGformative.WT","dRT.WT","AVGnaive.dKO","AVGformative.dKO","dRT.dKO")]
dat_sk <- dat_sk[dat_sk$dRT.WT > s | dat_sk$dRT.WT < -s,]
dat_sk$wt.state <- ifelse(dat_sk$dRT.WT > 0,"Earlier","Later")
dat_sk$transition.state <- ifelse(dat_sk$dRT.WT > s & dat_sk$dRT.dKO > s,"Earlier",
                           ifelse(dat_sk$dRT.WT < -s & dat_sk$dRT.dKO < -s,"Later","Unchanged"
                                                       ))

dat_sk <- dat_sk[,7:8] %>% make_long(wt.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#Sankey plots delta-RT WT vs DKO
s <- 0.5 #define log2 fold change cutoff
dat_sk <- dat[,c("AVGnaive.WT","AVGformative.WT","dRT.WT","AVGnaive.dCD","AVGformative.dCD","dRT.dCD")]
dat_sk <- dat_sk[dat_sk$dRT.WT > s | dat_sk$dRT.WT < -s,]
dat_sk$wt.state <- ifelse(dat_sk$dRT.WT > 0,"Earlier","Later")
dat_sk$transition.state <- ifelse(dat_sk$dRT.WT > s & dat_sk$dRT.dCD > s,"Earlier",
                                  ifelse(dat_sk$dRT.WT < -s & dat_sk$dRT.dCD < -s,"Later","Unchanged"
                                  ))

dat_sk <- dat_sk[,7:8] %>% make_long(wt.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

#EOF