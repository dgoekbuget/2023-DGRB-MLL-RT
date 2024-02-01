#Analysis replication Timing naive vs formative
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
setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/mapping/bg-20kb-windows")
plotdir <- "/blellochlab/data1/deniz/analysis/mll-rt-paper/exports" #where to export plots
#setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/mapping/bg-50kb-windows")

#General settings
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

#Average replicates
colnames(dat)
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

#write.table(dat_loess,"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-20kb-bins-loess.bed", row.names = F, quote = F, sep ="\t")

#Circular binary segmentation
library(DNAcopy)
set.seed(178)
##Empiric determination of optimized segmentation paramters
input <- dat_loess[1:1000,]

mypar(5,6)
dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
for (i in seq(12,13,by=0.1)) {
  for(j in c(1e-15,1e-200)){
  seg.cna  = segment(dat.cna, nperm = 1000, alpha = j, undo.splits = "sdundo", 
                     undo.SD = i, verbose = 0)
  plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
       pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
       )
  legend("topright",legend=paste0("alpha=",j," undo.SD=",i),cex=0.8,box.lty=0,bg="transparent")
  }
}

mypar(5,6)
dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "delta.WT")
for (i in seq(12,13,by=0.1)) {
  for(j in c(1e-15,1e-200)){
    seg.cna  = segment(dat.cna, nperm = 1000, alpha = j, undo.splits = "sdundo", 
                       undo.SD = i, verbose = 0)
    plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
         pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
    )
    legend("topright",legend=paste0("alpha=",j," undo.SD=",i),cex=0.8,box.lty=0,bg="transparent")
  }
}
#Plot selected parameters
mypar(1,2)
dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
seg.cna  = segment(dat.cna, nperm = 1000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 12, verbose = 2)
plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
     pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
)
legend("topright",legend=paste0("alpha=",1e-15," undo.SD=",12),cex=0.8,box.lty=0,bg="transparent")

dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "delta.WT")
seg.cna  = segment(dat.cna, nperm = 1000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 12, verbose = 0)
plot(subset(seg.cna,chromlist = "chr1"), pch = 19,
     pt.cols = c("gray","gray"),  xmaploc = T, ylim = c(-5,5),
)
legend("topright",legend=paste0("alpha=",1e-15," undo.SD=",12),cex=0.8,box.lty=0,bg="transparent")


##Segmentation of naive and delta(naive-to-form)
input <- dat_loess

dat.cna  = CNA(input$AVGnaive.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "naive.WT")
seg.cna.naive  = segment(dat.cna, nperm = 10000, alpha = 1e-15, undo.splits = "sdundo", 
                   undo.SD = 12, verbose = 2)

dat.cna  = CNA(input$dRT.WT, input$chr, input$start, data.type = "logratio", 
               sampleid  = "dRT.WT")
seg.cna.dRT  = segment(dat.cna, nperm = 10000, alpha = 1e-15, undo.splits = "sdundo", 
                         undo.SD = 12, verbose = 2)

# output <- seg.cna.naive$output
# output <- output[,c(2,3,4,6)]
# write.table(output,"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-domains-naive.bed", row.names = F, quote = F, sep ="\t")
# output <- seg.cna.dRT$output
# output <- output[,c(2,3,4,6)]
# write.table(output,"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-delta-domains.bed", row.names = F, quote = F, sep ="\t")


#Fig.1 
##PCA WT naive + formative
library(rafalib)
set.seed(178)
rt <- dat_loess[,grepl("WT.naive.n|WT.form.n",colnames(dat_loess))]
rt <- rt[complete.cases(rt),]
y <- rt - rowMeans(rt)
s <- svd(y)

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100

##Plot variance explained
plot(pvar,ylab="Percent variability explained",cex=2,bg="gray",pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

##PC1 only stripchart
mycols <- rep(c('#8491b4ff','#00a087ff'),each=3)
set.seed(111)
stripchart(PC1, method = "jitter",pch=21,bg=mycols,cex=2,vertical=T,
           ylab="PC1")

##Heatmap RT WT naive + formative
input <- dat_loess[,grep("WT.naive.n|WT.form.n|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, cluster_columns = T,km=10,
        column_labels = rep(c("Naive","Formative"),each=3),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 0) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

#Fig.3 - PCA WT, CKO, DKO
set.seed(178)
head(dat_loess)
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

##Heatmap RT WT,CKO,DKO,dCD naive
set.seed(178)
head(dat_loess)
input <- dat_loess[,grep("AVGnaive.WT|AVGnaive.KO|AVGnaive.dKO|AVGnaive.dCD|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, cluster_columns = T,km=10,
        column_labels = c("WT","3KO","dKO","dCD"),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

##Heatmap RT CKO,DKO,dCD vs WT naive
set.seed(178)
input <- dat_loess[,grep("AVGnaive.WT|AVGnaive.KO|AVGnaive.dKO|AVGnaive.dCD|chr",colnames(dat_loess))]
input <- data.frame(chr=input$chr,
                    KOvsWT=input$AVGnaive.KO-input$AVGnaive.WT,
                    dKOvsWT=input$AVGnaive.dKO-input$AVGnaive.WT,
                    dCDvsWT=input$AVGnaive.dCD-input$AVGnaive.WT)
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, raster_by_magick = T, cluster_columns = T,km=5,
        column_labels = c("3KO","dKO","dCD"),
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

##Heatmap delta RT WT,3KO,dKO,dCD
input <- dat_loess[,grep("dRT.WT|dRT.KO|dRT.dKO|dRT.dCD|chr",colnames(dat_loess))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
df <- df[order(df$dRT.WT,decreasing = T),]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T,raster_by_magick = T,cluster_rows = F,cluster_columns = T,
        name = paste0("log2(RT) N=",nrow(df)), #title of legend
        col = col_fun,
        row_names_gp = gpar(fontsize = 0) # Text size for row names
)

##PCA WT, CKO, dKO, dCD (naive)
set.seed(178)
rt <- dat_loess[,grepl("naive.n",colnames(dat_loess))]
y <- rt - rowMeans(rt)
s <- svd(y)

label <- factor(colnames(y))

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100
mycols <- brewer.pal(12,"Paired")
plot(pvar,ylab="Percent variability explained",cex=2,bg=mycols,pch=21,
     xlab="PC")
PC1 <- s$v[,1]*s$d[1]
PC2 <- s$v[,2]*s$d[2]
PC3 <- s$v[,3]*s$d[3]
PC4 <- s$v[,4]*s$d[4]

plot(PC1,PC2,xlab = paste0("PC1 (",round(pvar[1]),"% var)"), ylab = paste0("PC2 (",round(pvar[2]),"% var)"),
     bg = ifelse(grepl("WT",label),"black",ifelse(grepl("dKO",label),"#CC3333",ifelse(grepl("dCD",label),"#996699","#999999"))),
     pch=ifelse(grepl("naive",label),21,23),cex = 2)

##Heatmap RT 3KO, dKO, dCD naive
input <- dat_loess[,grep("chr|KO.naive.n|dKO.naive.n|dCD.naive.n",colnames(dat_loess))]
input[,-1] <- input[,-1]-dat_loess$AVGnaive.WT
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, cluster_columns = T,km=10,
        column_labels = rep(c("3KO","dKO","dCD"),each=3),
        name = "log2(RT)", #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

##Fold change PCA
library(rafalib)
library(sva)

rt <- dat_loess[,grepl("dRT.",colnames(dat_loess))]
rt <- rt[complete.cases(rt),]
y <- rt - rowMeans(rt)
y <- y[,c(1,4,3,2)] #reorder colnames by genotype
s <- svd(y)

mypar(1,1)
pvar <- s$d^2/sum(s$d^2)*100
mycols <- brewer.pal(4,"Paired")
plot(pvar,ylab="Percent variability explained",cex=2,bg=mycols,pch=21,
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

mypar(1,1)
plot(PC1,PC3,xlab=paste0("PC1 (",round(pvar[1]),"% var)"), 
     ylab = paste0("PC3 (",round(pvar[3]),"% var)"),
     bg = geno_colors, pch=21,cex = 2)
legend("topleft",colnames(y),col=geno_colors,
       pch=16,box.lwd=1,cex = 1.5)

##Heatmap deltaRT 3KO, dKO, dCD naive to formative
input <- dat_loess[abs(dat_loess$dRT.WT) > 1,]
input <- input[,grep("chr|dRT.KO|dRT.dKO|dRT.dCD",colnames(input))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
colnames(df)

col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, cluster_columns = T,km=10,
        column_labels = rep(c("dCD","dKO","3KO"),each=1),
        name = "log2(RT)", #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

#head(dat_loess[order(rowMeans(cbind(dat_loess$dRT.dKO-dat_loess$dRT.KO,dat_loess$dRT.dCD-dat_loess$dRT.KO))),c(1:3)])


#Correlation plots of RT in WT vs all genotypes in naive
# mypar(4,4)
# avgRT <- colnames(dat_loess)[grepl("AVGnaive",colnames(dat_loess))]
# avgRT <- avgRT[c(1,4,3,2)]
# for (i in 1:length(avgRT)) {
#   x <- dat_loess[,avgRT[i]]
#   for (j in 1:length(avgRT)) {
#     y <- dat_loess[,avgRT[j]]
#     density_est <- kde2d(x, y, n = 200)
#     mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
#     image(density_est, col = mycols, xlab = gsub("AVGnaive.","",avgRT[i]), 
#           ylab = gsub("AVGnaive.","",avgRT[j]),xlim=c(-5,5),ylim=c(-5,5),main=paste0("rho=",round(cor(x,y,method="spearman"),3)))
#     abline(0,1)
#   }
# }

#Histogram of RT distribution in naive
##compute densities
density1 <- density(dat_loess[,"AVGnaive.WT"])
density2 <- density(dat_loess[,"AVGnaive.KO"])
density3 <- density(dat_loess[,"AVGnaive.dKO"])
density4 <- density(dat_loess[,"AVGnaive.dCD"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-5,5),xlab="RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[1], lwd = 3)
lines(density2, col = geno_colors[2], lwd = 3)
lines(density3, col = geno_colors[3], lwd = 3)
lines(density4, col = geno_colors[4], lwd = 3)
legend("topright", legend = c("WT", "KO", "dKO", "dCD"), col = geno_colors, lty = 1, lwd = 2)

#Histogram of RT distribution in naive
##compute densities
density1 <- density(dat_loess[,"AVGnaive.KO"]-dat_loess[,"AVGnaive.WT"])
density2 <- density(dat_loess[,"AVGnaive.dKO"]-dat_loess[,"AVGnaive.WT"])
density3 <- density(dat_loess[,"AVGnaive.dCD"]-dat_loess[,"AVGnaive.WT"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-5,5),xlab="RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[2], lwd = 3)
lines(density2, col = geno_colors[3], lwd = 3)
lines(density3, col = geno_colors[4], lwd = 3)
legend("topright", legend = c("WT", "KO", "dKO", "dCD"), col = geno_colors, lty = 1, lwd = 2)

#Histogram of RT distribution in naive-to-form transition
##compute densities
density1 <- density(dat_loess[,"dRT.KO"])
density2 <- density(dat_loess[,"dRT.dKO"])
density3 <- density(dat_loess[,"dRT.dCD"])

# Create an empty plot to initialize the plotting area
plot(density1, type = "n",xlim=c(-5,5),ylim=c(0,1),xlab="RT")

# Plot density estimates as lines
lines(density1, col = geno_colors[2], lwd = 3)
lines(density2, col = geno_colors[3], lwd = 3)
lines(density3, col = geno_colors[4], lwd = 3)
legend("topright", legend = c("KO", "dKO", "dCD"), col = geno_colors[-1], lty = 1, lwd = 2)

#Correlation between changes in dKO and changes in dCD during nf transition
x <- dat_loess$dRT.dCD-dat_loess$dRT.KO
y <- dat_loess$dRT.dKO-dat_loess$dRT.KO

# Calculate the 2D density using kde2d()
density_est <- kde2d(x, y, n = 500)

# Create a contour plot
contour(density_est, xlab = "dCD", ylab = "dKO",xlim=c(-2,2),ylim=c(-2,2))

# Create a color-coded density heatmap
mycols <- colorRampPalette(c("white", "blue", "red", "orange"))(100)
image(density_est, col = mycols, xlab = "dCD", ylab = "dKO",xlim=c(-2,2),ylim=c(-2,2))

#Auto-correlation QC
mypar(1,8)
acf(dat_loess$AVGnaive.WT,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGformative.WT,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGnaive.KO,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGformative.KO,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGnaive.dKO,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGformative.dKO,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGnaive.dCD,lag = 1000,na.action = na.pass)$acf[2]
acf(dat_loess$AVGformative.dCD,lag = 1000,na.action = na.pass)$acf[2]

#Correlation plot
library(corrplot)
cc <- cor(dat_loess[,c(4:27)],method = "spearman")
cc <- cor(dat[,c(4:27)],method = "spearman")
mypar()
colnames(cc) <- gsub(".20kb.RT.Loess","",colnames(cc))
rownames(cc) <- gsub(".20kb.RT.Loess","",colnames(cc))
corrplot(cc,tl.col = "black",col.lim=c(0.8,1),is.corr = F,order="original")

#Generate 2D domain heatmap
input <- dat_loess[dat_loess$chr=="chr1",c("AVGform.WT")]
input <- dat_loess[dat_loess$chr=="chr1",c("AVGform.WT")]-dat_loess[dat_loess$chr=="chr1",c("AVGnaive.WT")]

# mat <- matrix(1,nrow=length(input),ncol=length(input))
# for (i in 1:(length(input)-2)) {
#   for (j in 1:(length(input)-2)) {
#     #mat[i,j] <- dist(c(input[i],input[j]),method="euclidian")
#     #mat[i,j] <- cosine_similarity(input[i],input[j])
#     #mat[i,j] <- signed_euclidean_distance(input[i],input[j])
#     mat[i,j] <- cor(
#       x=c(input[i],input[i+1],input[i+2]),
#                     y=c(input[j],input[j+1],input[j+2])^2
#     )
#   }
# }

# Create logical matrices for positive and negative values
positive_matrix <- outer(input > 0, input > 0, FUN = "&")
negative_matrix <- outer(input < 0, input < 0, FUN = "&")

# Initialize the result matrix with zeros
mat <- matrix(0, nrow = length(input), ncol = length(input))

# Fill in the result matrix
mat[positive_matrix] <- 1
mat[negative_matrix] <- -1

# Set diagonal elements to 0
diag(mat) <- 0

mycols <- colorRampPalette(c("dodgerblue", "white", "orange"))(100)
Heatmap(mat,cluster_rows = F,cluster_columns = F, use_raster = T)


# Create a single vector
my_vector <- c(1, 2, 4, 8, 10)

# Convert the vector into a matrix with one row
my_matrix <- matrix(my_vector, nrow = 1)

# Calculate the distance matrix using the dist() function
distance_matrix <- dist(my_matrix)

# Print the distance matrix
print(distance_matrix)

head(dist(input[1:100]))
head(input)

input <- dat_loess[abs(dat_loess$dRT.WT) > 1,]
input <- input[,grep("chr|dRT.KO|dRT.dKO|dRT.dCD",colnames(input))]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
colnames(df)

col_fun = circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue", "white", "orange"))
Heatmap(df, use_raster = T, cluster_columns = T,km=10,
        column_labels = rep(c("dCD","dKO","3KO"),each=1),
        name = "log2(RT)", #title of legend
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

#Changes in the RT
##Percentage change analysis
##WT naive to formative
##Scale naive to baseline and propagate errors to formative
control <- cbind(dat$WT.naive.n1.20kb.RT.Loess,dat$WT.naive.n2.20kb.RT.Loess,dat$WT.naive.n3.20kb.RT.Loess)
treat <- cbind(dat$WT.form.n1.20kb.RT.Loess,dat$WT.form.n2.20kb.RT.Loess,dat$WT.form.n3.20kb.RT.Loess)
f <- matrix(1,nrow=nrow(control),ncol=ncol(control))/2^control
control_scaled <- log2(2^control*f)
treat_scaled <- log2(2^treat*f)

mypar(mar = c(6, 2.5, 1.6, 1.1))
mycols <- brewer.pal(4,"Paired")
d <- data.frame("Pct.0.5.earlier" = 100*apply(treat_scaled > 0.5,2,sum)/nrow(treat_scaled),
                "Pct.1.earlier" = 100*apply(treat_scaled > 1,2,sum)/nrow(treat_scaled),
                "Pct.0.5.later" = 100*apply(treat_scaled < -0.5,2,sum)/nrow(treat_scaled),
                "Pct.1.later" = 100*apply(treat_scaled < -1,2,sum)/nrow(treat_scaled)
                )
stripchart(x=d, method = "jitter",pch=17,col=mycols,cex=2,vertical=T,
           ylab=c("% genome changed"),las=2)


##WT, KO, dKO naive
##Scale naive to baseline and propagate errors to formative
control <- cbind(dat$KO.naive.n1,dat$KO.naive.n2,dat$KO.naive.n3)
treat1 <- cbind(dat$dKO.naive.n1,dat$dKO.naive.n2,dat$dKO.naive.n3)
f <- matrix(1,nrow=nrow(control),ncol=ncol(control))/2^control
control_scaled <- log2(2^control*f)
treat1_scaled <- log2(2^treat1*f)

mypar(mar = c(6, 2.5, 1.6, 1.1))
mycols <- brewer.pal(4,"Paired")
d <- data.frame("Pct.0.5.earlier" = 100*apply(treat1_scaled > 0.5,2,sum)/nrow(treat1_scaled),
                "Pct.1.earlier" = 100*apply(treat1_scaled > 1,2,sum)/nrow(treat1_scaled),
                "Pct.0.5.later" = 100*apply(treat1_scaled < -0.5,2,sum)/nrow(treat1_scaled),
                "Pct.1.later" = 100*apply(treat1_scaled < -1,2,sum)/nrow(treat1_scaled)
)
stripchart(x=d, method = "jitter",pch=17,col=mycols,cex=2,vertical=T,
           ylab=c("% genome changed"),las=2)

#write.table(dat,"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/deltaRT-50kbBins.bed", row.names = F, quote = F, sep ="\t")

#Sankey plots delta-RT (naive to form), WT
##Simple
s <- 0.5 #define log2 fold change cutoff
dat_sk <- dat[,c("AVGnaive.WT","AVGformative.WT","dRT.WT")]
dat_sk$naive.state <- ifelse(dat_sk$AVGnaive.WT > 0,"early","late")
dat_sk$transition.state <- ifelse(dat_sk$dRT.WT > s,"earlier",
                                  ifelse(dat_sk$dRT.WT < -s,"later","unchanged"
                                         ))
dat_sk <- dat_sk[,4:5] %>% make_long(naive.state,transition.state)

dat_sk %>% ggplot(aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.75, node.color = "black",flow.color = "black") +
  scale_fill_viridis_d(option = "B", alpha = 0.95) +
  theme_sankey(base_size = 16)

##More categories
s <- 0.5 #define log2 fold change cutoff
dat_sk <- dat[,c("AVGnaive.WT","AVGformative.WT","dRT.WT")]
dat_sk$naive.state <- ifelse(dat_sk$AVGnaive.WT > 0,"early","late")
dat_sk$transition.state <- ifelse(dat_sk$dRT.WT > s & dat_sk$AVGnaive.WT < 0 & dat_sk$AVGformative.WT > 0,"late-to-early",
                                  ifelse(dat_sk$dRT.WT < -s & dat_sk$AVGnaive.WT > 0 & dat_sk$AVGformative.WT < 0,"early-to-late",
                                         ifelse(dat_sk$dRT.WT > s & dat_sk$AVGnaive.WT < 0 & dat_sk$AVGformative.WT < 0,"later-to-late",
                                                ifelse(dat_sk$dRT.WT > s & dat_sk$AVGnaive.WT > 0 & dat_sk$AVGformative.WT > 0,"early-to-earlier",
                                                       ifelse(dat_sk$dRT.WT < -s & dat_sk$AVGnaive.WT < 0 & dat_sk$AVGformative.WT < 0,"late-to-later",
                                                              ifelse(dat_sk$dRT.WT < -s & dat_sk$AVGnaive.WT > 0 & dat_sk$AVGformative.WT > 0,"earlier-to-early","unchanged"
                                                              ))))))

dat_sk <- dat_sk[,4:5] %>% make_long(naive.state,transition.state)

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

#Heatmaps

#Figure 1 addition: Genome wide replication timing changes

#Create input and filter for bins that change more than s in WT from naive to form
s=2
input <- dat[abs(dat$dRT.WT) > s,c("chr", "WT.naive.n1",       "WT.naive.n2" ,      "WT.naive.n3", "KO.naive.n1",       "KO.naive.n2",      "KO.naive.n3",
                  "dKO.naive.n1",      "dKO.naive.n2",      "dKO.naive.n3", "dCD.naive.n1",      "dCD.naive.n2",      "dCD.naive.n3",
                  "WT.formative.n1",   "WT.formative.n2",   "WT.formative.n3", "KO.formative.n1",   "KO.formative.n2",   "KO.formative.n3",        
                  "dKO.formative.n1",  "dKO.formative.n2",  "dKO.formative.n3",  "dCD.formative.n1",  "dCD.formative.n2",  "dCD.formative.n3")]
region <- input$chr %in% c(paste0("chr",1:19),"chrX")
chr <- factor(input$chr[region],levels = c(paste0("chr",1:19),"chrX"))
df <- input[region,-1]
Heatmap(df, use_raster = T, cluster_columns = T,km=5,
        #column_labels = c("Naive","Formative"),
        name = "log2(RT)", #title of legend
        show_row_dend = T,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
    Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))

##Z score version of heatmap
df_scaled <- t(scale(t(df)))
Heatmap(df_scaled, use_raster = T, cluster_columns = T,km=5,
        #column_labels = c("Naive","Formative"),
        name = "Z-score(RT)", #title of legend
        show_row_dend = T,
        row_names_gp = gpar(fontsize = 1) # Text size for row names
) +
  Heatmap(factor(chr), name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(factor(chr)))))



# #Custom segmentation based on sign changes for naive RT
# signal <- dat_loess[,c("AVGnaive.WT")]
# sign_changes <- which(diff(sign(signal)) != 0)
# bounds <- dat_loess[sign_changes,c(1:3)]
# doms <- data.frame()
# for (i in unique(bounds$chr)) {
#   x <- bounds[bounds$chr==i,]
#   for (j in 1:(nrow(x)-1)) {
#     doms <- rbind(doms,c(x[j,]$chr,x[j,]$stop+1,x[j+1,]$stop,abs(x[j+1,]$stop-x[j,]$stop+1)))
#   }
# }
# colnames(doms) <- c("chr","start","stop","dom.size")

# 
# #Calculate median RT per domain
# columnsToSummarize <- -c(1:3)
# ovlp <- findOverlaps(GRanges(doms),GRanges(dat_loess),ignore.strand=T)
# sumVals  <- apply(dat_loess[subjectHits(ovlp),columnsToSummarize], 2, function(x) tapply(x, factor(queryHits(ovlp)), median,na.rm=T))
# doms[,colnames(sumVals)] <- NA
# doms[unique(queryHits(ovlp)),5:ncol(doms)] <- sumVals
# colnames(doms)
# #write.table(doms[,c("chr","start","stop","dom.size","AVGnaive.WT")],"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-domains-naive.bed", row.names = F, quote = F, sep ="\t")
# 
# #Custom segmentation based on sign changes for delta RT
# signal <- dat_loess[,c("dRT.WT")]
# sign_changes <- which(diff(sign(signal)) != 0)
# bounds <- dat_loess[sign_changes,c(1:3)]
# deltadoms <- data.frame()
# for (i in unique(bounds$chr)) {
#   x <- bounds[bounds$chr==i,]
#   for (j in 1:(nrow(x)-1)) {
#     deltadoms <- rbind(deltadoms,c(x[j,]$chr,x[j,]$stop+1,x[j+1,]$stop,abs(x[j+1,]$stop-x[j,]$stop+1)))
#   }
# }
# head(deltadoms)
# colnames(deltadoms) <- c("chr","start","stop","dom.size")
# deltadoms <- deltadoms[deltadoms$dom.size > 100001,]

#Calculate median RT per domain
# columnsToSummarize <- -c(1:3)
# ovlp <- findOverlaps(GRanges(deltadoms),GRanges(dat_loess),ignore.strand=T)
# sumVals  <- apply(dat_loess[subjectHits(ovlp),columnsToSummarize], 2, function(x) tapply(x, factor(queryHits(ovlp)), median,na.rm=T))
# deltadoms[,colnames(sumVals)] <- NA
# deltadoms[unique(queryHits(ovlp)),5:ncol(deltadoms)] <- sumVals
# colnames(deltadoms)
#write.table(deltadoms[,c("chr","start","stop","dom.size","dRT.WT")],"/blellochlab/data1/deniz/analysis/mll-rt-paper/exports/RT-delta-domains.bed", row.names = F, quote = F, sep ="\t")

