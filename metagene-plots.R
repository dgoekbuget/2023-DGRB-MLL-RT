#This script takes tab separated output from plotProfile function of deepTools2 package as input

#Load packages
library(readr)
library(tidyr)
library(ggplot2)
library(gridExtra)

#Load and reformat input data
setwd("/Volumes/External_DG1/Box\ Sync/postdoc/data/bioinformatics/analysis/c4/analysis-k4m1fc1-sites/100kb-flanks")
infile <- dir(pattern="*.tab")[1]
dat <- read_table(infile,skip = 1)
dat <- data.frame(dat)
colnames(dat) <- c("sample","region",paste0("bin",1:(ncol(dat)-2)))
dat[,-c(1:2)] <- as.matrix(dat[,-c(1:2)])

#Format to long format
datL <- dat[,1:202] %>% pivot_longer(!c(sample,region),names_to="bins",values_to="counts")
datL <- data.frame(datL)
datL$counts <- as.numeric(datL$counts)

#Format sample names
samples <- unique(dat$sample)
regions <- unique(dat$region)
samples3 <- paste0(rep(samples,each=2),regions)
samples2 <- c("OK-seq","MLL3/4 ChIP","H3K4me1 CR CKO","H3K4me1 CR DKO","NAIL-seq","SNS-seq","H3.3 CT CKO",
              "H3.3 CT DKO","MCM2 CKO","pMCM2 CKO","MCM2 DKO","pMCM2 DKO")
samples <- data.frame("index"=samples3,"samples"=rep(samples,each=2),"region"=regions,"title"=rep(samples2,each=2))

#Perform baseline correction
##Generate vector of background values per sample based on average signal of peripheral bins
backgr <- sapply(1:nrow(samples), function(x){
  left_end <- datL[datL$sample == samples$samples[x] & 
                datL$region == samples$region[x] &
                datL$bins %in% paste0("bin",1),]$counts
  right_end <- datL[datL$sample == samples$samples[x] & 
                     datL$region == samples$region[x] &
                     datL$bins %in% paste0("bin",100),]$counts
  mean(c(left_end,right_end))
})

##Don't correct background for OK-seq data
backgr[c(1:4,9:12)] <- c(rep(0,8))

##Define y axis limits
ylim_dat2 <- data.frame(lower=rep(c(-0.04, 0, -0.2, -0.2, 0.5, 2.5, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2),each=2),
                       upper=rep(c(0.04, 0.04, 0.1, 0.1, 1, 5, 0.1, 0.1, 0.1 ,0.1, 0.1, 0.1),each=2))

#Subtract sample specific background from each bin of respective sample
plotlist <- list()
for (i in 1:nrow(samples)) {
  b <- backgr[i]
  dat_filtered <- datL[datL$sample == samples$samples[i] & 
                         datL$region == samples$region[i] &
                         datL$bins %in% paste0("bin",21:80),]
  dat_filtered$counts <- dat_filtered$counts - b
  plotlist[[i]] <- dat_filtered %>% ggplot(aes(1:60,counts)) + 
    ylim(as.numeric(ylim_dat2[i,])) + ggtitle(samples$title[i]) +
    geom_smooth(color="darkgray",span=0.2) + theme_classic() + theme(aspect.ratio = 0.5)
}
do.call("grid.arrange",c(list(ncol=4),plotlist))

#EOF