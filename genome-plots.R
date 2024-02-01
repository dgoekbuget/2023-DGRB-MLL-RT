library(trackplot)
library(GenomicRanges)
library(viridis)
library(RColorBrewer)
path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(path, "/c4/home/gdeniz/bwtool", sep = ":"))

#Define paths to bw files
folder="/blellochlab/data1/deniz/analysis/mll-rt-paper/mapping/bg-20kb-windows/bw"

naive=file.path(folder,"WT.naive.mean.RT.Loess.bw")
form=file.path(folder,"WT.form.mean.RT.Loess.bw")

mll3ko_n_merge=file.path(folder,"KO.naive.mean.RT.Loess.bw")
mll34dko_n_merge=file.path(folder,"dKO.naive.mean.RT.Loess.bw")
mll34dcd_n_merge=file.path(folder,"dCD.naive.mean.RT.Loess.bw")

mll3ko_f_merge=file.path(folder,"KO.form.mean.RT.Loess.bw")
mll34dko_f_merge=file.path(folder,"dKO.form.mean.RT.Loess.bw")
mll34dcd_f_merge=file.path(folder,"dCD.form.mean.RT.Loess.bw")

##Epigenetics bws
epifiles <- dir("/blellochlab/data1/deniz/ryan-naive-to-form/bigwigs",pattern="*.b*",full.names = T)
#pol2n=file.path("/blellochlab/data1/deniz/go-data/marks2012/naive/bw/SRR314986.sorted.dedup.CPM.bw")

#Define global parameters
gen <- "mm10"

#Figure 1 WT track at Fgf locus chr5:95,798,300-100,997,281
#Bw to plot
(bigWigs <- c(naive,form))

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10",sample_names = c("Naive","Formative"))

#Extract the signal for loci of interest
#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
loc = track_extract(colData = bigWigs, gene="Foxp1",padding=2*10^6,binsize = 10000,build=gen)
locus = attr(loc$data,"meta")$loci
t = track_extract(colData = bigWigs, loci=locus,binsize = 10000,build=gen)
#Define colors and plot
mycols <- c('#8491b4ff','#00a087ff')
track_plot(summary_list = t,groupAutoScale = T,col=mycols,gene_track_height = 4,
           track_names_to_left = F,show_axis = T, draw_gene_track = T)

####################################################################################################
####################################################################################################
####################################################################################################
#Figure 2 WT track at Fgf5, Klf4, Pou5f1
#Bw to plot
toPlot <- c("v65","H3K27a","H33","K4m1","RNApol","K4m3")

bigWigs <- c(naive,form)
for (i in toPlot) {
  bigWigs <- c(bigWigs,epifiles[grep(i,epifiles)])
}

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10",
                       sample_names = 1:12
                       )
bigWigs <- bigWigs[c(1,2,4,3,6,5,8,7,10,9,12,11),]
bigWigs$condition <- rep(c("RT","H3K27ac","H3.3","H3K4me1","P-Pol II","H3K4me3"),each=2)
bigWigs$bw_sample_names <- paste0(rep(c("RT","H3K27ac","H3.3","H3K4me1","P-Pol II","H3K4me3"),each=2),"_",rep(c("N","F")))

#Extract the signal for loci of interest
#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
loc1 <- track_extract(colData = bigWigs, gene="Fgf5",padding=3*10^5,binsize = 1000,build=gen,nthreads=12)
loc2 <- track_extract(colData = bigWigs, gene="Klf4",padding=2.5*10^5,binsize = 1000,build=gen,nthreads=12)
loc3 <- track_extract(colData = bigWigs, gene="Pou5f1",padding=2.5*10^5,binsize = 1000,build=gen,nthreads=12)
locus1 = attr(loc1$data,"meta")$loci
locus2 = attr(loc2$data,"meta")$loci
locus3 = attr(loc3$data,"meta")$loci

t1 <- track_extract(colData = bigWigs, loci=locus1,binsize = 1000,build=gen,nthreads=12)
t2 <- track_extract(colData = bigWigs, loci=locus2,binsize = 1000,build=gen,nthreads=12)
t3 <- track_extract(colData = bigWigs, loci=locus3,binsize = 1000,build=gen,nthreads=12)

#Define colors and plot
mycols <- rep(brewer.pal(6,"Paired"),each=2)
track_plot(summary_list = t1,col=mycols,show_ideogram = F,
           track_names_to_left = F,show_axis = T,gene_track_height = 4,
           y_max = c(2,2,2,2,17,17,5,5,3,3,8,8),
           y_min = c(-4,-4,rep(0,10))
           )

track_plot(summary_list = t2,col=mycols,show_ideogram = F,
           track_names_to_left = F,show_axis = T,gene_track_height = 4,
           y_max = c(3,3,11,11,5,5,4,4,8,8,9,9),
           y_min = c(-3,-3,rep(0,10))
)

track_plot(summary_list = t3,col=mycols,show_ideogram = F,
           track_names_to_left = F,show_axis = T,gene_track_height = 4,
           y_max = c(4,4,15,15,45,45,18,18,3,3,16,16),
           y_min = c(0,0,rep(0,10))
)


####################################################################################################
####################################################################################################
####################################################################################################
#Figure 5 Naive tracks affected region WT, 3KO, dKO
bigWigs <- c(naive, mll3ko_n_merge,mll34dko_n_merge)

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10",
                       sample_names = c("WT","MLL3KO","MLL3/4dKO")
)

#Extract the signal for loci of interest
#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
t1 <- track_extract(colData = bigWigs, loci="chr6:111000001-112000001",binsize = 10000,build=gen,nthreads=12)
t2 <- track_extract(colData = bigWigs, loci="chr8:96500001-98000001",binsize = 10000,build=gen,nthreads=12)

#Define colors and plot
mycols <- c("black","#999999","#CC3333","#996699")
track_plot(summary_list = t1,col=mycols,
           track_names_to_left = F,show_axis = T,
           y_max = c(2,2,2),
           y_min = c(-5,-5,-5)
)

track_plot(summary_list = t2,col=mycols,
           track_names_to_left = F,show_axis = T,
           y_max = c(2,2,2),
           y_min = c(-5,-5,-5)
)


####################################################################################################
####################################################################################################
####################################################################################################
#Figure 5 Naive tracks affected region 3KO, dKO, dCD
bigWigs <- c(mll3ko_n_merge,mll34dko_n_merge,mll34dcd_n_merge)

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10",
                       sample_names = c("MLL3KO","MLL3/4dKO","MLL3/4dCD")
)

#Extract the signal for loci of interest
#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
t1 <- track_extract(colData = bigWigs, loci="chr6:105000001-11250001",binsize = 10000,build=gen,nthreads=12)
t2 <- track_extract(colData = bigWigs, loci="chr8:96500001-98000001",binsize = 10000,build=gen,nthreads=12)

#Define colors and plot
mycols <- c("#999999","#CC3333","#996699")
track_plot(summary_list = t1,col=mycols,
           track_names_to_left = F,show_axis = T,
           y_max = c(2,2,2),
           y_min = c(-5,-5,-5)
)

track_plot(summary_list = t2,col=mycols,
           track_names_to_left = F,show_axis = T,
           y_max = c(2,2,2),
           y_min = c(-5,-5,-5)
)
