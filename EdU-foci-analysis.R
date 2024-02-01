library(ggplot2)
library(readr)
library(stringr)
library(viridis)

#Read data
setwd("/Volumes/External_DG1/Box\ Sync/postdoc/data/microscopy/phenix/20230228-cko-dko-10m-20m-edu-dapi__2023-02-28T16_21_23-Measurement\ 1/Evaluation1") #Path to Phenix output data
res <- read_delim("Objects_Population\ -\ Singlets.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE, 
                  skip = 9)
res <- data.frame(res)

spots <- read_delim("Objects_Population - Spots Alexa 647.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE, 
                  skip = 9)
spots <- data.frame(spots)

#Generate sample info data
unique(res$Column)
res$Treatment <- res$Column
res$Treatment[res$Column %in% c(1,3)] <- "MLL3KO"
res$Treatment[res$Column %in% c(2,4)] <- "MLL3/4dKO"
res$Treatment <- factor(res$Treatment)
res$Treatment <- relevel(res$Treatment,"MLL3KO")

spots$Treatment <- spots$Column
spots$Treatment[spots$Column %in% c(1,3)] <- "MLL3KO"
spots$Treatment[spots$Column %in% c(2,4)] <- "MLL3/4dKO"
spots$Treatment <- factor(spots$Treatment)
spots$Treatment <- relevel(spots$Treatment,"MLL3KO")

res$time <- res$Column
res$time[res$Column %in% c(1,2)] <- "10min"
res$time[res$Column %in% c(3,4)] <- "20min"
res$time <- factor(res$time)

spots$time <- spots$Column
spots$time[spots$Column %in% c(1,2)] <- "10min"
spots$time[spots$Column %in% c(3,4)] <- "20min"
spots$time <- factor(spots$time)

#Load packages
library(ggplot2)
library(hexbin)
library(LSD)
library(lattice)
library(dplyr)
library(tidyr)
library(gridExtra)

#what kind of intensity? ____________
#define column number to pull intensity for dapi and edu
colnames(res)

dna.sum <- 33
log.edu.sum <- 43
spots.1 <- 40
spots.2 <- 41

#subset data.frame by intensity values and convert to long format
dat <- res[,c(16,dna.sum,log.edu.sum,spots.1,spots.2,45:46)]
colnames(dat) <- c("area","dna_sum","log_edu_sum","spot.1","spot.2","treatment","time")

#Downsample each timepoint and group before EdU plot
(ctr_size <- nrow(dat[dat$treatment == "MLL3KO" & dat$time == "10min",]))
(dko_size <- nrow(dat[dat$treatment == "MLL3/4dKO" & dat$time == "10min",]))
ind <- sample(which(dat$treatment == "MLL3/4dKO" & dat$time == "10min"),dko_size-ctr_size)

datDS <- dat[-ind,]

(ctr_size <- nrow(datDS[datDS$treatment == "MLL3KO" & datDS$time == "20min",]))
(dko_size <- nrow(datDS[datDS$treatment == "MLL3/4dKO" & datDS$time == "20min",]))
ind2 <- sample(which(datDS$treatment == "MLL3/4dKO" & datDS$time == "20min"),dko_size-ctr_size)

datDS <- datDS[-ind2,]

library(MASS)
library(ggplot2)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
datDS$density <- 0

for (i in unique(datDS$treatment)){
  for(k in unique(datDS$time)){
      sub <- datDS[datDS$treatment==i & datDS$time==k & datDS$area < 500,]
      datDS$density[datDS$treatment==i & datDS$time==k & datDS$area < 500] <- get_density(sub$dna_sum, sub$log_edu_sum,n=500)
  }
}

#Gate Sub S phase populations
(quants <- quantile(datDS$dna_sum[datDS$log_edu_sum > 5.2 & datDS$dna_sum < 800000 & datDS$dna_sum > 250000]))
dat$subS <- NA
dat$subS[dat$log_edu_sum > 5 & dat$dna_sum >= quants[1] & dat$dna_sum < quants[2]] <- "S1"
dat$subS[dat$log_edu_sum > 5.2 & dat$dna_sum >= quants[2] & dat$dna_sum < quants[3]] <- "S2"
dat$subS[dat$log_edu_sum > 5.3 & dat$dna_sum >= quants[3] & dat$dna_sum < quants[4]] <- "S3"
dat$subS[dat$log_edu_sum > 5.3 & dat$dna_sum >= quants[4] & dat$dna_sum < quants[5]] <- "S4"

#Filter dead and >2N cells
keep <- dat$dna_sum > 250000 & dat$dna_sum < 800000
dat <- dat[keep,]

#Color-code selected sub S population
p1 <- tibble(dat) %>% group_by(treatment,subS) %>% sample_n(1000) %>% ggplot() + 
  geom_point(aes(dna_sum, log_edu_sum, color = subS),alpha=0.2) + 
  scale_color_manual(values = viridis(5)) + 
  xlim(0,1000000) +
  theme_minimal() +
  xlab("DNA Content") +
  facet_wrap(~time) +
  ylab("log(EdU Sum)")
p1

#Plot distribution of spots
#Gate Sub S phase populations (based on above DNA content distribution of Edu+ cells)
res$subS <- NA
res$subS[res$Singlets...logEdUSum > 5 & res$Singlets...DNA.Sum >= quants[1] & res$Singlets...DNA.Sum < quants[2]] <- "S1"
res$subS[res$Singlets...logEdUSum > 5.2 & res$Singlets...DNA.Sum >= quants[2] & res$Singlets...DNA.Sum < quants[3]] <- "S2"
res$subS[res$Singlets...logEdUSum > 5.3 & res$Singlets...DNA.Sum >= quants[3] & res$Singlets...DNA.Sum < quants[4]] <- "S3"
res$subS[res$Singlets...logEdUSum > 5.3 & res$Singlets...DNA.Sum >= quants[4] & res$Singlets...DNA.Sum < quants[5]] <- "S4"

#Filter dead and >2N cells
keep2 <- res$Singlets...DNA.Sum > 250000 & res$Singlets...DNA.Sum < 800000
res2 <- res[keep2,]

#Merge with spot data
library(digest)
nuc_meta <- res2[res2$time=="10min",c(1:5,9,10,16,47)]
spots_meta <- spots[res$time=="10min",c(1:4,23,1:22,25,26)]

# concatenate the values of the columns into a string
id_string <- apply(nuc_meta[,1:5], 1, function(x) paste(x, collapse = "_"))
id <- sapply(id_string,function(x) digest(x, algo = "md5"))
nuc_meta$id <- id

id_string <- apply(spots_meta[,1:5], 1, function(x) paste(x, collapse = "_"))
id <- sapply(id_string,function(x) digest(x, algo = "md5"))
spots_meta$id <- id

#merge
spots_merged <- merge(spots_meta,nuc_meta,by="id")

#Calculate eucledian distance for each focus
spots_merged$eucl.dist <- sqrt(
  (spots_merged$Position.X..µm..x-spots_merged$Position.X..µm..y)^2+
    (spots_merged$Position.Y..µm..x-spots_merged$Position.Y..µm..y)^2
)

#Plot eucledian distance histogram
library(plyr)
mu <- spots_merged %>% filter(subS != "<NA>") %>% ddply(c("Treatment","subS"), summarise, grp.median=median(eucl.dist/2*sqrt(Singlets...Nucleus.Area..µm../pi)))
mu$Treatment <- relevel(mu$Treatment,ref="MLL3KO")

#Normalized relative distance from nucleus center
p1 <-spots_merged %>% filter(subS != "<NA>") %>% 
  ggplot(aes(x=eucl.dist/2*sqrt(Singlets...Nucleus.Area..µm../pi),fill=Treatment)) + 
  geom_density(alpha=.4) +
  geom_vline(dat=mu,aes(xintercept=grp.median,color=Treatment), linetype="dashed", size=.75) +
  scale_fill_manual(values=c("#999999", "#FF0000")) +
  facet_wrap(~subS,ncol=1) +
  xlim(0,40) +
  xlab("Relative distance from nucleus center") +
  theme_minimal()
p1

#Test significance
wilcox.test(spots_merged$eucl.dist[spots_merged$subS == "S1" & spots_merged$Treatment =="MLL3KO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S1" & spots_merged$Treatment =="MLL3KO"]/pi),
            spots_merged$eucl.dist[spots_merged$subS == "S1" & spots_merged$Treatment =="MLL3/4dKO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S1" & spots_merged$Treatment =="MLL3/4dKO"]/pi))
wilcox.test(spots_merged$eucl.dist[spots_merged$subS == "S2" & spots_merged$Treatment =="MLL3KO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S2" & spots_merged$Treatment =="MLL3KO"]/pi),
            spots_merged$eucl.dist[spots_merged$subS == "S2" & spots_merged$Treatment =="MLL3/4dKO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S2" & spots_merged$Treatment =="MLL3/4dKO"]/pi))
wilcox.test(spots_merged$eucl.dist[spots_merged$subS == "S3" & spots_merged$Treatment =="MLL3KO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S3" & spots_merged$Treatment =="MLL3KO"]/pi),
            spots_merged$eucl.dist[spots_merged$subS == "S3" & spots_merged$Treatment =="MLL3/4dKO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S3" & spots_merged$Treatment =="MLL3/4dKO"]/pi))
wilcox.test(spots_merged$eucl.dist[spots_merged$subS == "S4" & spots_merged$Treatment =="MLL3KO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S4" & spots_merged$Treatment =="MLL3KO"]/pi),
            spots_merged$eucl.dist[spots_merged$subS == "S4" & spots_merged$Treatment =="MLL3/4dKO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$subS == "S4" & spots_merged$Treatment =="MLL3/4dKO"]/pi))

#sample sizes
length(spots_merged$eucl.dist[spots_merged$Treatment =="MLL3KO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$Treatment =="MLL3KO"]/pi))
length(spots_merged$eucl.dist[spots_merged$Treatment =="MLL3/4dKO"]/2*sqrt(spots_merged$Singlets...Nucleus.Area..µm..[spots_merged$Treatment =="MLL3/4dKO"]/pi))

#EOF