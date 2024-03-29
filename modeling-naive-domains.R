#This script is used for machine learning of features predictive of naive RT within naive RT segments

#Load packages
library(readr)
library(rafalib)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(caret)
library(glmnet)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(MASS)
library(ComplexHeatmap)


#Define location of background corrected peaks signal summed by RT segment for each chromatin feature
setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/ml/domains-naive")
(files <-dir(pattern="backgroundCorSum.bed"))

#Load signal data and merge into data.frame
for (i in 1:length(files)) {
  if(i == 1){
    rt <- read_delim(files[i],"\t", escape_double = FALSE, trim_ws = TRUE,col_names=F)[,1:6]
    colnames(rt) <- c("chr","start","end","RT",
                      gsub(".RTdomain.backgroundCorSum.bed",".N1",files[i]),
                      gsub(".RTdomain.backgroundCorSum.bed",".N2",files[i]))
  } else {
    rt <- cbind(rt,read_delim(files[i],"\t", escape_double = FALSE, trim_ws = TRUE,col_names=F)[,5:6])
    colnames(rt)[(ncol(rt)-1):ncol(rt)] <- c(gsub(".RTdomain.backgroundCorSum.bed",".N1",files[i]),gsub(".RTdomain.backgroundCorSum.bed",".N2",files[i]))
  }
}

#Create binary RT categories
rt$timing <- ifelse(rt$RT > 0, "E",ifelse(rt$RT < 0, "L",ifelse(rt$RT == 0,NA,NA)))

#Separate counts from meta
meta <- rt[,1:4]
counts <- rt[,-c(1:4)]
counts <- cbind(counts$timing,log2(counts[,-ncol(counts)]/(abs(meta$end-meta$start))+0.01))
colnames(counts)[1] <- "timing"

counts$timing <- factor(counts$timing)
dat <- counts[complete.cases(counts),]
meta <- meta[complete.cases(counts),]

#Simplify column names
colnames(dat) <- c("timing",paste0(rep(c("H2AX","H2AZ","H2BK5ac","H3.3","H3K14ac","H3K27ac",
                                         "H3K27me3","H3K36me1","H3K36me2","H3K4me2","H3K4me3",
                                         "H3K9ac","H3K9me2","H3K9me3","H4K14ac","H4K20me1","H4K8ac","H2AK119Ub",
                                         "H3K4me1","PPolII","Rad21"),each=2),c(".N1",".N2")))
dat <- dat[,c("timing",paste0(rep(c("PPolII","Rad21","H3K4me1","H3K4me2","H3K4me3",
         "H3K27ac","H3K36me1","H3K36me2","H2BK5ac","H3K14ac",
         "H2AX","H2AZ","H3.3", "H3K9ac", "H4K8ac","H4K14ac",
         "H4K20me1","H3K27me3","H3K9me2","H3K9me3","H2AK119Ub"),each=2),c(".N1",".N2")))
         ]

#Calculate average chromatin feature signal from replicates
mean_dat <- data.frame(dat$timing)
for (i in seq(1,42,by=2)) {
  m <- data.frame(rowMeans(dat[,c(i+1,i+2)]))
  colnames(m) <- gsub(".N1","",colnames(dat)[i+1])
  mean_dat <- cbind(mean_dat,m)
}

#Define color vectors
color21 <- colorRampPalette(c("#336699","#FF9966"))(100)

#Machine learning using Elastic Net Regression
set.seed(178)
X <- mean_dat[,-1]
X$RT <- meta$RT

control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = T)
fit <- train(RT ~ .,
             data = X,
             method = "glmnet",
             preProcess = c("center", "scale"),
             tuneLength = 25,
             trControl = control)

##Plot Observed over predicted
predictions <- extractPrediction(list(fit))

mypar()
plot(predictions$pred,predictions$obs,xlab = "Predicted",
     ylab ="Observed",main=paste0("RT naive"," rho=",round(cor(predictions$pred,predictions$obs,method="spearman"),3)),
     pch=21,bg="gray",cex=1,lwd=1.5,ylim=c(-5,5),xlim=c(-5,5))
abline(lm(predictions$obs~predictions$pred),col="dodgerblue",lty="dashed",lwd=3)

##Plot parameter weights (relative feature importance)
# Define a range of colors for the gradient
coefs <- coef(fit$finalModel,fit$bestTune$lambda)
rafalib::mypar(1,1,mar = c(2.5, 5, 1.6, 1.1))
barplot(coefs[-1][order(coefs[-1])],horiz=T,las=1,xlab="Parameter weight",col = "gray",
        cex.names=0.8,names=rownames(coefs)[-1][order(coefs[-1])],xlim=c(-0.5,0.5))

#Perform pairwise regression using linear regression or glm
residual_dispersion_results <- matrix(data=NA,nrow=21,ncol=21)

# Function to calculate residual dispersion for a glm model
calculate_residual_dispersion <- function(model) {
  residuals <- residuals(model, type = "deviance")
  degrees <- df.residual(model)
  dispersion <- sum(residuals^2) / degrees
  return(dispersion)
}

# Loop over all combinations of two columns
for (i in 1:(ncol(X) - 2)) {
  for (j in (i + 1):(ncol(X)-1)) {
    
    # Select the columns for the model
    columns <- c(i, j)
    
    # Fit the glm model
    model <- glm(RT ~ ., data = X[,c(columns,22)])
    adj.r2 <- summary(lm(RT ~ ., data = X[,c(columns,22)]))$adj.r.squared
    # Calculate residual dispersion
    dispersion <- calculate_residual_dispersion(model)
    
    #Save pairwise dispersion to matrix
    residual_dispersion_results[i,j] <- adj.r2
    residual_dispersion_results[j,i] <- adj.r2
  }
}
colnames(residual_dispersion_results) <- colnames(X)[-22]
rownames(residual_dispersion_results) <- colnames(X)[-22]

# Loop over each individual column
for (i in 1:(ncol(X)-1)) {
  
  # Fit the glm model
  model <- glm(RT ~ ., data = X[,c(i,22)])
  adj.r2 <- summary(lm(RT ~ ., data = X[,c(i,22)]))$adj.r.squared
  
  # Calculate residual dispersion
  dispersion <- calculate_residual_dispersion(model)
  
  # Save the results to diagonal of matrix
  diag(residual_dispersion_results)[i] <- adj.r2
}

#col_fun = circlize::colorRamp2(c(0,1,3.5), c("white", "blue", "white"))
col_fun = circlize::colorRamp2(c(0,0.7), c("white", "blue"))
Heatmap(residual_dispersion_results,col=col_fun,name="Adj. R2")

##Perform correlation analysis of selected features
set.seed(178)
selected <- rownames(coefs)[-1][order(abs(coefs[-1]),decreasing = T)][1:10]
X <- mean_dat[,-1][,selected]
X$RT <- meta$RT

##Separate signal of selected marks by kmeans into no signal - background/low - positive signal
k=2 #determine 2 clusters
dat_clustered <- X
for(i in 1:(ncol(X)-1)){
  input <- X[,c(i,ncol(X))]
  ind <- input[,1] <= min(X)
  input <- input[!ind,]
  dist_matrix <- dist(input,method="manhattan")  # Calculate the distance matrix
  hierarchical_result <- hclust(dist_matrix, method = "ward.D")
  cluster_assignment <- cutree(hierarchical_result, k = k)
  cluster_assignment[cluster_assignment == 1] <- "salmon"
  cluster_assignment[cluster_assignment == 2] <- "dodgerblue"
  dat_clustered$new <- "gray"
  dat_clustered$new[!ind] <- cluster_assignment
  colnames(dat_clustered)[ncol(dat_clustered)] <- paste0(colnames(X)[i],".cluster")
}
##Rename clusters as no signal, low, high
for (i in c(1:10)) {
  if(i %in% c(1:5,8,10)){
    dat_clustered[,11+i] <- gsub("dodgerblue","high",dat_clustered[,11+i])
    dat_clustered[,11+i] <- gsub("salmon","low",dat_clustered[,11+i])
    dat_clustered[,11+i] <- gsub("gray","negative",dat_clustered[,11+i]) 
  }else{
    dat_clustered[,11+i] <- gsub("salmon","high",dat_clustered[,11+i])
    dat_clustered[,11+i] <- gsub("dodgerblue","low",dat_clustered[,11+i])
    dat_clustered[,11+i] <- gsub("gray","negative",dat_clustered[,11+i]) 
  }
}

#Highlight clusters in scatterplot
mypar(5,2)
for(i in 1:(ncol(X)-1)){
    plot(y=dat_clustered[,i],x=dat_clustered$RT,bg=factor(dat_clustered[,ncol(X)+i]),
         ylab=paste0(colnames(X)[i]," Signal"),pch=21,cex=1,lwd=1.5,xlab="Naive RT",
         main=paste0("R2=",round(summary(lm(dat_clustered$RT ~ dat_clustered[,i]))$adj.r.squared,2)))
  
}

#Boxplot for signal containing clusters
mypar(2,5)
for(i in 1:(ncol(X)-1)){
  input <- dat_clustered$RT
  ind1 <- dat_clustered[,ncol(X)+i] == "low"
  ind2 <- dat_clustered[,ncol(X)+i] == "high"
  boxplot(input[ind1],input[ind2],main=colnames(dat_clustered)[i],
          col=c("lightblue","dodgerblue"),ylab="RT",names=c("Low","High"))
}

#Compute Jaccard similarity between high signal regions
jacc_res <- sapply(1:10,function(x){
  sapply(1:10, function(y){
    set1 <- colnames(dat_clustered)[11+x]
    set2 <- colnames(dat_clustered)[11+y]
    sets <- paste0(set1,".vs.",set2)
    set1 <- rownames(dat_clustered)[dat_clustered[,set1] == "high"]
    set2 <- rownames(dat_clustered)[dat_clustered[,set2] == "high"]
    intersection <- length(intersect(set1,set2))
    union <- length(union(set1,set2))
    jacc_similarity <- intersection/union
    jacc_similarity
  })
})
jacc_res <- data.frame(jacc_res)
colnames(jacc_res) <- colnames(dat_clustered)[12:ncol(dat_clustered)]
rownames(jacc_res) <- colnames(dat_clustered)[12:ncol(dat_clustered)]
col_fun = circlize::colorRamp2(c(0,0.5,1), c("red", "white", "blue"))
Heatmap(jacc_res,col=col_fun)

##Determine contribution of each selected feature to explain RT change
fit_results <- sapply(1:(ncol(X)-1),function(x){
  input <- X[,c(x,ncol(X))]
  ind <- dat_clustered[,c(11+x)] == "high"
  input <- input[ind,]
  fit <- lm(RT ~ .,data = input)
  fit_summary <- summary(fit)
  fstat <- fit_summary$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  out <- c(colnames(input)[1],fit_summary$adj.r.squared,fit_summary$coefficients[2,1],pval)
  names(out) <- c("Predictor","Adj.R2","Coefficient","P-value")
  out
})
fit_results <- data.frame(t(fit_results))
fit_results <- cbind(fit_results[,1],data.frame(apply(fit_results[,-1],2,as.numeric)))
colnames(fit_results)[1] <- "Feature"

rafalib::mypar(1,1,mar = c(2.5, 10, 1.6, 1.1))
barplot(fit_results$Adj.R2[order(fit_results$Adj.R2)],horiz=T,las=1,xlab=expression("Adjusted R" ^ 2),col = "gray",
        cex.names=0.8,names=fit_results$Feature[order(fit_results$Adj.R2)])

#EOF