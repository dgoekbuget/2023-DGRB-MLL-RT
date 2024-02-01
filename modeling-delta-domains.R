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

setwd("/blellochlab/data1/deniz/analysis/mll-rt-paper/ml/domains-delta")
files <-dir("naive",pattern="backgroundCorSum.bed")

#Load data
for (i in 1:length(files)) {
    if(i == 1){
        rt <- read_delim(paste0("naive/",files[i]),"\t", escape_double = FALSE, trim_ws = TRUE,col_names=F)[,1:6]
        colnames(rt) <- c("chr","start","end","dRT",
                          gsub(".RTdomain.backgroundCorSum.bed",".N1",files[i]),
                          gsub(".RTdomain.backgroundCorSum.bed",".N2",files[i]))
    } else {
        rt <- cbind(rt,read_delim(paste0("naive/",files[i]),"\t", escape_double = FALSE, trim_ws = TRUE,col_names=F)[,5:6])
        colnames(rt)[(ncol(rt)-1):ncol(rt)] <- c(gsub(".RTdomain.backgroundCorSum.bed",".naive.N1",files[i]),gsub(".RTdomain.backgroundCorSum.bed",".naive.N2",files[i]))
    }
}

for (i in 1:length(files)) {
        rt <- cbind(rt,read_delim(paste0("form/",files[i]),"\t", escape_double = FALSE, trim_ws = TRUE,col_names=F)[,5:6])
        colnames(rt)[(ncol(rt)-1):ncol(rt)] <- c(gsub(".RTdomain.backgroundCorSum.bed",".form.N1",files[i]),gsub(".RTdomain.backgroundCorSum.bed",".form.N2",files[i]))
}

#Create timing categories
rt$timing <- ifelse(rt$dRT > 0, "E",ifelse(rt$dRT < 0, "L",ifelse(rt$dRT == 0,NA,NA)))

#Separate counts from meta
colnames(rt)
meta <- rt[,1:4]
counts <- rt[,-c(1:4)]
counts <- cbind(counts$timing,log2(counts[,-ncol(counts)]/(abs(meta$end-meta$start))+0.01))
colnames(counts)[1] <- "timing"

counts$timing <- factor(counts$timing)
dat <- counts[complete.cases(counts),]
meta <- meta[complete.cases(counts),]

#Calculate differential signal
ncol(dat)
dat <- cbind(dat[,1],dat[,44:85]-dat[,2:43])
colnames(counts)[1] <- "timing"

colnames(dat) <- c("timing",paste0(rep(c("H2AX","H2AZ","H2BK5ac","H3.3","H3K14ac","H3K27ac",
                                         "H3K27me3","H3K36me1","H3K36me2","H3K4me2","H3K4me3",
                                         "H3K9ac","H3K9me2","H3K9me3","H4K14ac","H4K20me1","H4K8ac","H2AK119Ub",
                                         "H3K4me1","PPolII","Rad21"),each=2),c(".N1",".N2")))
dat <- dat[,c("timing",paste0(rep(c("PPolII","Rad21","H3K4me1","H3K4me2","H3K4me3",
                                    "H3K27ac","H3K36me1","H3K36me2","H2BK5ac","H3K14ac",
                                    "H2AX","H2AZ","H3.3", "H3K9ac", "H4K8ac","H4K14ac",
                                    "H4K20me1","H3K27me3","H3K9me2","H3K9me3","H2AK119Ub"),each=2),c(".N1",".N2")))
]

#Calculate average chromatin mark signal
mean_dat <- data.frame(dat$timing)
for (i in seq(1,42,by=2)) {
  m <- data.frame(rowMeans(dat[,c(i+1,i+2)]))
  colnames(m) <- gsub(".N1","",colnames(dat)[i+1])
  mean_dat <- cbind(mean_dat,m)
}

#Define color vectors
color21 <- colorRampPalette(c("#336699","#FF9966"))(100)

#Machine learning feature selection
set.seed(178)
X <- mean_dat[,-1]
X$RT <- meta$dRT

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
     ylab ="Observed",main=paste0("dRT"," rho=",round(cor(predictions$pred,predictions$obs,method="spearman"),3)),
     pch=21,bg="gray",cex=1,lwd=1.5,ylim=c(-5,5),xlim=c(-5,5))
abline(lm(predictions$obs~predictions$pred),col="dodgerblue",lty="dashed",lwd=3)

##Plot parameter weights (relative feature importance)
# Define a range of colors for the gradient
coefs <- coef(fit$finalModel,fit$bestTune$lambda)
rafalib::mypar(1,1,mar = c(2.5, 5, 1.6, 1.1))
barplot(coefs[-1][order(coefs[-1])],horiz=T,las=1,xlab="Parameter weight",col = "gray",
        cex.names=0.8,names=rownames(coefs)[-1][order(coefs[-1])],xlim=c(-0.5,0.5))

#Calculate pairwise prediction using glm
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

#col_fun = circlize::colorRamp2(c(0.5,1,1.5), c("white", "blue", "white"))
col_fun = circlize::colorRamp2(c(0,0.5), c("white", "blue"))
Heatmap(residual_dispersion_results,col=col_fun,name="Adj. R2")

##Perform correlation analysis of selected features
set.seed(178)
selected <- rownames(coefs)[-1][order(abs(coefs[-1]),decreasing = T)][1:10]
X <- mean_dat[,-1][,selected]
X$RT <- meta$dRT

##Separate signal of selected marks by kmeans into no signal - background/low - positive signal
k=3 #determine 3 clusters
dat_clustered <- X
for(i in 1:(ncol(X)-1)){
  input <- X[,c(i,ncol(X))]
  dist_matrix <- dist(input,method="manhattan")  # Calculate the distance matrix
  hierarchical_result <- hclust(dist_matrix, method = "ward.D")
  cluster_assignment <- cutree(hierarchical_result, k = k)
  dat_clustered$new <- cluster_assignment
  colnames(dat_clustered)[ncol(dat_clustered)] <- paste0(colnames(X)[i],".cluster")
}

#Highlight clusters in scatterplot
mypar(5,2)
for(i in 1:(ncol(X)-1)){
  plot(y=dat_clustered[,i],x=dat_clustered$RT,bg=factor(dat_clustered[,ncol(X)+i]),
       ylab=colnames(X)[i],pch=21,cex=1,lwd=1.5,xlab="delta RT",
       main=paste0("R2=",round(summary(lm(dat_clustered$RT ~ dat_clustered[,i]))$adj.r.squared,2)))
}

#Boxplot for signal containing clusters
mypar(1,6)
for(i in 1:(ncol(X)-1)){
  input <- dat_clustered$RT
  ind1 <- dat_clustered[,ncol(X)+i] == 1
  ind2 <- dat_clustered[,ncol(X)+i] == 2
  ind3 <- dat_clustered[,ncol(X)+i] == 3
  boxplot(input[ind1],input[ind2],input[ind3],main=colnames(dat_clustered)[i],
          col=c(1,2,3),ylab="RT")
}

#Predictive power of H3K4me1 vs PPol2
fit <- glm(RT ~ PPolII, data=X)
dev <- residuals(fit,type="deviance")
dfr <- df.residual(fit)
sum(dev^2) / dfr
mypar(1,3)

fit <- glm(RT ~ H3K4me1, data=X)
dev <- residuals(fit,type="deviance")
dfr <- df.residual(fit)
sum(dev^2) / dfr

fit <- glm(RT ~ PPolII+H3K4me1, data=X)
dev <- residuals(fit,type="deviance")
dfr <- df.residual(fit)
sum(dev^2) / dfr

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





























# Calculate Spearman correlations for each column
correlations <- lapply(dat[,-1], function(col) cor(meta$dRT, col, method = "spearman"))

# Combine the correlations and group information into a data frame
correlation_df <- data.frame(Group = gsub(".N1|.N2","",names((unlist(correlations)))), Replicate = rep(1:2, times = 21), Correlation = unlist(correlations))

# Calculate the mean correlation per group
mean_correlations <- correlation_df %>%
  group_by(Group) %>%
  dplyr::summarize(Mean_Correlation = mean(Correlation))

# Reorder the "Group" factor variable based on "Mean Correlation"
mean_correlations <- mean_correlations %>%
  arrange(Mean_Correlation) %>%
  mutate(Group = factor(Group, levels = Group))

# Create a bar plot
ggplot(mean_correlations, aes(x = Mean_Correlation, y = factor(Group), fill = factor(Group))) +
  geom_bar(stat = "identity", position = "identity") +
  geom_jitter(data = correlation_df, aes(x = Correlation, y = factor(Group)), shape=21,size=3, 
              position = position_jitterdodge(jitter.width = 6, dodge.width = 0.7)) +
  labs(x = "Mean Correlation", y = "Group") +
  theme_minimal() +
  theme(legend.position = "none") 

# Calculate Spearman correlations and p-values for each group
input <- dat[,-1]
stats <- sapply(seq(1,42,by=2), function(x){
  c(colnames(input[x]),cor.test(meta$dRT,input[,c(x:x+1)])$p.value)
})
nms <- stats[1,]
nms <- gsub(".N1","",nms)
stats <- as.numeric(stats[2,])
names(stats) <- nms
ind <- match(mean_correlations$Group,names(stats))
stats <- stats[ind]

# Define a function to map values to colors using a logarithmic scale
gradient_colors <- colorRampPalette(c("darkgray", "white"))(100)
map_values_to_colors_log <- function(values, gradient_colors) {
  # Ensure values are positive for the logarithmic scale
  values[values <= 0] <- 1e-300  # Avoid taking the logarithm of zero or negative values
  log_values <- log10(values)
  
  # Scale the logarithmic values to the range [0, 1]
  scaled_values <- (log_values - min(log_values)) / (max(log_values) - min(log_values))
  
  # Find the color indices
  color_indices <- floor(scaled_values * (length(gradient_colors) - 1)) + 1
  
  # Interpolate colors
  interpolated_colors <- gradient_colors[color_indices]
  
  return(interpolated_colors)
}

# Map values to colors with a logarithmic scale
mapped_colors_log <- map_values_to_colors_log(stats, gradient_colors)

ggplot(mean_correlations, aes(x = Mean_Correlation, y = factor(Group), fill=factor(Group))) +
  geom_bar(stat = "identity", position = "identity",color="black") +
  geom_jitter(data = correlation_df, aes(x = Correlation, y = factor(Group)), shape=21,size=3,
              position = position_jitterdodge(jitter.width = 6, dodge.width = 0.7)) +
  scale_fill_manual(values=mapped_colors_log) +
  labs(x = "Mean Correlation", y = "Group") +
  theme_minimal() +
  theme(legend.position = "none") 

#Repeat with caret
set.seed(178)
X <- dat[,-1]
X$RT <- meta$dRT

control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = F)

fit_results <- sapply(1:(ncol(X)-1),function(x){
  input <- X[,c(x,ncol(X))]
  fit <- train(RT ~ .,
               data = input,
               method = "glm",
               preProcess = c("center", "scale"),
               tuneLength = 25,
               trControl = control)
  out <- c(colnames(input)[1],fit$results$RMSE,fit$results$Rsquared,coef(summary(fit))[2,1])
  names(out) <- c("Predictor","RMSE","R2","Coefficient")
  out
})
fit_results <- data.frame(t(fit_results))
group <- gsub(".N1|.N2","",fit_results$Predictor)
fit_results <- data.frame(apply(fit_results[,-1],2,as.numeric))
fit_results$group <- group
fit_results$replicate <- rep(c(1,2),nrow(fit_results)/2)
n <- nrow(X)

mean_fit_results <- fit_results %>% group_by(group) %>%
  dplyr::summarize(
    mean_rmse = mean(RMSE),
    mean_r2 = mean(R2),
    mean_coef = mean(Coefficient),
    pvalue = t.test(c(Coefficient),mu=0)$p.value #2*pt(-abs(diff(Coefficient)/diff(SE_Coef)),n-2)
  ) %>%
  arrange(mean_coef) %>%mutate(group = factor(group, levels = group))

p1 <- ggplot(data= mean_fit_results, aes(y=group,x=mean_coef)) +
  geom_bar(stat = "identity", position = "identity",fill="gray",color="black") +
  geom_jitter(data = fit_results, aes(x = Coefficient, y = group),fill="gray", shape=21,size=3) +
  theme_minimal()

p2 <- ggplot(data= mean_fit_results, aes(y=group,x=mean_r2)) +
  geom_bar(stat = "identity", position = "identity",fill="gray",color="black") +
  geom_jitter(data = fit_results, aes(x = R2, y = group),fill="gray", shape=21,size=3) +
  theme_minimal()
grid.arrange(p1,p2,ncol=2)

#Fit a regularized model to obtain overall predictive power of epigenome on RT
set.seed(178)
X <- dat[,-1]
X$RT <- meta$dRT

control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = F)
fit <- train(RT ~ .,
             data = X,
             method = "glmnet",
             preProcess = c("center", "scale"),
             tuneLength = 25,
             trControl = control)

predictions <- extractPrediction(list(fit))
mypar()
plot(predictions$pred,predictions$obs,xlab = "Predicted",
     ylab ="Observed",main=paste0("dRT"," rho=",round(cor(predictions$pred,predictions$obs,method="spearman"),3)),
     pch=21,bg="gray",cex=1,lwd=1.5,xlim=c(-4.5,4.5),ylim=c(-4.5,4.5))
abline(lm(predictions$obs~predictions$pred),col="dodgerblue",lty="dashed",lwd=3)

