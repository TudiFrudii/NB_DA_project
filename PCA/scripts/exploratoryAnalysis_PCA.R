# initializing libraries
library("GEOquery")
library("useful")
library("rgl")

# import data from geo, extract data matrix
gse <- getGEO("GSE84046", destdir = ".", getGPL = FALSE)
gse <- gse[[1]]
ex <- exprs(gse)
dim(ex)

# exploratory analysis
boxplot(ex, main = "Boxplot of raw values")
# no need for normalization

# std way of dealing with missing values
# remove all lines that han one or more omits
ex2 <- na.omit(as.matrix(ex))
dim(ex2)
# no need, no null values


#PCA
pca <- prcomp(t(ex))
screeplot(pca)

#create array of colors that correspond to labels!
gse$characteristics_ch1.1
gse$characteristics_ch1.2

#set labels for pca
labels = c()
for (x in 1:44) {
  # BEFORE SAMPLING -> control group
  if (gse$characteristics_ch1.2[x] == "time of sampling (before/after): before"){
    labels[x] <- "green"
  } 
  # AFTER SAMPLING
  else {
    if (gse$characteristics_ch1.1[x] == "protein content restricted diet: high"){
      labels[x] <-  "red"
    } else {
      labels[x] <- "blue"
    }
  }
}

# 2D PCA PLOTS
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=labels)
plot(pca$x[,2], pca$x[,3], xlab="PCA2", ylab="PCA3", main="PCA for components 2&3", type="p", pch=10, col=labels)
plot(pca$x[,1], pca$x[,3], xlab="PCA1", ylab="PCA3", main="PCA for components 1&3", type="p", pch=10, col=labels)

text(pca$x[,1], pca$x[,2], rownames(pca$x), cex=0.75)

# 3D PCA PLOT
scores = as.data.frame(pca$x)
plot3d(scores[,1:3], size=5, col = labels)



