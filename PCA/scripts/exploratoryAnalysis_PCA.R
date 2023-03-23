# initializing libraries
library("GEOquery")
library("useful")
library("rgl")
library("plotly")


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
# ex2 <- na.omit(as.matrix(ex))
# dim(ex2)
# no need, no null values


#PCA
pca <- prcomp(t(ex))
screeplot(pca)

#create array of colors that correspond to labels!
gse$characteristics_ch1.1
gse$characteristics_ch1.2

#set labels for pca
labels = c()
group_labels = c()
complete_labels = c()
for (x in 1:44) {
  # BEFORE SAMPLING -> control group
  if (gse$characteristics_ch1.2[x] == "time of sampling (before/after): before"){
    labels[x] <- "green"
    group_labels[x] <- "control"
    complete_labels[x] <- paste(colnames(ex)[x], "(", "control", ")")
  } 
  # AFTER SAMPLING
  else {
    if (gse$characteristics_ch1.1[x] == "protein content restricted diet: high") {
      labels[x] <- "red"
      group_labels[x] <- "ER_HP"
      complete_labels[x] <- paste(colnames(ex)[x], "(", "ER_HP", ")")
    } else {
      labels[x] <- "blue"
      group_labels[x] <- "ER_SP"
      complete_labels[x] <- paste(colnames(ex)[x], "(", "ER_SP", ")")
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

t <- as.data.frame(pca$x)
plot_ly(t, x~PC1, y~PC2, z~PC3, col = labels)

