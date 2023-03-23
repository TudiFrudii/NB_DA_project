library("GEOquery")
library("useful")
library("rgl")

gse <- getGEO("GSE84046", destdir = ".", getGPL = FALSE)
gse <- gse[[1]]
ex <- exprs(gse)

# set labels 

labels = c()
group_labels = c()
complete_labels = c()
for (x in 1:44) {
  # BEFORE SAMPLING -> control group
  if (gse$characteristics_ch1.2[x] == "time of sampling (before/after): before"){
    labels[x] <- "green"
    group_labels[x] <- "control"
    complete_labels[x] <- paste(colnames(ex), "(", "control", ")")
  } 
  # AFTER SAMPLING
  else {
    if (gse$characteristics_ch1.1[x] == "protein content restricted diet: high") {
      labels[x] <- "red"
      group_labels[x] <- "ER_HP"
      complete_labels[x] <- paste(colnames(ex), "(", "ER_HP", ")")
    } else {
      labels[x] <- "blue"
      group_labels[x] <- "ER_SP"
      complete_labels[x] <- paste(colnames(ex), "(", "ER_SP", ")")
    }
  }
}

# K_MEANS

# number of clusters in my dataset
k <- 3
kmeans_results = kmeans(t(ex), k)
table(kmeans_results$cluster)
# plot(kmeans_results, data=t(ex)) + geom_text(aes(label = colnames(ex)), hjust =0, vjust=0)
plot(kmeans_results, data=t(ex)) + geom_text(aes(label = group_labels), hjust =0, vjust=0)


# HIERARCHICAL CLUSTERING

# square matrix n*n, distance values i-j
dist_matrix = dist(t(ex))
hc_results = hclust(dist_matrix, method = "ave")
hc_results = hclust(dist_matrix, method = "complete")
hc_results = hclust(dist_matrix, method = "single")

k <- 3
groups = cutree(hc_results, k=k)
# plot(hc_results, hang <-1, labels = groups)
plot(hc_results, hang <-1, labels = complete_labels, main = "hierarchical clustering")
rect.hclust(hc_results, k=3, which = NULL, x = NULL, h = NULL, cluster = NULL, border = 2)

