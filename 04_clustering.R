## cluster and visualize fcs file at high dimensionality
## author: alice yue
## input: raw fcs file
## output: preprocessed fcs file

# load packages
library("flowCore")
library("flowWorkspace") # also imports ggcyto and ggplot2 packages
library("Rphenograph") # also imports igraph package
library("FlowSOM")
library("Rtsne")
library("umap")

# directory to save results in
res_dir <- "/home/alice/res"

# set seed for randomness
set.seed(4)

# load gating set
gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))
cpops <- flowWorkspace::gs_get_leaf_nodes(gs)

# prepare colours for plots!
cpopcol <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
             "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")

# prepare data
starting_cpop <- "lymphocytes"
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(gs, starting_cpop), Class="list")[[1]] )
starting_inds <- flowWorkspace::gh_pop_get_indices(gs, starting_cpop)
markers <- flowCore::markernames(ff)
cpops <- flowWorkspace::gs_get_leaf_nodes(gs)

# prepare cell population labels from gating
cpops_matrix <- matrix(FALSE, ncol=length(cpops), nrow=sum(starting_inds))
colnames(cpops_matrix) <- cpops
for (cpop in cpops) {
    cpops_matrix[,cpop] <- 
        flowWorkspace::gh_pop_get_indices(gs, cpop)[starting_inds]
}


## cluster data
n_clust <- 30
n_sample <- 5000

# phenograph
# TRY CHANGING k (number of neighbours)
# TRY CHANGING names(markers) to only your markers of interest e.g. "PerCP-Cy5-5-A", "APC-A"
pc_ <- Rphenograph::Rphenograph(flowCore::exprs(ff)[, names(markers)], k=50) 
pc <- igraph::membership(pc_[[2]])

# flowsom
# TRY CHANGING nClus (# clusters)
# TRY CHANGING names(markers) to only your markers of interest e.g. "PerCP-Cy5-5-A", "APC-A"
fc_ <- FlowSOM::FlowSOM(ff, nClus=n_clust, colsToUse=names(markers))
fc <- as.numeric(FlowSOM::GetMetaclusters(fc_))


## sample data to reduce runtime/memory consumption
# ensure there are points from each cell population and cluster are sampled
ffsample <- sample(which(rowSums(cpops_matrix)>0), n_sample)

for (i in unique(pc)) {
    clusti <- which(pc==i)
    if (!any(clusti%in%ffsample)) {
        ffsample <- append(ffsample, sample(clusti, 3))
    }
}
for (i in unique(fc)) {
    clusti <- which(fc==i)
    if (!any(clusti%in%ffsample)) {
        ffsample <- append(ffsample, sample(clusti, 3))
    }
}
for (i in seq_len(ncol(cpops_matrix))) {
    clusti <- which(cpops_matrix[,i])
    if (!any(clusti%in%ffsample)) {
        ffsample <- append(ffsample, sample(clusti, 3))
    }
}
ffs <- flowCore::exprs(ff)[ffsample, names(markers)]

## make 2D scatterplots
print(names(markers))
# choose to markers to plot!
# TRY CHANGING fc to pc (flowsom vs phenograph)
df <- data.frame(marker1=ffs[,c("PerCP-Cy5-5-A")],
                 marker2=ffs[,c("APC-A")],
                 cluster=factor(fc[ffsample])) # prepare data frame
ggplot2::ggplot(df, ggplot2::aes(x=marker1, y=marker2, colour=cluster)) + 
    ggplot2::geom_point(size=.1)

## reduce dimensionality of data to 2D
t2 <- Rtsne::Rtsne(ffs)$Y
u2 <- umap::umap(ffs)$layout
# TRY adjusting "n_neighbors" and 
#               "metric" for distance metric used 
#               i.e. any of: "cosine", "manhattan", "hamming", "correlation"
# install.packages("uwot")
# u2 <- uwot::umap(ffs, n_neighbors=15, metric="euclidean")
colnames(t2) <- colnames(u2) <- c("x", "y")


## plot clusters and ground truth cell populations in 2D
tu2d <- data.frame(
    tx=t2[,1], ty=t2[,2], ux=u2[,1], uy=u2[,2], 
    phenograph=pc[ffsample], flowsom=fc[ffsample],
    cpops=sapply(apply(cpops_matrix[ffsample,], 1, function(x) 
        (cpops[x])), function(y) y[1]) )

# tsne
gptp <- ggplot2::ggplot(tu2d, ggplot2::aes(x=tx, y=ty, colour=factor(phenograph))) + ggplot2::geom_point(size=0.5)
gptf <- ggplot2::ggplot(tu2d, ggplot2::aes(x=tx, y=ty, colour=factor(flowsom))) + ggplot2::geom_point(size=0.5)
gptc <- ggplot2::ggplot(tu2d, ggplot2::aes(x=tx, y=ty, colour=factor(cpops))) + ggplot2::geom_point(size=0.5)

# umap
gpup <- ggplot2::ggplot(tu2d, ggplot2::aes(x=ux, y=uy, colour=factor(phenograph))) + ggplot2::geom_point(size=0.5)
gpuf <- ggplot2::ggplot(tu2d, ggplot2::aes(x=ux, y=uy, colour=factor(flowsom))) + ggplot2::geom_point(size=0.5)
gpuc <- ggplot2::ggplot(tu2d, ggplot2::aes(x=ux, y=uy, colour=factor(cpops))) + ggplot2::geom_point(size=0.5)

# display
gptp
gptf
gptc

gpup
gpuf
gpuc

# # saving a ggplot is a bit different from saving the default plot:
# ggplot2::ggsave(filename=paste0(res_dir, "/2D_tsne_phenograph.png"), plot=gpup)


## EXTRA: cluster evaluation ####
install.packages("TreeDist") # isntall package
install.packages("cluster")
install.packages("clValid")

data(iris) # load example data, a matrix with 4 columns + a 5th column indicating the true class
mydata <- as.matrix(iris[,1:4])
d <- dist(mydata)
mycol  <- as.factor(iris[,5])

# kmeans++ CLUSTERING automatically initializes clusters
klusters <- TreeDist::KMeansPP(mydata, k=5) # k = number of clusters
plot(mydata, col=klusters$cluster, main="k=5", pch=19, cex=0.5)
plot(mydata, col=mycol, main="truth", pch=19, cex=0.5)

# silhouette index; best plot looks like a triangle per cluster, not connected.
# finds dense round clusters based on cluster size
# minimize distance between events of the same cluster, normalized for cluster size.
si <- cluster::silhouette(x=klusters$cluster, dist=d)
plot(si, col = c("red", "green", "blue", "purple", "yellow"))

# dunn index; best if the value is small, bad if 1.
# finds dense round clusters based on cluster size
# minimize distance between events of the same cluster, normalized for cluster size.
clValid::dunn(d, klusters$cluster)

# other cluster evaluation statistics: 
# DB index, CS  index, Xie-Beni index, Dube's separation index, 
# Coulder and Odell's separation index, Scale-free weighted ratio, AIC/BIC
