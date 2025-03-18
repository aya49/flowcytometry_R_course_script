## batch normalization

## batch normalization (cytonorm) ####

## install and load package
devtools::install_github('saeyslab/CytoNorm')

# load example data
dir <- system.file("extdata", package = "CytoNorm")
files <- list.files(dir, pattern = "fcs$")
files

# create meta-data
data <- data.frame(
    File=files,
    Path=file.path(dir, files),
    Type=stringr::str_match(files, "_([12]).fcs")[, 2],
    Batch=stringr::str_match(files, "PTLG[0-9]*")[, 1],
    stringsAsFactors=FALSE
)
data$Type <- c("1"="Train", "2"="Validation")[data$Type]
data

train_data <- dplyr::filter(data, Type=="Train")
validation_data <- dplyr::filter(data, Type=="Validation")

# preprocess: transform data
ff <- flowCore::read.FCS(data$Path[1])
channels <- flowCore::colnames(ff)[
    c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
      40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
      49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30, 
      28, 29, 25, 35)]
transformList <- flowCore::transformList(channels, CytoNorm::cytofTransform)
transformList.reverse <- flowCore::transformList(channels, CytoNorm::cytofTransform.reverse)

# cluster data with flowSOM
fsom <- CytoNorm::prepareFlowSOM(
    train_data$Path, channels, nCells=6000,
    FlowSOM.params=list(xdim=5, ydim=5, nClus=10, scale=FALSE),
    transformList=transformList, seed=1)
# cvs <- CytoNorm::testCV(fsom, cluster_values = c(5, 10, 15)) # we take 10 clusters; variance < 1-2

# batch effect removal based on the mean of the training data
# normalized files are saved in folder getwd() > "Normalized"
model <- CytoNorm::CytoNorm.train(
    files=train_data$Path,
    labels=train_data$Batch,
    channels=channels,
    transformList=transformList,
    FlowSOM.params=list(nCells=6000, xdim=5, ydim=5, nClus=10, scale=FALSE),
    normMethod.train=CytoNorm::QuantileNorm.train,
    normParams=list(nQ=101, goal="mean"),
    seed=1,
    verbose=TRUE)
CytoNorm::CytoNorm.normalize(
    model=model,
    files=validation_data$Path,
    labels=validation_data$Batch,
    transformList=transformList,
    transformList.reverse=transformList.reverse,
    normMethod.normalize=CytoNorm::QuantileNorm.normalize,
    outputDir="Normalized", # output folder
    prefix="Norm_",
    clean=TRUE,
    verbose=TRUE)

## TRY: practice problems ####

# 1. Loop through the original "validation" files in validation_data$Path. Run their UMAP plots (see file "...clustering.R"). 

# 2. Compare the validation file's UMAP plots with their normalized version (see the "Normalized" folder). Do you see any changes?

# 3. Compare the normalized file's UMAP plots with those of the "training" files in train_data$Path; do the normalized files look more similar to the training files of their corresponding batch? (see batch in train_data$Batch, validation_data$Batch)
