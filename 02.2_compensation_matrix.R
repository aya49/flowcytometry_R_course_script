## creates the unmixing and compensation spill over matrix 
## for spectral and flow cytometry data.

## unmixing for spectral cytometry ####

## install and load package
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("flowSpecs")

# load example data
# beads = single fluorochrome staining, dead = dead marker, PBMC = no stain
library("flowSpecs")
data("fullPanel") # raw data
fullPanel
data("unmixCtrls") # control data
sampleNames(unmixCtrls)

# create unmixing matrix
spec_matrix <- flowSpecs::specMatCalc(
    unmixCtrls, 
    groupNames = c("Beads_", "Dead_"), # single-stained, unstained files
    autoFluoName = "PBMC_unstained.fcs" # control
)
spec_matrix
dim(spec_matrix)

# unmix and get unmixed fcs file
f <- flowSpecs::specUnmix(fullPanel, spec_matrix)
f


## compensation for flow cytometry (manual) ####

## install and load package
devtools::install_github("DillonHammill/CytoExploreRData")

# save and load example controls
library(CytoExploreRData)
CytoExploreR::cyto_save(Compensation, save_as="compensation_samples")
CytoExploreR::cyto_save(Activation, save_as="activation_samples")

# set compensation controls
gs <- CytoExploreR::cyto_setup(
    "compensation_samples", 
    gatingTemplate="compensation_gating_template.csv"
)

# transform fluorescent channels - logicle
gs <- CytoExploreR::cyto_transform(gs)

# gate cells / single cells
CytoExploreR::cyto_gate_draw(
    gs, parent="root", alias="Cells", channels=c("FSC-A", "SSC-A"))
CytoExploreR::cyto_gate_draw(
    gs, parent="Cells", alias="Single Cells", channels=c("FSC-A", "FSC-H"))

# calculate spillover matrix
# spill <- cyto_spillover_compute(
#     gs, parent="Single Cells", spillover="spillover_matrix.csv")
spill <- cyto_spillover_compute(
    gs, parent=NULL, spillover="spillover_matrix.csv")


## compensation for flow cytometry (auto) ####
# more reading: https://github.com/carlosproca/autospill
devtools::install_github("carlosproca/autospill")


## cleaning (flowcut, an alternative to peacoqc) ####
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("flowCut")

# load example data
library("flowCut")
data(flowCutData)

# clean a file using flowCut
res_fc <- flowCut::flowCut(flowCutData[[1]], FileID=1, Plot="All", Verbose=TRUE)
# for (i in seq_len(length(flowCutData))) {
#     res_fc <- flowCut::flowCut(flowCutData[[i]], FileID=i, Plot="All", Verbose=TRUE)
# }

# get metadata of results
res_fc$data

# get cleaned fcs file
f <- res_fc$frame
f
