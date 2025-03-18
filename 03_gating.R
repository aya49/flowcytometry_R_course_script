## gate fcs file
## author: alice yue
## input: raw fcs file
## output: gated flowworkspace

# load packages
library("flowCore")
library("PeacoQC")
library("flowWorkspace")
library("flowDensity")
# library("ggpointdensity") # for easier density scatterplots

# path to raw fcs file
# file from: http://flowrepository.org/id/FR-FCM-ZYXN
# download file: https://drive.google.com/file/d/1PpSM93GTj9zejVDZzD89_k3sx7Lc-TQl/view?usp=sharing
# we will be gollowing this gating strategy: https://docs.google.com/presentation/d/1dUamRDWtN6cuZXaUN1I1baDUopPhhufKTGXZAXOSP3o/edit?usp=sharing
fcs_path <- "/home/alice/sangerP2.fcs"

# directory to save results in
res_dir <- "/home/alice/res"
dir.create(res_dir)

gateplot_dir <- paste0(res_dir, "/gs_plots")
dir.create(gateplot_dir)

# load fcs file
f <- flowCore::read.FCS(fcs_path)

# explore fcs file
# flowCore::exprs(f) # raw cell x marker matrix
head(flowCore::exprs(f))
dim(flowCore::exprs(f))
flowCore::markernames(f) # marker names (excluding morphology columns)
f@parameters@data

## 1.1 compensate ####
spillover <- flowCore::keyword(f)$SPILL
f <- flowCore::compensate(f, spillover=spillover)

## 1.2 logicle transform ####
transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
f <- flowWorkspace::transform(f, transformList)

## 1.3 cleaning; see res_path for plot ####
fmr <- PeacoQC::RemoveMargins(f, channels=c(1,4), output="full")
pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=colnames(spillover),
                        plot=TRUE, save_fcs=FALSE, report=FALSE,
                        output_directory=res_dir)
f <- pQC[["FinalFF"]]

# clean memory by removing variables
rm("fmr")
rm("pQC")
gc()


## 2 GATING: cell population identification based on a gating strategy ####
# TRY: completing this gating strategy

## function to rotate 2D frame
## input: 2D matrix and angle
## output: rotated 2D matrix
rotate_data <- function(data, theta=pi/2 - atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])) {
    data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}

# initialize gating set (containing 1 file)
fs <- flowCore::flowSet(list(f))
# fs <- flowCore::flowSet(list(f1, f2, f3)) # can add more files if needed
gs <- flowWorkspace::GatingSet(fs)


## 2.1 gating all > singlets ####

# get threshold gates
gate_singlets <- flowDensity::deGate(
    f, channel="SSC-W", 
    use.upper=TRUE, upper=TRUE, tinypeak.removal=0.95)
gate_singlets_low <- flowDensity::deGate(
    f, channel="SSC-W", 
    use.percentile=TRUE, percentile=0.0001)

# gate
temp <- flowDensity::flowDensity(
    f, channels=c("FSC-A", "SSC-W"), 
    position=c(NA,FALSE), gates=c(NA, gate_singlets))
fd_singlets <- flowDensity::flowDensity(
    temp, channels=c("FSC-A", "SSC-W"), 
    position=c(NA,TRUE), gates=c(NA, gate_singlets_low))

# plot
png(paste0(gateplot_dir, "/01_all_singlets.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(f)[,"FSC-A"])
plot(d, ylab="", axes=FALSE, 
     main="all events (cells, particles, etc.) > singlets")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(f)[,"SSC-W"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=c(gate_singlets, gate_singlets_low), lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(f, channels=c("FSC-A", "SSC-W"), main="")
lines(fd_singlets@filter)
abline(h=c(gate_singlets, gate_singlets_low), lty="dashed", col="red")
graphics.off()


## 2.2 gating singlets > live ####

# get threshold gates
temp <- flowDensity::flowDensity(
    fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), 
    position=c(NA,F), gates=c(NA,50000))
gate_live <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="APC-Cy7-A")

# gate
fd_live <- flowDensity::flowDensity(
    fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), 
    position = c(FALSE,NA), gates=c(gate_live,NA))

# register gate into gs
gate <- flowCore::rectangleGate(
    filterId="live", 
    "APC-Cy7-A"=c(-1, gate_live), # include all
    "SSC-W"=c(gate_singlets_low, gate_singlets))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="root")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "live")

# plot and save as png
png(paste0(gateplot_dir, "/02_singlets_live.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(fd_singlets))[,"APC-Cy7-A"])
plot(d, ylab="", axes=FALSE, main="singlets > live")
abline(v=gate_live, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(fd_singlets))[,"SSC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")

par(mar=c(5,5,1,1))
flowDensity::plotDens(fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), main="")
abline(v=gate_live, lty="dashed", col="red")
graphics.off()

# clean memory by removing variables
rm(fd_singlets)
gc()


## 2.3 gating live > lymphocytes ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "live"), Class="list")[[1]] )

# get upper limit gates
gate_ssca_high <- flowDensity::deGate(
    ff, channel="SSC-A", 
    use.percentile=TRUE, percentile=0.9999999)
gate_fsca_high <- flowDensity::deGate(
    ff, channel="FSC-A", 
    use.percentile=TRUE, percentile=0.99999999)

# get threshold gate
gate_fsca <- flowDensity::deGate(
    ff, channel="FSC-A")

# register gate into gs
gate <- flowCore::rectangleGate(
    filterId="lymphocytes", 
    "FSC-A"=c(gate_fsca, gate_fsca_high), # include all
    "SSC-A"=c(0, gate_ssca_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="live")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "lymphocytes")

# plot
png(paste0(gateplot_dir, "/03_live_lymphocytes.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"FSC-A"])
plot(d, ylab="", axes=FALSE, main="live > lymphocytes")
abline(v=c(gate_fsca, gate_fsca_high), lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"SSC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=c(0, gate_ssca_high), lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("FSC-A", "SSC-A"), main="")
abline(v=c(gate_fsca, gate_fsca_high), 
       h=c(0, gate_ssca_high), lty="dashed", col="red")
graphics.off()


## 2.4 gating lymphocytes > not/granulocytes ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "lymphocytes"), Class="list")[[1]] )

# get upper limit gates
gate_cd11b_high <- flowDensity::deGate(
    ff, channel="BV510-A", 
    use.percentile=TRUE, percentile=0.9999999)
gate_cd5_high <- flowDensity::deGate(
    ff, channel="APC-A", 
    use.percentile=TRUE, percentile=0.9999999)

# get threshold gates
gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")
temp <- flowDensity::flowDensity(
    ff, channels=c("APC-A", "BV510-A"), 
    position=c(NA, TRUE), gates=c(NA, gate_cd11b)) #upper half
gate_cd5 <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="APC-A")

# register gates into gs
gate <- flowCore::rectangleGate(
    filterId="granulocytes", 
    "APC-A"=c(gate_cd5, gate_cd5_high), # include all
    "BV510-A"=c(gate_cd11b, gate_cd11b_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="lymphocytes")
flowWorkspace::recompute(gs)

fd_gran <- flowDensity::flowDensity(
    ff, channels=c("APC-A", "BV510-A"), 
    position=c(TRUE,TRUE), gates=c(gate_cd5, gate_cd11b))
fd_notgran <- flowDensity::notSubFrame(
    ff, channels=c("APC-A", "BV510-A"), 
    position="logical", gates="missing", fd_gran@filter)
gate <- flowCore::polygonGate(.gate=fd_notgran@filter)
node <- flowWorkspace::gs_pop_add(
    gs, gate, name="not_granulocytes", parent="lymphocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "granulocytes")
# flowWorkspace::gs_pop_remove(gs, "not_granulocytes")

# plot
png(paste0(gateplot_dir, "/04_lymphocytes_granulocytes.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(temp))[,"APC-A"], na.rm=TRUE)
plot(d, ylab="", axes=FALSE, main="lymphocytes > not/granuloctyes",
     xlim=range(flowCore::exprs(ff)[,"APC-A"], na.rm=TRUE))
abline(v=gate_cd5, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd11b, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("APC-A", "BV510-A"), main="")
abline(v=gate_cd5, h=gate_cd11b, lty="dashed", col="red")
graphics.off()


## 2.5 gating not granulocytes > not/monocytes ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_granulocytes"), Class="list")[[1]] )

# <your gates/register/plot here>
# populations you will gate: "monocytes", "not_monocytes"


## 2.6 gating not monocytes > not/eosoniphil ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_monocytes"), Class="list")[[1]] )

# <your gates/register/plot here>
# populations you will gate: "eosinophils", "not_not_eosinophils"


## 2.7 gating not eosoniphil > CD11b+/- ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_eosinophils"), Class="list")[[1]] )

# <your gates/register/plot here>
# cell populations you will gate: "CD161+", "CD161-"


## 2.8 gating CD161+ > NK/T ####
# NK: natural killer
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "CD161+"), Class="list")[[1]] )

# rotate cells
ffs <- ff
flowCore::exprs(ffs)[,c("APC-A","BV650-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("APC-A","BV650-A")], theta=-pi/6)

# get slanted threshold gates
gate_cd5_slant <- flowDensity::deGate(ffs, channel="APC-A")

# register gates into gs by manually making filter (convex hull; gate outline)
nkTF <- flowCore::exprs(ffs)[,"APC-A"] <= gate_cd5_slant

nkchull <- chull(flowCore::exprs(ff)[nkTF, c("APC-A","BV650-A"), drop=FALSE])
nkfilter <- flowCore::exprs(ff)[
    which(nkTF)[c(nkchull, nkchull[1])], c("APC-A","BV650-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=nkfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="NK", parent="CD161+")
flowWorkspace::recompute(gs)

nktchull <- chull(flowCore::exprs(ff)[!nkTF, c("APC-A","BV650-A"), drop=FALSE])
nktfilter <- flowCore::exprs(ff)[
    which(!nkTF)[c(nktchull, nktchull[1])], c("APC-A","BV650-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=nktfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="NKT", parent="CD161+")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NK")
# flowWorkspace::gs_pop_remove(gs, "NKT")

# plot
png(paste0(gateplot_dir, "/08_CD161p_NKT.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"APC-A"])
plot(d, ylab="", axes=FALSE, main="CD161+ > NK/T")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV650-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("APC-A", "BV650-A"), main="")
lines(nkfilter)
graphics.off()


## 2.9 gating NK > NK im/mature Ly6C+/- ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "NK"), Class="list")[[1]] )

# get loose threshold gates
gate_ly6c <- flowDensity::deGate(ff, channel="Alexa Fluor 700-A") 
gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")

# tighten threshold gates
temp <- flowDensity::flowDensity(
    ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
    position=c(FALSE,NA), gates=c(gate_ly6c, NA))
gate_ly6c <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A")

# register gates into gs
gate <- flowCore::quadGate(
    filterId=c("NK"),
    "Alexa Fluor 700-A"=gate_ly6c,
    "BV510-A"=gate_cd11b)
node <- flowWorkspace::gs_pop_add(gs, gate, parent="NK")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700-CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700+CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700-CD11b BV510+")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700+CD11b BV510+")

# plot
png(paste0(gateplot_dir, "/09_NK_immatureLy6C.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"Alexa Fluor 700-A"])
plot(d, ylab="", axes=FALSE, main="NK > im/mature Ly6C+/-")
abline(v=gate_ly6c, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd11b, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("Alexa Fluor 700-A", "BV510-A"), main="")
abline(v=gate_ly6c, h=gate_cd11b, lty="dashed", col="red")
graphics.off()


## 2.10 gating NKT > NKT im/mature Ly6C+/- ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "NKT"), Class="list")[[1]] )

# <your gates/register/plot here>
# cell populations you will gate: quad gate! names will be automatically given.


## 2.11 gating CD161- > not/tcells ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "CD161-"), Class="list")[[1]] )

# rotate cells
ffs <- ff
flowCore::exprs(ffs)[,c("FITC-A","APC-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("FITC-A","APC-A")], theta=-pi/4)

# get mhcii gates
gate_mhcii_low <- flowDensity::deGate(
    ff, channel="FITC-A", 
    use.percentile=TRUE, percentile=.001)
gate_mhcii <- flowDensity::deGate(ffs, channel="FITC-A")

# get cd5 gate
tTF <- flowCore::exprs(ffs)[,"FITC-A"] < gate_mhcii
fft <- ff
flowCore::exprs(fft) <- flowCore::exprs(ff)[tTF,, drop=FALSE]

gate_cd5 <- flowDensity::deGate(
    fft, channel="APC-A", 
    use.upper=TRUE, upper=FALSE)

# register gates into gs
fd_t <- flowDensity::flowDensity(
    fft, channels=c("FITC-A", "APC-A"), 
    position=c(TRUE,TRUE), gates=c(gate_mhcii_low, gate_cd5))
gate <- flowCore::polygonGate(.gate=fd_t@filter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="tcells", parent="CD161-")
flowWorkspace::recompute(gs)

fd_nott <- flowDensity::notSubFrame(
    ff, channels=c("FITC-A", "APC-A"), 
    position="logical", gates="missing", fd_t@filter)
gate <- flowCore::polygonGate(.gate=fd_nott@filter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="not_tcells", parent="CD161-")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "tcells")
# flowWorkspace::gs_pop_remove(gs, "not_tcells")

# plot
png(paste0(gateplot_dir, "/11_CD161n_tcells.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"FITC-A"])
plot(d, ylab="", axes=FALSE, main="CD161- > not/tcells")
abline(v=gate_mhcii_low, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"APC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd5, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("FITC-A", "APC-A"), main="")
lines(fd_t@filter)
abline(h=gate_cd5, v=gate_mhcii_low, lty="dashed", col="red")
graphics.off()


## 2.12 gating not tcells > cDC/bcell ####
#cDC: conventional dendritic cells
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_tcells"), Class="list")[[1]] )

# rotate cells
ffs <- ff
flowCore::exprs(ffs)[,c("PE-Cy7-A","BV786-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("PE-Cy7-A","BV786-A")], theta=-pi/4)

# <your gates/register/plot here>
# cell populations you will gate: "cDC", "bcells"


## 2.13 gating cDC > cDC CD11b+/- ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "cDC"), Class="list")[[1]] )

# get upper/lower limit gates
gate_cd11b_high <- flowDensity::deGate(
    ff, channel="BV510-A", 
    use.percentile=TRUE, percentile=0.999999)
gate_mhcii_high <- flowDensity::deGate(
    ff, channel="FITC-A", 
    use.percentile=TRUE, percentile=0.999999)
gate_cd11b_low <- flowDensity::deGate(
    ff, channel="BV510-A", 
    use.upper=TRUE, upper=FALSE)
gate_mhcii_low <- flowDensity::deGate(
    ff, channel="FITC-A", 
    use.upper=TRUE, upper=FALSE)

# get threshold gates
# <your gates/register/plot here>
# cell populations you will gate: "cDC_cd11b+", "cDC_cd11b-"


## 2.14 gating bcells > b1/b2 ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "bcells"), Class="list")[[1]] )

# rotate cells
ffs <- ff
flowCore::exprs(ffs)[,c("APC-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("APC-A","PE-A")], theta=-pi/12)

# <your gates/register/plot here>


## 2.15 gating b2bcells > preB/MZB/folB ####
# CD21: high > MZ (marginal zone) > fol (follicular) > pre > low
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "b2bcells"), Class="list")[[1]] )

# get threshold gates (frame)
gate_cd23_low <- flowDensity::deGate(
    ff, channel="BV421-A", use.upper=TRUE, upper=FALSE)
allTF <- flowCore::exprs(ff)[,"BV421-A"] >= gate_cd23_low
flowCore::exprs(ff) <- flowCore::exprs(ff)[allTF,, drop=FALSE]

# rotate, get middle gate, tighten
ffs <- ff
flowCore::exprs(ffs)[,c("BV421-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("BV421-A","PE-A")], theta=-pi/8)

# <your gates here>

# rotate, but bottom vertical gate, get preB
flowCore::exprs(ffs)[,c("BV421-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ffs)[,c("BV421-A","PE-A")], theta=-pi/8)

# <your gates/register/plot here>


## save gating set! ####
flowWorkspace::save_gs(gs, path=paste0(res_dir, "/gs"))
# gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))

# gating tree
pdf(paste0(gateplot_dir, "/tree.pdf"))
flowWorkspace::plot(gs)
graphics.off()

# you can also save the gatingset as a flowjo workspace
# you will need to install package CytoML:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("CytoML")
CytoML::gatingset_to_flowjo(gs, outFile=paste0(res_dir, "/gs.wsp"))



## TRY: practice problems ####

# 1. complete this gating script according to the given gating strategy. 
#    Search through the script for "# <your"... for where you willl need to 
#    fill the blanks. Read through the other gatings for ideas and guidance.
# Gating strategy: https://docs.google.com/presentation/d/1dUamRDWtN6cuZXaUN1I1baDUopPhhufKTGXZAXOSP3o/edit?usp=sharing

# 2. the sangerP2.fcs file was downloaded from 
#    http://flowrepository.org/id/FR-FCM-ZYPK
#    go to this link, click "download" and select any one random .fcs file to 
#    download. all files there follow the same gating strategy.
#    change the path to the .fcs file to your new file.
# File download link: https://drive.google.com/file/d/1PpSM93GTj9zejVDZzD89_k3sx7Lc-TQl/view?usp=sharing

# question: does your existing script work at gating your new file?
# question: can you modify your script such that it works on both files?
# note: if your file looks obnoxiously different from our file, it might be an
#       outlier, in that case, download another one :P!

# 3. can you make one plotting function that creates the 
# scatterplots/density plots to replace all previous plotting code?

# 4. in the solution script (gating_full.R), there are variables whose name 
# does not make sense, despite the script working with no problem. 
# Can you find them? Did they hinder your understanding? 

# 4. if you want more practice, try creating a gating script for any of the 
#    files in http://flowrepository.org/public_experiment_representations/1146
#    used in the article "An immune clock of pregnancy" (doi/10.1126/sciimmunol.aan2946),
#    whose gating strategy can be found on page 5 of its supplementary material:
#    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5701281/bin/STM-02-eaan2946-s001.pdf
