## gate fcs file
## author: alice yue
## input: raw fcs file
## output: gated flowworkspace

# load packages
library("flowCore")
library("PeacoQC")
library("flowWorkspace")
library("flowDensity")

# directory to save results in
res_dir <- "/home/alice/res"
dir.create(res_dir)

gateplot_dir <- paste0(res_dir, "/gs_plots")
dir.create(gateplot_dir)

# path to raw fcs file
fcs_path <- paste0(res_dir, "/sangerP2.fcs")

## load fcs file
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
    position=c(NA,TRUE), gates=c(NA, gate_cd11b)) #upper half
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

# get upper limit gates
gate_cd5_high <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A", 
    use.percentile=TRUE, percentile=0.9999999)

# get threshold gates
temp <- flowDensity::flowDensity(
    ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
    position=c(NA,TRUE), gates=c(NA, gate_cd11b)) #upper half
gate_cd5 <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A")

# register gates into gs
gate <- flowCore::rectangleGate(
    filterId="monocytes", 
    "Alexa Fluor 700-A"=c(gate_cd5, gate_cd5_high), # include all
    "BV510-A"=c(gate_cd11b, gate_cd11b_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="not_granulocytes")
flowWorkspace::recompute(gs)

fd_mono <- flowDensity::flowDensity(
    ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
    position=c(TRUE,TRUE), gates=c(gate_cd5, gate_cd11b))
fd_notmono <- flowDensity::notSubFrame(
    ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
    position="logical", gates="missing", fd_mono@filter)
gate <- flowCore::polygonGate(.gate=fd_notmono@filter)
node <- flowWorkspace::gs_pop_add(
    gs, gate, name="not_monocytes", parent="not_granulocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "monocytes")
# flowWorkspace::gs_pop_remove(gs, "not_monocytes")

# plot
png(paste0(gateplot_dir, "/05_notgranuloctyes_monocytes.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(temp))[,"Alexa Fluor 700-A"], na.rm=TRUE)
plot(d, ylab="", axes=FALSE, main="not granulocytes > not/monocytes",
     xlim=range(flowCore::exprs(ff)[,"Alexa Fluor 700-A"], na.rm=TRUE))
abline(v=gate_cd5, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd11b, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("Alexa Fluor 700-A", "BV510-A"), main="")
abline(v=gate_cd5, h=gate_cd11b, lty="dashed", col="red")
graphics.off()


## 2.6 gating not monocytes > not/eosoniphil ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_monocytes"), Class="list")[[1]] )

# get upper limit gates
gate_ssch_high <- flowDensity::deGate(
    ff, channel="SSC-H", 
    use.percentile=TRUE, percentile=0.999999) 

# get threshold gates
gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")

temp <- flowDensity::flowDensity(
    ff, channels=c("BV510-A", "SSC-H"), 
    position=c(TRUE,FALSE), gates=c(gate_cd11b, gate_ssch_high))
gate_ssch <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="SSC-H")

# register gates into gs
gate <- flowCore::rectangleGate(
    filterId="eosinophils", 
    "BV510-A"=c(gate_cd11b, gate_cd11b_high), # include all
    "SSC-H"=c(gate_ssch, gate_ssch_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="not_monocytes")
flowWorkspace::recompute(gs)

fd_eosi <- flowDensity::flowDensity(
    ff, channels=c("BV510-A", "SSC-H"), 
    position=c(TRUE,TRUE), gates=c(gate_cd11b, gate_ssch))
fd_noteosi <- flowDensity::notSubFrame(
    ff, channels=c("BV510-A", "SSC-H"), 
    position="logical", gates="missing", fd_eosi@filter)
gate <- flowCore::polygonGate(.gate=fd_noteosi@filter)
node <- flowWorkspace::gs_pop_add(
    gs, gate, name="not_eosinophils", parent="not_monocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "eosinophils")
# flowWorkspace::gs_pop_remove(gs, "not_eosinophils")

# plot
png(paste0(gateplot_dir, "/06_notmonocytes_eosinophils.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d, ylab="", axes=FALSE, main="not monocytes > not/eosinophils")
abline(v=gate_cd11b, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"SSC-H"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_ssch, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("BV510-A", "SSC-H"), main="")
abline(v=gate_cd11b, h=gate_ssch, lty="dashed", col="red")
graphics.off()


## 2.7 gating not eosoniphil > CD11b+/- ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "not_eosinophils"), Class="list")[[1]] )

# get upper/lower limit gates
gate_cd161_high <- flowDensity::deGate(
    ff, channel="BV650-A", 
    use.percentile=TRUE, percentile=0.999999) 
gate_cd19_low <- flowDensity::deGate(
    ff, channel="PE-Cy7-A", 
    use.percentile=TRUE, percentile=0.000001) 

# get loose threshold gates
gate_cd161 <- flowDensity::deGate(ff, channel="BV650-A")
gate_cd19 <- flowDensity::deGate(ff, channel="PE-Cy7-A")

# tighten threshold gates
temp <- flowDensity::flowDensity(
    ff, channels=c("BV650-A", "PE-Cy7-A"), 
    position=c(TRUE,FALSE), gates=c(gate_cd161, gate_cd19))
gate_cd161 <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="BV650-A")
gate_cd19 <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="PE-Cy7-A")

# register gates into gs
gate <- flowCore::rectangleGate(
    filterId="CD161+", 
    "BV650-A"=c(gate_cd161, gate_cd161_high), # include all
    "PE-Cy7-A"=c(gate_cd19_low, gate_cd19))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="not_eosinophils")
flowWorkspace::recompute(gs)

fd_cd161 <- flowDensity::flowDensity(
    ff, channels=c("BV650-A", "PE-Cy7-A"), 
    position=c(TRUE,FALSE), gates=c(gate_cd161, gate_cd19))
fd_notcd161 <- flowDensity::notSubFrame(
    ff, channels=c("BV650-A", "PE-Cy7-A"), 
    position="logical", gates="missing", fd_cd161@filter)
gate <- flowCore::polygonGate(.gate=fd_notcd161@filter)
node <- flowWorkspace::gs_pop_add(
    gs, gate, name="CD161-", parent="not_eosinophils")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "eosinophil")
# flowWorkspace::gs_pop_remove(gs, "not_eosinophils")

# plot
png(paste0(gateplot_dir, "/07_noteosinophils_cd161p.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"BV650-A"])
plot(d, ylab="", axes=FALSE, main="not eosinophils > CD161+/-")
abline(v=gate_cd161, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"PE-Cy7-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd19, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("BV650-A", "PE-Cy7-A"), main="")
abline(v=gate_cd161, h=gate_cd19, lty="dashed", col="red")
graphics.off()


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

# get threshold gates
gate_ly6c <- flowDensity::deGate(ff, channel="Alexa Fluor 700-A") 
gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")

# register gates into gs
gate <- flowCore::quadGate(
    filterId=c("NKT"),
    "Alexa Fluor 700-A"=gate_ly6c,
    "BV510-A"=gate_cd11b)
node <- flowWorkspace::gs_pop_add(gs, gate, parent="NKT")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700-CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700+CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700-CD11b BV510+")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700+CD11b BV510+")

# plot
png(paste0(gateplot_dir, "/10_NKT_immatureLy6C.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"Alexa Fluor 700-A"])
plot(d, ylab="", axes=FALSE, main="NKT > im/mature Ly6C+/-", 
     xlim=range(flowCore::exprs(ff)[,"Alexa Fluor 700-A"], na.rm=TRUE))
abline(v=gate_ly6c, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd11b, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("Alexa Fluor 700-A", "BV510-A"), main="")
abline(v=gate_ly6c, h=gate_cd11b, lty="dashed", col="red")
graphics.off()


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

# get threshold gates for bcells
gate_cd19 <- flowDensity::deGate(
    ff, channel="PE-Cy7-A")
gate_cd19_slant <- flowDensity::deGate(
    ffs, channel="PE-Cy7-A", 
    use.upper=TRUE, upper=FALSE)

# get threshold gates for cDC
gate_cd19_low <- flowDensity::deGate(
    ff, channel="PE-Cy7-A", 
    use.percentile=TRUE, percentile=0.0000001)
gate_cd11c_high <- flowDensity::deGate(
    ff, channel="BV786-A", 
    use.percentile=TRUE, percentile=0.9999999)

temp <- flowDensity::flowDensity(
    ff, channels=c("PE-Cy7-A", "BV786-A"), 
    position=c(TRUE,NA), gates=c(gate_cd19, NA))
gate_cd11c <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="BV786-A")

# register gates into gs
bTF <- flowCore::exprs(ffs)[,"PE-Cy7-A"] >= gate_cd19_slant &
    flowCore::exprs(ff)[,"PE-Cy7-A"] >= gate_cd19
bchull <- chull(flowCore::exprs(ff)[bTF, c("PE-Cy7-A","BV786-A"), drop=FALSE])
bfilter <- flowCore::exprs(ff)[
    which(bTF)[c(bchull, bchull[1])], c("PE-Cy7-A","BV786-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=bfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="bcells", parent="not_tcells")
flowWorkspace::recompute(gs)

gate <- flowCore::rectangleGate(
    filterId="cDC", 
    "PE-Cy7-A"=c(gate_cd19_low, gate_cd19), # include all
    "BV786-A"=c(gate_cd11c, gate_cd11c_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="not_tcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "monocytes")
# flowWorkspace::gs_pop_remove(gs, "not_monocytes")

# plot
png(paste0(gateplot_dir, "/12_nottcells_cDCbcells.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"PE-Cy7-A"], na.rm=TRUE)
plot(d, ylab="", axes=FALSE, main="not tcells > cDC/bcells")
abline(v=gate_cd19, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"BV786-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_cd11c, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("PE-Cy7-A", "BV786-A"), main="")
lines(bfilter)
abline(v=gate_cd19, h=gate_cd11c, lty="dashed", col="red")
graphics.off()


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
gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")

# register gates into gs
gate <- flowCore::rectangleGate(
    filterId="cDC_cd11b+", 
    "BV510-A"=c(gate_cd11b, gate_cd11b_high), # include all
    "FITC-A"=c(gate_mhcii_low, gate_mhcii_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="cDC")
flowWorkspace::recompute(gs)

gate <- flowCore::rectangleGate(
    filterId="cDC_cd11b-", 
    "BV510-A"=c(gate_cd11b_low, gate_cd11b), # include all
    "FITC-A"=c(gate_mhcii_low, gate_mhcii_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="cDC")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "cDC_cd11b+")
# flowWorkspace::gs_pop_remove(gs, "cDC_cd11b-")

# plot
png(paste0(gateplot_dir, "/13_cDC_CD11b.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"BV510-A"])
plot(d, ylab="", axes=FALSE, main="cDC > CD11b+/-")
abline(v=c(gate_cd11b_low, gate_cd11b), lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"FITC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=gate_mhcii_low, lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("BV510-A", "FITC-A"), main="")
abline(v=c(gate_cd11b_low, gate_cd11b), 
       h=gate_mhcii_low, lty="dashed", col="red")
graphics.off()


## 2.14 gating bcells > b1/b2 ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "bcells"), Class="list")[[1]] )

# rotate cells
ffs <- ff
flowCore::exprs(ffs)[,c("APC-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("APC-A","PE-A")], theta=-pi/12)

# get threshold gates
gate_cd5 <- flowDensity::deGate(
    ffs, channel="APC-A", 
    use.upper=TRUE, upper=TRUE)

# register gates into gs
b1TF <- flowCore::exprs(ffs)[,"APC-A"] >= gate_cd5

b1chull <- chull(flowCore::exprs(ff)[b1TF, c("APC-A","PE-A"), drop=FALSE])
b1filter <- flowCore::exprs(ff)[
    which(b1TF)[c(b1chull, b1chull[1])], c("APC-A","PE-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=b1filter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="b1bcells", parent="bcells")
flowWorkspace::recompute(gs)

b2chull <- chull(flowCore::exprs(ff)[!b1TF, c("APC-A","PE-A"), drop=FALSE])
b2filter <- flowCore::exprs(ff)[
    which(!b1TF)[c(b2chull, b2chull[1])], c("APC-A","PE-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=b2filter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="b2bcells", parent="bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "b1bcells")
# flowWorkspace::gs_pop_remove(gs, "b2bcells")

# plot
png(paste0(gateplot_dir, "/14_bcells_b1b2.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"APC-A"])
plot(d, ylab="", axes=FALSE, main="bcells > b1/2 bcells")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"PE-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("APC-A", "PE-A"), main="")
lines(b1filter)
graphics.off()


## 2.15 gating b2bcells > preB/MZB/folB ####
# CD21: high > MZ (marginal zone) > fol (follicular) > pre > low
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "b2bcells"), Class="list")[[1]] )

# get threshold gates (frame)
gate_cd23_low <- flowDensity::deGate(
    ff, channel="BV421-A", 
    use.upper=TRUE, upper=FALSE)
allTF <- flowCore::exprs(ff)[,"BV421-A"] >= gate_cd23_low
flowCore::exprs(ff) <- flowCore::exprs(ff)[allTF,, drop=FALSE]

# rotate, get middle gate, tighten
ffs <- ff
flowCore::exprs(ffs)[,c("BV421-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ff)[,c("BV421-A","PE-A")], theta=-pi/8)

gate_cd21_slanth_temp <- flowDensity::deGate(ffs, channel="PE-A")

ffs_ <- ffs
botTF <- flowCore::exprs(ffs)[,"PE-A"] >= gate_cd21_slanth_temp
flowCore::exprs(ffs_) <- flowCore::exprs(ffs)[botTF,c("BV421-A","PE-A")]
gate_cd21_slanth <- flowDensity::deGate(
    ffs_, channel="PE-A", use.upper=TRUE, upper=FALSE)

topTF <- flowCore::exprs(ffs)[,"PE-A"] >= gate_cd21_slanth

# rotate, but bottom vertical gate, get preB
flowCore::exprs(ffs)[,c("BV421-A","PE-A")] <- 
    rotate_data(flowCore::exprs(ffs)[,c("BV421-A","PE-A")], theta=-pi/8)
ffs_ <- ffs
flowCore::exprs(ffs_) <- flowCore::exprs(ffs)[!topTF,c("BV421-A","PE-A")]
gate_cd23_slantv <- flowDensity::deGate(
    ffs_, channel="BV421-A", 
    use.upper=TRUE, upper=TRUE)
gate_cd23_slanth <- flowDensity::deGate(ffs, channel="BV421-A")

# register gates into gs
mzbTF <- flowCore::exprs(ffs)[,"BV421-A"] <= gate_cd23_slanth
folbTF <- !mzbTF & topTF
prebTF <- !mzbTF & !topTF & flowCore::exprs(ffs)[,"BV421-A"] <= gate_cd23_slantv

mzbchull <- chull(flowCore::exprs(ff)[mzbTF, c("BV421-A","PE-A"), drop=FALSE])
mzbfilter <- flowCore::exprs(ff)[
    which(mzbTF)[c(mzbchull, mzbchull[1])], c("BV421-A","PE-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=mzbfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="MZbcells", parent="b2bcells")
flowWorkspace::recompute(gs)

folbchull <- chull(flowCore::exprs(ff)[folbTF, c("BV421-A","PE-A"), drop=FALSE])
folbfilter <- flowCore::exprs(ff)[
    which(folbTF)[c(folbchull, folbchull[1])], c("BV421-A","PE-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=folbfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="folbcells", parent="b2bcells")
flowWorkspace::recompute(gs)

prebchull <- chull(flowCore::exprs(ff)[prebTF, c("BV421-A","PE-A"), drop=FALSE])
prebfilter <- flowCore::exprs(ff)[
    which(prebTF)[c(prebchull, prebchull[1])], c("BV421-A","PE-A"), drop=FALSE]
gate <- flowCore::polygonGate(.gate=prebfilter)
node <- flowWorkspace::gs_pop_add(gs, gate, name="prebcells", parent="b2bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "MZbcells")
# flowWorkspace::gs_pop_remove(gs, "folbcells")
# flowWorkspace::gs_pop_remove(gs, "prebcells")

# plot
png(paste0(gateplot_dir, "/15_b2bcells_MZfolprebcells.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"BV421-A"])
plot(d, ylab="", axes=FALSE, main="b2bcells > MZ/fol/pre bcells")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"PE-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("BV421-A", "PE-A"), main="")
lines(mzbfilter)
lines(folbfilter)
lines(prebfilter)
graphics.off()


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
