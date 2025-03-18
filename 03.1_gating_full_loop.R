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

# path to all raw fcs files in directory
fcs_paths <- list.files(res_dir, pattern="fcs$", full.names=TRUE)

## load fcs file
fl <- lapply(fcs_paths, function(fp) {
    f <- flowCore::read.FCS(fp)

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
    return(f)
})
gc()
fn <- floor(log10(length(fl))) + 1 # number of digits in number of files
pad_fn <- function(x, fn) {
    stringr::str_pad(x, fn, side="left", pad=0)
}


## 2 GATING: cell population identification based on a gating strategy ####

## function to rotate 2D frame
## input: 2D matrix and angle
## output: rotated 2D matrix
rotate_data <- function(data, theta=pi/2 - atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])) {
    data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}

# initialize gating set (containing 1 file)
fs <- flowCore::flowSet(fl)
# fs <- flowCore::flowSet(list(f1, f2, f3)) # can add more files if needed
gs <- flowWorkspace::GatingSet(fs)


## 2.1 gating all > singlets ####

# get gate value
gates <- lapply(seq_len(length(gs)), function(fi) {
    # f <- flowWorkspace::cytoframe_to_flowFrame(
    #     as(flowWorkspace::gs_pop_get_data(gs, "live"), Class="list")[[fi]])
    f <- fl[[fi]]
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
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_01_all_singlets.png"))
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
    
    # plot and save as png
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_02_singlets_live.png"))
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
    
    return(c(gate_live, gate_singlets_low, gate_singlets))
})
# apply gates to make a list of gates
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="live", 
        "APC-Cy7-A"=c(-1, x[1]), # include all
        "SSC-W"=c(x[2], x[3]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="root")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "live")

# clean memory by removing variables
gc()


## 2.3 gating live > lymphocytes ####

fli <- as(flowWorkspace::gs_pop_get_data(gs, "live"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_03_live_lymphocytes.png"))
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
    
    return(c(gate_fsca, gate_fsca_high, gate_ssca_high))
})
rm(fli)
gc()
    
# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
    filterId="lymphocytes", 
    "FSC-A"=c(x[1], x[2]), # include all
    "SSC-A"=c(0, x[3]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="live")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "lymphocytes")


## 2.4 gating lymphocytes > not/granulocytes ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "lymphocytes"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    # get filter
    fd_gran <- flowDensity::flowDensity(
        ff, channels=c("APC-A", "BV510-A"), 
        position=c(TRUE,TRUE), gates=c(gate_cd5, gate_cd11b))
    fd_notgran <- flowDensity::notSubFrame(
        ff, channels=c("APC-A", "BV510-A"), 
        position="logical", gates="missing", fd_gran@filter)
    
    # get upper limit gates
    gate_ly6c_high <- flowDensity::deGate(
        flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A", 
        use.percentile=TRUE, percentile=0.9999999)
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_04_lymphocytes_granulocytes.png"))
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
    
    return(list(c(
        gate_cd5, gate_cd5_high, gate_cd11b, gate_cd11b_high, gate_ly6c_high),
        fd_notgran@filter))
})
rm(fli)
gc()

# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="granulocytes", 
        "APC-A"=c(x[[1]][1], x[[1]][2]), # include all
        "BV510-A"=c(x[[1]][3], x[[1]][4]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="lymphocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "granulocytes")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="not_granulocytes", parent="lymphocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "not_granulocytes")


## 2.5 gating not granulocytes > not/monocytes ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "not_granulocytes"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    gate_cd11b <- gates[[fi]][[1]][3]
    gate_cd11b_high <- gates[[fi]][[1]][4]
    gate_ly6c_high <- gates[[fi]][[1]][5]
    
    # get threshold gates
    temp <- flowDensity::flowDensity(
        ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
        position=c(NA,TRUE), gates=c(NA, gate_cd11b)) #upper half
    gate_ly6c <- flowDensity::deGate(
        flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A")
    
    fd_mono <- flowDensity::flowDensity(
        ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
        position=c(TRUE,TRUE), gates=c(gate_ly6c, gate_cd11b))
    fd_notmono <- flowDensity::notSubFrame(
        ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
        position="logical", gates="missing", fd_mono@filter)
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_05_notgranuloctyes_monocytes.png"))
    layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))
    
    par(mar=c(0,5,3,1))
    d <- density(flowCore::exprs(
        flowDensity::getflowFrame(temp))[,"Alexa Fluor 700-A"], na.rm=TRUE)
    plot(d, ylab="", axes=FALSE, main="not granulocytes > not/monocytes",
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

    return(list(
        c(gate_ly6c, gate_ly6c_high, gate_cd11b, gate_cd11b_high),
        fd_notmono@filter))
})
rm(fli)
gc()

# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="monocytes", 
        "Alexa Fluor 700-A"=c(x[[1]][1], x[[1]][2]), # include all
        "BV510-A"=c(x[[1]][3], x[[1]][4]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="not_granulocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "monocytes")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(
    gs, rg, name="not_monocytes", parent="not_granulocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "not_monocytes")


## 2.6 gating not monocytes > not/eosoniphil ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "not_monocytes"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    gate_cd11b_high <- gates[[fi]][[1]][4]
    
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
    
    # get filter
    fd_eosi <- flowDensity::flowDensity(
        ff, channels=c("BV510-A", "SSC-H"), 
        position=c(TRUE,TRUE), gates=c(gate_cd11b, gate_ssch))
    fd_noteosi <- flowDensity::notSubFrame(
        ff, channels=c("BV510-A", "SSC-H"), 
        position="logical", gates="missing", fd_eosi@filter)
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_06_notmonocytes_eosinophils.png"))
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

    return(list(c(gate_cd11b, gate_cd11b_high, gate_ssch, gate_ssch_high),
                fd_noteosi@filter))
})
rm(fli)
gc()

# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="eosinophils", 
        "BV510-A"=c(x[[1]][1], x[[1]][2]), # include all
        "SSC-H"=c(x[[1]][3], x[[1]][4]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="not_monocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "eosinophils")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(
    gs, rg, name="not_eosinophils", parent="not_monocytes")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "not_eosinophils")


## 2.7 gating not eosoniphil > CD11b+/- ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "not_eosinophils"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    fd_cd161 <- flowDensity::flowDensity(
        ff, channels=c("BV650-A", "PE-Cy7-A"), 
        position=c(TRUE,FALSE), gates=c(gate_cd161, gate_cd19))
    fd_notcd161 <- flowDensity::notSubFrame(
        ff, channels=c("BV650-A", "PE-Cy7-A"), 
        position="logical", gates="missing", fd_cd161@filter)
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_07_noteosinophils_cd161p.png"))
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
    
    return(list(c(gate_cd161, gate_cd161_high, gate_cd19_low, gate_cd19),
                fd_notcd161@filter))
})
rm(fli)
gc()

# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="CD161+", 
        "BV650-A"=c(x[[1]][1], x[[1]][2]), # include all
        "PE-Cy7-A"=c(x[[1]][3], x[[1]][4]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="not_eosinophils")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "CD161+")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(
    gs, rg, name="CD161-", parent="not_eosinophils")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "CD161-")


## 2.8 gating CD161+ > NK/T ####
# NK: natural killer
fli <- as(flowWorkspace::gs_pop_get_data(gs, "CD161+"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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

    nktchull <- chull(flowCore::exprs(ff)[!nkTF, c("APC-A","BV650-A"), drop=FALSE])
    nktfilter <- flowCore::exprs(ff)[
        which(!nkTF)[c(nktchull, nktchull[1])], c("APC-A","BV650-A"), drop=FALSE]
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_08_CD161p_NKT.png"))
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
    
    return(list(nkfilter, nktfilter))
})
rm(fli)
gc()

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[1]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="NK", parent="CD161+")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NK")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="NKT", parent="CD161+")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NKT")


## 2.9 gating NK > NK im/mature Ly6C+/- ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "NK"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
    # get loose threshold gates
    gate_ly6c <- flowDensity::deGate(ff, channel="Alexa Fluor 700-A") 
    gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")
    
    # tighten threshold gates
    temp <- flowDensity::flowDensity(
        ff, channels=c("Alexa Fluor 700-A", "BV510-A"), 
        position=c(FALSE,NA), gates=c(gate_ly6c, NA))
    gate_ly6c <- flowDensity::deGate(
        flowDensity::getflowFrame(temp), channel="Alexa Fluor 700-A")
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_09_NK_immatureLy6C.png"))
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
    
    return(c(gate_ly6c, gate_cd11b))
})
rm(fli)
gc()

# register gates into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::quadGate(
        filterId=c("NK"),
        "Alexa Fluor 700-A"=x[1],
        "BV510-A"=x[2])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="NK")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700-CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700+CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700-CD11b BV510+")
# flowWorkspace::gs_pop_remove(gs, "NK/Ly6C AF700+CD11b BV510+")


## 2.10 gating NKT > NKT im/mature Ly6C+/- ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "NKT"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
    # get threshold gates
    gate_ly6c <- flowDensity::deGate(ff, channel="Alexa Fluor 700-A") 
    gate_cd11b <- flowDensity::deGate(ff, channel="BV510-A")
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_10_NKT_immatureLy6C.png"))
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
    
    return(c(gate_ly6c, gate_cd11b))
})
rm(fli)
gc()

# register gates into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::quadGate(
        filterId=c("NKT"),
        "Alexa Fluor 700-A"=x[1],
        "BV510-A"=x[2])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="NKT")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700-CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700+CD11b BV510-")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700-CD11b BV510+")
# flowWorkspace::gs_pop_remove(gs, "NKT/Ly6C AF700+CD11b BV510+")


## 2.11 gating CD161- > not/tcells ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "CD161-"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    fd_t <- flowDensity::flowDensity(
        fft, channels=c("FITC-A", "APC-A"), 
        position=c(TRUE,TRUE), gates=c(gate_mhcii_low, gate_cd5))
    
    fd_nott <- flowDensity::notSubFrame(
        ff, channels=c("FITC-A", "APC-A"), 
        position="logical", gates="missing", fd_t@filter)
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_11_CD161n_tcells.png"))
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
    
    return(list(fd_t@filter, fd_nott@filter))
})
rm(fli)
gc()

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[1]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="tcells", parent="CD161-")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "tcells")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="not_tcells", parent="CD161-")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "not_tcells")


## 2.12 gating not tcells > cDC/bcell ####
#cDC: conventional dendritic cells
fli <- as(flowWorkspace::gs_pop_get_data(gs, "not_tcells"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_12_nottcells_cDCbcells.png"))
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

    return(list(c(gate_cd19_low, gate_cd19, gate_cd11c, gate_cd11c_high),
                bfilter))
})
rm(fli)
gc()

# register gate into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="cDC", 
        "PE-Cy7-A"=c(x[[1]][1], x[[1]][2]), # include all
        "BV786-A"=c(x[[1]][3], x[[1]][4]))
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="not_tcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "cDC")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="bcells", parent="not_tcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "bcells")


## 2.13 gating cDC > cDC CD11b+/- ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "cDC"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])
    
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
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_13_cDC_CD11b.png"))
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
    
    return(c(gate_cd11b, gate_cd11b_high, gate_cd11b_low,
             gate_mhcii_low, gate_mhcii_high))
})
rm(fli)
gc()


# register gates into gs
rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="cDC_cd11b+", 
        "BV510-A"=c(x[1], x[2]), # include all
        "FITC-A"=c(x[4], x[5]))
    return(gate)
})
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="cDC")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "cDC_cd11b+")

rg <- lapply(gates, function(x) {
    gate <- flowCore::rectangleGate(
        filterId="cDC_cd11b-", 
        "BV510-A"=c(x[3], x[1]), # include all
        "FITC-A"=c(x[4], x[5]))
    return(gate)
})
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, parent="cDC")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "cDC_cd11b-")


## 2.14 gating bcells > b1/b2 ####
fli <- as(flowWorkspace::gs_pop_get_data(gs, "bcells"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])

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

    b2chull <- chull(flowCore::exprs(ff)[!b1TF, c("APC-A","PE-A"), drop=FALSE])
    b2filter <- flowCore::exprs(ff)[
        which(!b1TF)[c(b2chull, b2chull[1])], c("APC-A","PE-A"), drop=FALSE]
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_14_bcells_b1b2.png"))
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

    return(list(b1filter, b2filter))
})
rm(fli)
gc()

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[1]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="b1bcells", parent="bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "b1bcells")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="b2bcells", parent="bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "b2bcells")


## 2.15 gating b2bcells > preB/MZB/folB ####
# CD21: high > MZ (marginal zone) > fol (follicular) > pre > low
fli <- as(flowWorkspace::gs_pop_get_data(gs, "b2bcells"), Class="list")
gates <- lapply(seq_len(length(gs)), function(fi) {
    ff <- flowWorkspace::cytoframe_to_flowFrame(fli[[fi]])

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
    
    mzbchull <- chull(flowCore::exprs(ff)[
        mzbTF, c("BV421-A","PE-A"), drop=FALSE])
    mzbfilter <- flowCore::exprs(ff)[
        which(mzbTF)[c(mzbchull, mzbchull[1])], 
        c("BV421-A","PE-A"), drop=FALSE]

    folbchull <- chull(flowCore::exprs(ff)[
        folbTF, c("BV421-A","PE-A"), drop=FALSE])
    folbfilter <- flowCore::exprs(ff)[
        which(folbTF)[c(folbchull, folbchull[1])], 
        c("BV421-A","PE-A"), drop=FALSE]

    prebchull <- chull(flowCore::exprs(ff)[
        prebTF, c("BV421-A","PE-A"), drop=FALSE])
    prebfilter <- flowCore::exprs(ff)[
        which(prebTF)[c(prebchull, prebchull[1])], 
        c("BV421-A","PE-A"), drop=FALSE]
    
    # plot
    png(paste0(gateplot_dir, "/", pad_fn(fi, fn), "_15_b2bcells_MZfolprebcells.png"))
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
    
    return(list(mzbfilter, folbfilter, prebfilter))
})
rm(fli)
gc()

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[1]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="MZbcells", parent="b2bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "MZbcells")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[2]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="folbcells", parent="b2bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "folbcells")

rg <- lapply(gates, function(x) {
    gate <- flowCore::polygonGate(.gate=x[[3]])
    return(gate)
})
# give the gates the same name as your files
names(rg) <- sampleNames(fs)
# register your gate
node <- flowWorkspace::add(gs, rg, name="prebcells", parent="b2bcells")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "prebcells")


## save gating set! ####
flowWorkspace::save_gs(gs, path=paste0(res_dir, "/gs"))
# gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))

# gating tree
pdf(paste0(gateplot_dir, "/tree.pdf"))
flowWorkspace::plot(gs)
graphics.off()

# ## plot everything ####
# # all gatings
# # BiocManager::install("ggcyto")
# gag <- ggcyto::autoplot(gs[[1]], bins=100)
# ggplot2::ggsave(filename=paste0(gateplot_dir, "/all_gates.png"), plot=gag)

# # one gating
# ggcyto::autoplot(gs, "live", bins=100) # replace "live" with the cell population of interest


# you can also save the gatingset as a flowjo workspace
# you will need to install package CytoML:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("CytoML")
# you will also need Docker!
CytoML::gatingset_to_flowjo(gs, outFile=paste0(res_dir, "/gs.wsp"))

