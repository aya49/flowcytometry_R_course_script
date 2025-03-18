## preprocessing fcs file
## author: alice yue
## input: raw fcs file
## output: preprocessed fcs file

# load packages
library("flowCore")
library("PeacoQC")

# change where you want to save PeacoQC cleaning plots!
res_dir <- "/home/alice/res"
dir.create(res_dir)

# path to fcs file 
fcs_path <- system.file("extdata", "111.fcs", package="PeacoQC")
# fcs_path <- "path/to/sangerP2.fcs"

# load fcs file
f <- flowCore::read.FCS(fcs_path)


## flowCore::exprs(f) is a matrix!
dim(flowCore::exprs(f))
head(flowCore::exprs(f))
tail(flowCore::exprs(f))
colnames(flowCore::exprs(f))
f@parameters@data # description for columns
# flowCore::markernames(f) # markers only, no morphology/time columns

## plotting a flow frame
channel <- "FSC-A"
flowDensity::plotDens(f, channels=c("Time", channel)) # this function doesn't take matrices as input, only flow frames
gate_margin <- 150000

# let's isolate only the cells with lower values
good_cells_ind <- flowCore::exprs(f)[,channel] < gate_margin
f1 <- f # replicate flow frame first, so we still have all original cells
flowCore::exprs(f1) <- flowCore::exprs(f1)[good_cells_ind, ]

# plot again
flowDensity::plotDens(f1, channels=c("Time", channel))

# kind of hard to compare because the range is different, let's fix that
channel_range <- range(flowCore::exprs(f)[,channel])
flowDensity::plotDens(f1, channels=c("Time", channel), ylim=channel_range)


## 1.1 compensate ####
# compensation not necessary for mass spectrometry
spillover <- flowCore::keyword(f)$SPILL
# spillover <- read.csv("path", header=TRUE, row.names=TRUE)
# flowCore::keyword(f)$SPILL <- spillover
f <- flowCore::compensate(f, spillover=spillover)

## 1.2 logicle transform ####

# let's look at the before transformation plot
flowDensity::plotDens(f, channels=c("G780-A", "G710-A"))

# transform
transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
f <- flowCore::transform(f, transformList)

transformList <- flowCore::estimateLogicle(f1, channels=colnames(spillover))
f1t <- flowCore::transform(f1, transformList)

# and the before compensation plots
flowDensity::plotDens(f1t, channels=c("G780-A", "G710-A"))

# let's look at the after compensation + transformation plot
flowDensity::plotDens(f, channels=c("G780-A", "G710-A"))


## 1.3 cleaning; see res_dir for plot ####
## parameter: Mean Absolute Deviation (MAD) distance (decrease = less strict)
## parameter: IsolationTree (IT) (decrease = more strict)

# let's look at the Time vs FSC-A plot to see the flow of cells
flowDensity::plotDens(f, channels=c("Time", "FSC-A"))

# change the channels based on what channels you want to base your cleaning!
channels <- c(1, 3, 5:14, 18, 21)
time_range <- range(flowCore::exprs(f)[,"Time"])
# RemoveMargins not necessary for mass spectrometry
fmr <- PeacoQC::RemoveMargins(f, channels=channels, output="full")
pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=channels,
                        plot=TRUE, save_fcs=FALSE, report=FALSE,
                        # IT_limit=0.6, remove_zeros=TRUE, time_units=50000 # change some parameters for mass spectrometry
                        output_directory=res_dir)


# final preprocessed fcs!
fc <- pQC[["FinalFF"]] 
flowDensity::plotDens(fc, channels=c("Time", channel), ylim=channel_range, xlim=time_range)


# clean memory by removing variables
rm("fmr") # remove variable
rm("pQC")
gc() # remove removed variables from memory
