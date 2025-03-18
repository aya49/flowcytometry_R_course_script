## obtain stats and plots from gated fcs (gating set)
## author: alice yue
## input: gated gating set
## output: stats and plots

# load packages
library("flowCore")
library("flowWorkspace")

# directory to save and obtain results in
res_dir <- "/home/alice/res"
gateplot_dir <- paste0(res_dir, "/gs_plots")

# load gating set
gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))
cpops <- flowWorkspace::gs_get_leaf_nodes(gs)

# function to get parent cell population
parent_cpop <- function(child_cpop) {
    a <- strsplit(child_cpop,"/")[[1]]
    return( paste0(a[-length(a)], collapse="/") )
}

# function to get leaf cell population
leaf_cpop <- function(full_cpop) {
    a <- strsplit(full_cpop,"/")[[1]]
    return( a[length(a)] )
}


## plots ####

# gating tree
pdf(paste0(gateplot_dir, "/tree.pdf"))
flowWorkspace::plot(gs)
graphics.off()

# # all gatings
# gag <- ggcyto::autoplot(gs[[1]], bins=100)
# ggplot2::ggsave(filename=paste0(gateplot_dir, "/all_gates.png"), plot=gag)
# 
# # one gating
# ggcyto::autoplot(gs, cpops[1], bins=100)


## stats ####

# get the cell count/percent
statsm <- data.frame(
    cell_population=cpops,
    cell_count=sapply(cpops, function(cp) 
        flowWorkspace::gh_pop_get_count(gs, cp)) )
total_count <- flowWorkspace::gh_pop_get_count(gs, "root") # root = all cells
statsm[["cell_percent"]] <- statsm$cell_count / total_count

# our fcs file
statsm

# get the median/mean/SD/CV FI for each cell population and marker
statsl <- list()
for (cp in cpops) {
    ff <- flowWorkspace::cytoset_to_flowSet(
        flowWorkspace::gs_pop_get_data(gs, cp))[[1]]
    marker <- flowWorkspace::markernames(ff)
    mefi <- apply(flowCore::exprs(ff)[,names(marker),drop=FALSE], 2, median)
    mfi <- apply(flowCore::exprs(ff)[,names(marker),drop=FALSE], 2, mean)
    sdfi <- apply(flowCore::exprs(ff)[,names(marker),drop=FALSE], 2, sd)
    # # for clusters subset-ed from the FULL flow frame, see clustering.R:
    # # replace for loop with: for (cp in unique(fc)) {
    # mefi <- apply(flowCore::exprs(ff)[fc==cluster_id,names(marker),drop=FALSE], 2, median)
    # mfi <- apply(flowCore::exprs(ff)[fc==cluster_id,names(marker),drop=FALSE], 2, mean)
    # sdfi <- apply(flowCore::exprs(ff)[fc==cluster_id,names(marker),drop=FALSE], 2, sd)
    
    
    cvfi <- sdfi/mfi
    
    statsl[[cp]] <- data.frame(mefi, mfi, sdfi, cvfi)
}
# you can save this list as a csv --- but note the column headers!
write.table(statsl, file=paste0(res_dir, "/cell_population_stats.csv"), sep=",")

# let's look at the mfi's for each cell population across each marker :3
mfi_covs <- data.frame()
for (cpop in cpops) {
    parent <- leaf_cpop(parent_cpop(cpop))
    mfi_covs <- rbind(mfi_covs, data.frame(
        parent=parent,
        cpop=paste0(parent, "/", leaf_cpop(cpop)), 
        marker=flowCore::markernames(ff), 
        median_fi=statsl[[cpop]][["mefi"]],
        variance_fi=statsl[[cpop]][["cvfi"]]))
}
cpmfi <- ggplot2::ggplot(mfi_covs, ggplot2::aes(
    x=marker, y=median_fi, colour=1/variance_fi)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group=cpop)) +
    ggplot2::facet_wrap(cpop~., scales="free_y") +
    ggplot2::scale_colour_continuous(type="viridis") + 
    ggplot2::theme(axis.text.x=ggplot2::element_text(
        angle=90, vjust = 0.5, hjust=1)) # turn x axis text sideways
cpmfi

## analyzing cell population stats across multiple files ####
# since we don't have multiple fcs files, we will generate data to emulate
# different types of experiments.
# typically, we analyze cell populations based on their cell count percent

# generate 100 fcs files (and their stats), (100 = 50 controls + 50 experiment)
statsm_multiple <- data.frame()
for (x in c(1:100)) {
    statsm_ <- statsm
    statsm_$file <- paste0("file_",x)
    if (x <= 50) {
        statsm_[["class"]] <- "control"
    } else {
        statsm_[["class"]] <- "experiment"
    }
    for (i in seq_len(nrow(statsm))) {
        ov <- statsm[i, "cell_count"]
        statsm_[i, "cell_count"] <- runif(1, min=ov*0.8, max=ov*1.2)
        if (i==12 & x > 50) {
            statsm_[i, "cell_count"] <- runif(1, min=ov*0.4, max=ov*0.6)
        }
    }
    statsm_multiple <- rbind(statsm_multiple, statsm_)
}
statsm_multiple[["cell_percent"]] <- statsm_multiple[["cell_count"]] / total_count
statsm_multiple[["cp"]] <- sapply(stringr::str_split(
    statsm_multiple[["cell_population"]], "/"), function(x) {
        if (length(x) < 2) {
            return(paste0(x, collapse="/"))
        }
        return(paste0(x[(length(x)-1):length(x)], collapse="/"))
    })
cps <- statsm_multiple[["cp"]][c(1:length(cpops))]

# plot boxplots
gp <- ggplot2::ggplot(statsm_multiple, ggplot2::aes(
    x=class, 
    y=cell_percent, 
    fill=class)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(cp~., scales="free_y") + 
    ggplot2::theme(axis.text.x=ggplot2::element_text(
        angle=90, vjust = 0.5, hjust=1)) # turn x axis text sideways
gp

# calculate p-values and log fold change between cell_percent of 
# the same cell population across files of different classes
p <- lfc <- c()
controli <- statsm_multiple[["class"]]=="control"
for (cp in cpops) {
    cpi <- statsm_multiple[["cell_population"]] == cp
    controlv <- statsm_multiple[controli & cpi,"cell_percent"]
    experimentv <- statsm_multiple[!controli & cpi,"cell_percent"]
    p[cp] <- t.test(controlv, experimentv)[["p.value"]]
    lfc[cp] <- log(experimentv / controlv)
}
plfc <- data.frame(cell_population=cps, p_value=p, log_fold_change=lfc, row.names=cps)
write.table(plfc, file=paste0(res_dir, "/differential_analysis.csv"), sep=",")

# volcano plot
plot(plfc[["log_fold_change"]], -log(plfc[["p_value"]]), xlab="ln fold change", ylab="-ln p-value")

# cell populations with lowest p-values
head(plfc[order(plfc[["p_value"]]),])


## TRY: practice problem ####

# 1. in 03_preprocess_loop.csv, there was a comment on how to save matrices
#    and data.frames as .csv files. select some of the data.frame or matrix
#    that you think is important to keep for record above and save them as .csv
#    fils.

# 2. another common approach to comparing cell populations is to get the 
#    percentage of cells in a cell population NOT over the total number of 
#    cells but over the number of cells in its PARENT cell population.
#    For example, cell population a/b/c has 10 cells, its
#    parent cell poulation a/b has 20 cells (and a has, let's say 30 cells), 
#    then a/b/c's percentage would be 10/20 = 0.5.
#    
# your task: replace the "cell_percent" column in statsm_multiple with the 
#            above said percentage for each cell poulation across all samples.
#            calculate p-values and log fold change again to see if there are
#            any changes to the results! (there shouldn't be, most of the time)
# 
# hint: use the parent_cpop function to identify parent cell populations.
parent_cpop("/live/lymphocytes/not_granulocytes/not_monocytes")

# <your code here>
# e.g. statsm_multiple[["cell_percent"]] <- ...

# 3. instead of comparing the cell count or the cell percentage, it is also
#    common to compare the "mfi" or the mean fluorescent intensity for each
#    marker across samples for the same cell population.
#    
# your task: for any one marker and any one cell population of your choosing, 
#            generate another table, like statsm_multiple, where
#            you simulate mfi values for 100 samples, 
