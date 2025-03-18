# Flow cytometry analysis on R/bioconductor

This course is a 3 hour course instructed by Alice Yue, hosted by Physalia.

- **Date and time**: 2025-03-20 @ 11am CEST
- **Instructor**: Alice Yue
- **Level**: Beginner

**Overview**: Flow cytometry is a gold standard in the analysis of immune cells for clinical diagnosis and research. This course introduces what flow cytometry is and why we use it to analyze cell population composition in biological samples. We will learn about best practices for how to analyze flow cytometry data using R/Bioconductor. Said practices include: preprocessing of the data (compensation, transformation, and quality control), multi-dimensional cell population identification via clustering, and visualizing the results in 2D. These tools can be applied to all types of cytometry including flow, mass, and spectral.

**Target audience**: This course is created for anyone interested in analyzing biological samples with single-cell flow cytometry. Background in flow cytometry and R/Bioconductor is necessary.

**Format**: Follow along with the numbered scripts in this repository to go through a full analysis of a single (and optionally multiple) flow cytometry sample.

## Learning outcomes

- Be able to set up the infrastructure for and write basic data analytic scripts in R.
- Describe and execute each step in the flow cytometry data analytics pipeline in R/Bioconductor.
- Be comfortable with interpreting and eliciting conclusions from the results of the flow cytometry data analytics pipeline.

## Pre-requisites

Before starting,

- Install [R](https://www.r-project.org/) and [Rstudio](https://www.rstudio.com/categories/rstudio-ide/): https://learnr-examples.shinyapps.io/ex-setup-r/
- Install required R packages using script: [01_packages](01_packages.R)

Download [sangerP2.fcs](https://drive.google.com/file/d/1PpSM93GTj9zejVDZzD89_k3sx7Lc-TQl/view?usp=sharing). We will be gating this file following [this gating strategy](https://docs.google.com/presentation/d/1dUamRDWtN6cuZXaUN1I1baDUopPhhufKTGXZAXOSP3o/edit?usp=sharing).