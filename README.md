## File Description

This repository is used for depositing the analysis scripts, derived data/results and corresponding figures for the paper "A Phylogeny-based Test of Mediation Effect in Microbiome". 

Here are brief descriptions for each folder, detailed explanation are commented within code and described in the paper:

* a folder named **Analyis**, which contains the real data anlysis R scripts. 
  - 0.prepareData.R
  - 1.runModel.R
  - 2.generatePlot.R
  - 3.generateTable.R
  - \*.utility.R

* a folder named **Simulation**, which contains the simulation shell scripts and R scripts.
  - **continuous_JC**: DAGman workflow for continuous outcome
  - **binary_JC**: DAGman workflow for binary outcome
  - By submitting the ./medtest.dag to High-Throughput Condor, you are expected to get a csv file in the ./\*_JC/RESULT folder.

* a folder named **Data**, which contains raw data and derived data.

* a folder named **Figs**, which contains the figures in the main text.

* a folder named **SuppFigs**, which contains the figures in the supporting information.

## R Package

The proposed method **PhyloMed** in the paper is implemented in the R package **miMediation**, which is publicly available at https://tangzheng1.github.io/tanglab/software

``` r
# Within R
install.packages("miMediation_0.1.tar.gz", repos=NULL, type = "source")
# Using RStudio and devtools
install_github("tangzheng1/miMediation")
# Using Terminal
R CMD INSTALL ./miMediation_0.1.tar.gz
```

## Contact

* Qilin (Kirin) Hong - qhong8@wisc.edu
