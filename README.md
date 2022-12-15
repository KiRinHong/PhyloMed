## File Description

This repository is used for depositing the analysis scripts, derived data/results and corresponding figures for the paper "PhyloMed: a phylogeny-based test of mediation effect in microbiome". 

Here are brief descriptions for each folder, detailed explanations are commented within code and described in the paper:

* a folder named **Analysis**, which contains the real data analysis R scripts. 
  - 0.prepareData.R
  - 1.runModel.R
  - 2.generatePlot.R
  - 3.runPairLeafAnalysis.R
  - 4.sensitivityAnalysis.R
  - 5.generateTab.R
  - \*.utility.R

* a folder named **Data**, which contains raw data and derived data.

* a folder named **Figs**, which contains the figures in the main text and supplementary information.

* a folder named **Tabs**, which contains the tables in the main text and supplmentary information.

* a folder named **Simulation**, which contains the simulation shell scripts and R scripts.
  - **ProposedModel_Main**: The mediating OTUs are clustered on the tree. The tree is constructed from top 100 OTUs.
  - **RareTaxa**: The mediating OTUs are clustered on the tree. The tree is constructed from all 819 OTUs.
  - **IncreaseNumOfCausalTaxa**: The mediating OTUs are randomly scattered with large or small effect size.
  - **DiffModels**: Evaluate empirical type I error of global mediation test when different subcomposition models were used.
  - **PropOfNull**: Evaluate bias and standard deviation of the estimated probabilities of three sub-nulls using different estimation approaches.
  - **PseudoCount**: Evaluate empirical type I error and power of global mediation test when different pseudo counts were used.
  - **IncreaseEffectSize**: Evaluate empirical FDR when the mediating OTUs are clustered on the tree with relatively small effect size.
  - Download a copy of R (R413.tar.gz) and create a portable copy of necessary R packages (See [tutorial](https://chtc.cs.wisc.edu/uw-research-computing/r-jobs.html)). By submitting the ./medtest.dag to High-Throughput Condor, you are expected to get a csv file in the corresponding *RESULT* folder.
  - generateData.R: Generate the zeeviD_\*.Rdata, which is used as a basis in the simulation.
  - generateTxt.R: Generate the input.txt and comb_sim.txt in *GETSIM* and *COMBSIM* folder, respectively.
  - generateSeqDepth.R: Generate the box plot of sequence depth between treatment and control group.
  - runDACT.R: Evaluate the performance of other mediating test method (DACT).
  - zeeviD_*.Rdata: The phyloseq-object used in the simulation.

## R Package

The proposed method **PhyloMed** in the paper is implemented in the R package **miMediation**, which is publicly available at https://github.com/KiRinHong/miMediation

## Contact

* Qilin (Kirin) Hong - qhong8@wisc.edu
