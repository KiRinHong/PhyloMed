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
  - **ProposedModel_Main**: The treatment-associated taxa set and the outcome-associated taxa set are completely overlapped. The taxa are clustered on the tree.
  - **PartialOverlap**: The treatment-associated taxa set and the outcome-associated taxa set are partially overlapped. The taxa are clustered on the tree.
  - **IncreaseNumOfCausalTaxa**: The treatment-associated taxa set and the outcome-associated taxa set are completely overlapped. The taxa are randomly scattered with large or small effect size.
  - **DiffModels**: Evaluate empirical type I error of global mediation test when different subcomposition models were used.
  - **PropOfNull**: Evaluate bias and standard deviation of the estimated probabilities of three sub-nulls using two differernt estimation approaches.
  - **SmallEffectSize**: Evaluate empirical type I error of different global mediation tests under two scenarios of mediation effect size.
  - Download a copy of R (R361.tar.gz) and create a portable copy of necessary R packages (See [tutorial](https://chtc.cs.wisc.edu/r-jobs)). By submitting the ./medtest.dag to High-Throughput Condor, you are expected to get a csv file in the corresponding *RESULT* folder.
  - generateTxt.R: Generate the inpuut.txt and comb_sim.txt in *GETSIM* and *COMBSIM* folder, respectively. 

* a folder named **Data**, which contains raw data and derived data.

* a folder named **Figs**, which contains the figures in the main text.

* a folder named **SuppFigs**, which contains the figures in the supporting information.

## R Package

The proposed method **PhyloMed** in the paper is implemented in the R package **miMediation**, which is publicly available at https://github.com/KiRinHong/miMediation

## Contact

* Qilin (Kirin) Hong - qhong8@wisc.edu
