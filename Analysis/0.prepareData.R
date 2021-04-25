rm(list = ls())
setwd("~/Project/QH0020/medtest/PhyloMed/Analysis/")
library(phyloseq)
library(TreeTools)
library(ape)
source("0a.utility.R")

phy_pheno = readRDS("../Data/Rawdata/Application_STAT.Rds")
# cecal full 100 tips
cecal_phy = subset_samples(phy_pheno, Location == "cecal")
cecal_phy = prune_samples(sample_sums(cecal_phy) >= 1000, cecal_phy)
tmp=otu_table(cecal_phy)
tmp=t(tmp@.Data)
tmp=tmp/rowSums(tmp)
sum(sort(colMeans(tmp), decreasing = TRUE)[1:100]) # 0.8015
cecal_phy.top100 = get_top_taxa(cecal_phy, n = 100, relative = T, discard_other = T)
cecal.tree = phy_tree(cecal_phy.top100)
tree = prepareTree(cecal.tree)
meta = sample_data(cecal_phy.top100)
Trt = ifelse(meta$Treatment == "C", 0, 1)
outcome = meta$pFat
M = t(otu_table(cecal_phy.top100))
M@.Data = M@.Data[,tree$tip.label]
data.cecal.full = list(treatment = Trt, mediators = M@.Data,
                  outcome = outcome, tree = tree)
save(data.cecal.full, file = "../Data/Deriveddata/cecal.full.rda")

# cecal merge 45 tips
cecal_phy.merge = tip_glom(cecal_phy.top100, h = 0.2)
cecal.tree = phy_tree(cecal_phy.merge)
tree = prepareTree(cecal.tree)
meta = sample_data(cecal_phy.merge)
Trt = ifelse(meta$Treatment == "C", 0, 1)
outcome = meta$pFat
M = t(otu_table(cecal_phy.merge))
M@.Data = M@.Data[,tree$tip.label]
data.cecal.merge = list(treatment = Trt, mediators = M@.Data,
                  outcome = outcome, tree = tree)
save(data.cecal.merge, file = "../Data/Deriveddata/cecal.merge.rda")

# fecal full 100 tips
fecal_phy = subset_samples(phy_pheno, Location == "fecal")
fecal_phy = prune_samples(sample_sums(fecal_phy) >= 1000, fecal_phy)
tmp=otu_table(fecal_phy)
tmp=t(tmp@.Data)
tmp=tmp/rowSums(tmp)
sum(sort(colMeans(tmp), decreasing = TRUE)[1:100]) # 0.8652
fecal_phy.top100 = get_top_taxa(fecal_phy, n = 100, relative = T, discard_other = T)
fecal.tree = phy_tree(fecal_phy.top100)
tree = prepareTree(fecal.tree)
meta = sample_data(fecal_phy.top100)
Trt = ifelse(meta$Treatment == "C", 0, 1)
outcome = meta$pFat
M = t(otu_table(fecal_phy.top100))
M@.Data = M@.Data[,tree$tip.label]
data.fecal.full = list(treatment = Trt, mediators = M@.Data,
                  outcome = outcome, tree = tree)
save(data.fecal.full, file = "../Data/Deriveddata/fecal.full.rda")

#### fecal merge 39 tips
fecal_phy.merge = tip_glom(fecal_phy.top100, h = 0.2)
fecal.tree = phy_tree(fecal_phy.merge)
tree = prepareTree(fecal.tree)
meta = sample_data(fecal_phy.merge)
Trt = ifelse(meta$Treatment == "C", 0, 1)
outcome = meta$pFat
M = t(otu_table(fecal_phy.merge))
M@.Data = M@.Data[,tree$tip.label]
data.fecal.merge = list(treatment = Trt, mediators = M@.Data,
                        outcome = outcome, tree = tree)
save(data.fecal.merge, file = "../Data/Deriveddata/fecal.merge.rda")

