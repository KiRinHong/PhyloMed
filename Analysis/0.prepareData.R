rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")
source("0a.utility.R")
library(phyloseq)
library(ape)
library(TreeTools)
library(dplyr)
library(tidyverse)

##### Mice Data, cecal microbiome #####
# https://www.nature.com/articles/nature11400.pdf
phy_pheno <- readRDS("../Data/Rawdata/Application_STAT.Rds")
phy <- subset_samples(phy_pheno, Location == "cecal")
phy <- prune_samples(sample_sums(phy) >= 1000, phy)
sum(taxa_sums(phy) != 0) # 4973
tax.clean <- as.data.frame(tax_table(phy))
tax.clean <- tax.clean[,c("Phylum","Class","Order","Family","Genus","Species")]
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="Other"] <- ""
tax.clean[] <- lapply(tax.clean, function(x) gsub("\"", "", x))
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    phylum <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:5] <- phylum
  } else if (tax.clean[i,3] == ""){
    class <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:5] <- class
  } else if (tax.clean[i,4] == ""){
    order <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:5] <- order
  } else if (tax.clean[i,5] == ""){
    genus <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:6] <- genus
  }else if (tax.clean[i,6] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}
tax_table(phy) <- tax_table(as.matrix(tax.clean))
tmp.full <- otu_table(phy)
tmp.full <- t(tmp.full@.Data)
tmp.full <- tmp.full/rowSums(tmp.full)
sum(colMeans(tmp.full>0) > 0.2) #757
phy.nonrare <- get_nonrare_taxa(phy, threshold = 0.2, discard_other = T)
tmp <- otu_table(phy.nonrare)
tmp <- t(tmp@.Data)
tmp <- tmp/rowSums(tmp)
tmp.id <- names(sort(colMeans(tmp), decreasing = TRUE)[1:100])
sum(colMeans(tmp.full)[tmp.id]) #0.8015
phy.top <- get_top_taxa(phy.nonrare, n = 100, relative = T, discard_other = T)
tree <- phy_tree(phy.top)
# data.cecal <- list()
# data.cecal$tree <- tree
# data.cecal$treatment <- ifelse(meta$Treatment == "C", 0, 1)
# data.cecal$mediators <- M@.Data
# data.cecal$outcome <- meta$pFat
# save(data.cecal, file = "../../miMediation/data/data.cecal.rda")
tree <- prepareTree(tree)
meta <- sample_data(phy.top)
meta$Treatment <- ifelse(meta$Treatment == "C", 0, 1)
M <- t(otu_table(phy.top))
M@.Data <- M@.Data[,tree$tip.label]
cecal.top <- list(phy = phy.top,
                  meta = meta, mediators = M@.Data, 
                  taxonomy = tax_table(phy.top),
                  tree = tree)
save(cecal.top, file = "../Data/Deriveddata/Cecal.filter20top100.rda")

# for (seed in 1:20) {
#   phy.rarefy <- rarefy_even_depth(phy, rngseed = seed+12, verbose = FALSE)
#   tmp.top <- prune_taxa(taxa_names(phy.top), phy.rarefy)
#   tree <- phy_tree(tmp.top)
#   tree <- prepareTree(tree)
#   meta <- sample_data(tmp.top)
#   meta$Treatment <- ifelse(meta$Treatment == "C", 0, 1)
#   M <- t(otu_table(tmp.top))
#   M@.Data <- M@.Data[,tree$tip.label]
#   cecal.rarefy <- list(phy = tmp.top,
#                        meta = meta, mediators = M@.Data, 
#                        taxonomy = tax_table(tmp.top),
#                        tree = tree)
#   save(cecal.rarefy, file = paste0("../Data/Deriveddata/Cecal.RRF", seed, ".filter20top100.rda"))
# }

##### COMBO data #####
load("../Data/Rawdata/combo.RData")
load("../Data/Rawdata/conf.Rdata")
id <- rownames(conf)
otu.tab <- otu_table(otu.tab[,id], taxa_are_rows = TRUE)
tax.clean <- as.data.frame(otu.names.mat)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="Other"] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    phylum <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:5] <- phylum
  } else if (tax.clean[i,3] == ""){
    class <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:5] <- class
  } else if (tax.clean[i,4] == ""){
    order <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:5] <- order
  } else if (tax.clean[i,5] == ""){
    tax.clean$Genus[i] <- paste("Unclassified ",tax.clean$Family[i], sep = " ")
  }
}
tax.tab <- tax_table(as.matrix(tax.clean))
meta.df <- data.frame(calor = ffq.adj[id, "calor"], fat = ffq.adj[id,"tfat"], bmi = bmi.c[id])
rownames(meta.df) <- id
meta.data <- sample_data(meta.df)
phy <- phyloseq(otu.tab, 
                tax.tab, 
                meta.data, 
                tree.rooted)
phy <- prune_samples(sample_sums(phy) >= 1000, phy)
mod1 <- lm(sample_sums(phy) ~ meta.data$fat)
summary(mod1)
cor.test(sample_sums(phy), meta.data$fat)
mod2 <- lm(sample_sums(phy) ~ meta.data$bmi)
summary(mod2)
tmp.full <- otu_table(phy)
tmp.full <- t(tmp.full@.Data)
tmp.full <- tmp.full/rowSums(tmp.full)
sum(colMeans(tmp.full>0) > 0.2) # 395 in total
tmp.id <- names(which(colMeans(tmp.full>0) > 0.2))
sum(colMeans(tmp.full)[tmp.id]) # 0.7921801
phy.nonrare <- get_nonrare_taxa(phy, threshold = 0.2, discard_other = T)
tree <- phy_tree(phy.nonrare)
tree <- prepareTree(tree)
meta <- sample_data(phy.nonrare)
cor.test(sample_sums(phy.nonrare), meta$fat)
cor.test(sample_sums(phy.nonrare), meta$calor)
cor.test(sample_sums(phy.nonrare), meta$bmi)

M <- t(otu_table(phy.nonrare))
M@.Data <- M@.Data[,tree$tip.label]
combo.filter <- list(phy = phy.nonrare,
                     meta = meta, mediators = M@.Data, 
                     taxonomy = tax_table(phy.nonrare),
                     tree = tree)
save(combo.filter, file = "../Data/Deriveddata/COMBO.filter20.rda")

# for (seed in 1:20) {
#   phy.rarefy <- rarefy_even_depth(phy, rngseed = seed+12, verbose = FALSE)
#   tmp.nonrare <- prune_taxa(taxa_names(phy.nonrare), phy.rarefy)
#   tree <- phy_tree(tmp.nonrare)
#   tree <- prepareTree(tree)
#   meta <- sample_data(tmp.nonrare)
#   M <- t(otu_table(tmp.nonrare))
#   M@.Data <- M@.Data[,tree$tip.label]
#   combo.rarefy <- list(phy = tmp.nonrare,
#                        meta = meta, mediators = M@.Data, 
#                        taxonomy = tax_table(tmp.nonrare),
#                        tree = tree)
#   save(combo.rarefy, file = paste0("../Data/Deriveddata/COMBO.RRF", seed, ".filter20.rda"))
# }
