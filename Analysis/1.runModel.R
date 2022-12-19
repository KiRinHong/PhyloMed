rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")
# load source code for PhyloMed
source("1a.PhyloMed_utility.R")
# load source code for MedTest
source("1b.MeTest_utility.R")
# load source code for MODIMA
source("1c.MODIMA_utility.R")
library(phyloseq)
library(vegan)
library(GUniFrac)
library(energy)
library(fdrtool)
library(SKAT)
library(harmonicmeanp)
# library(ccmm)
library(matrixStats)
library(MASS)
library(LDM)
library(HIMA)

runSeveralMeds <- function(input_data, exposure, covariate=NULL, outcome, 
                           fdr.alpha = 0.1, perm.prec = 0.05, seed = 84, runCMM = FALSE,
                           pseudoCount = 0.5){
  # input_data = cecal.top; exposure = "Treatment"; covariate = NULL; outcome = "pFat"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 84; runCMM = FALSE; pseudoCount = 0.5
  # input_data = combo.filter; exposure = "fat"; covariate = "calor"; outcome = "bmi"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 123; runCMM = TRUE
  Trt = input_data$meta[[exposure]]
  M = input_data$mediators; M2 = M + pseudoCount
  Y = input_data$meta[[outcome]]
  
  tree = input_data$tree
  gp = numeric(27)
  names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.clr.MedTest", "dist.omn.MedTest",
                "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.clr.MODIMA", "dist.bonf.MODIMA",
                "LDM.med.global", # cannot find Trt in the local environment
                "Prod.Simes.PhyloMed.asym", "Prod.Fisher.PhyloMed.asym", "Prod.HMP.PhyloMed.asym",
                "Minus.Simes.PhyloMed.asym", "Minus.Fisher.PhyloMed.asym", "Minus.HMP.PhyloMed.asym",
                "Prod.Simes.PhyloMed.perm", "Prod.Fisher.PhyloMed.perm", "Prod.HMP.PhyloMed.perm",
                "Minus.Simes.PhyloMed.perm", "Minus.Fisher.PhyloMed.perm", "Minus.HMP.PhyloMed.perm",
                "CMM.norm", "CMM.boot")
  nullprop.est = numeric(12)
  names(nullprop.est) = c("Prod.H00.asym","Prod.H10.asym","Prod.H01.asym",
                          "Minus.H00.asym","Minus.H10.asym","Minus.H01.asym",
                          "Prod.H00.perm","Prod.H10.perm","Prod.H01.perm",
                          "Minus.H00.perm","Minus.H10.perm","Minus.H01.perm")
  
  # perform distance-based MedTest
  set.seed(seed)
  unifracs = GUniFrac(M, tree)$unifracs 
  m.list = list(BC = as.matrix(vegdist(M, method = "bray")),
                JAC = as.matrix(vegdist(M, 'jaccard', binary = TRUE)),
                UniFrac=unifracs[, , c('d_UW')], # unweighted unifrac distance doesn't work because M2 doesn't contain zero
                WUniFrac=unifracs[, , c('d_1')],
                CLR = as.matrix(vegdist(M2, method = "aitchison")))
  rslt.dist = MedOmniTest(x = Trt, y = Y, m.list = m.list, z = NULL, nperm = 2000)
  gp[1] = rslt.dist$margPs[1]
  gp[2] = rslt.dist$margPs[2]
  gp[3] = rslt.dist$margPs[3]
  gp[4] = rslt.dist$margPs[4]
  gp[5] = rslt.dist$margPs[5]
  gp[6] = rslt.dist$permP
  cat("distance-based MedTest done\n")
  
  dist_Trt = dist(Trt)
  dist_outcome = dist(Y)
  gp[7] = modima(dist_Trt, m.list[[1]],  dist_outcome, nrep=2000)$p.value
  gp[8] = modima(dist_Trt, m.list[[2]],  dist_outcome, nrep=2000)$p.value
  gp[9] = modima(dist_Trt, m.list[[3]],  dist_outcome, nrep=2000)$p.value
  gp[10] = modima(dist_Trt, m.list[[4]],  dist_outcome, nrep=2000)$p.value
  gp[11] = modima(dist_Trt, m.list[[5]],  dist_outcome, nrep=2000)$p.value
  gp[12] = min(c(gp[7:11]*5, 1))
  cat("distance-based MODIMA done\n")
  
  return(gp[1:12])
  
  set.seed(seed)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      covariate = input_data$meta[[covariate]]
    }else{
      covariate = as.matrix(sample_data(input_data$meta[,covariate]))
    }
  }
  rslt.phylomed = phylomed(Trt, M, Y, tree, covariates = covariate, fdr.alpha = fdr.alpha, perm.prec = perm.prec, 
                           verbose = T, pseudoCount = pseudoCount)
  gp[14:19] = rslt.phylomed$PhyloMed.A$global.pval
  nullprop.est[1:6] = c(rslt.phylomed$PhyloMed.A$null.prop.prod,
                        rslt.phylomed$PhyloMed.A$null.prop.minus)
  node.pval.asym.jsmix.prod = rslt.phylomed$PhyloMed.A$node.pval.jsmix.prod
  node.pval.asym.jsmix.minus = rslt.phylomed$PhyloMed.A$node.pval.jsmix.minus
  node.pval.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js
  node.pval.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel
  gp[20:25] = rslt.phylomed$PhyloMed.P$global.pval
  nullprop.est[7:12] = c(rslt.phylomed$PhyloMed.P$null.prop.prod,
                         rslt.phylomed$PhyloMed.P$null.prop.minus)
  node.pval.perm.jsmix.prod = rslt.phylomed$PhyloMed.P$node.pval.jsmix.prod
  node.pval.perm.jsmix.minus = rslt.phylomed$PhyloMed.P$node.pval.jsmix.minus
  node.pval.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js
  
  pval.alphabeta.mat = rbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.P$pval.alpha,
                             rslt.phylomed$PhyloMed.A$pval.beta, rslt.phylomed$PhyloMed.P$pval.beta)
  rownames(pval.alphabeta.mat) = c("asym.alpha", "perm.alpha", "asym.beta", "perm.beta")
  
  node.pval.mat = rbind(node.pval.asym.jsmix.prod, node.pval.asym.jsmix.minus,
                        node.pval.asym.js, node.pval.sobel, 
                        node.pval.perm.jsmix.prod, node.pval.perm.jsmix.minus,
                        node.pval.perm.js)
  rownames(node.pval.mat) = c("asym.jsmix.prod", "asym.jsmix.minus","asym.js","sobel",
                              "perm.jsmix.prod", "perm.jsmix.minus", "perm.js")
  
  sig.nodeID.BH = apply(node.pval.mat, 1, .getSigNode.BH, alpha = fdr.alpha)
  sig.nodeID.HMP = apply(node.pval.mat, 1, .getSigNode.HMP, alpha = fdr.alpha, tree = tree)
  sig.nodeID = c(sig.nodeID.BH, sig.nodeID.HMP)
  names(sig.nodeID) = apply(expand.grid(c("asym.jsmix.prod", "asym.jsmix.minus","asym.js","sobel",
                                          "perm.jsmix.prod", "perm.jsmix.minus", "perm.js"), c("BH","HMP")), 1, paste, collapse=".")
  cat("tree-based PhyloMed done\n")
  
  if(runCMM){
    M2 = M + pseudoCount
    M2.P = M2/rowSums(M2)
    ## here perform ccmm
    rslt.ccmm.boot <- ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "bootstrap", n.boot = 2000, sig.level = fdr.alpha)
    rslt.ccmm.norm <- ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "normal")
    gp[26] = 2*(1 - pnorm(abs(rslt.ccmm.norm$TIDE)/sqrt(rslt.ccmm.norm$Var.TIDE)))
    gp[27] = ifelse(rslt.ccmm.boot$TIDE.CI[1] > 0 | rslt.ccmm.boot$TIDE.CI[2] < 0, 0, 1)
    cat("ccmm done\n")
    rslt = list(gp = gp, nullprop.est = nullprop.est, 
                pval.alphabeta = pval.alphabeta.mat, 
                node.pval = node.pval.mat,
                sig.node = sig.nodeID,
                rslt.ccmm.norm = rslt.ccmm.norm,
                rlst.ccmm.boot = rslt.ccmm.boot)
    return(rslt)
  }else{
    gp[26] = NA
    gp[27] = NA
    rslt = list(gp = gp, nullprop.est = nullprop.est, 
                pval.alphabeta = pval.alphabeta.mat, 
                node.pval = node.pval.mat,
                sig.node = sig.nodeID)
    return(rslt)
  }
}
# runMedTestwInt <- function(input_data, exposure, covariate=NULL, outcome, interaction = interaction,
#                            fdr.alpha = 0.1, perm.prec = 0.05, seed = 84){
#   # input_data = cecal.top; exposure = "Treatment"; covariate = NULL; outcome = "pFat"; fdr.alpha = 0.1;
#   # perm.prec = 0.05; seed = 123
#   Trt = input_data$meta[[exposure]]
#   M = input_data$mediators
#   Y = input_data$meta[[outcome]]
#   
#   tree = input_data$tree
#   gp = numeric(6)
#   nullprop.est = numeric(6)
#   names(nullprop.est) = c("H00.asym","H10.asym","H01.asym","H00.perm","H10.perm","H01.perm")
#   
#   set.seed(seed)
#   if(!is.null(covariate)){
#     if(length(covariate) == 1){
#       covariate = input_data$meta[[covariate]]
#     }else{
#       covariate = as.matrix(sample_data(input_data$meta[,covariate]))
#     }
#   }
#   rslt.phylomed = phylomed(Trt, M, Y, tree, covariates = covariate, interaction = interaction,
#                            fdr.alpha = fdr.alpha, perm.prec = perm.prec, verbose = F)
#   
#   gp[1:3] = rslt.phylomed$PhyloMed.A$global.pval
#   nullprop.est[1:3] = rslt.phylomed$PhyloMed.A$null.prop
#   node.pval.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix
#   node.pval.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js
#   node.pval.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel
#   gp[4:6] = rslt.phylomed$PhyloMed.P$global.pval
#   nullprop.est[4:6] = rslt.phylomed$PhyloMed.P$null.prop
#   node.pval.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix
#   node.pval.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js
#   
#   pval.alphabeta.mat = rbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.P$pval.alpha,
#                              rslt.phylomed$PhyloMed.A$pval.beta, rslt.phylomed$PhyloMed.P$pval.beta)
#   rownames(pval.alphabeta.mat) = c("asym.alpha", "perm.alpha", "asym.beta", "perm.beta")
#   
#   node.pval.mat = rbind(node.pval.asym.jsmix, node.pval.asym.js, node.pval.sobel, 
#                         node.pval.perm.jsmix, node.pval.perm.js)
#   rownames(node.pval.mat) = c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js")
#   
#   sig.nodeID.BH = apply(node.pval.mat, 1, .getSigNode.BH, alpha = fdr.alpha)
#   sig.nodeID.HMP = apply(node.pval.mat, 1, .getSigNode.HMP, alpha = fdr.alpha, tree = tree)
#   sig.nodeID = c(sig.nodeID.BH, sig.nodeID.HMP)
#   names(sig.nodeID) = apply(expand.grid(c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js"), c("BH","HMP")), 1, paste, collapse=".")
#   cat("tree-based PhyloMed done\n")
#   
#   rslt = list(gp = gp, nullprop.est = nullprop.est, 
#               pval.alphabeta = pval.alphabeta.mat, 
#               node.pval = node.pval.mat,
#               sig.node = sig.nodeID)
#   return(rslt)
#   
# }
load("../Data/Deriveddata/rslt.runModel.rda")
gp.mat = matrix(NA, nrow=39, ncol=6, dimnames = list(c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.clr.MedTest", "dist.omn.MedTest",
                                                       "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.clr.MODIMA", "dist.bonf.MODIMA",
                                                       "LDM.med.global",
                                                       "Prod.Simes.PhyloMed.asym", "Prod.Fisher.PhyloMed.asym", "Prod.HMP.PhyloMed.asym",
                                                       "Minus.Simes.PhyloMed.asym", "Minus.Fisher.PhyloMed.asym", "Minus.HMP.PhyloMed.asym",
                                                       "Prod.Simes.PhyloMed.perm", "Prod.Fisher.PhyloMed.perm", "Prod.HMP.PhyloMed.perm",
                                                       "Minus.Simes.PhyloMed.perm", "Minus.Fisher.PhyloMed.perm", "Minus.HMP.PhyloMed.perm",
                                                       "CMM.norm", "CMM.boot",
                                                       "Prod.H00.asym","Prod.H10.asym","Prod.H01.asym","Prod.H00.perm","Prod.H10.perm","Prod.H01.perm",
                                                       "Minus.H00.asym","Minus.H10.asym","Minus.H01.asym","Minus.H00.perm","Minus.H10.perm","Minus.H01.perm"),
                                                     c("Cecal.top100.psc01", "Cecal.top100.psc05", "Cecal.top100.psc10", 
                                                       "COMBO.filter20.psc01", "COMBO.filter20.psc05", "COMBO.filter20.psc10")))

load("../Data/Deriveddata/Cecal.filter20top100.rda")
rslt.cecal.top100.psc01 = runSeveralMeds(cecal.top, exposure = "Treatment", outcome = "pFat", runCMM = FALSE, pseudoCount = 0.1)
rslt.cecal.top100.psc05 = runSeveralMeds(cecal.top, exposure = "Treatment", outcome = "pFat", runCMM = FALSE, pseudoCount = 0.5)
rslt.cecal.top100.psc10 = runSeveralMeds(cecal.top, exposure = "Treatment", outcome = "pFat", runCMM = FALSE, pseudoCount = 1)

# rslt.cecal.top100.Int = runMedTestwInt(cecal.top, exposure = "Treatment", outcome = "pFat", interaction = TRUE)
# gp.mat[c(11:16, 19:24),2] = c(rslt.cecal.top100.Int$gp, rslt.cecal.top100.Int$nullprop.est)

load("../Data/Deriveddata/COMBO.filter20.rda")
rslt.combo.filter20.psc01 = runSeveralMeds(combo.filter, exposure = "fat", covariate = "calor", outcome = "bmi", runCMM = FALSE, pseudoCount = 0.1, seed = 123)
rslt.combo.filter20.psc05 = runSeveralMeds(combo.filter, exposure = "fat", covariate = "calor", outcome = "bmi", runCMM = FALSE, pseudoCount = 0.5, seed = 123)
rslt.combo.filter20.psc10 = runSeveralMeds(combo.filter, exposure = "fat", covariate = "calor", outcome = "bmi", runCMM = FALSE, pseudoCount = 1, seed = 123)

set.seed(123)
Trt = cecal.top$meta[["Treatment"]]; M = cecal.top$mediators; Y = cecal.top$meta[["pFat"]]
rslt.ldm = ldm(M ~ Trt + Y, test.mediation = TRUE)
rslt.cecal.top100.psc01$gp[13] = rslt.ldm$med.p.global.omni
rslt.cecal.top100.psc05$gp[13] = rslt.ldm$med.p.global.omni
rslt.cecal.top100.psc10$gp[13] = rslt.ldm$med.p.global.omni

# no mediator selected
M2 = M + 0.1
M2.p = M2/rowSums(M2)
rslt.hima.psc01 = microHIMA(Trt, Y, M2.p, NULL, 0.1)
rslt.hima.psc01
M2 = M + 0.5
M2.p = M2/rowSums(M2)
rslt.hima.psc05 = microHIMA(Trt, Y, M2.p, NULL, 0.1)
rslt.hima.psc05
M2 = M + 1
M2.p = M2/rowSums(M2)
rslt.hima.psc10 = microHIMA(Trt, Y, M2.p, NULL, 0.1)
rslt.hima.psc10

Trt = combo.filter$meta[["fat"]]; M = combo.filter$mediators; Y = combo.filter$meta[["bmi"]]; covariate = combo.filter$meta[["calor"]]
rslt.ldm = ldm(M | covariate ~ Trt + Y, test.mediation = TRUE)
rslt.combo.filter20.psc01$gp[13] = rslt.ldm$med.p.global.omni
rslt.combo.filter20.psc05$gp[13] = rslt.ldm$med.p.global.omni
rslt.combo.filter20.psc10$gp[13] = rslt.ldm$med.p.global.omni

# M2 = M + 0.1
# M2.p = M2/rowSums(M2)
# rslt.hima.psc01 = microHIMA(Trt, Y, M2.p, covariate, 0.1)
# rslt.hima.psc01
# Running debiased Lasso...     (2022-12-01 15:39:00)
# Error in hommel::hommel(P_adj_DLASSO, simes = FALSE) : 
#   missing values in input p-values
# In addition: Warning messages:
#   1: In sqrt(sigma_e2) : NaNs produced
# 2: In sqrt(sigma_e2) : NaNs produced
# M2 = M + 0.5
# M2.p = M2/rowSums(M2)
# rslt.hima.psc05 = microHIMA(Trt, Y, M2.p, covariate, 0.1)
# rslt.hima.psc05
# M2 = M + 1
# M2.p = M2/rowSums(M2)
# rslt.hima.psc10 = microHIMA(Trt, Y, M2.p, covariate, 0.1)
# rslt.hima.psc10


gp.mat[,1] = c(rslt.cecal.top100.psc01$gp, rslt.cecal.top100.psc01$nullprop.est)
gp.mat[,2] = c(rslt.cecal.top100.psc05$gp, rslt.cecal.top100.psc05$nullprop.est)
gp.mat[,3] = c(rslt.cecal.top100.psc10$gp, rslt.cecal.top100.psc10$nullprop.est)
gp.mat[,4] = c(rslt.combo.filter20.psc01$gp, rslt.combo.filter20.psc01$nullprop.est)
gp.mat[,5] = c(rslt.combo.filter20.psc05$gp, rslt.combo.filter20.psc05$nullprop.est)
gp.mat[,6] = c(rslt.combo.filter20.psc10$gp, rslt.combo.filter20.psc10$nullprop.est)
View(round(gp.mat,3))
# rslt.combo.filter20.Int = runMedTestwInt(combo.filter, exposure = "fat", covariate = "calor", outcome = "bmi")
# gp.mat[c(11:16, 19:24),4] = c(rslt.combo.filter20.Int$gp, rslt.combo.filter20.Int$nullprop.est)

# #### Rarefy
# cecal.rrf.tmp = combo.rrf.tmp = vector(mode = "list", length = 20)
# tmp = matrix(NA, nrow = 24, ncol = 40)
# for (seed in 1:20) {
#   cat(sprintf("======Processing rarefy version #%d======\n", seed))
#   #load(paste0("../Data/Deriveddata/Cecal.RRF", seed, ".filter20top100.rda"))
#   #cecal.rrf.tmp[[seed]] = runSeveralMeds(cecal.rarefy, exposure = "Treatment", outcome = "pFat")
#   tmp[,2*seed-1] = c(cecal.rrf.tmp[[seed]]$gp, cecal.rrf.tmp[[seed]]$nullprop.est)
#   
#   #load(paste0("../Data/Deriveddata/COMBO.RRF", seed, ".filter20.rda"))
#   #combo.rrf.tmp[[seed]] = runSeveralMeds(combo.rarefy, exposure = "fat", covariate = "calor", outcome = "bmi")
#   tmp[,2*seed] = c(combo.rrf.tmp[[seed]]$gp, combo.rrf.tmp[[seed]]$nullprop.est)
# }
# 
# gp.mat[,5] = rowMeans(tmp[,seq(1,40,2)])
# gp.mat[,6] = rowMeans(tmp[,seq(2,40,2)])
# 
# 
# View(t(gp.mat))
# floor((which(tmp[16,] < 0.1) + 1)/2)
# 1 5 6 9 15 16 20
# cecal.rrf.tmp[[16]]$sig.node["perm.jsmix.BH"]
# 15 19
# combo.rrf.tmp[[19]]$sig.node["perm.jsmix.HMP"]
# load(paste0("../Data/Deriveddata/COMBO.RRF", 19, ".filter20.rda"))
# tree = combo.rarefy$tree
# tax = combo.rarefy$taxonomy
# treestructure = .phylostructure(tree)
# tmp = tax[tree$tip.label[(treestructure$descendant[.ntaxa(tree)+336,] == TRUE)],]
# View(tmp)
save(rslt.cecal.top100.psc01, rslt.cecal.top100.psc05, rslt.cecal.top100.psc10,
     rslt.combo.filter20.psc01, rslt.combo.filter20.psc05, rslt.combo.filter20.psc10, 
     gp.mat,
     file = "../Data/Deriveddata/rslt.runModel.rda")
