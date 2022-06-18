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
runSeveralMeds <- function(input_data, exposure, covariate=NULL, outcome, 
                           fdr.alpha = 0.1, perm.prec = 0.05, seed = 84, runCMM = FALSE){
  # input_data = cecal.top; exposure = "Treatment"; covariate = NULL; outcome = "pFat"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 123; runCMM = FALSE
  # input_data = combo.filter; exposure = "fat"; covariate = "calor"; outcome = "bmi"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 123; runCMM = TRUE
  Trt = input_data$meta[[exposure]]
  M = input_data$mediators
  Y = input_data$meta[[outcome]]
  
  tree = input_data$tree
  gp = numeric(18)
  names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.bonf.MODIMA",
                "Simes.PhyloMed.asym", "Fisher.PhyloMed.asym", "HMP.PhyloMed.asym",
                "Simes.PhyloMed.perm", "Fisher.PhyloMed.perm", "HMP.PhyloMed.perm",
                "CMM.norm", "CMM.boot")
  nullprop.est = numeric(6)
  names(nullprop.est) = c("H00.asym","H10.asym","H01.asym","H00.perm","H10.perm","H01.perm")
  
  # perform distance-based MedTest
  unifracs = GUniFrac(M, tree)$unifracs 
  m.list = list(BC = as.matrix(vegdist(M, method = "bray")),
                JAC = as.matrix(vegdist(M, 'jaccard', binary = TRUE)),
                UniFrac=unifracs[, , c('d_UW')], # unweighted unifrac distance doesn't work because M2 doesn't contain zero
                WUniFrac=unifracs[, , c('d_1')])
  rslt.dist = MedOmniTest(x = Trt, y = Y, m.list = m.list, z = NULL, nperm = 2000)
  gp[1] = rslt.dist$margPs[1]
  gp[2] = rslt.dist$margPs[2]
  gp[3] = rslt.dist$margPs[3]
  gp[4] = rslt.dist$margPs[4]
  gp[5] = rslt.dist$permP
  cat("distance-based MedTest done\n")
  
  dist_Trt = dist(Trt)
  dist_outcome = dist(Y)
  gp[6] = modima(dist_Trt, m.list[[1]],  dist_outcome, nrep=2000)$p.value
  gp[7] = modima(dist_Trt, m.list[[2]],  dist_outcome, nrep=2000)$p.value
  gp[8] = modima(dist_Trt, m.list[[3]],  dist_outcome, nrep=2000)$p.value
  gp[9] = modima(dist_Trt, m.list[[4]],  dist_outcome, nrep=2000)$p.value
  gp[10] = min(c(gp[6:9]*4, 1))
  cat("distance-based MODIMA done\n")
  
  set.seed(seed)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      covariate = input_data$meta[[covariate]]
    }else{
      covariate = as.matrix(sample_data(input_data$meta[,covariate]))
    }
  }
  rslt.phylomed = phylomed(Trt, M, Y, tree, covariates = covariate, fdr.alpha = fdr.alpha, perm.prec = perm.prec, verbose = T)
  
  gp[11:13] = rslt.phylomed$PhyloMed.A$global.pval
  nullprop.est[1:3] = rslt.phylomed$PhyloMed.A$null.prop
  node.pval.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix
  node.pval.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js
  node.pval.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel
  gp[14:16] = rslt.phylomed$PhyloMed.P$global.pval
  nullprop.est[4:6] = rslt.phylomed$PhyloMed.P$null.prop
  node.pval.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix
  node.pval.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js
  
  pval.alphabeta.mat = rbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.P$pval.alpha,
                             rslt.phylomed$PhyloMed.A$pval.beta, rslt.phylomed$PhyloMed.P$pval.beta)
  rownames(pval.alphabeta.mat) = c("asym.alpha", "perm.alpha", "asym.beta", "perm.beta")
  
  node.pval.mat = rbind(node.pval.asym.jsmix, node.pval.asym.js, node.pval.sobel, 
                        node.pval.perm.jsmix, node.pval.perm.js)
  rownames(node.pval.mat) = c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js")
  
  sig.nodeID.BH = apply(node.pval.mat, 1, .getSigNode.BH, alpha = fdr.alpha)
  sig.nodeID.HMP = apply(node.pval.mat, 1, .getSigNode.HMP, alpha = fdr.alpha, tree = tree)
  sig.nodeID = c(sig.nodeID.BH, sig.nodeID.HMP)
  names(sig.nodeID) = apply(expand.grid(c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js"), c("BH","HMP")), 1, paste, collapse=".")
  cat("tree-based PhyloMed done\n")
  
  if(runCMM){
    M2 = M + 0.5
    M2.P = M2/rowSums(M2)
    ## here perform ccmm
    rslt.ccmm.boot <- ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "bootstrap", n.boot = 2000, sig.level = fdr.alpha)
    rslt.ccmm.norm <- ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "normal")
    gp[17] = 2*(1 - pnorm(abs(rslt.ccmm.norm$TIDE)/sqrt(rslt.ccmm.norm$Var.TIDE)))
    gp[18] = ifelse(rslt.ccmm.boot$TIDE.CI[1] > 0 | rslt.ccmm.boot$TIDE.CI[2] < 0, 0, 1)
    cat("ccmm done\n")
    rslt = list(gp = gp, nullprop.est = nullprop.est, 
                pval.alphabeta = pval.alphabeta.mat, 
                node.pval = node.pval.mat,
                sig.node = sig.nodeID,
                rslt.ccmm.norm = rslt.ccmm.norm,
                rlst.ccmm.boot = rslt.ccmm.boot)
    return(rslt)
  }else{
    gp[17] = NA
    gp[18] = NA
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

gp.mat = matrix(NA, nrow=24, ncol=2, dimnames = list(c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                       "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.bonf.MODIMA",
                                                       "Simes.PhyloMed.asym", "Fisher.PhyloMed.asym", "HMP.PhyloMed.asym",
                                                       "Simes.PhyloMed.perm", "Fisher.PhyloMed.perm", "HMP.PhyloMed.perm",
                                                       "CMM.norm", "CMM.boot",
                                                       "H00.asym","H10.asym","H01.asym","H00.perm","H10.perm","H01.perm"),
                                                     c("Cecal.top100", "COMBO.filter20")))

load("../Data/Deriveddata/Cecal.filter20top100.rda")
rslt.cecal.top100 = runSeveralMeds(cecal.top, exposure = "Treatment", outcome = "pFat", runCMM = FALSE)
gp.mat[,1] = c(rslt.cecal.top100$gp, rslt.cecal.top100$nullprop.est)
# rslt.cecal.top100.Int = runMedTestwInt(cecal.top, exposure = "Treatment", outcome = "pFat", interaction = TRUE)
# gp.mat[c(11:16, 19:24),2] = c(rslt.cecal.top100.Int$gp, rslt.cecal.top100.Int$nullprop.est)

load("../Data/Deriveddata/COMBO.filter20.rda")
rslt.combo.filter20 = runSeveralMeds(combo.filter, exposure = "fat", covariate = "calor", outcome = "bmi", runCMM = FALSE)
gp.mat[,2] = c(rslt.combo.filter20$gp, rslt.combo.filter20$nullprop.est)
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
save(rslt.cecal.top100, rslt.combo.filter20, 
     gp.mat,
     file = "../Data/Deriveddata/rslt.runModel.rda")
