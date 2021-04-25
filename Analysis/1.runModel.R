rm(list = ls())
setwd("~/Project/QH0020/medtest/PhyloMed/Analysis/")
# load source code for PhyloMed
source("1a.PhyloMed_utility.R")
# load source code for MedTest
source("1b.MeTest_utility.R")
# load source code for MODIMA
source("1c.MODIMA_utility.R")

library(vegan)
library(GUniFrac)
library(energy)
library(fdrtool)
library(SKAT)
library(harmonicmeanp)
library(ccmm)

runSeveralMeds <- function(input_data, fdr.alpha = 0.1, perm.prec = 0.1, seed = 123){
  Trt = input_data$treatment
  M = input_data$mediators
  Y = input_data$outcome
  tree = input_data$tree
  gp = numeric(17)
  names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
                "Simes.PhyloMed.asym", "Fisher.PhyloMed.asym", "HMP.PhyloMed.asym",
                "Simes.PhyloMed.perm", "Fisher.PhyloMed.perm", "HMP.PhyloMed.perm",
                "ccmm.boot.tide", "ccmm.norm.tide")
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
  cat("distance-based MODIMA done\n")
  
  if(grepl("merge", deparse(substitute(input_data)))){
    M2 = M + 0.5 
    M2.P = M2/rowSums(M2)
    ## here perform ccmm
    rslt.ccmm.boot = ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "bootstrap", n.boot = 2000)
    rslt.ccmm.norm = ccmm(y = Y, M = M2.P, tr = Trt, method.est.cov = "normal")
    gp[16] = ifelse(rslt.ccmm.boot$TIDE.CI[1] > 0 | rslt.ccmm.boot$TIDE.CI[2] < 0, 0, 1)
    gp[17] = 2*(1 - pnorm(abs(rslt.ccmm.norm$TIDE)/sqrt(rslt.ccmm.norm$Var.TIDE)))
    cat("ccmm done\n")
  }else{
    gp[16:17] = rep(NA, 2)
  }
  
  set.seed(seed)
  rslt.phylomed = phylomed(Trt, M, Y, tree, fdr.alpha = fdr.alpha, perm.prec = perm.prec, verbose = T)
  
  gp[10:12] = rslt.phylomed$PhyloMed.A$global.pval
  nullprop.est[1:3] = rslt.phylomed$PhyloMed.A$null.prop
  node.pval.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix
  node.pval.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js
  node.pval.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel
  gp[13:15] = rslt.phylomed$PhyloMed.P$global.pval
  nullprop.est[4:6] = rslt.phylomed$PhyloMed.P$null.prop
  node.pval.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix
  node.pval.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js
  
  pval.alphabeta.mat = rbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.P$pval.alpha,
                             rslt.phylomed$PhyloMed.A$pval.beta, rslt.phylomed$PhyloMed.P$pval.beta)
  rownames(pval.alphabeta.mat) = c("asym.alpha", "perm.alpha", "asym.beta", "perm.beta")
  
  node.pval.mat = rbind(node.pval.asym.jsmix, node.pval.asym.js, node.pval.sobel, 
                        node.pval.perm.jsmix, node.pval.perm.js)
  rownames(node.pval.mat) = c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js")
  
  sig.nodeID = apply(node.pval.mat, 1, .getSigNode, alpha = fdr.alpha)
  names(sig.nodeID) = c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js")
  
  rslt = list(gp = gp, nullprop.est = nullprop.est, 
              pval.alphabeta = pval.alphabeta.mat, 
              node.pval = node.pval.mat,
              sig.node = sig.nodeID)
  cat("tree-based PhyloMed done\n")
  
  return(rslt)
}

load("../Data/Deriveddata/cecal.full.rda")
load("../Data/Deriveddata/fecal.full.rda")
load("../Data/Deriveddata/cecal.merge.rda")
load("../Data/Deriveddata/fecal.merge.rda")
gp.mat = matrix(NA, nrow=23, ncol=4, dimnames = list(c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                       "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
                                                       "Simes.PhyloMed.asym", "Fisher.PhyloMed.asym", "HMP.PhyloMed.asym",
                                                       "Simes.PhyloMed.perm", "Fisher.PhyloMed.perm", "HMP.PhyloMed.perm",
                                                       "ccmm.boot.tide", "ccmm.norm.tide",
                                                       "H00.asym","H10.asym","H01.asym","H00.perm","H10.perm","H01.perm"),
                                                     c("cecal.full", "fecal.full", "cecal.merge", "fecal.merge")))

rslt.cecal.full = runSeveralMeds(data.cecal.full, perm.prec = 0.05)
rslt.fecal.full = runSeveralMeds(data.fecal.full, perm.prec = 0.05)
rslt.cecal.merge = runSeveralMeds(data.cecal.merge, perm.prec = 0.05)
rslt.fecal.merge = runSeveralMeds(data.fecal.merge, perm.prec = 0.05)

gp.mat[,1] = c(rslt.cecal.full$gp, rslt.cecal.full$nullprop.est)
gp.mat[,2] = c(rslt.fecal.full$gp, rslt.fecal.full$nullprop.est)
gp.mat[,3] = c(rslt.cecal.merge$gp, rslt.cecal.merge$nullprop.est)
gp.mat[,4] = c(rslt.fecal.merge$gp, rslt.fecal.merge$nullprop.est)

save(rslt.cecal.full, rslt.fecal.full, rslt.cecal.merge, rslt.fecal.merge, gp.mat,
     file = "../Data/Deriveddata/rslt.runModel.rda")
