# compare the global p-value of true proportion, JC.prod proportion and Storey.minus proportion
# N=200 continuous outcome, causalType=1/2, small/large signal
# output alpha00/10/01 estimate and global p-value
# Type1 error and power

rm(list=ls())

########### input arguments
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 3){
  causal.type = as.numeric(args[1]) # can take values 1,2 corresponding to two approaches to simulate causal species
  A = as.numeric(args[2]) # correspond to A in the document, control the perturbation level for species S_ab in treatment group
  B = as.numeric(args[3]) # correspond to B in the document, control the mediation effect in the outcome model for species S_ab 
} else {
  cat('usage: Rscript evalPropNull.R <causalType> <A> <B>\n', file = stderr())
  stop()
}

# rm(list = ls())
# setwd("~/Documents/Project/PhyloMed/Simulation/PropOfNull/eval_asym/")
# ########## input arguments
# causal.type = 2  # can take values 1,2 corresponding to two approaches to simulate causal species
# A = 0.1  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 0.5 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
n.sample = 200 # sample size could be 50 or 200
seed = 12 # overall simulation seed for this job 
C = 1 # correspond to C in the document, control direct effect
########### input arguments END

source("sim_utility.R")
source("AVGcodes.R")
source("SLIMcodes.R")
library(SKAT)
library(MASS)
library(phyloseq)
library(matrixStats)
library(vegan)
library(fdrtool)
library(harmonicmeanp)

load("zeeviD.Rdata")
count = zeeviD$count  # rare species (over 90% obs are zero) have been removed
n.total = nrow(count) # 900 subjects
K = ncol(count) 

P = count/rowSums(count)
P.avg = colMeans(P)

tree = zeeviD$phy.tree

#### sample treatment and confounder (if present)
n1 = n.sample/2
n2 = n1
Trt.ori = c(rep(1,n1),rep(0,n2))  # binary treatment

treestructure = .phylostructure(tree)
# JC estimate pi0. pi.0 and Storey estimate pi00
gp.mat = alpha00.mat = alpha10.mat = alpha01.mat = alpha11.mat = alphastar.mat =
  alpha00.ori.mat = alpha01.ori.mat = alpha10.ori.mat = 
  alpha1.mat = alpha2.mat = matrix(NA, nrow = 2000, ncol = 10, 
                dimnames = list(c(), c("true", "jc.prod", "avg.prod", 
                                       "jc.minus.pmin", "avg.minus.pmin", # pmin to estimate pi00
                                       "avg.minus.revised", # revised storey to estimate pi00
                                       "jc.minus.pmax", "avg.minus.pmax",
                                       "jc.minus.multi", "avg.minus.multi")))
# rawp.jc.prod.mat = rawp.avg.prod.mat = 
#   rawp.jc.minus.pmin.mat = rawp.avg.minus.pmin.mat = 
#   rawp.avg.minus.rev.mat = 
#   rawp.jc.minus.pmax.mat = rawp.avg.minus.pmax.mat = 
#   
causalNode.mat = matrix(NA, nrow = 2000, ncol = tree$Nnode)

for (sim in 1:nrow(gp.mat)) {
  causalNode = numeric(tree$Nnode); causalLeaf = numeric(K)
  if(sim %% 100 == 0){
    print(paste0("Now process #", sim))
  }
  seed.sim = seed + sim
  set.seed(seed.sim)
  #### randomly select n subject
  M = count[sample(1:n.total, n.sample),]
  N = rowSums(M)
  
  #### determine species in sets: S_ab, S_a, S_b, S
  S_ab = NULL
  causalNodeWLeaf = which(rowSums(treestructure$descendant, na.rm = T) == 3)

  if(causal.type==1){
    K_ab = 3
    random.id = sample(causalNodeWLeaf, 1)
    S_ab = which(treestructure$descendant[random.id,] == 1)
  }else if(causal.type==2){
    K_ab = 6
    random.id = sample(causalNodeWLeaf, 2)
    S_ab = which(colSums(treestructure$descendant[random.id,]) == 1)
  }
  
  # S_ab = sample(1:K, 3*causal.type)
  # K_ab = length(S_ab)
  # 1 denotes the causal leaves
  causalLeaf[S_ab] = 1
  # 1 denotes the ancestor of causal species
  causalNode[which(rowSums(treestructure$descendant[(K+1) : (K+tree$Nnode), S_ab]) > 0)] = 1 # the index starts from 1, not 101
  
  #### perturb species in S_ab for subjects in treatment group
  if(any(A != 0, B!= 0)){
    tmp.id = sample(0:1, K_ab, replace = TRUE)
    mul.id = S_ab[tmp.id == 1]
    div.id = S_ab[tmp.id == 0]
    if(length(mul.id) > 0){
      for (j in mul.id) {
        M[1:n1,j] = M[1:n1,j] + rbinom(n1, N[1:n1], A * P.avg[j])
      }
    }
    if(length(div.id) > 0){
      for (j in div.id) {
        M[n1+1:n2,j] = M[n1+1:n2,j] + rbinom(n2, N[n1+1:n2], A * P.avg[j])
      }
    }
  }
  
  #### simulate outcome using species in S_ab and S_b
  M2 = M + 0.5 # outcome model use M2 (no zero in it)
  M2.P = M2/rowSums(M2)
  beta = runif(K_ab, 0, B)
  beta = beta - sum(beta)/length(beta)
  beta.t = runif(1,0,C)
  
  outcome.ori = rnorm(n.sample, beta.t*Trt.ori + log(M2.P[,S_ab]) %*% beta, sd=1)
  conf.ori = matrix(1, nrow = n.sample, ncol = 1)
  
  ## here perform tree-based method
  chi.stat.alpha = chi.stat.beta = z.stat.alpha = z.stat.beta = pval.alpha = pval.beta = numeric(tree$Nnode)

  #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
  for (i in (K + 1):(K + tree$Nnode)) {
    child.left = treestructure$phylochildren[i,1]
    child.right = treestructure$phylochildren[i,2]
    
    Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
    Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
    if(mean(Mc.left) < mean(Mc.right)){
      Mc = cbind(Mc.right, Mc.left)
    }else{
      Mc = cbind(Mc.left, Mc.right)
    }
    
    Mc2 = Mc + 0.5
    idx = which(rowSums(Mc) != 0)
    if(length(idx) != 0){
      Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
      G = log(Mc2[idx,-2]/as.numeric(Mc2[idx,2]))
    }else{
      Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
      G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
    }
    
    # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
    condition1 = any(colSums(Mc) == 0)
    condition2 = (qr(cbind(conf,Trt))$rank < ncol(conf)+1)
    if(any(condition1,condition2)){
      chi.stat.alpha[i-K] = NA
      z.stat.alpha[i-K] = NA
      pval.alpha[i-K] = NA
      chi.stat.beta[i-K] = NA
      z.stat.beta[i-K] = NA
      pval.beta[i-K] = NA
      next
    }
    
    mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
    est = mod.full$coefficients["Trt","Estimate"]
    mod.reduce = summary(glm(G~0+conf,family = quasi()))
    mod.reduce.disp = mod.reduce$dispersion
    mod.reduce.resid = mod.reduce$deviance.resid
    TestAlpha = .test_alpha(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
    stat = TestAlpha$stat
    pval = TestAlpha$pval
    
    chi.stat.alpha[i-K] = stat
    z.stat.alpha[i-K] = sqrt(stat)*sign(est)
    pval.alpha[i-K] = pval
    var.alpha = mod.full$cov.scaled["Trt", "Trt"]
    # fit = lm(G~0+conf)
    # r = fit$residuals
    # dispersion = sum(r^2)/fit$df.residual
    # fit.full = lm(G~0+conf+Trt)
    # est = coef(fit.full)["Trt"]
    # var.alpha = vcov(fit.full)["Trt","Trt"]
    # x2.1 = qr.resid(fit$qr, Trt)
    # stat = (sum(x2.1*r))^2/sum(x2.1*x2.1)/dispersion
    # pval = 1 - pchisq(stat, 1)
    
    # chi.stat.alpha[i-K] = stat
    # z.stat.alpha[i-K] = sqrt(stat)*sign(est)
    # pval.alpha[i-K] = pval
    
    if(length(unique(outcome)) > 2) {
      # continuous traits
      mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
      mod.s2 = mod$sigma^2
      mod.resid = mod$residuals
      TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
      
    } else {
      # binary traits
      mod = glm(outcome~cbind(conf[,-1], Trt), family = "binomial")
      mod.est = mod$coefficients
      TestBeta = .test_beta(outcome, G, Trt, conf, est.obs = mod.est, test.type="mv")
    }
    var.beta = TestBeta$var.beta
    chi.stat.beta[i-K] = TestBeta$stat
    z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(TestBeta$est)
    pval.beta[i-K] = TestBeta$pval
  }
  
  alpha1.true = ifelse(A==0, 1, mean(causalNode == 0))
  alpha2.true = ifelse(B==0, 1, mean(causalNode == 0))
  if(A==0){
    if(B==0){
      alphastar.true = 1
      alpha00.true = alpha00.ori = 1; alpha01.true = alpha10.true = alpha01.ori = alpha10.ori = 0; alpha11.true = 0
    }else{
      alphastar.true = 1
      alpha10.true = alpha10.ori = 0; alpha00.true = alpha00.ori = alpha2.true; alpha01.true = alpha01.ori = 1 - alpha00.true; alpha11.true = 0
    }
  }else{
    if(B==0){
      alphastar.true = 1
      alpha01.true = alpha01.ori = 0; alpha00.true = alpha00.ori = alpha1.true; alpha10.true = alpha10.ori = 1 - alpha00.true; alpha11.true = 0
    }else{
      alphastar.true = alpha1.true
      alpha00.ori = alpha1.true; alpha01.true = alpha10.true =alpha01.ori = alpha10.ori = 0; alpha11.true = 1-alpha1.true
      alpha00.true = 1
    }
  }
  tmp.true = .nullEstimation_true(pval.alpha, pval.beta, alpha1.true, alpha2.true, alpha00.true, alpha01.true, alpha10.true, alpha11.true)
  rawp.true = tmp.true$rawp
  null.prop.true = c(tmp.true$alpha00, tmp.true$alpha10, tmp.true$alpha01, tmp.true$alpha11)
  rawp.true.rm = na.omit(rawp.true)
  nullprop.ori.true = c(alpha00.ori, alpha10.true, alpha01.true)
  
  # JC estimate pi0. pi.0
  alpha1.est.jc = .pi0_JC(na.omit(z.stat.alpha))
  alpha2.est.jc = .pi0_JC(na.omit(z.stat.beta))
  tmp.jc.prod = .nullEstimation_prod(pval.alpha, pval.beta, alpha1.est.jc, alpha2.est.jc, z.stat.alpha, z.stat.beta, est.method = "prod")
  rawp.jc.prod = tmp.jc.prod$rawp
  null.prop.jc.prod = c(tmp.jc.prod$alpha00, tmp.jc.prod$alpha10, tmp.jc.prod$alpha01, tmp.jc.prod$alpha11)
  alphastar.jc.prod = tmp.jc.prod$alphastar
  rawp.jc.prod.rm = na.omit(rawp.jc.prod)
  nullprop.ori.jc.prod = c(tmp.jc.prod$alpha00, tmp.jc.prod$alpha10, tmp.jc.prod$alpha01)
  
  tmp.jc.minus.pmin = .nullEstimation_prod(pval.alpha, pval.beta, alpha1.est.jc, alpha2.est.jc, z.stat.alpha, z.stat.beta, est.method = "pmin")
  rawp.jc.minus.pmin = tmp.jc.minus.pmin$rawp
  null.prop.jc.minus.pmin = c(tmp.jc.minus.pmin$alpha00, tmp.jc.minus.pmin$alpha10, tmp.jc.minus.pmin$alpha01, tmp.jc.minus.pmin$alpha11)
  alphastar.jc.minus.pmin = tmp.jc.minus.pmin$alphastar
  rawp.jc.minus.pmin.rm = na.omit(rawp.jc.minus.pmin)
  nullprop.ori.jc.minus.pmin = c(tmp.jc.minus.pmin$alpha00.r, tmp.jc.minus.pmin$alpha10.r, tmp.jc.minus.pmin$alpha01.r)
  
  tmp.jc.minus.pmax = .nullEstimation_prod(pval.alpha, pval.beta, alpha1.est.jc, alpha2.est.jc, z.stat.alpha, z.stat.beta, est.method = "pmax")
  rawp.jc.minus.pmax = tmp.jc.minus.pmax$rawp
  null.prop.jc.minus.pmax = c(tmp.jc.minus.pmax$alpha00, tmp.jc.minus.pmax$alpha10, tmp.jc.minus.pmax$alpha01, tmp.jc.minus.pmax$alpha11)
  alphastar.jc.minus.pmax = tmp.jc.minus.pmax$alphastar
  rawp.jc.minus.pmax.rm = na.omit(rawp.jc.minus.pmax)
  nullprop.ori.jc.minus.pmax = c(tmp.jc.minus.pmax$alpha00.r, tmp.jc.minus.pmax$alpha10.r, tmp.jc.minus.pmax$alpha01.r)
  
  tmp.jc.minus.multi = .nullEstimation_pmax_multistage(pval.alpha, pval.beta, alpha1.est.jc, alpha2.est.jc, z.stat.alpha, z.stat.beta)
  rawp.jc.minus.multi = tmp.jc.minus.multi$rawp
  null.prop.jc.minus.multi = c(tmp.jc.minus.multi$alpha00, tmp.jc.minus.multi$alpha10, tmp.jc.minus.multi$alpha01, tmp.jc.minus.multi$alpha11)
  alphastar.jc.minus.multi = tmp.jc.minus.multi$alphastar
  rawp.jc.minus.multi.rm = na.omit(rawp.jc.minus.multi)
  nullprop.ori.jc.minus.multi = c(tmp.jc.minus.multi$alpha00.r, tmp.jc.minus.multi$alpha10.r, tmp.jc.minus.multi$alpha01.r)
  
  # rawp.jc.prod.mat[sim,] = rawp.jc.prod
  # rawp.jc.minus.pmin.mat[sim,] = rawp.jc.minus.pmin
  # rawp.jc.minus.pmax.mat[sim,] = rawp.jc.minus.pmax
  
  # AVG estimate pi0. and pi.0
  alpha1.est.avg = AverageEstimate(na.omit(pval.alpha))
  alpha2.est.avg = AverageEstimate(na.omit(pval.beta))
  tmp.avg.prod = .nullEstimation_prod(pval.alpha, pval.beta, alpha1.est.avg, alpha2.est.avg, z.stat.alpha, z.stat.beta, est.method = "prod")
  rawp.avg.prod = tmp.avg.prod$rawp
  null.prop.avg.prod = c(tmp.avg.prod$alpha00, tmp.avg.prod$alpha10, tmp.avg.prod$alpha01, tmp.avg.prod$alpha11)
  alphastar.avg.prod = tmp.avg.prod$alphastar
  rawp.avg.prod.rm = na.omit(rawp.avg.prod)
  nullprop.ori.avg.prod = c(tmp.avg.prod$alpha00.r, tmp.avg.prod$alpha10.r, tmp.avg.prod$alpha01.r)
  
  tmp.avg.minus.pmin = .nullEstimation_minus(pval.alpha, pval.beta, alpha1.est.avg, alpha2.est.avg, est.method = "pmin")
  rawp.avg.minus.pmin = tmp.avg.minus.pmin$rawp
  null.prop.avg.minus.pmin = c(tmp.avg.minus.pmin$alpha00, tmp.avg.minus.pmin$alpha10, tmp.avg.minus.pmin$alpha01, tmp.avg.minus.pmin$alpha11)
  alphastar.avg.minus.pmin = tmp.avg.minus.pmin$alphastar
  rawp.avg.minus.pmin.rm = na.omit(rawp.avg.minus.pmin)
  nullprop.ori.avg.minus.pmin = c(tmp.avg.minus.pmin$alpha00.r, tmp.avg.minus.pmin$alpha10.r, tmp.avg.minus.pmin$alpha01.r)
  
  tmp.avg.minus.pmax = .nullEstimation_minus(pval.alpha, pval.beta, alpha1.est.avg, alpha2.est.avg, est.method = "pmax")
  rawp.avg.minus.pmax = tmp.avg.minus.pmax$rawp
  null.prop.avg.minus.pmax = c(tmp.avg.minus.pmax$alpha00, tmp.avg.minus.pmax$alpha10, tmp.avg.minus.pmax$alpha01, tmp.avg.minus.pmax$alpha11)
  alphastar.avg.minus.pmax = tmp.avg.minus.pmax$alphastar
  rawp.avg.minus.pmax.rm = na.omit(rawp.avg.minus.pmax)
  nullprop.ori.avg.minus.pmax = c(tmp.avg.minus.pmax$alpha00.r, tmp.avg.minus.pmax$alpha10.r, tmp.avg.minus.pmax$alpha01.r)
  
  tmp.avg.minus.rev = .nullEstimation_minus(pval.alpha, pval.beta, alpha1.est.avg, alpha2.est.avg, est.method = "revised")
  rawp.avg.minus.rev = tmp.avg.minus.rev$rawp
  null.prop.avg.minus.rev = c(tmp.avg.minus.rev$alpha00, tmp.avg.minus.rev$alpha10, tmp.avg.minus.rev$alpha01, tmp.avg.minus.rev$alpha11)
  alphastar.avg.minus.rev = tmp.avg.minus.rev$alphastar
  rawp.avg.minus.rev.rm = na.omit(rawp.avg.minus.rev)
  nullprop.ori.avg.minus.rev = c(tmp.avg.minus.rev$alpha00.r, tmp.avg.minus.rev$alpha10.r, tmp.avg.minus.rev$alpha01.r)
  
  tmp.avg.minus.multi = .nullEstimation_pmax_multistage(pval.alpha, pval.beta, alpha1.est.avg, alpha2.est.avg, z.stat.alpha, z.stat.beta)
  rawp.avg.minus.multi = tmp.avg.minus.multi$rawp
  null.prop.avg.minus.multi = c(tmp.avg.minus.multi$alpha00, tmp.avg.minus.multi$alpha10, tmp.avg.minus.multi$alpha01, tmp.avg.minus.multi$alpha11)
  alphastar.avg.minus.multi = tmp.avg.minus.multi$alphastar
  rawp.avg.minus.multi.rm = na.omit(rawp.avg.minus.multi)
  nullprop.ori.avg.minus.multi = c(tmp.avg.minus.multi$alpha00.r, tmp.avg.minus.multi$alpha10.r, tmp.avg.minus.multi$alpha01.r)
  # rawp.avg.prod.mat[sim,] = rawp.avg.prod
  # rawp.avg.minus.pmin.mat[sim,] = rawp.avg.minus.pmin
  # rawp.avg.minus.pmax.mat[sim,] = rawp.avg.minus.pmax
  # rawp.avg.minus.rev.mat[sim,] = rawp.avg.minus.rev

  causalNode.mat[sim,] = causalNode
  
  for (i in 1:4) {
    tmp = c(null.prop.true[i], null.prop.jc.prod[i], null.prop.avg.prod[i],
            null.prop.jc.minus.pmin[i], null.prop.avg.minus.pmin[i], null.prop.avg.minus.rev[i],
            null.prop.jc.minus.pmax[i], null.prop.avg.minus.pmax[i],
            null.prop.jc.minus.multi[i], null.prop.avg.minus.multi[i])
    tmp.ori = c(nullprop.ori.true[i], nullprop.ori.jc.prod[i], nullprop.ori.avg.prod[i], 
                nullprop.ori.jc.minus.pmin[i], nullprop.ori.avg.minus.pmin[i], nullprop.ori.avg.minus.rev[i],
                nullprop.ori.jc.minus.pmax[i], nullprop.ori.avg.minus.pmax[i],
                nullprop.ori.jc.minus.multi[i], nullprop.ori.avg.minus.multi[i])
    if(i == 1){
      alpha00.mat[sim,] = tmp
      alpha00.ori.mat[sim,] = tmp.ori
    }else if(i == 2){
      alpha10.mat[sim,] = tmp
      alpha10.ori.mat[sim,] = tmp.ori
    }else if(i == 3){
      alpha01.mat[sim,] = tmp
      alpha01.ori.mat[sim,] = tmp.ori
    }else if(i == 4){
      alpha11.mat[sim,] = tmp
    }
  }
  alpha1.mat[sim,] = c(alpha1.true, alpha1.est.jc, alpha1.est.avg, alpha1.est.jc, alpha1.est.avg, alpha1.est.avg, alpha1.est.jc, alpha1.est.avg, alpha1.est.jc, alpha1.est.avg)
  alpha2.mat[sim,] = c(alpha2.true, alpha2.est.jc, alpha2.est.avg, alpha2.est.jc, alpha2.est.avg, alpha2.est.avg, alpha2.est.jc, alpha2.est.avg, alpha2.est.jc, alpha2.est.avg)
  alphastar.mat[sim,] = c(alphastar.true, 
                          alphastar.jc.prod, alphastar.avg.prod, alphastar.jc.minus.pmin, alphastar.avg.minus.pmin, 
                          alphastar.avg.minus.rev, alphastar.jc.minus.pmax, alphastar.avg.minus.pmax,
                          alphastar.jc.minus.multi, alphastar.avg.minus.multi)
  
  if(length(which(rawp.true.rm==0))>0)
    rawp.true.rm[rawp.true.rm == 0] = runif(length(which(rawp.true.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.jc.prod.rm==0))>0)  
    rawp.jc.prod.rm[rawp.jc.prod.rm == 0] = runif(length(which(rawp.jc.prod.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.jc.minus.pmax.rm==0))>0)  
    rawp.jc.minus.pmax.rm[rawp.jc.minus.pmax.rm == 0] = runif(length(which(rawp.jc.minus.pmax.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.jc.minus.multi.rm==0))>0)  
    rawp.jc.minus.multi.rm[rawp.jc.minus.multi.rm == 0] = runif(length(which(rawp.jc.minus.multi.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.jc.minus.pmin.rm==0))>0)  
    rawp.jc.minus.pmin.rm[rawp.jc.minus.pmin.rm == 0] = runif(length(which(rawp.jc.minus.pmin.rm==0)), min = 0, max = 1e-7)

  if(length(which(rawp.avg.prod.rm==0))>0)  
    rawp.avg.prod.rm[rawp.avg.prod.rm == 0] = runif(length(which(rawp.avg.prod.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.avg.minus.pmax.rm==0))>0)  
    rawp.avg.minus.pmax.rm[rawp.avg.minus.pmax.rm == 0] = runif(length(which(rawp.avg.minus.pmax.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.avg.minus.multi.rm==0))>0)  
    rawp.avg.minus.multi.rm[rawp.avg.minus.multi.rm == 0] = runif(length(which(rawp.avg.minus.multi.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.avg.minus.pmin.rm==0))>0)  
    rawp.avg.minus.pmin.rm[rawp.avg.minus.pmin.rm == 0] = runif(length(which(rawp.avg.minus.pmin.rm==0)), min = 0, max = 1e-7)
  if(length(which(rawp.avg.minus.rev.rm==0))>0)  
    rawp.avg.minus.rev.rm[rawp.avg.minus.rev.rm == 0] = runif(length(which(rawp.avg.minus.rev.rm==0)), min = 0, max = 1e-7)
  
  L = length(rawp.true.rm)
  gp.mat[sim,] = c(p.hmp(rawp.true.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.jc.prod.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.avg.prod.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.jc.minus.pmin.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.avg.minus.pmin.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.avg.minus.rev.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.jc.minus.pmax.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.avg.minus.pmax.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.jc.minus.multi.rm, w = rep(1/L, L), L = L),
                   p.hmp(rawp.avg.minus.multi.rm, w = rep(1/L, L), L = L))
}

rslt = list(gp.mat = gp.mat, 
            alpha1.mat = alpha1.mat,
            alpha2.mat = alpha2.mat,
            alpha00.mat = alpha00.mat, alpha10.mat = alpha10.mat, 
            alpha01.mat = alpha01.mat, alpha11.mat = alpha11.mat,
            alpha00.ori.mat = alpha00.ori.mat,
            alpha10.ori.mat = alpha10.ori.mat,
            alpha01.ori.mat = alpha01.ori.mat,
            alphastar.mat = alphastar.mat,
            # rawp.jc.prod.mat = rawp.jc.prod.mat, 
            # rawp.jc.minus.pmax.mat = rawp.jc.minus.pmax.mat, 
            # rawp.jc.minus.pmin.mat = rawp.jc.minus.pmin.mat,
            # rawp.avg.prod.mat = rawp.avg.prod.mat, 
            # rawp.avg.minus.pmax.mat = rawp.avg.minus.pmax.mat, 
            # rawp.avg.minus.pmin.mat = rawp.avg.minus.pmin.mat, 
            # rawp.avg.minus.rev.mat = rawp.avg.minus.rev.mat,
            causalNode.mat = causalNode.mat)
filename = paste0("Result/finalType", causal.type, "Nsample", n.sample, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)
