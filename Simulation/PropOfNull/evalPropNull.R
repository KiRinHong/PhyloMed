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
# setwd("~/Project/QH0020/medtest/PhyloMed/Simulation/PropOfNull/")
# ########## input arguments
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# A = 0  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 0 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
n.sample = 200 # sample size could be 50 or 200
seed = 12 # overall simulation seed for this job 
C = 1 # correspond to C in the document, control direct effect
########### input arguments END

source("sim_utility.R")
library(SKAT)
library(MASS)
library(phyloseq)
library(matrixStats)
library(vegan)
library(fdrtool)

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
Trt = c(rep(1,n1),rep(0,n2))  # binary treatment

treestructure = .phylostructure(tree)

gp.mat = matrix(NA, nrow = 2000, ncol = 3, 
                dimnames = list(c(), c("Simes.true", "Simes.prod", "Simes.minus")))
alpha00.mat = alpha10.mat = alpha01.mat = matrix(NA, nrow = nrow(gp.mat), ncol = 3, 
                     dimnames = list(c(), c("true","prod","minus")))
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
  
  # 1 denotes the causal leaves
  causalLeaf[S_ab] = 1
  # 1 denotes the ancestor of causal species
  causalNode[which(rowSums(treestructure$descendant[(K+1) : (K+tree$Nnode), S_ab]) > 0)] = 1 # the index starts from 1, not 101
  
  #### perturb species in S_ab for subjects in treatment group
  for(j in S_ab ){
    M[1:n1,j] = M[1:n1,j] + rbinom(n1, N[1:n1], A * P.avg[j])
  }
  
  #### simulate outcome using species in S_ab and S_b
  M2 = M + 0.5 # outcome model use M2 (no zero in it)
  M2.P = M2/rowSums(M2)
  beta = runif(K_ab, 0, B)
  beta = beta - sum(beta)/length(beta)
  beta.t = runif(1,0,C)
  
  outcome = rnorm(n.sample, beta.t*Trt + log(M2.P[,S_ab]) %*% beta, sd=1)
  conf = matrix(1, nrow = n.sample, ncol = 1)
  
  ## here perform tree-based method
  chi.stat.alpha = chi.stat.beta = z.stat.alpha = z.stat.beta = pval.alpha = pval.beta = numeric(tree$Nnode)

  #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
  for (i in (K + 1):(K + tree$Nnode)) {
    child.left = treestructure$phylochildren[i,1]
    child.right = treestructure$phylochildren[i,2]
    
    Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
    Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
    Mc2.left = rowSums(M2[,treestructure$descendant[child.left,], drop = FALSE])
    Mc2.right = rowSums(M2[,treestructure$descendant[child.right,], drop = FALSE])
    if(mean(Mc.left) < mean(Mc.right)){
      Mc = cbind(Mc.right, Mc.left); Mc2 = cbind(Mc2.right, Mc2.left);
    }else{
      Mc = cbind(Mc.left, Mc.right); Mc2 = cbind(Mc2.left, Mc2.right);
    }
    
    G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
    
    # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
    condition1 = any(colSums(Mc) == 0)
    
    if(condition1){
      chi.stat.alpha[i-K] = NA
      z.stat.alpha[i-K] = NA
      pval.alpha[i-K] = NA
      chi.stat.beta[i-K] = NA
      z.stat.beta[i-K] = NA
      pval.beta[i-K] = NA
      next
    }
    
    fit = lm(G~0+conf)
    r = fit$residuals
    dispersion = sum(r^2)/fit$df.residual
    fit.full = lm(G~0+conf+Trt)
    est = coef(fit.full)["Trt"]
    var.alpha = vcov(fit.full)["Trt","Trt"]
    x2.1 = qr.resid(fit$qr, Trt)
    stat = (sum(x2.1*r))^2/sum(x2.1*x2.1)/dispersion
    pval = 1 - pchisq(stat, 1)
    
    chi.stat.alpha[i-K] = stat
    z.stat.alpha[i-K] = sqrt(stat)*sign(est)
    pval.alpha[i-K] = pval
    
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
      alpha00.true = 1; alpha01.true = alpha10.true = 0
    }else{
      alpha10.true = 0; alpha00.true = alpha2.true; alpha01.true = 1 - alpha00.true
    }
  }else{
    if(B==0){
      alpha01.true = 0; alpha00.true = alpha1.true; alpha10.true = 1 - alpha00.true
    }else{
      alpha00.true = 1; alpha01.true = alpha10.true = 0
    }
  }
  tmp.true = .nullEstimation_true(pval.alpha, pval.beta, alpha1.true, alpha2.true,
                                  alpha00.true, alpha01.true, alpha10.true)
  rawp.true = tmp.true$rawp
  null.prop.true = c(tmp.true$alpha00, tmp.true$alpha10, tmp.true$alpha01)
  rawp.true.rm = na.omit(rawp.true)
  L.true = length(rawp.true.rm)
  
  alpha1.est = .pi0_JC(na.omit(z.stat.alpha))
  alpha2.est = .pi0_JC(na.omit(z.stat.beta))
  tmp.prod = .nullEstimation_prod(pval.alpha, pval.beta, alpha1.est, alpha2.est)
  rawp.prod = tmp.prod$rawp
  null.prop.prod = c(tmp.prod$alpha00, tmp.prod$alpha10, tmp.prod$alpha01)
  rawp.prod.rm = na.omit(rawp.prod)
  L.prod = length(rawp.prod.rm)
  
  # alpha00.minus = .pi0_JC(na.omit(.p_which_max(list(z.stat.alpha, z.stat.beta),
  #                                              list(abs(z.stat.alpha), abs(z.stat.beta)))))
  tmp.minus = .nullEstimation_minus(pval.alpha, pval.beta)
  rawp.minus = tmp.minus$rawp
  null.prop.minus = c(tmp.minus$alpha00, tmp.minus$alpha10, tmp.minus$alpha01)
  rawp.minus.rm = na.omit(rawp.minus)
  L.minus = length(rawp.minus.rm)
  
  alpha00.mat[sim,] = c(null.prop.true[1], null.prop.prod[1], null.prop.minus[1])
  alpha10.mat[sim,] = c(null.prop.true[2], null.prop.prod[2], null.prop.minus[2])
  alpha01.mat[sim,] = c(null.prop.true[3], null.prop.prod[3], null.prop.minus[3])
  
  gp.mat[sim,] = c(min(L.true * rawp.true.rm/rank(rawp.true.rm)),
                   min(L.prod * rawp.prod.rm/rank(rawp.prod.rm)),
                   min(L.minus * rawp.minus.rm/rank(rawp.minus.rm)))
}

rslt = list(gp.mat = gp.mat, alpha00.mat = alpha00.mat, alpha10.mat = alpha10.mat, alpha01.mat = alpha01.mat)
filename = paste0("Result/allType", causal.type, "Nsample", n.sample, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)

