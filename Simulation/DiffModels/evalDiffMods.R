# compare Type 1 error of log ratio model, dm-mean model and quasi-likelihood model
# N=200 continuous outcome, causalType=1/2
rm(list=ls())

########### input arguments
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 3){
  causal.type = as.numeric(args[1]) # can take values 1,2 corresponding to two approaches to simulate causal species
  A = as.numeric(args[2]) # correspond to A in the document, control the perturbation level for species S_ab in treatment group
  B = as.numeric(args[3]) # correspond to B in the document, control the mediation effect in the outcome model for species S_ab 
} else {
  cat('usage: Rscript evalDiffMods.R <causalType> <A> <B>\n', file = stderr())
  stop()
}

# rm(list = ls())
# setwd("~/Project/QH0020/medtest/PhyloMed/Revised/DiffModels/")
# ######### input arguments
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# A = 1  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 0 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
n.sample = 200 # sample size 
seed = 12 # overall simulation seed for this job 
C = 1 # correspond to C in the document, control direct effect
########### input arguments END

source("sim_utility.R")
library(SKAT)
library(MASS)
library(phyloseq)
library(matrixStats)
library(vegan)
library(nnet)
library(fdrtool)
library(HMP)

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

gp.mat = matrix(NA, nrow = 2000, ncol = 3, dimnames = list(c(), c("Simes.LR", "Simes.DM", "Simes.QS")))
rawp.lr.mat = rawp.dm.mat = rawp.qs.mat = matrix(NA, nrow = nrow(gp.mat), ncol = tree$Nnode)

for (sim in 1:nrow(gp.mat)) {
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
  
  chi.stat.alpha.lr = z.stat.alpha.lr = pval.alpha.lr = 
    chi.stat.alpha.dm = z.stat.alpha.dm = pval.alpha.dm = 
    chi.stat.alpha.qs = z.stat.alpha.qs = pval.alpha.qs = 
    chi.stat.beta = z.stat.beta = pval.beta = numeric(tree$Nnode)
  
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
      chi.stat.alpha.lr[i-K] = NA
      z.stat.alpha.lr[i-K] = NA
      pval.alpha.lr[i-K] = NA
      chi.stat.alpha.dm[i-K] = NA
      z.stat.alpha.dm[i-K] = NA
      pval.alpha.dm[i-K] = NA
      chi.stat.alpha.qs[i-K] = NA
      z.stat.alpha.qs[i-K] = NA
      pval.alpha.qs[i-K] = NA
      chi.stat.beta[i-K] = NA
      z.stat.beta[i-K] = NA
      pval.beta[i-K] = NA
      next
    }
    
    # log-ratio linear regression, score test
    fit = lm(G~0+conf)
    r = fit$residuals
    dispersion = sum(r^2)/fit$df.residual
    fit.full = lm(G~0+conf+Trt)
    est = coef(fit.full)["Trt"]
    x2.1 = qr.resid(fit$qr, Trt)
    stat.lr = (sum(x2.1*r))^2/sum(x2.1*x2.1)/dispersion
    pval.lr = 1 - pchisq(stat.lr, 1)
    chi.stat.alpha.lr[i-K] = stat.lr
    z.stat.alpha.lr[i-K] = sqrt(stat.lr)*sign(est)
    pval.alpha.lr[i-K] = pval.lr
    
    idx = which(rowSums(Mc) != 0)
    Mc.obs = Mc[idx,,drop=FALSE]; conf.obs = conf[idx,,drop=FALSE]
    # # MGLM wald test
    # mod.mglm = tryCatch(MGLMreg(Mc.obs~Trt[idx], dist = "DM"), # only one non-zero in right, output error
    #                     error=function(e) NA)
    # if(is.na(mod.mglm)){
    #   chi.stat.alpha.dm1[i-K] = NA
    #   z.stat.alpha.dm1[i-K] = NA
    #   pval.alpha.dm1[i-K] = NA
    # }else{
    #   stat.dm1 = mod.mglm@wald.value[2] # treatment term
    #   pval.dm1 = mod.mglm@wald.p[2]
    #   chi.stat.alpha.dm1[i-K] = stat.dm1
    #   z.stat.alpha.dm1[i-K] = sqrt(stat.dm1)*sign(est)
    #   pval.alpha.dm1[i-K] = pval.dm1
    # }

    # HMP DM-mean test
    xdc = tryCatch(Xmcupo.sevsample(list(Mc[which(Trt==1),], Mc[which(Trt==0),])), 
                   error=function(e){list(`Xmcupo statistics`=0, `p value`=1)})
    stat.dm = ifelse(xdc$`Xmcupo statistics` < 0, 0, xdc$`Xmcupo statistics`)
    pval.dm = xdc$`p value`
    chi.stat.alpha.dm[i-K] = stat.dm
    z.stat.alpha.dm[i-K] = sqrt(stat.dm)*sign(est)
    pval.alpha.dm[i-K] = pval.dm
    
    # quasi-likelihood test
    n.alpha = ncol(conf.obs) + 1
    # multinom function requires to put the reference taxon as the first one and doesn't need to include a vector of 1 for intercept in X
    input.data = list(Y=cbind(Mc.obs[,2], Mc.obs[,-2]))
    multinom.out = multinom(Y ~ 1, input.data, maxit=10000, abstol = 1.0e-10, reltol = 1.0e-15, trace=FALSE)
    
    est.alpha.obs = rep(0, n.alpha) # last item is for treatment estimate
    est.alpha.obs[-n.alpha] = as.numeric(coef(multinom.out))
    TestAlpha = .test_alpha_qs(Mc[idx,,drop=FALSE], Trt[idx], conf[idx,,drop=FALSE], est.reduce.beta = est.alpha.obs, test.type="mv")
    chi.stat.alpha.qs[i-K] = TestAlpha$stat
    z.stat.alpha.qs[i-K] = sqrt(TestAlpha$stat)*sign(TestAlpha$est)
    pval.alpha.qs[i-K] = TestAlpha$pval
    

    mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
    mod.s2 = mod$sigma^2
    mod.resid = mod$residuals
    TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
    
    var.beta = TestBeta$var.beta
    chi.stat.beta[i-K] = TestBeta$stat
    z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(TestBeta$est)
    pval.beta[i-K] = TestBeta$pval
  }
  
  alpha1.lr = .pi0_JC(na.omit(z.stat.alpha.lr))
  alpha1.dm = .pi0_JC(na.omit(z.stat.alpha.dm))  
  alpha1.qs = .pi0_JC(na.omit(z.stat.alpha.qs))
  alpha2 = .pi0_JC(na.omit(z.stat.beta))

  tmp.lr = .nullEstimation(pval.alpha.lr, pval.beta, alpha1.lr, alpha2)
  tmp.dm = .nullEstimation(pval.alpha.dm, pval.beta, alpha1.dm, alpha2)  
  tmp.qs = .nullEstimation(pval.alpha.qs, pval.beta, alpha1.qs, alpha2)
  rawp.lr = tmp.lr$rawp
  rawp.dm = tmp.dm$rawp
  rawp.qs = tmp.qs$rawp
    
  rawp.lr.rm = na.omit(rawp.lr)
  rawp.dm.rm = na.omit(rawp.dm)
  rawp.qs.rm = na.omit(rawp.qs)
  L.lr = length(rawp.lr.rm)
  L.dm = length(rawp.dm.rm)
  L.qs = length(rawp.qs.rm)
  
  rawp.lr.mat[sim,] = rawp.lr
  rawp.dm.mat[sim,] = rawp.dm
  rawp.qs.mat[sim,] = rawp.qs
  gp.mat[sim,] = c(min(L.lr * rawp.lr.rm/rank(rawp.lr.rm)), 
                   min(L.dm * rawp.dm.rm/rank(rawp.dm.rm)), 
                   min(L.qs * rawp.qs.rm/rank(rawp.qs.rm)))
}

rslt = list(gp.mat = gp.mat,
            rawp.lr.mat = rawp.lr.mat,
            rawp.dm.mat = rawp.dm.mat,
            rawp.qs.mat = rawp.qs.mat)
filename = paste0("Result/allType", causal.type, "Nsample", n.sample, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)


