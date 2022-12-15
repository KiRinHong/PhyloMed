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

######### input arguments
# setwd("~/Documents/Project/PhyloMed/Simulation/DiffModels")
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# A = 1  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 0 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
n.sample = 200 # sample size 
seed = 12 # overall simulation seed for this job 
C = 1 # correspond to C in the document, control direct effect
# ########### input arguments END

source("sim_utility.R")
library(SKAT)
library(MASS)
library(phyloseq)
library(matrixStats)
library(vegan)
library(nnet)
library(fdrtool)
library(HMP)
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

gp.mat = matrix(NA, nrow = 2000, ncol = 4, 
                dimnames = list(c(), c("HMP.LR", "HMP.DM", "HMP.QSC", "HMP.QSG")))
rawp.lr.mat = rawp.dm.mat = rawp.qsc.mat = rawp.qsg.mat = matrix(NA, nrow = nrow(gp.mat), ncol = tree$Nnode)

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
  
  chi.stat.alpha.lr = z.stat.alpha.lr = pval.alpha.lr = 
    chi.stat.alpha.dm = z.stat.alpha.dm = pval.alpha.dm = 
    chi.stat.alpha.qsc = z.stat.alpha.qsc = pval.alpha.qsc = 
    chi.stat.alpha.qsg = z.stat.alpha.qsg = pval.alpha.qsg =
    chi.stat.beta = z.stat.beta = pval.beta = numeric(tree$Nnode)
  
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
    if(any(condition1, condition2)){
      chi.stat.alpha.lr[i-K] = NA
      z.stat.alpha.lr[i-K] = NA
      pval.alpha.lr[i-K] = NA
      chi.stat.alpha.dm[i-K] = NA
      z.stat.alpha.dm[i-K] = NA
      pval.alpha.dm[i-K] = NA
      chi.stat.alpha.qsc[i-K] = NA
      z.stat.alpha.qsc[i-K] = NA
      pval.alpha.qsc[i-K] = NA
      chi.stat.alpha.qsg[i-K] = NA
      z.stat.alpha.qsg[i-K] = NA
      pval.alpha.qsg[i-K] = NA
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
    
    
    # HMP DM-mean test
    xdc = tryCatch(Xmcupo.sevsample(list(Mc[which(Trt==1),], Mc[which(Trt==0),])), 
                   error=function(e){list(`Xmcupo statistics`=0, `p value`=1)})
    stat.dm = ifelse(xdc$`Xmcupo statistics` < 0, 0, xdc$`Xmcupo statistics`)
    pval.dm = ifelse(xdc$`p value` == 0, 1e-7, xdc$`p value`)
    chi.stat.alpha.dm[i-K] = stat.dm
    z.stat.alpha.dm[i-K] = sqrt(stat.dm)*sign(est)
    pval.alpha.dm[i-K] = pval.dm
    
    # quasi-likelihood test, model count
    n.alpha = ncol(conf) + 1
    # multinom function requires to put the reference taxon as the first one and doesn't need to include a vector of 1 for intercept in X
    input.data = list(Y=cbind(Mc[idx,2,drop=FALSE], Mc[idx,-2,drop=FALSE]))
    multinom.out = multinom(Y ~ 1, input.data, maxit=10000, abstol = 1.0e-10, reltol = 1.0e-15, trace=FALSE)
    
    est.alpha.obs = rep(0, n.alpha) # last item is for treatment estimate
    est.alpha.obs[-n.alpha] = as.numeric(coef(multinom.out))
    TestAlpha = .test_alpha_qsc(Mc[idx,,drop=FALSE], Trt, conf, est.reduce.beta = est.alpha.obs, test.type="mv")
    chi.stat.alpha.qsc[i-K] = TestAlpha$stat
    z.stat.alpha.qsc[i-K] = sqrt(TestAlpha$stat)*sign(TestAlpha$est)
    pval.alpha.qsc[i-K] = TestAlpha$pval
    
    # # quasi-likelihood test, model log-ratio, wald test
    # fit.full.robust = geepack::geeglm(G~0+conf+Trt, id = 1:length(G))
    # stat.qsgw = summary(fit.full.robust)$coefficients["Trt", "Wald"]
    # pval.qsgw = 1 - pchisq(stat.qsgw, 1)
    # chi.stat.alpha.qsgw[i-K] = stat.qsgw
    # z.stat.alpha.qsgw[i-K] = sqrt(stat.qsgw)*sign(est)
    # pval.alpha.qsgw[i-K] = pval.qsgw

    # quasi-likelihood test, model log-ratio, score test
    mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
    est = mod.full$coefficients["Trt","Estimate"]
    mod.reduce = summary(glm(G~0+conf,family = quasi()))
    mod.reduce.disp = mod.reduce$dispersion
    mod.reduce.resid = mod.reduce$deviance.resid
    TestAlpha = .test_alpha_qsg(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
    chi.stat.alpha.qsg[i-K] = TestAlpha$stat
    z.stat.alpha.qsg[i-K] = sqrt(TestAlpha$stat)*sign(est)
    pval.alpha.qsg[i-K] = TestAlpha$pval
    
    # beta regression
    # fit.full.br =  tryCatch(betareg(G.p~0+conf+Trt, link = "logit"),
    #                         error=function(e){NA})
    # z.stat.alpha.br[i-K] = ifelse(is.na(fit.full.br), 0, summary(fit.full.br)$coefficients$mean["Trt", "z value"])
    # chi.stat.alpha.br[i-K] = ifelse(is.na(fit.full.br), 0, z.stat.alpha.br[i-K]^2)
    # pval.alpha.br[i-K] = ifelse(is.na(fit.full.br), 0, summary(fit.full.br)$coefficients$mean["Trt", "Pr(>|z|)"])
    # fit.full.br =  betareg(G.p~0+conf+Trt, link = "logit")
    # z.stat.alpha.br[i-K] = summary(fit.full.br)$coefficients$mean["Trt", "z value"]
    # chi.stat.alpha.br[i-K] = z.stat.alpha.br[i-K]^2
    # pval.alpha.br[i-K] = summary(fit.full.br)$coefficients$mean["Trt", "Pr(>|z|)"]
    
    mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
    mod.s2 = mod$sigma^2
    mod.resid = mod$residuals
    TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
    
    chi.stat.beta[i-K] = TestBeta$stat
    z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(TestBeta$est)
    pval.beta[i-K] = TestBeta$pval
  }
  
  alpha1.lr = .pi0_JC(na.omit(z.stat.alpha.lr))
  alpha1.dm = .pi0_JC(na.omit(z.stat.alpha.dm))  
  alpha1.qsc = .pi0_JC(na.omit(z.stat.alpha.qsc))
  alpha1.qsg = .pi0_JC(na.omit(z.stat.alpha.qsg))
  alpha2 = .pi0_JC(na.omit(z.stat.beta))

  tmp.lr = .nullEstimation(pval.alpha.lr, pval.beta, alpha1.lr, alpha2)
  tmp.dm = .nullEstimation(pval.alpha.dm, pval.beta, alpha1.dm, alpha2)  
  tmp.qsc = .nullEstimation(pval.alpha.qsc, pval.beta, alpha1.qsc, alpha2)
  tmp.qsg = .nullEstimation(pval.alpha.qsg, pval.beta, alpha1.qsg, alpha2)
  rawp.lr = tmp.lr$rawp
  rawp.dm = tmp.dm$rawp
  rawp.qsc = tmp.qsc$rawp
  rawp.qsg = tmp.qsg$rawp
  rawp.lr.rm = na.omit(rawp.lr)
  rawp.dm.rm = na.omit(rawp.dm)
  rawp.qsc.rm = na.omit(rawp.qsc)
  rawp.qsg.rm = na.omit(rawp.qsg)
  L.lr = length(rawp.lr.rm)
  L.dm = length(rawp.dm.rm)
  L.qsc = length(rawp.qsc.rm)
  L.qsg = length(rawp.qsg.rm)
  rawp.lr.mat[sim,] = rawp.lr
  rawp.dm.mat[sim,] = rawp.dm
  rawp.qsc.mat[sim,] = rawp.qsc
  rawp.qsg.mat[sim,] = rawp.qsg
  
  gp.mat[sim,] = c(p.hmp(rawp.lr.rm, w = rep(1/L.lr, L.lr), L = L.lr), 
                   p.hmp(rawp.dm.rm, w = rep(1/L.dm, L.dm), L = L.dm), 
                   p.hmp(rawp.qsc.rm, w = rep(1/L.qsc, L.qsc), L = L.qsc),
                   p.hmp(rawp.qsg.rm, w = rep(1/L.qsg, L.qsg), L = L.qsg))
}

rslt = list(gp.mat = gp.mat,
            rawp.lr.mat = rawp.lr.mat,
            rawp.dm.mat = rawp.dm.mat,
            rawp.qsc.mat = rawp.qsc.mat,
            rawp.qsg.mat = rawp.qsg.mat)
filename = paste0("Result/allType", causal.type, "Nsample", n.sample, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)

