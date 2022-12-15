# focus on JC to estimate alpha1/alpha2

# Run from the command line via, e.g.,
#   Rscript getRslt_sim.R causal.type sim n.sample A B
rm(list=ls())

########### input arguments
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 5){
  causal.type = as.numeric(args[1]) # can take values 1,2 corresponding to two approaches to simulate causal species
  n.sample = as.numeric(args[2]) # number of samples
  sim = as.numeric(args[3]) # replicate
  A = as.numeric(args[4]) # correspond to A in the document, control the perturbation level for species S_ab in treatment group
  B = as.numeric(args[5]) # correspond to B in the document, control the mediation effect in the outcome model for species S_ab 
} else {
  cat('usage: Rscript getRslt_sim.R <causalType> <n.sample> <sim> <A> <B>\n', file = stderr())
  stop()
}

# rm(list = ls())
# setwd("~/Documents/Project/PhyloMed/Simulation/ProposedModel_Main/continuous_JC/GETSIM/shared/")
# ########## input arguments
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# sim = 1000 # simulation replicate for this job (use different sim seed for different parallel jobs on CHTC)
# n.sample = 200 # sample size could be 50 or 200
# A = 0  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 0 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
# The values of A and B can be changed such that the power of the mediation test is not too small/large
seed = 12 # overall simulation seed for this job 
C = 1 # correspond to C in the document, control direct effect
########### input arguments END

source("sim_utility.R")
library(SKAT)
library(MASS)
library(phyloseq)
library(matrixStats)
library(vegan)
library(GUniFrac)
library(ccmm)
library(energy)
library(fdrtool)
library(harmonicmeanp)
library(LDM)
library(qvalue) # install_github("jdstorey/qvalue")
library(HDMT)
library(HIMA)

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
causalNode = numeric(tree$Nnode); causalLeaf = numeric(K)

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

# # rescale to original sequencing depth
# M = round(M*N/rowSums(M))
#### simulate outcome using species in S_ab and S_b
M2 = M + 0.5 # outcome model use M2 (no zero in it)
M2.P = M2/rowSums(M2)
beta = runif(K_ab, 0, B)
beta = beta - sum(beta)/length(beta)
beta.t = runif(1,0,C)

outcome = rnorm(n.sample, beta.t*Trt + log(M2.P[,S_ab]) %*% beta, sd=1)

gp = numeric(20)
names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.clr.MedTest", "dist.omn.MedTest",
              "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.clr.MODIMA",
              "ldm.med.global",
              "simes.asym", "fisher.asym", "hmp.asym", 
              "simes.perm", "fisher.perm", "hmp.perm", 
              "ccmm.boot.tide", "ccmm.norm.tide")
   
## here perform distance-based
source("MedOmniTest.R")
unifracs = GUniFrac(M, tree)$unifracs 
m.list = list(BC = as.matrix(vegdist(M, method = "bray")),
              JAC = as.matrix(vegdist(M, 'jaccard', binary = TRUE)),
              UniFrac=unifracs[, , c('d_UW')], # unweighted unifrac distance doesn't work because M2 doesn't contain zero
              WUniFrac=unifracs[, , c('d_1')],
              CLR = as.matrix(vegdist(M2, method = "aitchison")))
rslt.dist = MedOmniTest(x = Trt, y = outcome, m.list = m.list, z = NULL, nperm = 2000)
gp[1] = rslt.dist$margPs[1]
gp[2] = rslt.dist$margPs[2]
gp[3] = rslt.dist$margPs[3]
gp[4] = rslt.dist$margPs[4]
gp[5] = rslt.dist$margPs[5]
gp[6] = rslt.dist$permP
print("distance-based MedTest done")

## here perform MODIMA
source("modima.R")
dist_Trt = dist(Trt)
dist_outcome = dist(outcome)
gp[7] = modima(dist_Trt, m.list[[1]], dist_outcome, nrep=2000)$p.value
gp[8] = modima(dist_Trt, m.list[[2]], dist_outcome, nrep=2000)$p.value
gp[9] = modima(dist_Trt, m.list[[3]], dist_outcome, nrep=2000)$p.value
gp[10] = modima(dist_Trt, m.list[[4]], dist_outcome, nrep=2000)$p.value
gp[11] = modima(dist_Trt, m.list[[5]], dist_outcome, nrep=2000)$p.value
print("distance-based MODIMA done")

# here perform LDM
rslt.ldm = tryCatch(ldm(M ~ Trt + outcome, test.mediation = TRUE), 
                    error = function(cond) {
                      message(cond)
                      return(NA)
                    })
if(is.na(rslt.ldm)){
  gp[12] = NA
}else{
  gp[12] = rslt.ldm$med.p.global.omni
}

# here perform PhyloMed
rslt.phylomed = phylomed(Trt, M, outcome, tree, fdr.alpha = 0.05, verbose = TRUE)
gp[13:15] = rslt.phylomed$PhyloMed.A$global.pval
gp[16:18] = rslt.phylomed$PhyloMed.P$global.pval

if(all(A>0, B>0)){
  tmp.js.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)
  
  tmp.js.asym.BY = .calEFDRnDR.BY(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.asym.BY = .calEFDRnDR.BY(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.BY = .calEFDRnDR.BY(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.BY = .calEFDRnDR.BY(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.BY = .calEFDRnDR.BY(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)
  
  tmp.js.asym.HIMA = .calEFDRnDR.HIMA(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.asym.HIMA = .calEFDRnDR.HIMA(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.HIMA = .calEFDRnDR.HIMA(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.HIMA = .calEFDRnDR.HIMA(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.HIMA = .calEFDRnDR.HIMA(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)
  
  tmp.js.asym.StoreyQ = .calEFDRnDR.StoreyQ(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.asym.StoreyQ = .calEFDRnDR.StoreyQ(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.StoreyQ = .calEFDRnDR.StoreyQ(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.StoreyQ = .calEFDRnDR.StoreyQ(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.StoreyQ = .calEFDRnDR.StoreyQ(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)
  
  tmp.js.asym.SABHA = .calEFDRnDR.SABHA(rslt.phylomed$PhyloMed.A$node.pval.js, tree, causalNode, random.id-K)
  tmp.jsmix.asym.SABHA = .calEFDRnDR.SABHA(rslt.phylomed$PhyloMed.A$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.js.perm.SABHA = .calEFDRnDR.SABHA(rslt.phylomed$PhyloMed.P$node.pval.js, tree, causalNode, random.id-K)
  tmp.jsmix.perm.SABHA = .calEFDRnDR.SABHA(rslt.phylomed$PhyloMed.P$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.sobel.SABHA = .calEFDRnDR.SABHA(rslt.phylomed$PhyloMed.A$node.pval.sobel, tree, causalNode, random.id-K)

  node.pval.js.asym = rslt.phylomed$PhyloMed.A$node.pval.js
  node.pval.js.asym[node.pval.js.asym == 0 & !is.na(node.pval.js.asym)] = runif(sum(node.pval.js.asym == 0, na.rm = TRUE), min = 0, max = 1e-7)
  node.pval.sobel.asym = rslt.phylomed$PhyloMed.A$node.pval.sobel
  node.pval.sobel.asym[node.pval.sobel.asym == 0 & !is.na(node.pval.sobel.asym)] = runif(sum(node.pval.sobel.asym == 0, na.rm = TRUE), min = 0, max = 1e-7)
  
  tmp.js.asym.HMP = .calEFDRnDR.HMP(node.pval.js.asym, tree, causalNode, random.id-K)
  tmp.jsmix.asym.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.A$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.js.perm.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.P$node.pval.js, tree, causalNode, random.id-K)
  tmp.jsmix.perm.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.P$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.sobel.HMP = .calEFDRnDR.HMP(node.pval.sobel.asym, tree, causalNode, random.id-K)
  
  input.pvals.asym = cbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.A$pval.beta)
  tmp.jsmix.asym.HDMT = .calEFDRnDR.HDMT(input.pvals.asym, rslt.phylomed$PhyloMed.A$null.prop, causalNode, random.id-K)
  input.pvals.perm = cbind(rslt.phylomed$PhyloMed.P$pval.alpha, rslt.phylomed$PhyloMed.P$pval.beta)
  tmp.jsmix.perm.HDMT = .calEFDRnDR.HDMT(input.pvals.perm, rslt.phylomed$PhyloMed.P$null.prop, causalNode, random.id-K)
  
  tmp = rbind(tmp.js.asym.BH, tmp.jsmix.asym.BH, tmp.js.perm.BH, tmp.jsmix.perm.BH, tmp.sobel.BH,
              tmp.js.asym.BY, tmp.jsmix.asym.BY, tmp.js.perm.BY, tmp.jsmix.perm.BY, tmp.sobel.BY,
              tmp.js.asym.HIMA, tmp.jsmix.asym.HIMA, tmp.js.perm.HIMA, tmp.jsmix.perm.HIMA, tmp.sobel.HIMA,
              tmp.js.asym.StoreyQ, tmp.jsmix.asym.StoreyQ, tmp.js.perm.StoreyQ, tmp.jsmix.perm.StoreyQ, tmp.sobel.StoreyQ,
              tmp.js.asym.SABHA, tmp.jsmix.asym.SABHA, tmp.js.perm.SABHA, tmp.jsmix.perm.SABHA, tmp.sobel.SABHA,
              tmp.js.asym.HMP, tmp.jsmix.asym.HMP, tmp.js.perm.HMP, tmp.jsmix.perm.HMP, tmp.sobel.HMP,
              tmp.jsmix.asym.HDMT, tmp.jsmix.perm.HDMT)
  efdr = tmp[,1] # empirical FDR
  dr = tmp[,2] # discover rate on the most common ancestor
  names(efdr) = names(dr) = c("js.asym.BH", "jsmix.asym.BH", "js.perm.BH", "jsmix.perm.BH", "sobel.BH",
                              "js.asym.BY", "jsmix.asym.BY", "js.perm.BY", "jsmix.perm.BY", "sobel.BY",
                              "js.asym.HIMA", "jsmix.asym.HIMA", "js.perm.HIMA", "jsmix.perm.HIMA", "sobel.HIMA",
                              "js.asym.StoreyQ", "jsmix.asym.StoreyQ", "js.perm.StoreyQ", "jsmix.perm.StoreyQ", "sobel.StoreyQ",
                              "js.asym.SABHA", "jsmix.asym.SABHA", "js.perm.SABHA", "jsmix.perm.SABHA", "sobel.SABHA",
                              "js.asym.HMP", "jsmix.asym.HMP", "js.perm.HMP", "jsmix.perm.HMP", "sobel.HMP",
                              "jsmix.asym.HDMT", "jsmix.perm.HDMT")
}else{
  efdr = dr = NULL
}

if(n.sample == 200){
  ## here perform ccmm
  rslt.ccmm.boot = ccmm(y = outcome, M = M2.P, tr = Trt, method.est.cov = "bootstrap", n.boot = 2000)
  rslt.ccmm.norm = ccmm(y = outcome, M = M2.P, tr = Trt, method.est.cov = "normal")
  gp[19] = ifelse(rslt.ccmm.boot$TIDE.CI[1] > 0 | rslt.ccmm.boot$TIDE.CI[2] < 0, 0, 1)
  gp[20] = 2*(1 - pnorm(abs(rslt.ccmm.norm$TIDE)/sqrt(rslt.ccmm.norm$Var.TIDE)))
  print("ccmm done")
}else{
  gp[19:20] = rep(NA, 2)
}

rslt = list(gp = gp, efdr = efdr, dr = dr,
            pval.alpha.asym = rslt.phylomed$PhyloMed.A$pval.alpha, 
            pval.beta.asym = rslt.phylomed$PhyloMed.A$pval.beta, 
            pval.alpha.perm = rslt.phylomed$PhyloMed.P$pval.alpha, 
            pval.beta.perm = rslt.phylomed$PhyloMed.P$pval.beta, 
            rawp.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js,
            rawp.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix,
            rawp.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js,
            rawp.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix,
            rawp.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel,
            causalLeaf = causalLeaf, causalNode = causalNode)

filename = paste0("compType", causal.type, "Nsample", n.sample, "Sim", sim, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)
