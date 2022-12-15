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
# setwd("~/Documents/Project/PhyloMed/Simulation/IncreaseEffectSize/continuous/GETSIM/shared/")
# ########## input arguments
# causal.type = 1  # can take values 1,2,3 corresponding to three approaches to simulate causal species
# sim = 1 # simulation replicate for this job (use different sim seed for different parallel jobs on CHTC)
# n.sample = 200 
# A = 1  # correspond to A in the document, control the perturbation level for species S_ab in treatment group
# B = 1 # correspond to B in the document, control the mediation effect in the outcome model for species S_ab
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
}else if(causal.type==3){
  K_ab = 15
  random.id = sample(causalNodeWLeaf, 5)
  S_ab = which(colSums(treestructure$descendant[random.id,]) == 1)
}

# 1 denotes the causal leaves
causalLeaf[S_ab] = 1
# 1 denotes the ancestor of causal species
causalNode[which(rowSums(treestructure$descendant[(K+1) : (K+tree$Nnode), S_ab]) > 0)] = 1 # the index starts from 1, not 101

#### perturb species in S_ab for subjects in treatment group
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

#### simulate outcome using species in S_ab and S_b
M2 = M + 0.5 # outcome model use M2 (no zero in it)
M2.P = M2/rowSums(M2)
beta = runif(K_ab, 0, B)
beta = beta - sum(beta)/length(beta)
beta.t = runif(1,0,C)

outcome = rnorm(n.sample, beta.t*Trt + log(M2.P[,S_ab]) %*% beta, sd=1)

gp = numeric(6)
names(gp) = c("simes.asym", "fisher.asym", "hmp.asym", "simes.perm", "fisher.perm", "hmp.perm")

# here perform PhyloMed
rslt.phylomed = phylomed(Trt, M, outcome, tree, fdr.alpha = 0.05, verbose = TRUE)
gp[1:3] = rslt.phylomed$PhyloMed.A$global.pval
gp[4:6] = rslt.phylomed$PhyloMed.P$global.pval

# no identified --> 0
if(any(A==0,B==0)){
  tmp.js.asym.BH = .calEFDRnDR.BH.H0(rslt.phylomed$PhyloMed.A$node.pval.js) 
  tmp.jsmix.asym.BH = .calEFDRnDR.BH.H0(rslt.phylomed$PhyloMed.A$node.pval.jsmix)
  tmp.js.perm.BH = .calEFDRnDR.BH.H0(rslt.phylomed$PhyloMed.P$node.pval.js)
  tmp.jsmix.perm.BH = .calEFDRnDR.BH.H0(rslt.phylomed$PhyloMed.P$node.pval.jsmix)
  tmp.sobel.BH = .calEFDRnDR.BH.H0(rslt.phylomed$PhyloMed.A$node.pval.sobel)
}else{
  tmp.js.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K) 
  tmp.jsmix.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)
}


tmp = rbind(tmp.js.asym.BH, tmp.jsmix.asym.BH, tmp.js.perm.BH, tmp.jsmix.perm.BH, tmp.sobel.BH)
efdr = tmp[,1] # empirical FDR
dr = tmp[,2] # discover rate on the most common ancestor
names(efdr) = names(dr) = c("js.asym.BH", "jsmix.asym.BH", "js.perm.BH", "jsmix.perm.BH", "sobel.BH")

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
