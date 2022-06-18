# focus on JC to estimate alpha1/alpha2

# Run from the command line via, e.g.,
#   Rscript getRslt_sim.R causal.type n.sample sim A B
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
# setwd("~/Documents/Project/PhyloMed/Simulation/ProposedModel_Main/binary_JC/GETSIM/shared/")
# ########## input arguments
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# sim = 100 # simulation replicate for this job (use different sim seed for different parallel jobs on CHTC)
# n.sample = 200 # sample size could be 50 or 200
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
library(energy)
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

outcome = rbinom(n.sample, 1, 1/(1+ exp(-(beta.t*Trt + log(M2.P[,S_ab]) %*% beta))))

gp = numeric(17)
names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
              "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
              "simes.asym", "fisher.asym", "hmp.asym", 
              "simes.perm", "fisher.perm", "hmp.perm", 
              "ccmm.boot.tide", "ccmm.norm.tide")
   
## here perform distance-based
source("MedOmniTest.R")
unifracs = GUniFrac(M, tree)$unifracs 
m.list = list(BC = as.matrix(vegdist(M, method = "bray")),
              JAC = as.matrix(vegdist(M, 'jaccard', binary = TRUE)),
              UniFrac=unifracs[, , c('d_UW')], # unweighted unifrac distance doesn't work because M2 doesn't contain zero
              WUniFrac=unifracs[, , c('d_1')])
rslt.dist = MedOmniTest(x = Trt, y = outcome, m.list = m.list, z = NULL, nperm = 2000)
gp[1] = rslt.dist$margPs[1]
gp[2] = rslt.dist$margPs[2]
gp[3] = rslt.dist$margPs[3]
gp[4] = rslt.dist$margPs[4]
gp[5] = rslt.dist$permP
print("distance-based MedTest done")

## here perform MODIMA
source("modima.R")
dist_Trt = dist(Trt)
dist_outcome = dist(outcome)
gp[6] = modima(dist_Trt, m.list[[1]], dist_outcome, nrep=2000)$p.value
gp[7] = modima(dist_Trt, m.list[[2]], dist_outcome, nrep=2000)$p.value
gp[8] = modima(dist_Trt, m.list[[3]], dist_outcome, nrep=2000)$p.value
gp[9] = modima(dist_Trt, m.list[[4]], dist_outcome, nrep=2000)$p.value
print("distance-based MODIMA done")

# here perform PhyloMed
rslt.phylomed = phylomed(Trt, M, outcome, tree, fdr.alpha = 0.05, verbose = TRUE)
gp[10:12] = rslt.phylomed$PhyloMed.A$global.pval
gp[13:15] = rslt.phylomed$PhyloMed.P$global.pval

if(all(A==1, B==1)){
  tmp.js.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.asym.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.jsmix, causalNode, random.id-K)
  tmp.js.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.js, causalNode, random.id-K)
  tmp.jsmix.perm.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.P$node.pval.jsmix, causalNode, random.id-K)
  tmp.sobel.BH = .calEFDRnDR.BH(rslt.phylomed$PhyloMed.A$node.pval.sobel, causalNode, random.id-K)

  node.pval.js.asym = rslt.phylomed$PhyloMed.A$node.pval.js
  node.pval.js.asym[node.pval.js.asym == 0 & !is.na(node.pval.js.asym)] = runif(sum(node.pval.js.asym == 0, na.rm = TRUE), min = 0, max = 1e-7)
  node.pval.sobel.asym = rslt.phylomed$PhyloMed.A$node.pval.sobel
  node.pval.sobel.asym[node.pval.sobel.asym == 0 & !is.na(node.pval.sobel.asym)] = runif(sum(node.pval.sobel.asym == 0, na.rm = TRUE), min = 0, max = 1e-7)

  tmp.js.asym.HMP = .calEFDRnDR.HMP(node.pval.js.asym, tree, causalNode, random.id-K)
  tmp.jsmix.asym.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.A$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.js.perm.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.P$node.pval.js, tree, causalNode, random.id-K)
  tmp.jsmix.perm.HMP = .calEFDRnDR.HMP(rslt.phylomed$PhyloMed.P$node.pval.jsmix, tree, causalNode, random.id-K)
  tmp.sobel.HMP = .calEFDRnDR.HMP(node.pval.sobel.asym, tree, causalNode, random.id-K)
  tmp = rbind(tmp.js.asym.BH, tmp.jsmix.asym.BH, tmp.js.perm.BH, tmp.jsmix.perm.BH, tmp.sobel.BH,
              tmp.js.asym.HMP, tmp.jsmix.asym.HMP, tmp.js.perm.HMP, tmp.jsmix.perm.HMP, tmp.sobel.HMP)
  efdr = tmp[,1] # empirical FDR
  dr = tmp[,2] # discover rate on the most common ancestor
  names(efdr) = names(dr) = c("js.asym.BH", "jsmix.asym.BH", "js.perm.BH", "jsmix.perm.BH", "sobel.BH",
                              "js.asym.HMP", "jsmix.asym.HMP", "js.perm.HMP", "jsmix.perm.HMP", "sobel.HMP")
}else{
  efdr = dr = NULL
}

## here perform cmmb
# rslt.cmmb.boot = tryCatch(cmmb(Y = outcome, M = M2.P@.Data, tr = Trt, X = NULL),
#                           error=function(e){return(NA)}) 
# gp[16] = ifelse(is.na(rslt.cmmb.boot), NA, ifelse(rslt.cmmb.boot$total["TIDE",2] > 0 | rslt.cmmb.boot$total["TIDE",3] < 0, 0, 1))
gp[16] = NA
gp[17] = NA
print("cmmb done")

rslt = list(gp = gp, efdr = efdr, dr = dr,
            rawp.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js,
            rawp.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix,
            rawp.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js,
            rawp.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix,
            rawp.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel,
            causalLeaf = causalLeaf, causalNode = causalNode)

filename = paste0("compType", causal.type, "Nsample", n.sample, "Sim", sim, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)
