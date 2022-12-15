# Run from the command line via, e.g.,
#   Rscript getRslt_sim.R causal.type sim n.sample A B
rm(list=ls())

########### input arguments
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 5){
  causal.type = as.numeric(args[1]) # can take values 1,2,3 corresponding to, respectively 
  n.sample = as.numeric(args[2]) # number of samples
  sim = as.numeric(args[3]) # replicate
  A = as.numeric(args[4]) # correspond to A in the document, control the perturbation level for species S_ab in treatment group
  B = as.numeric(args[5]) # correspond to B in the document, control the mediation effect in the outcome model for species S_ab 
} else {
  cat('usage: Rscript getRslt_sim.R <causalType> <n.sample> <sim> <A> <B>\n', file = stderr())
  stop()
}

# rm(list = ls())
# setwd("~/Documents/Project/PhyloMed/Simulation/RareTaxa/continuous/GETSIM/shared/")
# ########## input arguments
# causal.type = 1  # can take values 1,2 corresponding to two approaches to simulate causal species
# sim = 1 # simulation replicate for this job (use different sim seed for different parallel jobs on CHTC)
# n.sample = 50 # sample size could be 50 or 200
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
library(LDM)

load("zeeviD_pseq_bacteria.Rdata")
seed.sim = seed + sim
set.seed(seed.sim)
n.total = nsamples(zeeviD.pseq.bact)
id.sel = sample_data(zeeviD.pseq.bact)$subjectID[sample(1:n.total, n.sample)]
pseq = subset_samples(zeeviD.pseq.bact, subjectID %in% id.sel)
pseq = filter_taxa(pseq, function(x) sum(x != 0) > 0, prune = TRUE)
count = t(otu_table(pseq))
tax = tax_table(pseq)
phy.tree = phy_tree(pseq)
zeeviD = list(count=count, tax=tax, phy.tree=phy.tree)

count = zeeviD$count
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
# # rescale to original sequencing depth
# M = round(M*N/rowSums(M))
#### simulate outcome using species in S_ab and S_b
M2 = M + 0.5 # outcome model use M2 (no zero in it)
M2.P = M2/rowSums(M2)
beta = runif(K_ab, 0, B)
beta = beta - sum(beta)/length(beta)
beta.t = runif(1, 0, C)

outcome = rbinom(n.sample, 1, 1/(1+ exp(-(beta.t*Trt + log(M2.P[,S_ab]) %*% beta))))

gp = numeric(18)
names(gp) = c("dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.clr.MedTest", "dist.omn.MedTest",
              "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.clr.MODIMA",
              "ldm.med.global",
              "simes.asym", "fisher.asym", "hmp.asym", 
              "simes.perm", "fisher.perm", "hmp.perm")

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

rslt = list(gp = gp, effNumOfTaxa = ntaxa(pseq))

filename = paste0("compType", causal.type, "Nsample", n.sample, "Sim", sim, "A", A, "B", B, ".Rdata")
save(rslt, file = filename)
