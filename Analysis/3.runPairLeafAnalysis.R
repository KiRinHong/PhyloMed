# Suppose for an identified node, there are K1 and K2 leave descendants 
# under its left and right child node, respectively. 
# For each pair of left-side leave taxon and right-side leave taxon, take the log-ratio, 
# use it as the mediator, 
# and produce the mediation test p-value (JS test). 
# You can organize all mediation p-values into a K1 x K2 upper triangular matrix. 

rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")
source("1a.PhyloMed_utility.R")
generatePairLeafMat <- function(data.lst, rslt.lst, treatment, covariate, outcome){
  signode.id = as.numeric(unlist(strsplit(rslt.lst$sig.node["perm.jsmix.BH"], split = ",")))
  pval.js.mat.lst = vector("list", length(signode.id))
  names(pval.js.mat.lst) = paste0("Node", as.character(signode.id))
  tree = data.lst$tree
  M = data.lst$mediators
  M2 = M + 0.5
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  Trt.ori = data.lst$meta[[treatment]]
  outcome.ori = data.lst$meta[[outcome]]
  id.na = which(is.na(outcome.ori))
  if(length(id.na) > 0){
    Trt.ori = Trt.ori[-id.na]
    M = M[-id.na,]
    M2 = M2[-id.na,]
    outcome.ori = outcome.ori[-id.na]
  }
  n.sample = length(outcome.ori)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      conf.ori = cbind(1, data.lst$meta[[covariate]])
    }else{
      conf.ori = cbind(1, as.matrix(sample_data(data.lst$meta[,covariate])))
    }
    if(length(id.na) > 0){
      conf.ori = conf.ori[-id.na,]
    }
  }else{
    conf.ori = matrix(1, nrow = n.sample, ncol = 1)
  }
  
  for (id in 1:length(signode.id)) {
    
    child.left = treestructure$phylochildren[K+signode.id[id],1]
    child.right = treestructure$phylochildren[K+signode.id[id],2]
    Mc.left = M[,treestructure$descendant[child.left,], drop = FALSE]
    Mc.right = M[,treestructure$descendant[child.right,], drop = FALSE]
    Mc2.left = M2[,treestructure$descendant[child.left,], drop = FALSE]
    Mc2.right = M2[,treestructure$descendant[child.right,], drop = FALSE]
    
    pval.alpha.mat = pval.beta.mat = matrix(NA, nrow = ncol(Mc2.left), ncol = ncol(Mc2.right))
    for(i in 1:ncol(Mc2.left)){
      for (j in 1:ncol(Mc2.right)) {
        Mc = cbind(Mc.left[,i], Mc.right[,j])
        idx = which(rowSums(Mc) != 0)
        if(length(idx) != 0){
          Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
          G = log(Mc2.left[idx,i]/as.numeric(Mc2.right[idx,j]))
        }else{
          Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
          G = log(Mc2.left[,i]/as.numeric(Mc2.right[,j]))
        }
        # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
        condition1 = any(colSums(Mc) == 0)
        if(condition1){
          warnings("Some children have all observations being 0, skip internal node")
        }
        
        # remove the subjects with all subcomponents being zero may result in collinearity
        # the outcome may be unique after removing the subjects when is binary
        condition2 = (qr(cbind(conf,Trt))$rank < (ncol(conf)+1)) | (length(unique(outcome)) == 1)
        if(any(condition1, condition2)){
          pval.alpha.mat[i,j] = NA
          pval.beta.mat[i,j] = NA
          next
        }
        mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
        est = mod.full$coefficients["Trt","Estimate"]
        var.alpha = mod.full$cov.scaled["Trt", "Trt"]
        mod.reduce = summary(glm(G~0+conf,family = quasi()))
        mod.reduce.disp = mod.reduce$dispersion
        mod.reduce.resid = mod.reduce$deviance.resid
        TestAlpha = .test_alpha(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
        pval.alpha.mat[i,j] = TestAlpha$pval
        
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
        pval.beta.mat[i,j] = TestBeta$pval
      }
    }
    pval.js.mat = round(pmax(pval.alpha.mat, pval.beta.mat),4)
    rownames(pval.js.mat) = tree$tip.label[which(treestructure$descendant[child.left,])]
    colnames(pval.js.mat) = tree$tip.label[which(treestructure$descendant[child.right,])]
    pval.js.mat.lst[[id]] = pval.js.mat
  }
  return(pval.js.mat.lst)
}

load("../Data/Deriveddata/rslt.runModel.rda")
load("../Data/Deriveddata/Cecal.filter20top100.rda")
tmp = generatePairLeafMat(cecal.top, rslt.cecal.top100, treatment = "Treatment", covariate = NULL, outcome = "pFat")
View(tmp$Node45)
View(tmp$Node75)

tree = cecal.top$tree
tax = cecal.top$taxonomy
treestructure = .phylostructure(tree)
tmp = tax[tree$tip.label[(treestructure$descendant[.ntaxa(tree)+45,] == TRUE)],]
View(tmp)
tmp = tax[tree$tip.label[(treestructure$descendant[.ntaxa(tree)+75,] == TRUE)],]
View(tmp)

rslt.cecal.top100$node.pval["perm.jsmix",45]
rslt.cecal.top100$node.pval["perm.jsmix",75]

rslt.cecal.top100$pval.alphabeta[c("perm.alpha","perm.beta"),45]
rslt.cecal.top100$pval.alphabeta[c("perm.alpha","perm.beta"),75]

load("../Data/Deriveddata/COMBO.filter20.rda")
rslt.combo.filter20$sig.node["perm.jsmix.BH"]
tmp = generatePairLeafMat(combo.filter, rslt.combo.filter20, treatment = "fat", covariate = "calor", outcome = "bmi")
View(tmp$Node297)

tree = combo.filter$tree
tax = combo.filter$taxonomy
treestructure = .phylostructure(tree)
tmp = tax[tree$tip.label[(treestructure$descendant[.ntaxa(tree)+297,] == TRUE)],]
View(tmp)

rslt.combo.filter20$node.pval["perm.jsmix",297]
# perm.jsmix 
# 0.0001264831 
rslt.combo.filter20$pval.alphabeta[c("perm.alpha","perm.beta"),297]
# perm.alpha  perm.beta 
# 0.00603787 0.00197000 
