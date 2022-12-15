rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")
source("1a.PhyloMed_utility.R")
library(phyloseq)
library(vegan)
library(GUniFrac)
library(fdrtool)
library(SKAT)
library(harmonicmeanp)
library(matrixStats)
library(MASS)
library(mediation)
generateSensSigNode <- function(data.lst, rslt.lst, treatment, covariate, outcome){
  #data.lst=cecal.top;rslt.lst=rslt.cecal.top100.psc05;treatment = "Treatment";covariate=NULL;outcome="pFat"
  signode.id = as.numeric(unlist(strsplit(rslt.lst$sig.node["perm.jsmix.prod.BH"], split = ",")))
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
  
  #cor.sig = numeric(length(signode.id))
  res = vector("list", length = length(signode.id))
  names(res) = paste0("NodeID", signode.id)
  for (id in 1:length(signode.id)) {
    child.left = treestructure$phylochildren[K+signode.id[id],1]
    child.right = treestructure$phylochildren[K+signode.id[id],2]
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
    mod1 = lm(G~0+conf+Trt)
    alpha.resid = summary(mod1)$residuals
    alpha.est = summary(mod1)$coefficients["Trt",1]
    mod2 = lm(outcome~0+conf+Trt+G)
    beta.resid = summary(mod2)$residuals
    beta.est = summary(mod2)$coefficients["G",1]
    #cor.sig[id] = cor(alpha.resid, beta.resid)
    
    med.res = mediate(mod1, mod2, treat = "Trt", mediator = "G", sims = 1000, robustSE = TRUE, conf.level = 0.95)
    # med.res = mediate(mod1, mod2, treat = "Trt", mediator = "G", sims = 1000, boot = TRUE, conf.level = 0.95)
    #summary(med.res)
    ## sensitivity
    sens.out = medsens(med.res, rho.by = 0.01, effect.type = "indirect")
    figure.out = paste0(signode.id[id],".pdf")
    pdf(file = paste0("../Figs/", figure.out))
    plot(sens.out, sens.par = "rho", main = "Sensitivity analysis", ylim = c(-2, 2)) 
    ## check empirical rho 
    # plot(alpha.resid, beta.resid, 
    #      xlab = "Error term in mediator regression", 
    #      ylab = "Error term in outcome regression",
    #      main = "Correlation between two error terms")
    dev.off()
    res[[id]] = list(alpha.resid = alpha.resid, beta.resid = beta.resid,
                     alpha.est = alpha.est, beta.est = beta.est,
                     cor.est = cor(alpha.resid, beta.resid), 
                     med.res = med.res, sens.out = sens.out)
  }
  return(res)
}

set.seed(123)
load("../Data/Deriveddata/rslt.runModel.rda")
load("../Data/Deriveddata/Cecal.filter20top100.rda")
Cecal.Sens = generateSensSigNode(cecal.top, rslt.cecal.top100.psc05, treatment = "Treatment", covariate = NULL, outcome = "pFat")
# Node75: Rho at which ACME = 0: 0.39,
# 0.07 ~ 0.61 
# -1.71e-17
load("../Data/Deriveddata/COMBO.filter20.rda")
Combo.Sens = generateSensSigNode(combo.filter, rslt.combo.filter20.psc05, treatment = "fat", covariate = "calor", outcome = "bmi")
# Node 297: Rho at which ACME = 0: -0.44
# -0.61 ~ -0.20
# 5.32e-17
save(Cecal.Sens, Combo.Sens,  file = "../Data/Deriveddata/rslt.sens.rda")

load("../Data/Deriveddata/rslt.sens.rda")
pdf(file = "../Figs/sens.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plot(Cecal.Sens$NodeID75$sens.out, sens.par = "rho", main = "Mouse cecal study", ylim = c(-2, 2), cex.lab=1.5, cex.main = 2)
plot(Combo.Sens$NodeID297$sens.out, sens.par = "rho", main = "Human gut study", ylim = c(-2, 2), cex.lab=1.5, cex.main = 2)
dev.off()
## check empirical rho 
# plot(alpha.resid, beta.resid, 
#      xlab = "Error term in mediator regression", 
#      ylab = "Error term in outcome regression",
#      main = "Correlation between two error terms")

# runPhyloMedWPCoA <- function(input_data, exposure, covariate=NULL, outcome, PCoA.dist = "BC",
#                              fdr.alpha = 0.1, perm.prec = 0.05, seed = 84){
#   # input_data = cecal.top; exposure = "Treatment"; covariate = NULL; outcome = "pFat"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 84
#   # input_data = combo.filter; exposure = "fat"; covariate = "calor"; outcome = "bmi"; fdr.alpha = 0.1; perm.prec = 0.05; seed = 123
#   Trt = input_data$meta[[exposure]]
#   M = input_data$mediators; M2 = M + 0.5
#   Y = input_data$meta[[outcome]]
#   tree = input_data$tree
#   
#   # PCoA
#   # bray_curtis_pcoa <- ecodist::pco(m.list[[1]])
#   # bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
#   #                                   pcoa2 = bray_curtis_pcoa$vectors[,2])
#   # bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
#   #   geom_point() +
#   #   labs(x = "PC1",
#   #        y = "PC2", 
#   #        title = "Bray-Curtis PCoA") +
#   #   theme(title = element_text(size = 10)) # makes titles smaller
#   # bray_curtis_plot
#   
#   # euclidean_pcoa <- ecodist::pco(m.list[[5]])
#   # euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
#   #                                 pcoa2 = euclidean_pcoa$vectors[,2])
#   #   euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
#   #   geom_point() +
#   #   labs(x = "PC1",
#   #        y = "PC2",
#   #        title = "Euclidean PCoA with CLR transformation") +
#   #   theme(title = element_text(size = 12)) # makes titles smaller
#   # euclidean_plot
#   if(PCoA.dist == "BC"){
#     pcoa = ecodist::pco(vegdist(M, method = "bray"))
#   }else if(PCoA.dist == "JAC"){
#     pcoa = ecodist::pco(vegdist(M, 'jaccard', binary = TRUE))
#   }else if(PCoA.dist == "CLR"){
#     pcoa = ecodist::pco(vegdist(M2, method = "aitchison"))
#   }else if(PCoA.dist == "UW"){
#     unifracs = GUniFrac(M, tree)$unifracs 
#     pcoa = ecodist::pco(as.dist(unifracs[, , c('d_UW')]))
#   }else if(PCoA.dist == "W"){
#     unifracs = GUniFrac(M, tree)$unifracs 
#     pcoa = ecodist::pco(as.dist(unifracs[, , c('d_1')]))
#   }
#   
#   pcoa.df = cbind(pcoa$vectors[,1], pcoa$vectors[,2])
#   print(cumsum(pcoa$values/sum(pcoa$values))[1:2])
#   
#   gp = numeric(6)
#   names(gp) = c("Simes.PhyloMed.asym", "Fisher.PhyloMed.asym", "HMP.PhyloMed.asym",
#                 "Simes.PhyloMed.perm", "Fisher.PhyloMed.perm", "HMP.PhyloMed.perm")
#   nullprop.est = numeric(6)
#   names(nullprop.est) = c("H00.asym","H10.asym","H01.asym","H00.perm","H10.perm","H01.perm")
#   
#   set.seed(seed)
#   if(!is.null(covariate)){
#     if(length(covariate) == 1){
#       covariate = cbind(input_data$meta[[covariate]], pcoa.df)
#     }else{
#       covariate = cbind(as.matrix(sample_data(input_data$meta[,covariate])), pcoa.df)
#     }
#   }else{
#     covariate = pcoa.df
#   }
#   rslt.phylomed = phylomed(Trt, M, Y, tree, covariates = covariate, fdr.alpha = fdr.alpha, perm.prec = perm.prec, verbose = T, pseudoCount = 0.5)
#   gp[1:3] = rslt.phylomed$PhyloMed.A$global.pval[1:3]
#   nullprop.est[1:3] = rslt.phylomed$PhyloMed.A$null.prop.prod
#   node.pval.asym.jsmix = rslt.phylomed$PhyloMed.A$node.pval.jsmix
#   node.pval.asym.js = rslt.phylomed$PhyloMed.A$node.pval.js
#   node.pval.sobel = rslt.phylomed$PhyloMed.A$node.pval.sobel
#   gp[4:6] = rslt.phylomed$PhyloMed.P$global.pval[1:3]
#   nullprop.est[4:6] = rslt.phylomed$PhyloMed.P$null.prop.prod
#   node.pval.perm.jsmix = rslt.phylomed$PhyloMed.P$node.pval.jsmix
#   node.pval.perm.js = rslt.phylomed$PhyloMed.P$node.pval.js
#   
#   pval.alphabeta.mat = rbind(rslt.phylomed$PhyloMed.A$pval.alpha, rslt.phylomed$PhyloMed.P$pval.alpha,
#                              rslt.phylomed$PhyloMed.A$pval.beta, rslt.phylomed$PhyloMed.P$pval.beta)
#   rownames(pval.alphabeta.mat) = c("asym.alpha", "perm.alpha", "asym.beta", "perm.beta")
#   
#   node.pval.mat = rbind(node.pval.asym.jsmix, node.pval.asym.js, node.pval.sobel, 
#                         node.pval.perm.jsmix, node.pval.perm.js)
#   rownames(node.pval.mat) = c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js")
#   
#   sig.nodeID.BH = apply(node.pval.mat, 1, .getSigNode.BH, alpha = fdr.alpha)
#   sig.nodeID.HMP = apply(node.pval.mat, 1, .getSigNode.HMP, alpha = fdr.alpha, tree = tree)
#   sig.nodeID = c(sig.nodeID.BH, sig.nodeID.HMP)
#   names(sig.nodeID) = apply(expand.grid(c("asym.jsmix","asym.js","sobel","perm.jsmix","perm.js"), c("BH","HMP")), 1, paste, collapse=".")
#   cat("tree-based PhyloMed done\n")
#   
#   rslt = list(gp = gp, nullprop.est = nullprop.est, 
#               pval.alphabeta = pval.alphabeta.mat, 
#               node.pval = node.pval.mat,
#               sig.node = sig.nodeID, 
#               explainVar = (pcoa$values/sum(pcoa$values))[1:2])
#   return(rslt)
#   
# }
# generateSensSigNode.CondProb <- function(data.lst, rslt.lst, treatment, covariate, outcome){
#   #data.lst=cecal.top;rslt.lst=rslt.cecal.top100.psc05;treatment = "Treatment";covariate=NULL;outcome="pFat"
# 
#   signode.id = as.numeric(unlist(strsplit(rslt.lst$sig.node["perm.jsmix.prod.BH"], split = ",")))
#   signode.id = signode.id[length(signode.id)]
#   tree = data.lst$tree
#   M = data.lst$mediators
#   M2 = M + 0.5
#   K = .ntaxa(tree)
#   
#   B.max = 1e5
#   R.sel = .choose_r(0.1/K, 0.05)
#   
#   treestructure = .phylostructure(tree)
#   Trt.ori = data.lst$meta[[treatment]]
#   outcome.ori = data.lst$meta[[outcome]]
#   id.na = which(is.na(outcome.ori))
#   if(length(id.na) > 0){
#     Trt.ori = Trt.ori[-id.na]
#     M = M[-id.na,]
#     M2 = M2[-id.na,]
#     outcome.ori = outcome.ori[-id.na]
#   }
#   n.sample = length(outcome.ori)
#   if(!is.null(covariate)){
#     if(length(covariate) == 1){
#       conf.ori = cbind(1, data.lst$meta[[covariate]])
#     }else{
#       conf.ori = cbind(1, as.matrix(sample_data(data.lst$meta[,covariate])))
#     }
#     if(length(id.na) > 0){
#       conf.ori = conf.ori[-id.na,]
#     }
#   }else{
#     conf.ori = matrix(1, nrow = n.sample, ncol = 1)
#   }
#   
#   child.left = treestructure$phylochildren[K+signode.id,1]
#   child.right = treestructure$phylochildren[K+signode.id,2]
#   Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
#   Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
#   if(mean(Mc.left) < mean(Mc.right)){
#     Mc = cbind(Mc.right, Mc.left)
#   }else{
#     Mc = cbind(Mc.left, Mc.right)
#   }
#   Mc2 = Mc + 0.5
#   idx = which(rowSums(Mc) != 0)
#   if(length(idx) != 0){
#     Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
#     G = log(Mc2[idx,-2]/as.numeric(Mc2[idx,2]))
#   }else{
#     Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
#     G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
#   }
#   mod1 = glm(G~0+conf+Trt, family = quasi())
#   r1 = summary(mod1)$deviance.resid
#   sd1 = sqrt(summary(mod1)$dispersion)
#   G.pred = predict(mod1)
#   alpha.est = summary(mod1)$coefficients["Trt",1]
#   
#   mod2 = lm(outcome~0+conf+Trt+G)
#   r2 = summary(mod2)$residuals
#   sd2 = sd(r2)
#   beta.est = summary(mod2)$coefficients["G",1]
#   #cor.sig[id] = cor(alpha.resid, beta.resid)
#   
#   # compute residual forming matrix Rconf
#   Rconf = diag(nrow(conf)) - conf %*% ginv(t(conf) %*% conf) %*% t(conf)
#   
#   rho.seq = seq(-0.9, 0.9, 0.1)
#   chi.stat.alpha.mat = z.stat.alpha.mat = pval.alpha.asym.mat = pval.alpha.perm.mat = ide.est.mat = matrix(NA, length(rho.seq), 100)
#   set.seed(123)
#   for (i in 1:length(rho.seq)) {
#     print(i)
#     rho = rho.seq[i]
#     for (j in 1:100) {
#       tmp = rnorm(length(G), mean(r1)+rho*sd1/sd2*(r2-mean(r2)), sd1^2*(1-rho^2))
#       G.new = G.pred + tmp
#       
#       mod.full = summary(glm(G.new~0+conf+Trt,family = quasi()))
#       est = mod.full$coefficients["Trt","Estimate"]
#       mod.reduce = summary(glm(G.new~0+conf,family = quasi()))
#       mod.reduce.disp = mod.reduce$dispersion
#       mod.reduce.resid = mod.reduce$deviance.resid
#       TestAlpha = .test_alpha(G.new, Trt, conf, mod.reduce.resid, mod.reduce.disp)
#       stat = TestAlpha$stat
#       pval = TestAlpha$pval
#       chi.stat.alpha.mat[i,j] = stat
#       z.stat.alpha.mat[i,j] = sqrt(stat)*sign(est)
#       pval.alpha.asym.mat[i,j] = pval
#       ide.est.mat[i,j] = est*beta.est
#       #### permutation for alpha
#       if(is.infinite(stat)){
#         pval.alpha.perm[i-K] = 1/(B.max+1)
#         cat(sprintf("# of permutations for alpha: %d, Test statistic = Inf\n", B.max))
#         cat("Use Pseudo-ECDF approximation p-value\n")
#       }else{
#         m = 1
#         Nexc = 0
#         alpha.stat.perm = numeric(B.max)
#         while (Nexc < R.sel & m < B.max) {
#           TestAlpha_perm = .test_alpha(G.new, Rconf[sample(nrow(Rconf)),] %*% Trt, conf, mod.reduce.resid, mod.reduce.disp)
#           stat_perm = TestAlpha_perm$stat
#           # x2.1_perm = qr.resid(fit$qr, Rconf[sample(nrow(Rconf)),] %*% Trt)
#           # stat_perm = (sum(x2.1_perm*r))^2/sum(x2.1_perm*x2.1_perm)/dispersion
#           if(stat_perm >= stat) Nexc = Nexc + 1
#           alpha.stat.perm[m] = stat_perm
#           m = m + 1
#         }
#         if(m < B.max){
#           pval.alpha.perm.mat[i,j] = Nexc/m
#           cat(sprintf("# of permutations for alpha: %g\n", m))
#           cat("Use ECDF approximation p-value\n")
#         }else{
#           if(Nexc <= 10){
#             pval.alpha.perm.mat[i,j] = tryCatch(.gpd_approx(alpha.stat.perm, 250, stat), error=function(err) NA)
#            
#             cat(sprintf("# of permutations for alpha: %g\n", B.max))
#             cat("Use GPD approximation p-value\n")
#             
#             if(is.na(pval.alpha.perm.mat[i,j])){
#               pval.alpha.perm.mat[i,j] = (Nexc+1)/(B.max+1)
#               cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
#             }
#           }else{
#             pval.alpha.perm.mat[i,j] = Nexc/B.max
#             
#             cat(sprintf("# of permutations for alpha: %g\n", B.max))
#             cat("Use ECDF approximation p-value\n")
#           }
#         }
#       }
#     }
#   }
#   
#   return(list(chi.stat.alpha.mat = chi.stat.alpha.mat,
#               z.stat.alpha.mat = z.stat.alpha.mat,
#               pval.alpha.asym.mat = pval.alpha.asym.mat,
#               pval.alpha.perm.mat = pval.alpha.perm.mat, 
#               ide.est.mat = ide.est.mat))
# }
# 
# generateSensSelectProb <- function(data.lst, rslt.lst, sens.lst, treatment, covariate, outcome, fdr.alpha = 0.1){
#   #data.lst=combo.filter;rslt.lst=rslt.combo.filter20.psc05;treatment = "fat";covariate="calor";outcome="bmi"
#   signode.id = as.numeric(unlist(strsplit(rslt.lst$sig.node["perm.jsmix.prod.BH"], split = ",")))
#   signode.id = signode.id[length(signode.id)]
#   tree = data.lst$tree
#   M = data.lst$mediators
#   K = .ntaxa(tree)
#   treestructure = .phylostructure(tree)
#   Trt.ori = data.lst$meta[[treatment]]
#   outcome.ori = data.lst$meta[[outcome]]
#   id.na = which(is.na(outcome.ori))
#   if(length(id.na) > 0){
#     Trt.ori = Trt.ori[-id.na]
#     M = M[-id.na,]
#     M2 = M2[-id.na,]
#     outcome.ori = outcome.ori[-id.na]
#   }
#   n.sample = length(outcome.ori)
#   if(!is.null(covariate)){
#     if(length(covariate) == 1){
#       conf.ori = cbind(1, data.lst$meta[[covariate]])
#     }else{
#       conf.ori = cbind(1, as.matrix(sample_data(data.lst$meta[,covariate])))
#     }
#     if(length(id.na) > 0){
#       conf.ori = conf.ori[-id.na,]
#     }
#   }else{
#     conf.ori = matrix(1, nrow = n.sample, ncol = 1)
#   }
#   
#   
#   ## here perform tree-based method
#   chi.stat.alpha = chi.stat.beta = z.stat.alpha = z.stat.beta = pval.alpha.asym = pval.beta.asym = numeric(tree$Nnode)
#   
#   #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
#   for (i in (K + 1):(K + tree$Nnode)) {
#     child.left = treestructure$phylochildren[i,1]
#     child.right = treestructure$phylochildren[i,2]
#     
#     Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
#     Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
#     if(mean(Mc.left) < mean(Mc.right)){
#       Mc = cbind(Mc.right, Mc.left)
#     }else{
#       Mc = cbind(Mc.left, Mc.right)
#     }
#     Mc2 = Mc + 0.5
#     idx = which(rowSums(Mc) != 0)
#     if(length(idx) != 0){
#       Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
#       G = log(Mc2[idx,-2]/as.numeric(Mc2[idx,2]))
#     }else{
#       Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
#       G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
#     }
#     
#     
#     # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
#     condition1 = any(colSums(Mc) == 0)
#     if(condition1){
#       warnings("Some children have all observations being 0, skip internal node")
#     }
#     
#     # rank < 3
#     condition2 = FALSE
#     
#     # remove the subjects with all subcomponents being zero may result in collinearity
#     # the outcome may be unique after removing the subjects when is binary
#     condition3 = (qr(cbind(conf,Trt))$rank < (ncol(conf)+1)) | (length(unique(outcome)) == 1)
#     
#     if(any(condition1,condition2,condition3)){
#       chi.stat.alpha[i-K] = NA
#       z.stat.alpha[i-K] = NA
#       pval.alpha.asym[i-K] = NA
#       chi.stat.beta[i-K] = NA
#       z.stat.beta[i-K] = NA
#       pval.beta.asym[i-K] = NA
#       next
#     }
#     
#     mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
#     est = mod.full$coefficients["Trt","Estimate"]
#     var.alpha = mod.full$cov.scaled["Trt", "Trt"]
#     mod.reduce = summary(glm(G~0+conf,family = quasi()))
#     mod.reduce.disp = mod.reduce$dispersion
#     mod.reduce.resid = mod.reduce$deviance.resid
#     TestAlpha = .test_alpha(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
#     stat = TestAlpha$stat
#     pval = TestAlpha$pval
#     
#     chi.stat.alpha[i-K] = stat
#     z.stat.alpha[i-K] = sqrt(stat)*sign(est)
#     pval.alpha.asym[i-K] = pval
#     
#     mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
#     mod.s2 = mod$sigma^2
#     mod.resid = mod$residuals
#     TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
#     
#     var.beta = TestBeta$var.beta
#     chi.stat.beta[i-K] = TestBeta$stat
#     z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(sum(TestBeta$est))
#     pval.beta.asym[i-K] = TestBeta$pval
#   }
#   
#   pval.alpha.perm = rslt.lst$pval.alphabeta["perm.alpha",]
#   pval.beta.perm = rslt.lst$pval.alphabeta["perm.beta",]
#   alpha2.perm = .pi0_JC(na.omit(abs(qnorm(pval.beta.perm/2, lower.tail = FALSE))*sign(z.stat.beta)))
#   
#   globalp.perm.mat = sigNode.perm.mat = matrix(NA, nrow(sens.lst$pval.alpha.perm.mat), 100)
#   for (i in 1:nrow(globalp.perm.mat)) {
#     for (j in 1:100) {
#       pval.alpha.perm[signode.id] = sens.lst$pval.alpha.perm.mat[i,j]
#       z.stat.beta[signode.id] = sens.lst$z.stat.alpha.mat[i,j]
#       alpha1.perm = .pi0_JC(na.omit(abs(qnorm(pval.alpha.perm/2, lower.tail = FALSE))*sign(z.stat.alpha)))
#       tmp.perm = .nullEstimation_prod(pval.alpha.perm, pval.beta.perm, alpha1.perm, alpha2.perm)
#       
#       rawp.perm = tmp.perm$rawp
#       if(length(which(rawp.perm==0))>0)  rawp.perm[rawp.perm == 0] = runif(length(which(rawp.perm==0)), min = 0, max = 1e-8)
#       rawp.perm.rm = na.omit(rawp.perm)
#       L = length(rawp.perm.rm)
#       globalp.perm.mat[i,j] = p.hmp(rawp.perm.rm, w = rep(1/L, L), L = L)
#       rawp.perm.adj = p.adjust(rawp.perm, method = "BH")
#       id.sig = which(rawp.perm.adj < fdr.alpha)
#       sigNode.perm.mat[i,j] = ifelse(signode.id %in% id.sig, 1, 0)
#     }
#   }
#   return(list(globalp.perm.mat = globalp.perm.mat, sigNode.perm.mat = sigNode.perm.mat))
# }
# 
# Cecal.Sens.CondProb = generateSensSigNode.CondProb(cecal.top, rslt.cecal.top100.psc05, treatment = "Treatment", covariate = NULL, outcome = "pFat")
# Cecal.Sens.SelMat = generateSensSelectProb(cecal.top, rslt.cecal.top100.psc05, Cecal.Sens.CondProb, treatment = "Treatment", covariate = NULL, outcome = "pFat")
# Combo.Sens.CondProb = generateSensSigNode.CondProb(combo.filter, rslt.combo.filter20.psc05, treatment = "fat", covariate = "calor", outcome = "bmi")
# Combo.Sens.SelMat = generateSensSelectProb(combo.filter, rslt.combo.filter20.psc05, Combo.Sens.CondProb, treatment = "fat", covariate = "calor", outcome = "bmi")
# View(Combo.Sens.SelMat$sigNode.perm.mat)
# save(Cecal.Sens, Combo.Sens, Cecal.Sens.CondProb, Cecal.Sens.SelMat,
#      Combo.Sens.CondProb, Combo.Sens.SelMat, file = "../Data/Deriveddata/rslt.sens.rda")
