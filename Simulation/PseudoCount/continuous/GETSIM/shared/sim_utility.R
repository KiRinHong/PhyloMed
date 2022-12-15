phylomed <- function(treatment, mediators, outcome, tree, covariates = NULL, interaction = FALSE,  
                     fdr.alpha = 0.05, n.perm = 1e5, perm.prec = 0.1, verbose = FALSE, pseudoCount = 0.5){
  
  M = mediators
  
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  
  Trt.ori = treatment
  outcome.ori = outcome
  n.sample = length(Trt.ori)
  if(is.null(covariates)){
    conf.ori = matrix(1, nrow = n.sample, ncol = 1)
  }else{
    conf.ori = cbind(1, covariates)
  }
  
  ## here perform tree-based method
  chi.stat.alpha = chi.stat.beta = z.stat.alpha = z.stat.beta = pval.alpha.asym = pval.beta.asym = pval.alpha.perm = pval.beta.perm = numeric(tree$Nnode)
  rawp.sobel = numeric(tree$Nnode)
  
  #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
  B.max = n.perm
  rsel = .choose_r(fdr.alpha/K, perm.prec)
  for (i in (K + 1):(K + tree$Nnode)) {
    if(verbose){
      cat(sprintf("Processing internal node #%d\n", i))
    }
    child.left = treestructure$phylochildren[i,1]
    child.right = treestructure$phylochildren[i,2]
    
    Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
    Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
    if(mean(Mc.left) < mean(Mc.right)){
      Mc = cbind(Mc.right, Mc.left)
    }else{
      Mc = cbind(Mc.left, Mc.right)
    }
    
    Mc2 = Mc + pseudoCount
    idx = which(rowSums(Mc) != 0)
    if(length(idx) != 0){
      Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
      G = log(Mc2[idx,-2]/as.numeric(Mc2[idx,2]))
    }else{
      Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
      G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
    }
    
    if(interaction){
      G = cbind(G, Trt*G)
    }
    
    # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
    condition1 = any(colSums(Mc) == 0)
    if(condition1){
      warnings("Some children have all observations being 0, skip internal node")
    }
    
    # rank < 3
    condition2 = FALSE
    if(interaction){
      condition2 = (qr(cbind(Trt,G))$rank < 3)
      if(condition2) warnings(sprintf("Matrix (T, G, Trt*G) is not full rank, skip internal node #%g", i))
    }
    
    # remove the subjects with all subcomponents being zero may result in collinearity 
    condition3 = (qr(cbind(conf,Trt))$rank < ncol(conf)+1)
    
    if(any(condition1,condition2,condition3)){
      chi.stat.alpha[i-K] = NA
      z.stat.alpha[i-K] = NA
      pval.alpha.asym[i-K] = NA
      pval.alpha.perm[i-K] = NA
      chi.stat.beta[i-K] = NA
      z.stat.beta[i-K] = NA
      pval.beta.asym[i-K] = NA
      pval.beta.perm[i-K] = NA
      rawp.sobel[i-K] = NA
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
    pval.alpha.asym[i-K] = pval
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
    # pval.alpha.asym[i-K] = pval
    
    if(length(unique(outcome)) > 2) {
      # continuous traits
      if(interaction){
        obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type="C")
        TestBeta = .test_beta(outcome, G, Trt, conf, obj = obj, test.type="vc") # est[1] ~ mediator est[2] ~ exposure * mediator
      }else{
        mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
        mod.s2 = mod$sigma^2
        mod.resid = mod$residuals
        TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
      }
    } else {
      # binary traits
      if(interaction){
        obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type="D")
        TestBeta = .test_beta(outcome, G, Trt, conf, obj = obj, test.type="vc") # est[1] ~ mediator est[2] ~ exposure * mediator
      }else{
        mod = glm(outcome~cbind(conf[,-1], Trt), family = "binomial")
        mod.est = mod$coefficients
        TestBeta = .test_beta(outcome, G, Trt, conf, est.obs = mod.est, test.type="mv")
      }
    }
    var.beta = TestBeta$var.beta
    chi.stat.beta[i-K] = TestBeta$stat
    z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(TestBeta$est)
    pval.beta.asym[i-K] = TestBeta$pval
    
    # calculate the Sobel test
    ts.sobel = (est*TestBeta$est)/sqrt(est^2*var.beta+TestBeta$est^2*var.alpha)
    rawp.sobel[i-K] = 1 - pchisq(ts.sobel^2, 1)
    
    #### permutation for alpha
    if(is.infinite(stat)){
      pval.alpha.perm[i-K] = 1/(B.max+1)
      if(verbose){
        cat(sprintf("# of permutations for alpha: %d, Test statistic = Inf\n", B.max))
        cat("Use Pseudo-ECDF approximation p-value\n")
      }
    }else{
      m = 1
      Nexc = 0
      alpha.stat.perm = numeric(B.max)
      while (Nexc < rsel & m < B.max) {
        TestAlpha_perm = .test_alpha(G, sample(Trt), conf, mod.reduce.resid, mod.reduce.disp)
        stat_perm = TestAlpha_perm$stat
        if(stat_perm >= stat) Nexc = Nexc + 1
        alpha.stat.perm[m] = stat_perm
        m = m + 1
      }
      if(m < B.max){
        pval.alpha.perm[i-K] = Nexc/(m-1)
        if(verbose){
          cat(sprintf("# of permutations for alpha: %g\n", m-1))
          cat("Use ECDF approximation p-value\n")
        }
      }else{
        if(Nexc <= 10){
          pval.alpha.perm[i-K] = .gpd_approx(alpha.stat.perm, 250, stat)
          if(verbose){
            cat(sprintf("# of permutations for alpha: %g\n", B.max))
            cat("Use GPD approximation p-value\n")
          }
          if(is.na(pval.alpha.perm[i-K])){
            pval.alpha.perm[i-K] = (Nexc+1)/(B.max+1)
            if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
          }
        }else{
          pval.alpha.perm[i-K] = Nexc/B.max
          if(verbose){
            cat(sprintf("# of permutations for alpha: %g\n", B.max))
            cat("Use ECDF approximation p-value\n")
          }
        }
      }
    }
    
    #### permutation for beta
    if(is.infinite(TestBeta$stat)){
      pval.beta.perm[i-K] = 1/(B.max+1)
      if(verbose){
        cat(sprintf("# of permutations for beta: %d, Test statistic = Inf\n", B.max))
        cat("Use Pseudo-ECDF approximation p-value\n")
      }
    }else{
      m = 1
      Nexc = 0
      beta.stat.perm = numeric(B.max)
      while (Nexc < rsel & m < B.max) {
        if(interaction){
          tmp_beta = .test_beta(outcome, G[sample(nrow(G)),], Trt, conf, obj = obj, test.type="vc") # est[1] ~ mediator est[2] ~ exposure * mediator
        }else{
          if(length(unique(outcome)) > 2){
            tmp_beta = .test_beta(outcome, sample(G), Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type="mv")
          }else{
            tmp_beta = .test_beta(outcome, sample(G), Trt, conf, est.obs = mod.est, test.type="mv")
          }
        }
        if(tmp_beta$stat >= TestBeta$stat) Nexc = Nexc + 1
        beta.stat.perm[m] = tmp_beta$stat
        m = m + 1
      }
      if(m < B.max){
        pval.beta.perm[i-K] = Nexc/(m-1)
        if(verbose){
          cat(sprintf("# of permutations for beta: %g\n", m-1))
          cat("Use ECDF approximation p-value\n")
        }
      }else{
        if(Nexc <= 10){
          pval.beta.perm[i-K] = .gpd_approx(beta.stat.perm, 250, TestBeta$stat)
          if(verbose){
            cat(sprintf("# of permutations for beta: %g\n", B.max))
            cat("Use GPD approximation p-value\n")
          }
          if(is.na(pval.beta.perm[i-K])){
            pval.beta.perm[i-K] = (Nexc+1)/(B.max+1)
            if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
          }
        }else{
          pval.beta.perm[i-K] = Nexc/B.max
          if(verbose){
            cat(sprintf("# of permutations for beta: %g\n", B.max))
            cat("Use ECDF approximation p-value\n")
          }
        }
      }
    }
  }
  
  alpha1.asym = .pi0_JC(na.omit(z.stat.alpha))
  alpha2.asym = .pi0_JC(na.omit(z.stat.beta))
  alpha1.perm = .pi0_JC(na.omit(abs(qnorm(pval.alpha.perm/2, lower.tail = FALSE))*sign(z.stat.alpha)))
  alpha2.perm = .pi0_JC(na.omit(abs(qnorm(pval.beta.perm/2, lower.tail = FALSE))*sign(z.stat.beta)))
  
  rawp.asym.js = pmax(pval.alpha.asym, pval.beta.asym)
  tmp.asym = .nullEstimation(pval.alpha.asym, pval.beta.asym, alpha1.asym, alpha2.asym)
  rawp.asym = tmp.asym$rawp
  if(length(which(rawp.asym==0))>0)  rawp.asym[rawp.asym == 0] = runif(length(which(rawp.asym==0)), min = 0, max = 1e-7)
  null.prop.est.asym = c(tmp.asym$alpha00, tmp.asym$alpha10, tmp.asym$alpha01)
  names(null.prop.est.asym) = c("H00","H10","H01")
  p.asym.adj = p.adjust(rawp.asym, method = "BH")
  sig.nodeID.asym = which(p.asym.adj < fdr.alpha)
  rawp.asym.rm = na.omit(rawp.asym)
  L = length(rawp.asym.rm)
  globalp.asym = c(min(L * rawp.asym.rm/rank(rawp.asym.rm)), 
                   1 - pchisq(-2 * sum(log(rawp.asym.rm)), df = 2 * L), 
                   p.hmp(rawp.asym.rm, w = rep(1/L, L), L = L))
  names(globalp.asym) = c("Simes", "Fisher", "HMP")
  
  rawp.perm.js = pmax(pval.alpha.perm, pval.beta.perm)
  tmp.perm = .nullEstimation(pval.alpha.perm, pval.beta.perm, alpha1.perm, alpha2.perm)
  rawp.perm = tmp.perm$rawp
  if(length(which(rawp.perm==0))>0)  rawp.perm[rawp.perm == 0] = runif(length(which(rawp.perm==0)), min = 0, max = 1e-7)
  null.prop.est.perm = c(tmp.perm$alpha00, tmp.perm$alpha10, tmp.perm$alpha01)
  names(null.prop.est.perm) = c("H00","H10","H01")
  p.perm.adj = p.adjust(rawp.perm, method = "BH")
  sig.nodeID.perm = which(p.perm.adj < fdr.alpha)
  
  rawp.perm.rm = na.omit(rawp.perm)
  globalp.perm = c(min(L * rawp.perm.rm/rank(rawp.perm.rm)), 
                   1 - pchisq(-2 * sum(log(rawp.perm.rm)), df = 2 * L), 
                   p.hmp(rawp.perm.rm, w = rep(1/L, L), L = L))
  names(globalp.perm) = c("Simes", "Fisher", "HMP")
  
  rslt = list(PhyloMed.A = list(pval.alpha = pval.alpha.asym, pval.beta = pval.beta.asym,
                                node.pval.jsmix = rawp.asym, node.pval.js = rawp.asym.js, node.pval.sobel = rawp.sobel,
                                sig.node = sig.nodeID.asym, null.prop = null.prop.est.asym, global.pval = globalp.asym),
              PhyloMed.P = list(pval.alpha = pval.alpha.perm, pval.beta = pval.beta.perm,
                                node.pval.jsmix = rawp.perm, node.pval.js = rawp.perm.js, 
                                sig.node = sig.nodeID.perm, null.prop = null.prop.est.perm, global.pval = globalp.perm))
  return(rslt)
}


#### INTERNAL FUNCTIONS ####

# .ntaxa
.ntaxa <- function(tree){
  length(tree$tip.label)
}
# .is_rooted
.is_rooted <- function(tree, K){
  if(!is.null(tree$root.edge)) return(TRUE)
  if(tabulate(tree$edge[,1])[K+1]>2) FALSE else TRUE
}
# .is_binary
.is_binary <- function(tree){
  .ntaxa(tree)-tree$Nnode+.is_rooted(tree,.ntaxa(tree)) == 2
}

.phylostructure <- function (tree) {
  K = .ntaxa(tree)
  phyloparent <- numeric(tree$Nnode + K)
  phylochildren <- matrix(0, tree$Nnode + K, 2)
  for (i in 1:nrow(tree$edge)) {
    i1 <- tree$edge[i, 1]
    i2 <- tree$edge[i, 2]
    if (i1 <= K) 
      stop(sprintf("Internal node label is not larger than %d", K))
    if (i2 > K && i1 > i2) 
      stop("Parent node label is larger than child internal node label")
    
    phyloparent[i2] <- i1
    if (phylochildren[i1, 1] == 0) 
      phylochildren[i1, 1] <- i2
    else phylochildren[i1, 2] <- i2
  }
  descendant <- matrix(FALSE, tree$Nnode + K, K)
  for (i in 1:K) descendant[i, i] <- TRUE
  processed <- logical(tree$Nnode + K)
  processed[1:K] <- TRUE
  while (!all(processed)) {
    for (i in (K + 1):(K + tree$Nnode)) {
      if (all(processed[phylochildren[i, ]])) {
        descendant[i, descendant[phylochildren[i, 1], ]] <- TRUE
        descendant[i, descendant[phylochildren[i, 2], ]] <- TRUE
        processed[i] <- TRUE
      }
    }
  }
  list(phylotree = tree, phylochildren = phylochildren, 
       phyloparent = phyloparent, descendant = descendant)
}
# use Nexc most exterme test statistic to approximate GPD 
.gpd_approx <- function(yperm, nexc, obsts){
  yperm = na.omit(yperm)
  y = sort(yperm, decreasing = TRUE)
  tmp = .gpd_params_est(nexc, y)
  a_hat = tmp[1]
  k_hat = tmp[2]
  t = tmp[3]
  ### goodness-fit test
  nexc_re = .gpd_goft(nexc, y[1:nexc]-t)
  if(nexc_re[2] == nexc){
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc)
    return(p)
  }else{
    tmp = .gpd_params_est(nexc_re[2], y)
    a_hat = tmp[1]
    k_hat = tmp[2]
    t = tmp[3]
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc_re[2])
    return(p)
  }
}

# Parameter estimators for the generalized Pareto distribution
.gpd_params_est <- function(nexc, y){
  t = (y[nexc] + y[nexc+1])/2
  z = y[1:nexc] - t
  z2 = z^2
  m = mean(z)
  m2 = mean(z2)
  a_hat = m*m2/(m2-m^2)/2
  k_hat = (m^2/(m2-m^2)-1)/2
  return(c(a_hat, k_hat, t))
}

# p-value for the generalized Pareto distribution approximation
.gpd_pval <- function(a_hat, k_hat, z0, nperm, nexc){
  p = nexc/nperm*((1-k_hat*z0/a_hat)**(1/k_hat))
  return(p)
}

# iteratively reduce Nexc by 10 until a goodness-of-fit satisfy
.gpd_goft <- function(nexc, y){
  # y: the sorted test statistic - threshold t
  nexc = length(y)
  p = .gp_test(y)
  
  nexc.c = seq(0,nexc-10,10)
  z = y
  i = 0
  
  re = c()
  for(i in nexc.c) {
    z = y[1:(nexc-i)]
    p = .gp_test(z)
    re = rbind(re, c(i, p))
    if(!is.na(p) & p > 0.05) break
    i = i + 10
  }
  
  if(nrow(re) >=2) {
    nexc.c2 = seq(re[(nrow(re)-1),1]+1, re[nrow(re),1],1)
    re = c()
    for(i in nexc.c2) {
      z = y[1:(nexc-i)]
      p = .gp_test(z)
      re = rbind(re, c(i, p))
      if(!is.na(p) & p > 0.05) break
      i = i + 10
    }
  }
  
  p = re[nrow(re),2]
  len = nexc-re[nrow(re),1]
  
  return(c(p, len))
}

# Bootstrap goodness-of-fit test for the Generalized Pareto distribution
.gp_test <- function(x, B = 2999){
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   #  sample size without NA values
  samplerange <- max(x) - min(x)
  gammap <- .amle_method(x, k = ceiling(.2 * n))[1]    
  gamman <- .combined_method(x)[1]
  r1 <- .R1(x)     # observed value of R^-
  r2 <- .R2(x)     # observed value of R^+
  p.value1 <- sum(replicate(B, .R1(.rgp(n, shape = gamman))) < r1) / B  # bootstrap p-value for H_0^- 
  p.value2 <- sum(replicate(B, .R2(.rgp(n, shape = gammap))) < r2) / B  # bootstrap p-value for H_0^+ 
  p.value  <- max(p.value1, p.value2)    # p-value of the intersection-union test
  return(p.value)
}

# Asymptotic maximum likelihood estimators 
.amle_method <- function(x, k){
  x  <- sort(x)
  n  <- length(x)
  nk <- n - k
  x1 <- x[(nk+1):n]
  w  <- log(x1)
  g  <-  - (w[1] - sum(w) / k)  
  sigma <- g * exp(w[1] + g * log(k / n))
  return(c(g, sigma))
}

# Combined estimators 
.combined_method <- function(x){
  m     <- mean(x)
  maxi  <- max(x)
  g     <- m / (m - maxi) 
  sigma <- - g * maxi 
  return(c(g, sigma))
}

# Test statistic for H_0^-
.R1 <- function(x){
  gamma_neg <- .combined_method(x)[1]
  Fn        <- ecdf(x)
  x1        <- x[x != max(x)]
  z1        <- (1 - Fn(x1))^( - gamma_neg) 
  return(abs(cor(x1, z1)))
}

# Test statistic for H_0^+
.R2  <- function(x){
  n              <- length(x)
  Fn             <- ecdf(x)
  gamma_positive <- .amle_method(x, ceiling(.2 * n))[1]
  x1             <- x[x != max(x)]
  y1             <- (1 - Fn(x1))^( - gamma_positive) 
  x.star         <- log(x1)
  y.star         <- log( y1 -1 )
  if (gamma_positive <= 0.5)	return(cor(x1, y1))  
  if (gamma_positive >  0.5)  return((cor(x.star, y.star)))
}

# Simulation of random numbers from the gPd
.rgp  <- function (n, shape){
  if (shape != 0) 
    return((1 / shape) * (runif(n)^(-shape) - 1))
  else return(rexp(n, 1))
}

# Test alpha in mediator model
.test_alpha <- function(G, Trt, covariates, resid, disp){
  Trt = matrix(Trt, nrow = length(Trt), ncol = 1)
  X1 = model.matrix(G~0+covariates)
  D.resid2 = diag(resid^2)/disp
  A11 = t(X1) %*% X1; B11 = t(X1) %*% D.resid2 %*% X1
  A12 = t(X1) %*% Trt; B12 = t(X1) %*% D.resid2 %*% Trt
  A21 = t(Trt) %*% X1; B21 = t(Trt) %*% D.resid2 %*% X1
  A22 = t(Trt) %*% Trt; B22 = t(Trt) %*% D.resid2 %*% Trt
  # continuous traits
  W = B22 - A21 %*% solve(A11) %*% B12 - B21 %*% solve(A11) %*% A12 + 
    A21 %*% solve(A11) %*% B11 %*% solve(A11) %*% A12
  U = as.vector(t(Trt) %*% resid)/disp
  V = W/disp
  stat = as.numeric(U %*% solve(V) %*% U)
  pval = 1 - pchisq(stat, 1) 
  return(list(stat=stat, pval=pval))
}
# Test beta in outcome model
.test_beta <- function(outcome, G, Trt, conf, resid.obs=NULL, s2.obs=NULL, est.obs=NULL, obj=NULL, test.type="mv"){
  m = ncol(G)
  n = nrow(G)
  if(is.vector(G)){
    m = 1
    n = length(G)
  }
  ## get U and V
  if(test.type=="mv"){
    tmp = .cal_sumstats(outcome, G, cbind(conf[,-1], Trt), resid.obs, s2.obs, est.obs)
    U = tmp$U
    V = tmp$V
    
    stat = as.numeric(U %*% solve(V) %*% U)
    pval = 1 - pchisq(stat, m) 
  }
  
  if(test.type=="vc"){
    pval = SKAT(G, obj, is_check_genotype=FALSE, kernel="linear")$p.value
    stat = qchisq(1-pval, df=1)
  }
  ## estimate parameters
  if(length(unique(outcome)) > 2 ){
    # continuous trait
    tmp.mod = lm(outcome ~ G + cbind(conf[,-1], Trt))
    est = summary(tmp.mod)$coefficients[1+1:m,1]
    var.beta = vcov(tmp.mod)[1+1:m,1+1:m]
  }else{
    # binary trait
    tmp.mod = glm(outcome ~ G + cbind(conf[,-1], Trt), family = "binomial")
    est = summary(tmp.mod)$coefficients[1+1:m,1]
    var.beta = vcov(tmp.mod)[1+1:m,1+1:m]
  }
  
  return(list(stat=stat, est=est, pval=pval, var.beta=var.beta))
}

# modify Lan's code to Bowen: cal_score_sumstats.R 07/09/2019
# continuous: resid/s2 binary: est
.cal_sumstats <- function(Y, G, covariates, resid, s2, est){
  m = ncol(G)
  n = nrow(G)
  if(is.vector(G)){
    m = 1
    n = length(G)
  }
  if(!is.matrix(G)){
    G = matrix(G, nrow = n, ncol = m)
  }
  
  X1 <- model.matrix(Y~covariates)
  if(length(unique(Y)) > 2 ) {
    # continuous traits
    W <- t(G) %*% G-
      (t(G) %*% X1) %*%
      solve(t(X1) %*% X1) %*% (t(X1) %*% G)
    
    obs.stat <- list(U = as.vector(t(G) %*% resid) / s2, V = W/s2)
  } else {
    # binary traits
    mod <- glm(Y~covariates, family = "binomial")
    sum.U <- numeric(m)
    sum1 <- matrix(0, m, m)
    sum2 <- matrix(0, m, ncol(X1))
    sum3 <- matrix(0, ncol(X1), ncol(X1))
    sum4 <- matrix(0, ncol(X1), m)
    for(i in 1:n){
      lambdaX <- as.numeric(est %*% X1[i, ])
      b1 <- exp(lambdaX) / (1 + exp(lambdaX))
      b2 <- exp(lambdaX) / ((1 + exp(lambdaX)))^2
      
      U.part1 <- (Y[i] - b1) * G[i, ]
      U.part1[is.na(U.part1)] <- 0
      V.part1 <- b2 * G[i, ] %*% t(G[i, ])
      V.part1[is.na(V.part1)] <- 0
      V.part2 <- b2 * G[i, ] %*% t(X1[i, ])
      V.part2[is.na(V.part2)] <- 0
      V.part3 <- b2 * X1[i, ] %*% t(X1[i, ])
      V.part3[is.na(V.part3)] <- 0
      V.part4 <- b2 * X1[i, ] %*% t(G[i, ])
      V.part4[is.na(V.part4)] <- 0
      
      sum.U <- sum.U + U.part1
      sum1 <- sum1 + V.part1
      sum2 <- sum2 + V.part2
      sum3 <- sum3 + V.part3
      sum4 <- sum4 + V.part4
    }
    
    obs.stat <- list(U = sum.U, V = sum1 - sum2 %*% solve(sum3) %*% sum4)
  }
  
  return(obs.stat)
}

.pi0_JC <- function(z){
  xi = c(0:100)/100
  tmax = sqrt(log(length(z)))
  tt = seq(0, tmax, 0.05)
  
  epsest = NULL
  
  for (j in 1:length(tt)) {
    t = tt[j]
    f = t*xi
    f = exp(f^2/2)
    w = (1 - abs(xi))
    co = 0*xi
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z))
    }
    epshat = sum(w*f*co)/sum(w)
    epsest = c(epsest,epshat)
  }
  tmp = min(epsest)
  if(tmp > 1) tmp = 1
  return(tmp)
}

.nullEstimation <- function (pval.alpha, pval.beta, alpha1, alpha2) {
  # alpha00: a=0 and b=0
  # alpha01: a=0 and b!=0
  # alpha10: a!=0 and b=0
  input.pvals = cbind(pval.alpha, pval.beta)
  idx.na = which(!complete.cases(input.pvals))
  input.pvals = input.pvals[complete.cases(input.pvals), ]
  
  alpha00 = alpha1 * alpha2
  alpha01 = alpha1 * (1-alpha2)
  alpha10 = (1-alpha1) * alpha2
  w = alpha00+alpha01+alpha10
  alpha00 = alpha00/w
  alpha01 = alpha01/w
  alpha10 = alpha10/w
  
  tmp = pmax(input.pvals[,1], input.pvals[,2])
  nmed = length(tmp)
  cdf12 = input.pvals
  input.pvals = input.pvals + runif(tmp, min = 0, max = 1e-7)
  xx1 = c(0, input.pvals[order(input.pvals[, 1]), 1])
  yy1 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit1 = gcmlcm(xx1, yy1, type = "lcm")
  xknots1 = gfit1$x.knots[-1]
  Fknots1 = cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
  xx2 = c(0, input.pvals[order(input.pvals[, 2]), 2])
  yy2 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit2 = gcmlcm(xx2, yy2, type = "lcm")
  xknots2 = gfit2$x.knots[-1]
  Fknots2 = cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
  if (alpha1 != 1) 
    Fknots1 = (Fknots1 - alpha1 * xknots1)/(1 - alpha1)
  else Fknots1 = rep(0, length(xknots1))
  if (alpha2 != 1) 
    Fknots2 = (Fknots2 - alpha2 * xknots2)/(1 - alpha2)
  else Fknots2 = rep(0, length(xknots2))
  orderq1 = orderq2 = gcdf1 = gcdf2 = tmp
  for (i in 1:length(xknots1)) {
    if (i == 1) {
      gcdf1[orderq1 <= xknots1[i]] = (Fknots1[i]/xknots1[i]) * orderq1[orderq1 <= xknots1[i]]
    }
    else {
      if (sum(orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]) > 0) {
        temp = orderq1[orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]]
        gcdf1[orderq1 > xknots1[i - 1] & orderq1 <= 
                xknots1[i]] = Fknots1[i - 1] + (Fknots1[i] - 
                                                  Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                                                                                          1]) * (temp - xknots1[i - 1])
      }
    }
  }
  for (i in 1:length(xknots2)) {
    if (i == 1) {
      gcdf2[orderq2 <= xknots2[i]] = (Fknots2[i]/xknots2[i]) * orderq2[orderq2 <= xknots2[i]]
    }
    else {
      if (sum(orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]) > 0) {
        temp = orderq2[orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]]
        gcdf2[orderq2 > xknots2[i - 1] & orderq2 <= 
                xknots2[i]] = Fknots2[i - 1] + (Fknots2[i] - 
                                                  Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                                                                                          1]) * (temp - xknots2[i - 1])
      }
    }
  }
  gcdf1 = ifelse(gcdf1 > 1, 1, gcdf1)
  gcdf2 = ifelse(gcdf2 > 1, 1, gcdf2)
  cdf12[, 1] = gcdf1
  cdf12[, 2] = gcdf2
  rawp = (tmp * cdf12[, 2] * alpha01) + (tmp * cdf12[, 1] * alpha10) + (tmp^2 * alpha00)
  
  rawp.wNA = numeric(length(pval.alpha))
  if(length(idx.na) > 0){
    rawp.wNA[idx.na] = NA
    rawp.wNA[-idx.na] = rawp
  }else{
    rawp.wNA = rawp
  }
  rslt = list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2, rawp = rawp.wNA)
  return(rslt)
}


# generate empirical FDR and discover rate of most common ancestor of causalLeaves
.calEFDRnDR.BH <- function(pval, causalNode, ancestorNode, fdr.alpha = 0.05){
  pval.adj = p.adjust(pval, method = "BH")
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}

.calEFDRnDR.BY <- function(pval, causalNode, ancestorNode, fdr.alpha = 0.05){
  pval.adj = p.adjust(pval, method = "BY")
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}

# Mediation effect selection in high-dimensional and compositional microbiome data
.calEFDRnDR.HIMA <- function(pval, causalNode, ancestorNode, fdr.alpha = 0.05){
  # pval = rslt.phylomed$PhyloMed.A$node.pval.sobel; ancestorNode = random.id-K; fdr.alpha = 0.05
  pval.rm = na.omit(pval)
  set = which(pval < fdr.alpha)
  set.rm = which(pval.rm < fdr.alpha)
  hom = hommel::hommel(pval.rm, simes = FALSE)
  
  if (length(set.rm) > 0){
    N1 = hommel::discoveries(hom, set.rm, incremental = TRUE, alpha=0.05)
    L = length(set.rm)
    N2 = matrix(0,1,L)
    N2[2:L] = N1[1:(L-1)]
    N0 = N1 - N2
    ID_FDR = set[which(N0 > 0)]
    if(length(ID_FDR) == 0){
      fdr = 0
    }else{
      tp = length(intersect(which(causalNode == T), ID_FDR))
      fp = length(intersect(which(causalNode == F), ID_FDR))
      fdr = fp / (tp + fp)
    }
  }else{
    ID_FDR = NULL
    fdr = 0
  }
  
  power.mca = mean(ancestorNode %in% ID_FDR)
  return(c(fdr, power.mca))
}

.calEFDRnDR.StoreyQ <- function(pval, causalNode, ancestorNode, fdr.alpha = 0.05){
  # pval = rslt.phylomed$PhyloMed.A$node.pval.jsmix
  pval.adj = qvalue(pval, lambda = 0.5)$qvalues
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}

# A Multiple Testing Procedure for High-Dimensional Mediation Hypotheses, JASA
.calEFDRnDR.HDMT <- function(input.pvals, nullprop, causalNode, ancestorNode, fdr.alpha = 0.05){
  # input.pvals=input.pvals.perm; nullprop=rslt.phylomed$PhyloMed.P$null.prop; ancestorNode=random.id-K
  id.na = which(!complete.cases(input.pvals))
  input.pvals = input.pvals[complete.cases(input.pvals),]
  id.tmp1 = which(input.pvals != 1); id.tmp2 = which(input.pvals == 1)
  if(length(id.tmp1) > 0)
    input.pvals[id.tmp1] = input.pvals[id.tmp1] + runif(length(id.tmp1), min = 0, max = 1e-7)
  if(length(id.tmp2) > 0)
    input.pvals[id.tmp2] = input.pvals[id.tmp2] - runif(length(id.tmp2), min = 0, max = 1e-7)
  alpha00 = nullprop["H00"]; alpha01 = nullprop["H01"]; alpha10 = nullprop["H10"]
  alpha1 = alpha00+alpha01; alpha2 = alpha00+alpha10
  pval.adj = .fdr_est(alpha00, alpha01, alpha10, alpha1, alpha2, input.pvals, exact=1)
  for (k in id.na) {
    pval.adj = append(pval.adj, NA, k-1)
  }
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}

# adapted from HDMT, change the q-value function with lambda = 0.5 when exact = 1 (finite sample adjustment)
# otherwise, Error in smooth.spline(lambda, pi0, df = smooth.df) : 
# missing or infinite values in inputs are not allowed, which means cannot estimate pi0 in this case
.fdr_est <- function (alpha00, alpha01, alpha10, alpha1, alpha2, input_pvalues, 
                      exact = 0) 
{
  if (is.null(ncol(input_pvalues))) 
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2) 
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues)) 
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues), 
  ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) < 
      1) 
    stop("input_pvalues doesn't have valid p-values")
  pmax <- apply(input_pvalues, 1, max)
  nmed <- length(pmax)
  efdr1 <- rep(0, nmed)
  if (exact == 0) {
    for (i in 1:nmed) {
      fdr11 <- (pmax[i] * alpha01)/mean(pmax <= pmax[i])
      fdr12 <- (pmax[i] * alpha10)/mean(pmax <= pmax[i])
      fdr2 <- (pmax[i] * pmax[i] * alpha00)/mean(pmax <= pmax[i])
      efdr1[i] <- fdr11 + fdr12 + fdr2
    }
  }
  if (exact == 1) {
    ish11 <- qvalue(input_pvalues[, 1], lambda = 0.5)$qvalue < 0.25 & qvalue(input_pvalues[, 2], lambda = 0.5)$qvalue < 0.25
    out1 <- input_pvalues[!ish11, ]
    nmed1 <- nrow(out1)
    nmed <- nrow(input_pvalues)
    cdf12 <- input_pvalues
    xx1 <- c(0, out1[order(out1[, 1]), 1])
    yy1 <- c(0, seq(1, nmed1, by = 1)/nmed1)
    gfit1 <- gcmlcm(xx1, yy1, type = "lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
    xx2 <- c(0, out1[order(out1[, 2]), 2])
    yy2 <- c(0, seq(1, nmed1, by = 1)/nmed1)
    gfit2 <- gcmlcm(xx2, yy2, type = "lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
    alpha1 <- (alpha00 + alpha01)/(alpha00 + alpha01 + alpha10)
    alpha2 <- (alpha00 + alpha10)/(alpha00 + alpha01 + alpha10)
    if (alpha1 != 1) 
      Fknots1 <- (Fknots1 - alpha1 * xknots1)/(1 - alpha1)
    else Fknots1 <- rep(0, length(xknots1))
    if (alpha2 != 1) 
      Fknots2 <- (Fknots2 - alpha2 * xknots2)/(1 - alpha2)
    else Fknots2 <- rep(0, length(xknots2))
    orderq1 <- pmax
    orderq2 <- pmax
    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i == 1) {
        gcdf1[orderq1 <= xknots1[i]] <- (Fknots1[i]/xknots1[i]) * 
          orderq1[orderq1 <= xknots1[i]]
      }
      else {
        if (sum(orderq1 > xknots1[i - 1] & orderq1 <= 
                xknots1[i]) > 0) {
          temp <- orderq1[orderq1 > xknots1[i - 1] & 
                            orderq1 <= xknots1[i]]
          gcdf1[orderq1 > xknots1[i - 1] & orderq1 <= 
                  xknots1[i]] <- Fknots1[i - 1] + (Fknots1[i] - 
                                                     Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                                                                                             1]) * (temp - xknots1[i - 1])
        }
      }
    }
    for (i in 1:length(xknots2)) {
      if (i == 1) {
        gcdf2[orderq2 <= xknots2[i]] <- (Fknots2[i]/xknots2[i]) * 
          orderq2[orderq2 <= xknots2[i]]
      }
      else {
        if (sum(orderq2 > xknots2[i - 1] & orderq2 <= 
                xknots2[i]) > 0) {
          temp <- orderq2[orderq2 > xknots2[i - 1] & 
                            orderq2 <= xknots2[i]]
          gcdf2[orderq2 > xknots2[i - 1] & orderq2 <= 
                  xknots2[i]] <- Fknots2[i - 1] + (Fknots2[i] - 
                                                     Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                                                                                             1]) * (temp - xknots2[i - 1])
        }
      }
    }
    gcdf1 <- ifelse(gcdf1 > 1, 1, gcdf1)
    gcdf2 <- ifelse(gcdf2 > 1, 1, gcdf2)
    cdf12[, 1] <- gcdf1
    cdf12[, 2] <- gcdf2
    for (i in 1:nmed) {
      fdr11 <- (pmax[i] * cdf12[i, 2] * alpha01)/mean(pmax <= pmax[i])
      fdr12 <- (pmax[i] * cdf12[i, 1] * alpha10)/mean(pmax <= pmax[i])
      fdr2 <- (pmax[i] * pmax[i] * alpha00)/mean(pmax <= pmax[i])
      efdr1[i] <- fdr11 + fdr12 + fdr2
    }
  }
  efdr1.order <- efdr1[order(pmax, decreasing = T)]
  for (i in 2:nmed) {
    efdr1.order[i] <- min(efdr1.order[i], efdr1.order[i - 1])
  }
  efdr1 <- efdr1.order[rank(-pmax)]
  return(efdr = efdr1)
}
# structure-adaptive Benjamini–Hochberg algorithm
.calEFDRnDR.SABHA <- function(pval, tree, causalNode, ancestorNode, fdr.alpha = 0.05){
  edges = tree$edge - .ntaxa(tree) 
  edges = edges[apply(edges, 1, function(x) all(x>0)),]
  id.na = which(is.na(pval))
  if(length(id.na) == 0){
    edges.new = edges
  }else{
    edges.new = edges
    for (k in 1:length(id.na)) {
      id.na.sel = id.na[k]-k+1
      row.rm = which(edges.new == id.na.sel) %% nrow(edges.new)
      row.rm = ifelse(row.rm == 0, nrow(edges.new), row.rm)
      if(id.na.sel == 1){
        edges.rm = edges.new[row.rm,,drop=FALSE]
        edges.new = rbind(edges.new[-row.rm,],c(edges.rm[,2])) 
      }else{
        if(length(row.rm) == 1){
          edges.new = edges.new[-row.rm,]
        }else{
          edges.rm = edges.new[row.rm,,drop=FALSE]
          edges.reorg = edges.rm[-nrow(edges.rm),,drop=FALSE]
          edges.reorg[,1] = edges.rm[nrow(edges.rm),1]
          edges.new = rbind(edges.new[-row.rm,], edges.reorg)
        }
      }
      edges.new[which(edges.new > id.na.sel)] = edges.new[which(edges.new > id.na.sel)]-1
    }
  }
  edges.new = edges.new[order(edges.new[,2]),]
  
  ADMM_params = c(10^2, 10^3, 5, 5000, 1e-4) # alpha_ADMM, beta, eta, max_iters, converge_thr (parameters for ADMM)
  tau = 0.5; eps = 0.1 # parameters for SABHA total variation method
  TV_bd_seq = 0.5 # value of "m" for the SABHA total variation norm bound
  
  pval.rm = na.omit(pval)
  est_q = .Solve_q_TV(pval.rm, tau, eps, edges.new, TV_bd_seq, ADMM_params)
  tmp = .SABHA_method(pval.rm, est_q, fdr.alpha, tau)
  for (k in id.na) {
    tmp = append(tmp, NA, k-1)
  }
  pval.adj = tmp
  
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}
.calEFDRnDR.HMP <- function(pval, tree, causalNode, ancestorNode, fdr.alpha = 0.05){
  edge = tree$edge
  tips = setdiff(edge[,2], edge[,1])
  K = .ntaxa(tree)
  L = length(pval)
  W = rep(1/L, L)
  pval.adj = numeric(L)
  for(node in (K+1):(K+tree$Nnode)){
    x = edge[edge[,1] == node, 2]
    repeat{
      xx = x
      x = sort(unique(c(x, edge[,2][edge[,1] %in% x])))
      if(identical(x, xx)) break
    }
    x = setdiff(x, tips)
    x = c(node, x)
    if(any(is.na(pval[x-K]))){
      na.id = which(is.na(pval[x-K]))
      pval.rm = pval[x-K][-na.id]
      W.rm = W[x-K][-na.id]
      pval.adj[node-K] = p.hmp(pval.rm, W.rm, length(x)-length(na.id))/sum(W.rm)
    }else{
      pval.adj[node-K] = p.hmp(pval[x-K], W[x-K], length(x))/sum(W[x-K])
    }
    
  }
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = 0
  }else{
    fdr = fp / (tp + fp)
  }
  power.mca = mean(pval.adj[ancestorNode] < fdr.alpha)# most common ancestor
  return(c(fdr, power.mca))
}

.choose_r <- function(alpha, c) {
  error <- alpha * c
  R <- 0
  foundR <- FALSE
  while(!foundR) {
    R <- R + 1
    brange <- qnbinom(c(0.1586553, 0.8413447), R, alpha)
    pvalRange <- R / (R + brange)
    diff <- max(abs(pvalRange - alpha))
    if(diff < error) {
      foundR <- TRUE
    }
  }
  return(R)
}

.SABHA_method = function(pvals, qhat, alpha, tau){
  pvals[pvals>tau] = Inf
  khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
  qhat*pvals*length(pvals)/khat
  # which(qhat*pvals<=alpha*khat/length(pvals))
}
.Solve_q_TV = function(Pvals, tau, eps, edges, TV_bd, ADMM_params){
  # edges is a e-by-2 matrix giving the edges of the adjacency graph
  # edges[i,1:2] gives the indices of the nodes on the i-th edge
  # constraint: sum_{i=1,..,e} |q[edges[i,1]] - q[edges[i,2]]| <= TV_bd
  L1_proj = .create_L1_function(TV_bd)
  nedge = dim(edges)[1]; n = length(Pvals)
  M = matrix(0,nedge,n); for(i in 1:nedge){M[i,edges[i,1]]=1; M[i,edges[i,2]]=-1}
  q = .Solve_q_ADMM(Pvals, tau, eps, M, L1_proj, ADMM_params)
  q
}  
.Solve_q_ADMM = function(Pvals, tau, eps, M, projection, ADMM_params){
  # min -sum_i (B[i]*log((1-tau) q[i]) + (1-B[i])*log(1-(1-tau) q[i]))
  # subject to (1) q \in Qset (characterized by M*q \in Mset)
  # and (2) sum_i B[i]/q[i] <= gamma and (3) eps<=q<=1
  # introduce auxiliary variables x, y under the constraint Mq = x, q = y
  # ADMM optimization:
  # minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
  # where qt is the previous iteration's q value
  
  # ADMM_params are: alpha, beta, eta, max_iters, converge_thr
  alpha_ADMM = ADMM_params[1]
  beta = ADMM_params[2]
  eta = ADMM_params[3]
  max_iters = ADMM_params[4]
  converge_thr = ADMM_params[5]
  
  n = length(Pvals)
  B = (Pvals > tau) 
  gamma = n*(1-tau) # bound on sum_i (Pvals[i]>tau) / q[i]*(1-tau)
  q = y = rep(1,n)
  v = rep(0,n)
  u = x = rep(0,dim(M)[1])
  
  .converge_check = function(q,x,y,u,v,q_old,x_old,y_old,u_old,v_old){
    max(c(sqrt(sum((q-q_old)^2))/sqrt(1+sum(q_old^2)),
          sqrt(sum((x-x_old)^2))/sqrt(1+sum(x_old^2)),
          sqrt(sum((y-y_old)^2))/sqrt(1+sum(y_old^2)),
          sqrt(sum((u-u_old)^2))/sqrt(1+sum(u_old^2)),
          sqrt(sum((v-v_old)^2))/sqrt(1+sum(v_old^2))))
  }
  
  stop = FALSE
  iter = 0
  while(!stop){
    iter = iter+1
    q_old = q; x_old = x; y_old = y; u_old = u; v_old = v
    q = .q_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    x = .x_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta, projection)
    y = .y_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    u = .u_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    v = .v_update(B, M, tau,eps,q,x,y,u,v,alpha_ADMM,gamma,beta, eta)
    if(.converge_check(q,x,y,u,v,q_old,x_old,y_old,u_old,v_old)<=converge_thr){stop=TRUE}
    if(iter>=max_iters){stop=TRUE}
  }
  
  return(q)
  
}

# inverse_sum_prox solves: min{1/2 ||x-y||^2 : x_i>0, sum_i 1/x_i <= bound}
# Used in y-update step of ADMM
inverse_sum_prox = function(y,bound){
  
  y = pmax(0,y) # the solution will have all positive x_i's now
  # and we can now ignore the constraint x_i>0
  
  if(sum(1/y)<= bound){
    x=y
  }else{ # use Lagrange multipliers
    
    # we should have - lambda * d/dx_j (sum_i 1/x_i) = d/dx_j (1/2 ||x-y||^2)
    # for all j, for some single lambda>0
    # in other words, lambda / x^2 = x-y (this holds elementwise)
    # rearranging, lambda = x^3 - x^2*y
    # let c = log(lambda) so that it's real-valued
    # we need to solve x^3 - x^2*y - exp(c) = 0 (elementwise)
    
    .cuberoot = function(c){ # this solves the cubic equation x^3-x^2*y-exp(c)=0
      temp1 = ((y/3)^3 + exp(c)/2 + (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
      temp2 = ((y/3)^3 + exp(c)/2 - (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
      x = sign(temp1)*abs(temp1)^(1/3) + sign(temp2)*abs(temp2)^(1/3) + (y/3)
      x
    }
    
    # now we need to choose c, i.e. choose the lagrange multiplier lambda=exp(c)
    # the right value of c is the one that produces an x satisfying sum_i 1/x_i = bound
    
    c = uniroot(function(c){sum(1/.cuberoot(c))-bound},c(-100,100))$root
    x = .cuberoot(c)
  }
  x
}

.q_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
  # minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
  # where qt is the previous iteration's q value
  # equivalently, -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + (alpha eta + beta)/2 * ||q-w||_2^2
  # where w = - (M'(ut + alpha (M qt - xt)) + (vt - beta yt - alpha eta qt))/(alpha eta + beta)
  
  w = - ( t(M)%*%(u + alpha*(M%*%q - x)) + (v - beta*y - alpha*eta*q) )/(alpha*eta + beta)
  
  q[B==1] = (w[which(B==1)]+sqrt(w[which(B==1)]^2+4/(alpha*eta + beta)))/2
  q[B==0] = ((w[which(B==0)]+1/(1-tau))-sqrt((w[which(B==0)]-1/(1-tau))^2+4/(alpha*eta+beta)))/2
  q[q<eps] = eps
  q[q>1] = 1
  q
}

.x_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta, projection){
  # Proj_Mset (M q + u/alpha)
  x = projection(M%*%q + u/alpha) 
}

.y_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
  # Prof_B (q + v/beta)
  # where B = {sum_i B[i]/y[i]<= gamma}
  y = q + v/beta
  y[which(B==1)] = inverse_sum_prox((q+v/beta)[which(B==1)], gamma)
  y
}

.u_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
  u = u + alpha * (M%*%q -x)
  u
}

.v_update = function(B, M, tau,eps,q,x,y,u,v,alpha,gamma,beta, eta){
  v = v + beta * (q-y)
  v
}

.create_L1_function = function(bound){
  function(y){
    # solving: min{1/2 ||x-y||^2_2 : ||x||_1 <= bound}
    if(sum(abs(y))<=bound){x=y} else{
      mu = sort(abs(y), decreasing = TRUE)
      xi = max(which(mu - (cumsum(mu)-bound)/(1:length(mu))>0))
      theta = (sum(mu[1:xi])-bound)/xi
      tmp = abs(y)-theta
      x = rep(0, length(tmp))
      x[which(tmp>0)] = tmp[which(tmp>0)]
      x[which(tmp<=0)] = 0
      x = x*sign(y)
    }
    x
  }
}
