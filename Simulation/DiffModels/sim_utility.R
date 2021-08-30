#### INTERNAL FUNCTIONS ####

# .ntaxa
.ntaxa <- function(tree){
  length(tree$tip.label)
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

library(nnet)

# Test alpha in quasi-likelihood model
.test_alpha_qs <- function(Mc, Trt, conf, est.reduce.beta, test.type="mv"){
  
  X = cbind(conf, Trt)
  p = ncol(X)
  X.par.index = p
  
  Y = Mc
  nY = rowSums(Y)
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  X.reduce = X[,-X.par.index, drop=FALSE]
  par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
  
  n.par.interest.beta = length(par.interest.index.beta) 
  
  data.beta = list(Y=Y, X=X)
  Score.reduce.beta = .fun.score.i.beta(est.reduce.beta, data.beta)
  
  # for resampling: S.beta.list, I.beta.list
  # S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
  tmp = .fun.hessian.beta(est.reduce.beta, data.beta, save.list=TRUE)
  # I.beta.list = tmp$I.beta.list  
  
  Hess.reduce.beta = tmp$Hessian.beta
  
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                            cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  ## estimate parameters
  input.data = list(Y=cbind(Y[,m], Y[,-m]), X = X[,-1])
  multinom.out = multinom(Y ~ X , input.data, maxit=10000, abstol = 1.0e-10, reltol = 1.0e-15, trace=FALSE)
  est = coef(multinom.out)[,X.par.index]
  
  if(test.type=="mv"){
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    B2 =  matrix(0, n.beta, n.beta)
    
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
    }
    
    B = B1 %*% B2 %*% t(B1)
    
    stat = as.numeric(A %*% ginv(B) %*% A)
    pval = 1 - pchisq(stat, n.par.interest.beta) 
  }
  
  return(list(stat=stat, est=est, pval=pval))
}

.Ei.beta <- function(m, p, beta, X.i, Y.i){
  Ei.out = rep(NA,m)
  for(j in 1:(m-1)){
    Ei.out[j] = exp(beta[((j-1)*p+1):(j*p)] %*% X.i)
  }
  Ei.out[m] = 1
  return (Ei.out)
}

.fun.score.i.beta <- function(beta, data){
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    warning("Dim of initial beta does not match the dim of covariates")
  }else{
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      P.i[is.na(P.i)] = 1
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )
    }
    return (Score.beta.i)
  }
}

.fun.hessian.beta <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; X = data$X
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      P.i[is.na(P.i)] = 1
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) ) 
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
      }
    }
    
    if(save.list){
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
    }else{
      return (Hessian.beta)
    }
  }
}
