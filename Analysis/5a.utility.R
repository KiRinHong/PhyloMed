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
.inverse_sum_prox = function(y,bound){
  
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
  y[which(B==1)] = .inverse_sum_prox((q+v/beta)[which(B==1)], gamma)
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