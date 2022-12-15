library(patchwork)
library(lattice)

.plotTree <- function(input_data, rslt_data, layout = "rectangular"){
  # input_data = combo.filter; rslt_data = rslt.combo.filter20; layout = "circular"
  tree = input_data$tree
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  
  tree.vis = tree
  rawp = rslt_data$node.pval["perm.jsmix.prod",] # p.adjust(rslt_data$node.pval["perm.jsmix",], method = "BH")
  tree.vis$node.label = rawp
  
  #sig.nodeID.HMP = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.HMP"], ",")))
  sig.nodeID.BH = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.prod.BH"], ",")))
  min.id = which.min(rslt_data$node.pval["perm.jsmix.prod",sig.nodeID.BH])
  #mrca.id = getMRCA(tree, which(colSums(treestructure$descendant[sig.nodeID.BH+K,,drop=FALSE]) == 1)) - K
  
  # if(mrca.id %in% sig.nodeID.HMP){
  #   p <- ggtree(tree.vis, layout = layout, branch.length = "none") +
  #     geom_point2(aes(subset = !isTip), shape = 21, size= -log10(as.numeric(rawp))*3, fill = "red") +
  #     geom_hilight(node = mrca.id+K, fill = "darkgreen", alpha = .4)
  # }else{
  #   p <- ggtree(tree.vis, layout = layout, branch.length = "none") +
  #     geom_point2(aes(subset = !isTip), shape = 21, size = -log10(as.numeric(rawp))*3, fill = "red") +
  #     geom_hilight(node = sig.nodeID.HMP[length(sig.nodeID.HMP)]+K, fill = "darkgreen", alpha = .4)
  # }

  # lables = rep(NA, 2*K-1); lables[sig.nodeID.BH+K] = LETTERS[1:length(sig.nodeID.BH)]
  # p <- p + geom_text2(aes(subset=!isTip, label=lables), vjust=-.5, hjust=-.5, angle = 0, size=5)
  # for (i in 1:length(sig.nodeID.BH)) {
  #   p <- p + geom_hilight(node = sig.nodeID.BH[i]+K, fill = "steelblue", alpha = .4)
  # }
  if(layout == "rectangular"){
    p <- ggtree(tree.vis, layout = layout, branch.length = "none") +
      geom_point2(aes(subset = !isTip), shape = 21, size= -log10(as.numeric(rawp))*4, fill = "red") +
      geom_hilight(node = sig.nodeID.BH[min.id]+K, fill = "steelblue", alpha = .6) +
      theme_tree(plot.margin=margin(5,5,5,5)) +
      scale_x_reverse() + coord_flip()
    p
  }else if(layout == "circular"){
    p <- ggtree(groupClade(tree.vis, sig.nodeID.BH[min.id]+K), aes(color = group), layout = layout, branch.length = "none") + 
      scale_color_manual(values=c("black", "steelblue")) + 
      geom_point2(aes(subset = !isTip), shape = 21, size= -log10(as.numeric(rawp))*2, fill = "red") +
      theme_tree(plot.margin=margin(5,5,5,5), legend.position = "none")
    p <- rotate_tree(scaleClade(p, node = sig.nodeID.BH[min.id]+K, scale = 3), angle = 180)
    p <- p + geom_hilight(node = sig.nodeID.BH[min.id]+K, fill = "steelblue", alpha = .6)
    p
  }

}
.plotSubTree <- function(input_data, rslt_data){
  # input_data = cecal.top; rslt_data = rslt.cecal.top100
  tree = input_data$tree
  sub_tree = subtrees(tree)
  sig.nodeID.BH = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.prod.BH"], ",")))
  min.id = which.min(rslt_data$node.pval["perm.jsmix.prod",sig.nodeID.BH])
  sub_tree_node = sub_tree[[sig.nodeID.BH[min.id]]]
  
  p <- ggtree(sub_tree_node, layout = "rectangular", branch.length = "none") + #  geom_rootpoint(color="red", size=5) +
    geom_tippoint(shape=23, size=5, fill="black") +
    theme_tree(plot.margin=margin(5,5,5,5)) +
    scale_x_reverse() + coord_flip()
  p
}
.plotScatterPie <- function(input_data, rslt_data, treatment, covariate, outcome){
  Trt = input_data$meta[[treatment]]
  M = input_data$mediators
  Y = input_data$meta[[outcome]]
  id.na = which(is.na(Y))
  if(length(id.na) > 0){
    Trt = Trt[-id.na]
    M = M[-id.na,]
    Y = Y[-id.na]
  }
  n.sample = length(Y)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      conf = cbind(1, input_data$meta[[covariate]])
    }else{
      conf = cbind(1, as.matrix(sample_data(input_data$meta[,covariate])))
    }
    if(length(id.na) > 0){
      conf = conf[-id.na,]
    }
  }else{
    conf = matrix(1, nrow = n.sample, ncol = 1)
  }
  sig.nodeID.BH = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.prod.BH"], ",")))
  tree = input_data$tree
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  for (i in 1:length(sig.nodeID.BH)) {
    sig.nodeID = sig.nodeID.BH[i]
    child.left = treestructure$phylochildren[sig.nodeID+K,1]
    child.right = treestructure$phylochildren[sig.nodeID+K,2]
    Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
    Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
    if(mean(Mc.left) < mean(Mc.right)){
      Mc = cbind(Mc.right, Mc.left)
    }else{
      Mc = cbind(Mc.left, Mc.right)
    }
    idx = which(rowSums(Mc) != 0)
    mod = summary(lm(Y[idx]~0+conf[idx,]+Trt[idx]))
    mod.resid = mod$residuals
    df_tmp = data.frame(Trt = scale(Trt[idx]), Y = scale(mod.resid), Left = Mc[idx,1], Right = Mc[idx,2])
    
    # pval = rslt_data$node.pval["perm.jsmix", sig.nodeID]
    p = ggplot() + 
      geom_scatterpie(data = df_tmp, aes(x=Trt, y=Y, r=.08), cols=c("Left", "Right")) +
      # annotate(geom = 'text', x = -Inf, y = max(scale(mod.resid)), label = labels[i], parse = F,
      #          hjust = -0.2, size = 8, col = "#f03b20") +
      # annotate(geom = 'text', x = -Inf, y = max(scale(Y)), label = paste(" =", format(pval, digits = 2, scientific = T)),
      #          hjust = -0.3, size = 8, col = "#f03b20") +
      labs(x = "Fat intake", y = "BMI") + 
      theme_bw(base_size = 14) + 
      theme(axis.title = element_text(face = "bold",size = rel(2)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(2)),
            axis.text.x = element_text(size = rel(2)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.position="none",
            plot.margin = unit(c(5, 10, 5, 10), "mm"))
    print(p)
  }
}

.plotScatterScatter <- function(input_data, rslt_data, treatment, covariate, outcome){
  #input_data = combo.top; rslt_data = combo.top.calor.adj2; treatment = "calor"; covariates = c("tfat", "aofib")
  Trt = input_data$meta[[treatment]]
  M = input_data$mediators
  Y = input_data$meta[[outcome]]
  id.na = which(is.na(Y))
  if(length(id.na) > 0){
    Trt = Trt[-id.na]
    M = M[-id.na,]
    Y = Y[-id.na]
  }
  n.sample = length(Y)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      conf = cbind(1, input_data$meta[[covariate]])
    }else{
      conf = cbind(1, as.matrix(sample_data(input_data$meta[,covariate])))
    }
    if(length(id.na) > 0){
      conf = conf[-id.na,]
    }
  }else{
    conf = matrix(1, nrow = n.sample, ncol = 1)
  }
  sig.nodeID.BH = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.prod.BH"], ",")))
  tree = input_data$tree
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  labels = paste0("Node ", LETTERS[1:length(sig.nodeID.BH)])
  for (i in 1:length(sig.nodeID.BH)) {
    sig.nodeID = sig.nodeID.BH[i]
    pval.alpha = rslt_data$pval.alphabeta["perm.alpha", sig.nodeID]
    pval.beta = rslt_data$pval.alphabeta["perm.beta", sig.nodeID]
    child.left = treestructure$phylochildren[sig.nodeID+K,1]
    child.right = treestructure$phylochildren[sig.nodeID+K,2]
    Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
    Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
    if(mean(Mc.left) < mean(Mc.right)){
      Mc = cbind(Mc.right, Mc.left)
    }else{
      Mc = cbind(Mc.left, Mc.right)
    }
    idx = which(rowSums(Mc) != 0)
    Mc2 = Mc + 0.5
    G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
    mod1 = summary(glm(G[idx]~0+conf[idx,],family = quasi()))
    mod1.resid = mod1$deviance.resid
    mod2 = summary(lm(Y[idx]~0+conf[idx,]+Trt[idx]))
    mod2.resid = mod2$residuals
    df_tmp = data.frame(Trt = Trt[idx], G = G[idx], G.res = mod1.resid, outcome.res = mod2.resid)
    p1 <- ggplot(df_tmp, aes(x = Trt, y = G.res)) +
      geom_point() + 
      annotate(geom = 'text', x = -Inf, y = max(mod1.resid), label = 'pval', parse = T,
               hjust = -0.2, size = 8, col = "#f03b20") +
      annotate(geom = 'text', x = -Inf, y = max(mod1.resid), label = paste(" =", format(pval.alpha, digits = 2, scientific = T)),
               hjust = -0.3, size = 8, col = "#f03b20") +
      geom_smooth(method = lm, se = FALSE, aes(group=1), col = "red") + 
      labs(x = "Fat intake", y = "Residual of adjusted Subcompostion in log-ratio") +
      theme_bw(base_size=14) +
      theme(axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.position="none",
            plot.margin = unit(c(5, 10, 5, 10), "mm"))
    p2 <- ggplot(df_tmp, aes(x = G, y = outcome.res)) +
      geom_point() + 
      annotate(geom = 'text', x = -Inf, y = max(mod2.resid), label = 'pval', parse = T,
               hjust = -0.2, size = 8, col = "#f03b20") +
      annotate(geom = 'text', x = -Inf, y = max(mod2.resid), label = paste(" =", format(pval.beta, digits = 2, scientific = T)),
               hjust = -0.3, size = 8, col = "#f03b20") +
      geom_smooth(method = lm, se = FALSE, aes(group=1), col = "red") + 
      labs(x = "Subcomposition in log-ratio", y = "Residual of adjusted BMI") +
      theme_bw(base_size=14) +
      theme(axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.position="none",
            plot.margin = unit(c(5, 10, 5, 10), "mm"))
    p <- p1|p2
    print(p)
  }
}

.plotBoxScatter <- function(input_data, rslt_data, treatment, covariate, outcome){
  Trt = input_data$meta[[treatment]]
  M = input_data$mediators
  Y = input_data$meta[[outcome]]
  id.na = which(is.na(Y))
  if(length(id.na) > 0){
    Trt = Trt[-id.na]
    M = M[-id.na,]
    Y = Y[-id.na]
  }
  n.sample = length(Y)
  if(!is.null(covariate)){
    if(length(covariate) == 1){
      conf = cbind(1, input_data$meta[[covariate]])
    }else{
      conf = cbind(1, as.matrix(sample_data(input_data$meta[,covariate])))
    }
    if(length(id.na) > 0){
      conf = conf[-id.na,]
    }
  }else{
    conf = matrix(1, nrow = n.sample, ncol = 1)
  }
  sig.nodeID.BH = as.numeric(unlist(strsplit(rslt_data$sig.node["perm.jsmix.prod.BH"], ",")))
  min.id = which.min(rslt_data$node.pval["perm.jsmix.prod",sig.nodeID.BH])
  
  tree = input_data$tree
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  sig.nodeID = sig.nodeID.BH[min.id]
  pval.alpha = rslt_data$pval.alphabeta["perm.alpha", sig.nodeID]
  pval.beta = rslt_data$pval.alphabeta["perm.beta", sig.nodeID]
  child.left = treestructure$phylochildren[sig.nodeID+K,1]
  child.right = treestructure$phylochildren[sig.nodeID+K,2]
  Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
  Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
  if(mean(Mc.left) < mean(Mc.right)){
    Mc = cbind(Mc.right, Mc.left)
  }else{
    Mc = cbind(Mc.left, Mc.right)
  }
  idx = which(rowSums(Mc) != 0)
  Mc2 = Mc + 0.5
  G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
  mod1 = summary(glm(G[idx]~0+conf[idx,],family = quasi()))
  mod1.resid = mod1$deviance.resid
  mod2 = summary(lm(Y[idx]~0+conf[idx,]+Trt[idx]))
  mod2.resid = mod2$residuals
  df_tmp = data.frame(Trt = Trt[idx], G = G[idx], G.res = mod1.resid, outcome.res = mod2.resid)
  df_tmp$Trt = factor(df_tmp$Trt, levels = c(0,1), labels = c("Control", "Antibiotics"))
  p1 <- ggplot(df_tmp, aes(x = Trt, y = G.res, shape = Trt, color = Trt)) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(0.1), size = 4) + 
    scale_color_brewer(palette="Set2") + 
    scale_shape_manual(values=c(17,19)) + 
    annotate(geom = 'text', x = -Inf, y = max(mod1.resid), label = 'pval', parse = F,
             hjust = -0.2, size = 10) +
    annotate(geom = 'text', x = -Inf, y = max(mod1.resid), label = paste(" =", format(pval.alpha, digits = 2, scientific = T)),
             hjust = -0.5, size = 10) +
    labs(x = "Treatment", y = "Subcompostion in log-ratio") +
    theme_bw(base_size=14) +
    theme(axis.title = element_text(face = "bold",size = rel(2)),
          axis.title.y = element_text(angle = 90, vjust = 2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(2)),
          axis.text.x = element_text(size = rel(2)),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          legend.position="none",
          plot.margin = unit(c(5, 10, 5, 10), "mm"))
  p2 <- ggplot(df_tmp, aes(x = G, y = outcome.res, shape = Trt, color = Trt)) +
    geom_point(size = 4) + 
    scale_color_brewer(palette="Set2") + 
    scale_shape_manual(values=c(17,19)) + 
    annotate(geom = 'text', x = -Inf, y = max(mod2.resid), label = 'pval', parse = F,
             hjust = -0.2, size = 10) +
    annotate(geom = 'text', x = -Inf, y = max(mod2.resid), label = paste(" =", format(pval.beta, digits = 2, scientific = T)),
             hjust = -0.5, size = 10) +
    geom_smooth(method = lm, se = FALSE, aes(group=1), col = "black") + 
    labs(x = "Subcomposition in log-ratio", y = "Body fat (%)") +
    theme_bw(base_size=14) +
    theme(axis.title = element_text(face = "bold",size = rel(2)),
          axis.title.y = element_text(angle = 90,vjust = 2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(2)),
          axis.text.x = element_text(size = rel(2)),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          legend.position="none",
          plot.margin = unit(c(5, 10, 5, 10), "mm"))
  p <- p1|p2
  print(p)
  
}

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

# credit: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
.qqunifPlot<-function(pvalues,
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=FALSE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=19, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  #error checking
  if(is.list(pvalues)){
    pvalues = lapply(pvalues, .filterpval)
  }else{
    pvalues = .filterpval(pvalues)
  }
  
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  print(xyplot(pvalues~exp.x, groups=grp, xlab=list(xlab, cex=1.5), ylab=list(ylab,cex=1.5), aspect=aspect, 
               prepanel=prepanel, scales=list(axs="i", cex=1.5), pch=pch,
               panel = function(x, y, ...) {
                 if (draw.conf) {
                   panel.qqconf(n, conf.points=conf.points, 
                                conf.col=conf.col, conf.alpha=conf.alpha)
                 };
                 panel.xyplot(x,y, ...);
                 panel.abline(0,1);
               }, par.settings=par.settings, ...
  ))
}

theme_Publication <- function(base_size=14, base_family="") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5), hjust = 0.5),
            text = element_text(size = 20),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA), # key background
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1.25, "cm"),
            legend.key.height = NULL,                # key height (unit)
            legend.key.width = NULL,                 # key width (unit)
            legend.text = element_text(size=rel(1.5)),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", size = rel(1.5))
    ))
}

.filterpval <- function(x){
  return(na.omit(x[x>0 & x<=1]))
}

.calEFDR.BH <- function(pval, causalNode, fdr.alpha = 0.05){
  pval.adj = p.adjust(pval, method = "BH")
  tp = length(which(pval.adj < fdr.alpha & causalNode == T))
  fp = length(which(pval.adj < fdr.alpha & causalNode == F))
  if(length(which(pval.adj < fdr.alpha)) == 0){
    fdr = NA
  }else{
    fdr = fp / (tp + fp)
  }
  return(fdr)
}
