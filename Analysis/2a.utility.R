library(patchwork)
library(lattice)
.plotTreeBoxScatter <- function(input_data, rslt_data){
  Trt = input_data$treatment
  n.sample = length(Trt)
  conf = matrix(1, nrow = n.sample, ncol = 1)
  
  M = input_data$mediators
  M2 = M + 0.5
  
  Y = input_data$outcome
  
  tree = input_data$tree
  K = .ntaxa(tree)
  treestructure = .phylostructure(tree)
  
  sig.nodeID = as.numeric(rslt_data$sig.node["perm.jsmix"])
  pval.alpha = rslt_data$pval.alphabeta["perm.alpha", sig.nodeID]
  pval.beta = rslt_data$pval.alphabeta["perm.beta", sig.nodeID]
  
  tree.vis = tree
  rawp = rslt_data$node.pval["perm.jsmix",]
  tree.vis$node.label = rawp
  p1 <- ggtree(tree.vis, branch.length = "none") + 
    geom_point2(aes(subset=!isTip), shape=21, size=-log10(as.numeric(rawp))*3, fill = "red") +
    geom_hilight(node = sig.nodeID+K, fill = "steelblue", alpha = 0.5) + 
    theme_tree(plot.margin=margin(5,5,5,5)) +
    scale_x_reverse() + coord_flip()

  child.left = treestructure$phylochildren[sig.nodeID+K,1]
  child.right = treestructure$phylochildren[sig.nodeID+K,2]
  
  Mc2.left = rowSums(M2[,treestructure$descendant[child.left,], drop = FALSE])
  Mc2.right = rowSums(M2[,treestructure$descendant[child.right,], drop = FALSE])
  if(mean(Mc2.left) < mean(Mc2.right)){
    Mc2 = cbind(Mc2.right, Mc2.left);
  }else{
    Mc2 = cbind(Mc2.left, Mc2.right);
  }
  
  G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
  mod = summary(lm(Y~cbind(conf[,-1], Trt)))
  mod.resid = mod$residuals
  
  df_tmp = data.frame(Trt = factor(Trt, levels = c(0,1), labels = c("Control", "Antibiotics")), 
                      G = G, outcome.Residual = mod.resid)
  
  p2 <- ggplot(df_tmp, aes(x = Trt, y = G, shape = Trt, color = Trt)) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(0.1)) + 
    scale_color_brewer(palette="Set2") + 
    scale_shape_manual(values=c(17,19)) + 
    annotate(geom = 'text', x = -Inf, y = max(G), label = 'P[alpha]', parse = T,
             hjust = -0.2, size = 8, col = "#f03b20") +
    annotate(geom = 'text', x = -Inf, y = max(G), label = paste(" =", format(pval.alpha, digits = 2, scientific = T)),
             hjust = -0.3, size = 8, col = "#f03b20") +
    labs(title = "", x = "Treatment", y = "Subcomposition in log-ratio") +
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

  p3 <- ggplot(df_tmp, aes(x = G, y = mod.resid, shape = Trt, color = Trt)) +
    geom_point() + 
    scale_color_brewer(palette="Set2") + 
    scale_shape_manual(values=c(17,19)) + 
    annotate(geom = 'text', x = -Inf, y = max(mod.resid), label = 'P[beta]', parse = T,
             hjust = -0.2, size = 8, col = "#f03b20") +
    annotate(geom = 'text', x = -Inf, y = max(mod.resid), label = paste(" =", format(pval.beta, digits = 2, scientific = T)),
             hjust = -0.3, size = 8, col = "#f03b20") +
    geom_smooth(method = lm, se = FALSE, aes(group=1), col = "black") + 
    labs(x = "Subcomposition in log-ratio", y = "Residual of body fat ~ treatment") +
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

  p <- p1/(p2|p3)
  p
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
  pvalues = lapply(pvalues, .filterpval)
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
