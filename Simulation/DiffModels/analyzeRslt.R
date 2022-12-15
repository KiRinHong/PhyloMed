rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Simulation/DiffModels/Result/")
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 8, 
              dimnames = list(c(),
                              c("CausalType", "N", "A", "B", "LR", "DM", "QSC", "QSG")))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:4] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab[f, 5:8] <- round(colMeans(rslt$gp.mat < 0.05), 3)
}
tab <- as.data.frame(tab)

library(lattice)
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
.filterpval <- function(x){
  return(na.omit(x[x>0 & x<=1]))
}

load(fileNames[1])
range(colMeans(rslt$rawp.lr.mat < 0.05, na.rm = T))
range(colMeans(rslt$rawp.dm.mat < 0.05, na.rm = T))

# N=200 S_alpha=0 S_beta=0
rawp.list <- list("LR"=c(rslt$rawp.lr.mat),"DM"=c(rslt$rawp.dm.mat),"QS"=c(rslt$rawp.qs.mat))
.qqunifPlot(rawp.list, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1, cex.title = 1.5, col = "red")))
