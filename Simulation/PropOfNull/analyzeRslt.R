rm(list = ls())
library(matrixStats)
library(dplyr)
library(kableExtra)
setwd("~/Documents/Project/PhyloMed/Simulation/PropOfNull/eval_asym/Result/")
fileNames <- list.files("./", pattern = "^final.*.Rdata")
sink("EvalPropOfNull_asym_latex.txt")
for (f in 1:length(fileNames)) {
  tab <- matrix(NA, nrow = 11, ncol = 10, 
                dimnames = list(c("pi.0", "pi0.", "pi*", "pi00", "pi01", "pi10", "pi11", 
                                  "pi00.r", "pi01.r", "pi10.r", "Type1/Power"),
                                c("true", "jc.prod", "avg.prod",
                                  "jc.minus.pmin", "avg.minus.pmin",
                                  "avg.minus.pmin.storey",
                                  "jc.minus.pmax", "avg.minus.pmax",
                                  "jc.minus.multi", "avg.minus.multi")))
  tmp <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab["pi0.",] <- colMeans(rslt$alpha1.mat)
  tab["pi.0",] <- colMeans(rslt$alpha2.mat)
  tab["pi*",] <- colMeans(rslt$alphastar.mat)
  tab["pi00",] <- colMeans(rslt$alpha00.ori.mat)
  tab["pi01",] <- colMeans(rslt$alpha01.ori.mat)
  tab["pi10",] <- colMeans(rslt$alpha10.ori.mat)
  tab["pi11",] <- colMeans(rslt$alpha11.mat)
  tab["pi00.r",] <- colMeans(rslt$alpha00.mat)
  tab["pi01.r",] <- colMeans(rslt$alpha01.mat)
  tab["pi10.r",] <- colMeans(rslt$alpha10.mat)
  tab["Type1/Power",] <- colMeans(rslt$gp.mat < 0.05)
  tab <- round(tab, 3)
  cat(tab %>%
        kbl(caption= paste0("CausalType",tmp[1], "NSample", tmp[2], "A",tmp[3], "B",tmp[4]),
            format= "latex",
            align = "c") %>%
        kable_paper(full_width = F,latex_options = "scale_down"))
  cat("\n")
}
sink()


fileNames <- list.files("./", pattern = "^final.*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 9, 
              dimnames = list(c(),
                              c("CausalType", "A", "B",
                                "Pi00.prod-true", "Pi00.multi-true", 
                                "Pi10.prod-true", "Pi10.multi-true", 
                                "Pi01.prod-true", "Pi01.multi-true")))

for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:3] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))[c(1,3,4)]
  load(fileNames[f])
  tmp <- sweep(rslt$alpha00.mat[,-1], 1, rslt$alpha00.mat[,1], "-")[,c(1,6)]
  tab[f, 4:5] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
  tmp <- sweep(rslt$alpha10.mat[,-1], 1, rslt$alpha10.mat[,1], "-")[,c(1,6)]
  tab[f, 6:7] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
  tmp <- sweep(rslt$alpha01.mat[,-1], 1, rslt$alpha01.mat[,1], "-")[,c(1,6)]
  tab[f, 8:9] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
}
tab <- as.data.frame(tab)
tab <- tab[order(tab$A, tab$B, tab$CausalType),]
cat(tab %>%
      kbl(caption= "",
          format= "latex",
          align = "c") %>%
      kable_paper(full_width = F,latex_options = "scale_down"))

library(patchwork)
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

load("allType1Nsample200A0.1B0.5.Rdata")
rawp.prod <- rslt$rawp.prod.mat
rawp.minus.jcmin <- rslt$rawp.minus.jcmin.mat
rawp.minus.jcmax <- rslt$rawp.minus.jcmax.mat
rawp.minus.slim <- rslt$rawp.minus.slim.mat
rawp.minus.avg <- rslt$rawp.minus.avg.mat
globalp.mat = matrix(NA, nrow = nrow(rslt$rawp.prod.mat), 5)
for (i in 1:nrow(globalp.mat)) {
  rawp.prod.rm = na.omit(rawp.prod[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.jcmin.rm = na.omit(rawp.minus.jcmin[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.jcmax.rm = na.omit(rawp.minus.jcmax[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.slim.rm = na.omit(rawp.minus.slim[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.avg.rm = na.omit(rawp.minus.avg[i,which(rslt$causalNode.mat[i,] == 0)])
  L = length(rawp.prod.rm)
  globalp.mat[i,] = c(p.hmp(rawp.prod.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.jcmin.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.jcmax.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.slim.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.avg.rm, w = rep(1/L, L), L = L))
  }
my.pvalue.list <-list("Prod"=globalp.mat[,1],
                      "Minus.jcmin"=globalp.mat[,2],
                      "Minus.jcmax"=globalp.mat[,3],
                      "Minus.slim"=globalp.mat[,4],
                      "Minus.avg"=globalp.mat[,5])
postscript("../Figs/qq_nullTaxaCausalType1_globalp.eps")  
.qqunifPlot(my.pvalue.list, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("blue","red","green","gray","orange"))))
dev.off()

load("allType2Nsample200A0.1B0.5.Rdata")
rawp.prod <- rslt$rawp.prod.mat
rawp.minus.jcmin <- rslt$rawp.minus.jcmin.mat
rawp.minus.jcmax <- rslt$rawp.minus.jcmax.mat
rawp.minus.slim <- rslt$rawp.minus.slim.mat
rawp.minus.avg <- rslt$rawp.minus.avg.mat
globalp.mat = matrix(NA, nrow = nrow(rslt$rawp.prod.mat), 5)
for (i in 1:nrow(globalp.mat)) {
  rawp.prod.rm = na.omit(rawp.prod[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.jcmin.rm = na.omit(rawp.minus.jcmin[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.jcmax.rm = na.omit(rawp.minus.jcmax[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.slim.rm = na.omit(rawp.minus.slim[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.avg.rm = na.omit(rawp.minus.avg[i,which(rslt$causalNode.mat[i,] == 0)])
  L = length(rawp.prod.rm)
  globalp.mat[i,] = c(p.hmp(rawp.prod.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.jcmin.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.jcmax.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.slim.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.avg.rm, w = rep(1/L, L), L = L))
}
my.pvalue.list <-list("Prod"=globalp.mat[,1],
                      "Minus.jcmin"=globalp.mat[,2],
                      "Minus.jcmax"=globalp.mat[,3],
                      "Minus.slim"=globalp.mat[,4],
                      "Minus.avg"=globalp.mat[,5])
postscript("../Figs/qq_nullTaxaCausalType2_globalp.eps")  
.qqunifPlot(my.pvalue.list, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("blue","red","green","gray","orange"))))
dev.off()

load("Result_rawp/rawp.allType1Nsample200A0.1B0.Rdata")
rawp.prod <- rslt$rawp.prod.mat
rawp.minus.slim <- rslt$rawp.minus.slim.mat
rawp.minus.avg <- rslt$rawp.minus.avg.mat
globalp.mat = matrix(NA, nrow = nrow(rslt$rawp.prod.mat), 3)
for (i in 1:nrow(globalp.mat)) {
  rawp.prod.rm = na.omit(rawp.prod[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.slim.rm = na.omit(rawp.minus.slim[i,which(rslt$causalNode.mat[i,] == 0)])
  rawp.minus.avg.rm = na.omit(rawp.minus.avg[i,which(rslt$causalNode.mat[i,] == 0)])
  L = length(rawp.prod.rm)
  globalp.mat[i,] = c(p.hmp(rawp.prod.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.slim.rm, w = rep(1/L, L), L = L),
                      p.hmp(rawp.minus.avg.rm, w = rep(1/L, L), L = L))
}
my.pvalue.list <-list("Prod"=globalp.mat[,1],
                      "Minus.slim"=globalp.mat[,2],
                      "Minus.avg"=globalp.mat[,3])
postscript("Figs/qq_nullTaxaCausalType1A0.1B0_globalp.eps")
.qqunifPlot(my.pvalue.list, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("blue","red","green"))))
dev.off()
