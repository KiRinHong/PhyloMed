rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")

source("2a.utility.R")

library(ggplot2)
library(ggtree)
library(ggpubr)
library(scatterpie)
library(ape)
library(matrixStats)
library(phyloseq)
load("../Data/Deriveddata/Cecal.filter20top100.rda")
load("../Data/Deriveddata/rslt.runModel.rda")
#### Cecal: Trt ~ pFat ####
dat_tmp = data.frame(Trt = factor(cecal.top$meta[["Treatment"]], levels = c(0,1), labels = c("Control", "Antibiotics")), 
                     outcome = cecal.top$meta[["pFat"]])

p <- ggplot(dat_tmp, aes(x = Trt, y = outcome, shape = Trt, color = Trt)) +
  # facet_grid(cols = vars(type)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1), size = 4) + 
  stat_compare_means(label.x = 1.3, label.y = 30, size = 10) + 
  scale_color_brewer(palette="Set2") + 
  scale_shape_manual(values=c(17,19)) + 
  theme_Publication() + 
  theme(legend.position="none") +
  labs(title = "", x = "Treatment", y = "Body fat (%)")
p
ggsave("../Figs/Cecal.boxplot.antibiotics_pFat.pdf", plot = p, width = 15, height = 15)

pdf("../Figs/Cecal.tree.pdf", width = 20, height = 10)
.plotTree(cecal.top, rslt.cecal.top100, "rectangular")
dev.off()
pdf("../Figs/Cecal.subtree.pdf", width = 15, height = 15)
.plotSubTree(cecal.top, rslt.cecal.top100)
dev.off()

# rslt_data$node.pval["perm.jsmix",c(45,75)]
# p.adjust(rslt_data$node.pval["perm.jsmix",], method = "BH")[c(45,75)]
# p <- rslt_data$node.pval["perm.jsmix",]
# n <- length(p)
# nm <- names(p)
# p0 <- setNames(p, nm)
# nna <- !is.na(p)
# p <- p[nna]
# lp <- length(p)
# i <- lp:1L
# o <- order(p, decreasing = TRUE)
# ro <- order(o)
# p0[nna] <- pmin(1, cummin(98/i * p[o]))[ro]
# p0[c(45,75)]
# pmin(1, 98/i * p[o])[ro]
# https://www.science.org/doi/full/10.1126/sciadv.abd6989

pdf("../Figs/Cecal.boxscatter.pdf", width = 20, height = 10)
.plotBoxScatter(cecal.top, rslt.cecal.top100, "Treatment", NULL, "pFat")
dev.off()

rawp.list.cecal <- list("PhyloMed"=rslt.cecal.top100$node.pval["perm.jsmix",],
                        "JS"=rslt.cecal.top100$node.pval["perm.js",],
                        "Sobel"=rslt.cecal.top100$node.pval["sobel",])
pdf("../Figs/Cecal.qqplot.pdf", width = 10, height = 10)
.qqunifPlot(rawp.list.cecal, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
dev.off()

#### COMBO: fat ~ bmi, adjusted by calor ####
load("../Data/Deriveddata/COMBO.filter20.rda")
dat_tmp = data.frame(Trt = combo.filter$meta[["fat"]], 
                     outcome = combo.filter$meta[["bmi"]])

p <- ggplot(dat_tmp, aes(x = Trt, y = outcome)) +
  geom_point(size = 4)+
  stat_smooth(method = "lm",
              col = "red",
              se = FALSE,
              size = 1) +
  stat_cor(method = "pearson", label.x = -2, label.y = 40, size = 10, ) + 
  theme_Publication() + 
  theme(legend.position="none") +
  labs(title = "", x = "Fat intake", y = "BMI")
p
ggsave("../Figs/COMBO.scatterplot.fat_bmi.pdf", plot = p, width = 15, height = 15)

pdf("../Figs/COMBO.tree.pdf", width = 15, height = 15)
.plotTree(combo.filter, rslt.combo.filter20, "circular")
dev.off()
pdf("../Figs/COMBO.scatterpie.pdf", width = 15, height = 15)
.plotScatterPie(combo.filter, rslt.combo.filter20, "fat", "calor", "bmi")
dev.off()

rawp.list.combo <- list("PhyloMed"=rslt.combo.filter20$node.pval["perm.jsmix",],
                        "JS"=rslt.combo.filter20$node.pval["perm.js",],
                        "Sobel"=rslt.combo.filter20$node.pval["sobel",])

pdf("../Figs/COMBO.qqplot.pdf", width = 10, height = 10)
.qqunifPlot(rawp.list.combo, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
dev.off()

##### Fig2: qqplot continuous outcome #####
# qq-plot, continuous outcome, n.sample = 200
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)),
                              "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)),
                              "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)),
                              "CMM"=gp[,"ccmm.norm.tide"])
# qq-plot, continuous outcome, n.sample = 50
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))

setEPS()# Set postscript arguments
postscript("../Figs/qq_lscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/qq_lscon10.eps")  
.qqunifPlot(my.pvalue.list.con10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/qq_lscon01.eps")  
.qqunifPlot(my.pvalue.list.con01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/qq_sscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_sscon10.eps")  
.qqunifPlot(my.pvalue.list.con10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_sscon01.eps")  
.qqunifPlot(my.pvalue.list.con01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()

##### FigS1: qqplot binary outcome #####
# qq-plot, binary outcome, n.sample = 200
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample200A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample200A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
# qq-plot, binary outcome, n.sample = 50
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample50A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1)))
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample50A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*4,1), na.rm = T)) # 26 replicates produce NA

setEPS()# Set postscript arguments
postscript("../Figs/qq_lsbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_lsbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_lsbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_ssbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_ssbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/qq_ssbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()


#### power bar plot, continuous outcome and binary outcome, large effect size ####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
power.con.tab <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
power.con.tab <- power.con.tab[which(power.con.tab$A == 1 & power.con.tab$B == 1), c("CausalType", "N", method.sel)]
rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.bin.tab <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
power.bin.tab <- power.bin.tab[order(power.bin.tab$CausalType, power.bin.tab$N),]
power.bin.tab <- power.bin.tab[which(power.bin.tab$A == 1 & power.bin.tab$B == 1), c("CausalType", "N", method.sel)]
rownames(power.bin.tab) <-  paste0("bin_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.tab <- rbind(power.con.tab, power.bin.tab)
power.tab <- power.tab[,method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
power.tab$outcome_type <- factor(rep(c("con", "bin"), each = 4), levels = c("con", "bin"), labels = c("Continuous outcome", "Binary outcome"))
power.tab$sample_size <- factor(rep(c(50, 200), 4), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(rep(c(3, 6), each = 2), 2), levels = c(3, 6), labels = c("3", "6"))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("outcome_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating taxa",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
  facet_grid(sample_size ~ outcome_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt")))
p
ggsave("../Figs/power_barplot.pdf", plot = p, width = 20, height = 15)

#### power bar plot, supp, continuous outcome, N=200 and small effect size ####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
power.tab <- read.csv("../Simulation/SmallEffectSize/continuous/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.tab <- power.tab[order(power.tab$CausalType),]
power.tab <- power.tab[which(power.tab$A == 0.001 & power.tab$B == 0.5), method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
power.tab$numberoftaxa <- factor(c(3,6), levels = c(3,6), labels = c("3", "6"))
power.tab.long <- reshape2::melt(power.tab, id.vars = "numberoftaxa")
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating taxa",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt")))
p
ggsave("../Figs/power_barplot_Con200SmlEffSize.pdf", plot = p, width = 20, height = 15)

#### power bar plot, supp, continuous and binary outcome, partially overlapped ####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
power.con.tab <- read.csv("../Simulation/PartialOverlap/continuous/RESULT/PowerAllCont.csv", row.names = 1)
power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
power.con.tab <- power.con.tab[which(power.con.tab$A == 1 & power.con.tab$B == 1), c("CausalType", "N", method.sel)]
power.con.tab[c(1,3),"ccmm.norm.tide"] <- NA
rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.bin.tab <- read.csv("../Simulation/PartialOverlap/binary/RESULT/PowerAllCont.csv", row.names = 1)
power.bin.tab <- power.bin.tab[order(power.bin.tab$CausalType, power.bin.tab$N),]
power.bin.tab <- power.bin.tab[which(power.bin.tab$A == 1 & power.bin.tab$B == 1), c("CausalType", "N", method.sel)]
rownames(power.bin.tab) <-  paste0("bin_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.tab <- rbind(power.con.tab, power.bin.tab)
power.tab <- power.tab[,method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
power.tab$outcome_type <- factor(rep(c("con", "bin"), each = 4), levels = c("con", "bin"), labels = c("Continuous outcome", "Binary outcome"))
power.tab$sample_size <- factor(rep(c(50, 200), 4), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(rep(c(1, 2), each = 2), 2), levels = c(1, 2), labels = c("1", "2"))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("outcome_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating taxa",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
  facet_grid(sample_size ~ outcome_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt")))
p
ggsave("../Figs/power_barplot_PartialOverlap.pdf", plot = p, width = 20, height = 15)


#### power bar plot, supp, continuous outcome, N=50/200, increase number of causal taxa ####
# randomly select and fully overlapped, fix the effect size A=10e-3 B=0.5
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
power.tab.A103 <- read.csv("../Simulation/IncreaseNumOfCausalTaxa/continuousA103B5/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.tab.A103 <- power.tab.A103[order(power.tab.A103$CausalType, power.tab.A103$N), method.sel]
rownames(power.tab.A103) <-  paste0("Seffsize_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2", "ss_type5", "ls_type5"))
power.tab.A1 <- read.csv("../Simulation/IncreaseNumOfCausalTaxa/continuousA1B1/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.tab.A1 <- power.tab.A1[order(power.tab.A1$CausalType, power.tab.A1$N), method.sel]
rownames(power.tab.A1) <-  paste0("Leffsize_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2", "ss_type5", "ls_type5"))
power.tab <- rbind(power.tab.A1, power.tab.A103)
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
power.tab$size_type <- factor(rep(c("l", "s"), each = 6), levels = c("l", "s"), labels = c("Large effect size", "Small effect size"))
power.tab$sample_size <- factor(rep(c(50, 200), 6), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(rep(c(3,6,15), each = 2), 2), levels = c(3,6,15), labels = as.character(c(3,6,15)))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("size_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating taxa",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
  facet_grid(sample_size ~ size_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt")))
p
ggsave("../Figs/power_barplot_IncreaseNumOfTaxa.pdf", plot = p, width = 20, height = 15)
