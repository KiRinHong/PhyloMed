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
load("../Data/Deriveddata/rslt.runModel.rda")
##### RealData: Cecal: Trt ~ pFat #####
load("../Data/Deriveddata/Cecal.filter20top100.rda")
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
ggsave("../Figs/S5.Cecal.boxplot.antibiotics_pFat.pdf", plot = p, width = 15, height = 15)

pdf("../Figs/4a.Cecal.tree.pdf", width = 20, height = 10)
.plotTree(cecal.top, rslt.cecal.top100.psc05, "rectangular")
dev.off()
pdf("../Figs/4d.Cecal.subtree.pdf", width = 15, height = 15)
.plotSubTree(cecal.top, rslt.cecal.top100.psc05)
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

pdf("../Figs/4c.Cecal.boxscatter.pdf", width = 20, height = 10)
.plotBoxScatter(cecal.top, rslt.cecal.top100.psc05, "Treatment", NULL, "pFat")
dev.off()

rawp.list.cecal <- list("PhyloMed"=rslt.cecal.top100.psc05$node.pval["perm.jsmix.prod",],
                        "JS"=rslt.cecal.top100.psc05$node.pval["perm.js",],
                        "Sobel"=rslt.cecal.top100.psc05$node.pval["sobel",])
pdf("../Figs/4b.Cecal.qqplot.pdf", width = 10, height = 10)
.qqunifPlot(rawp.list.cecal, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
dev.off()

##### RealData: COMBO: fat ~ bmi, adjusted by calor #####
load("../Data/Deriveddata/COMBO.filter20.rda")
dat_tmp = data.frame(Trt = combo.filter$meta[["fat"]], 
                     outcome = combo.filter$meta[["bmi"]])

p <- ggplot(dat_tmp, aes(x = Trt, y = outcome)) +
  geom_point(size = 4)+
  stat_smooth(method = "lm",
              col = "red",
              se = FALSE,
              linewidth = 1) +
  stat_cor(method = "pearson", label.x = -2, label.y = 40, size = 10, ) + 
  theme_Publication() + 
  theme(legend.position="none") +
  labs(title = "", x = "Fat intake", y = "BMI")
p
ggsave("../Figs/S6.COMBO.scatterplot.fat_bmi.pdf", plot = p, width = 15, height = 15)

pdf("../Figs/S7.COMBO.tree.pdf", width = 15, height = 15)
.plotTree(combo.filter, rslt.combo.filter20.psc05, "circular")
dev.off()
pdf("../Figs/5.COMBO.scatterpie.pdf", width = 15, height = 15)
.plotScatterPie(combo.filter, rslt.combo.filter20.psc05, "fat", "calor", "bmi")
dev.off()

# rawp.list.combo <- list("PhyloMed"=rslt.combo.filter20.psc05$node.pval["perm.jsmix.prod",],
#                         "JS"=rslt.combo.filter20.psc05$node.pval["perm.js",],
#                         "Sobel"=rslt.combo.filter20.psc05$node.pval["sobel",])
# 
# pdf("../Figs/COMBO.qqplot.pdf", width = 10, height = 10)
# .qqunifPlot(rawp.list.combo, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
#             par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
# dev.off()

##### Simulation: qqplot continuous outcome #####
# qq-plot, continuous outcome, n.sample = 200
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"],
                              "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0.5B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A0.5B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"],
                              "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.5.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A0B0.5.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"],
                              "CMM"=gp[,"ccmm.norm.tide"])
# qq-plot, continuous outcome, n.sample = 50
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0.5B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A0.5B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.5.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A0B0.5.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])

setEPS()# Set postscript arguments
postscript("../Figs/2.qq_lscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey","orange"))))
dev.off()
postscript("../Figs/2.qq_lscon10.eps")  
.qqunifPlot(my.pvalue.list.con10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey","orange"))))
dev.off()
postscript("../Figs/2.qq_lscon01.eps")
.qqunifPlot(my.pvalue.list.con01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey","orange"))))
dev.off()
postscript("../Figs/2.qq_sscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/2.qq_sscon10.eps")
.qqunifPlot(my.pvalue.list.con10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/2.qq_sscon01.eps")  
.qqunifPlot(my.pvalue.list.con01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()

##### Simulation: qqplot binary outcome #####
# qq-plot, binary outcome, n.sample = 200
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0.5B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample200A0.5B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.5.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample200A0B0.5.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
# qq-plot, binary outcome, n.sample = 50
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0.5B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample50A0.5B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1)),
                              "LDM-med"=gp[,"ldm.med.global"])
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.5.Rdata")
gp <- rslt_all$gp
load("../Simulation/ProposedModel_Main/binary_JC/OUTPUT/shared/Type/allType2Nsample50A0B0.5.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.omn.MedTest"], 
                              "MODIMA"=rowMins(cbind(gp[,which(endsWith(colnames(gp), "MODIMA"))]*5,1), na.rm = T),
                              "LDM-med"=na.omit(gp[,"ldm.med.global"]))

setEPS()# Set postscript arguments
postscript("../Figs/S1.qq_lsbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/S1.qq_lsbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/S1.qq_lsbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/S1.qq_ssbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/S1.qq_ssbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()
postscript("../Figs/S1.qq_ssbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","grey"))))
dev.off()


##### Simulation: power bar plot, continuous and binary outcome, top100 taxa, A=0.5 B=0.5####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "hmp.asym", "hmp.perm", "ldm.med.global", "ccmm.norm.tide")
power.con.tab <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
power.con.tab <- power.con.tab[which(power.con.tab$A == 0.5 & power.con.tab$B == 0.5), c("CausalType", "N", method.sel)]
rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.bin.tab <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
power.bin.tab <- power.bin.tab[order(power.bin.tab$CausalType, power.bin.tab$N),]
power.bin.tab <- power.bin.tab[which(power.bin.tab$A == 0.5 & power.bin.tab$B == 0.5), c("CausalType", "N", method.sel)]
rownames(power.bin.tab) <-  paste0("bin_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.tab <- rbind(power.con.tab, power.bin.tab)
power.tab <- power.tab[,method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "LDM-med", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "LDM-med", "CMM")]
power.tab$outcome_type <- factor(rep(c("con", "bin"), each = 4), levels = c("con", "bin"), labels = c("Continuous outcome", "Binary outcome"))
power.tab$sample_size <- factor(rep(c(50, 200), 4), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(rep(c(3, 6), each = 2), 2), levels = c(3, 6), labels = c("3", "6"))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("outcome_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating OTUs",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "gray", "orange")) + 
  facet_grid(sample_size ~ outcome_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt"))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
p
ggsave("../Figs/3.power_barplot_ClusterdTaxa.pdf", plot = p, width = 20, height = 15)

##### Simulation: power bar plot, continuous and binary outcome, all taxa, A=0.5 B=0.5 #####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "hmp.asym", "hmp.perm", "ldm.med.global")
power.con.tab <- read.csv("../Simulation/RareTaxa/continuous/RESULT/PowerAllCont.csv", row.names = 1)
power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
power.con.tab <- power.con.tab[which(power.con.tab$A == 0.5 & power.con.tab$B == 0.5), c("CausalType", "N", method.sel)]
rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.bin.tab <- read.csv("../Simulation/RareTaxa/binary/RESULT/PowerAllBin.csv", row.names = 1)
power.bin.tab <- power.bin.tab[order(power.bin.tab$CausalType, power.bin.tab$N),]
power.bin.tab <- power.bin.tab[which(power.bin.tab$A == 0.5 & power.bin.tab$B == 0.5), c("CausalType", "N", method.sel)]
rownames(power.bin.tab) <-  paste0("bin_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.tab <- rbind(power.con.tab, power.bin.tab)
power.tab <- power.tab[,method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "LDM-med")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "LDM-med")]
power.tab$outcome_type <- factor(rep(c("con", "bin"), each = 4), levels = c("con", "bin"), labels = c("Continuous outcome", "Binary outcome"))
power.tab$sample_size <- factor(rep(c(50, 200), 4), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(rep(c(3, 6), each = 2), 2), levels = c(3, 6), labels = c("3", "6"))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("outcome_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating OTUs",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "gray")) + 
  facet_grid(sample_size ~ outcome_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt"))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
p
ggsave("../Figs/S2.power_barplot_allTaxa.pdf", plot = p, width = 20, height = 15)


# #### Archive: power bar plot, supp, continuous outcome, N=200 and small effect size ####
# method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
# power.tab <- read.csv("../Simulation/SmallEffectSize/continuous/RESULT/Type1PowerAllCont.csv", row.names = 1)
# power.tab <- power.tab[order(power.tab$CausalType),]
# power.tab <- power.tab[which(power.tab$A == 0.001 & power.tab$B == 0.5), method.sel]
# colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
# power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
# power.tab$numberoftaxa <- factor(c(3,6), levels = c(3,6), labels = c("3", "6"))
# power.tab.long <- reshape2::melt(power.tab, id.vars = "numberoftaxa")
# p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
#   geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
#   labs(title = "",
#        subtitle = "",
#        x = "Number of mediating OTUs",
#        y = "Power") +
#   scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
#   theme_Publication(base_size = 20) +
#   theme(#legend.spacing.x = unit(0.5, 'cm'),
#     legend.text = element_text(margin = margin(r = 15, unit = "pt")))
# p
# ggsave("../Figs/power_barplot_Con200SmlEffSize.pdf", plot = p, width = 20, height = 15)
# 
# #### Archive: power bar plot, supp, continuous and binary outcome, partially overlapped ####
# method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
# power.con.tab <- read.csv("../Simulation/PartialOverlap/continuous/RESULT/PowerAllCont.csv", row.names = 1)
# power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
# power.con.tab <- power.con.tab[which(power.con.tab$A == 1 & power.con.tab$B == 1), c("CausalType", "N", method.sel)]
# power.con.tab[c(1,3),"ccmm.norm.tide"] <- NA
# rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
# power.bin.tab <- read.csv("../Simulation/PartialOverlap/binary/RESULT/PowerAllCont.csv", row.names = 1)
# power.bin.tab <- power.bin.tab[order(power.bin.tab$CausalType, power.bin.tab$N),]
# power.bin.tab <- power.bin.tab[which(power.bin.tab$A == 1 & power.bin.tab$B == 1), c("CausalType", "N", method.sel)]
# rownames(power.bin.tab) <-  paste0("bin_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
# power.tab <- rbind(power.con.tab, power.bin.tab)
# power.tab <- power.tab[,method.sel]
# colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "CMM")
# power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "CMM")]
# power.tab$outcome_type <- factor(rep(c("con", "bin"), each = 4), levels = c("con", "bin"), labels = c("Continuous outcome", "Binary outcome"))
# power.tab$sample_size <- factor(rep(c(50, 200), 4), levels = c(50, 200), labels = c("n = 50", "n = 200"))
# power.tab$numberoftaxa <- factor(rep(rep(c(1, 2), each = 2), 2), levels = c(1, 2), labels = c("1", "2"))
# power.tab.long <- reshape2::melt(power.tab, id.vars = c("outcome_type", "sample_size", "numberoftaxa"))
# p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
#   geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
#   labs(title = "",
#        subtitle = "",
#        x = "Number of mediating OTUs",
#        y = "Power") +
#   scale_fill_manual(values = c("pink", "red", "green", "blue", "orange")) + 
#   facet_grid(sample_size ~ outcome_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
#   theme_Publication(base_size = 20) +
#   theme(#legend.spacing.x = unit(0.5, 'cm'),
#     legend.text = element_text(margin = margin(r = 15, unit = "pt")))
# p
# ggsave("../Figs/power_barplot_PartialOverlap.pdf", plot = p, width = 20, height = 15)

##### Simulation: power bar plot, continuous outcome, N=50/200, increase number of mediating taxa 3/6/15 ####
# randomly select and fully overlapped, fix the effect size A=1e-3 or 0.5 B=0.5
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "hmp.asym", "hmp.perm", "ldm.med.global", "ccmm.norm.tide")
power.tab <- read.csv("../Simulation/IncreaseNumOfCausalTaxa/continuous/RESULT/PowerAllCont.csv", row.names = 1)
power.tab <- power.tab[order(power.tab$CausalType),]
power.tab <- power.tab[which(power.tab$A %in% c(0.001, 0.5)), method.sel]
colnames(power.tab) <- c("MedTest", "MODIMA", "PhyloMed.A", "PhyloMed.P", "LDM-med", "CMM")
power.tab <- power.tab[,c("PhyloMed.A", "PhyloMed.P", "MedTest", "MODIMA", "LDM-med", "CMM")]
power.tab$size_type <- factor(rep(c("s", "l"), 6), levels = c("l", "s"), labels = c("Large mediation effect", "Small mediation effect"))
power.tab$sample_size <- factor(rep(rep(c(50, 200), each = 2), 3), levels = c(50, 200), labels = c("n = 50", "n = 200"))
power.tab$numberoftaxa <- factor(rep(c(3,6,15), each = 4), levels = c(3,6,15), labels = as.character(c(3,6,15)))
power.tab.long <- reshape2::melt(power.tab, id.vars = c("size_type", "sample_size", "numberoftaxa"))
p <- ggplot(na.omit(power.tab.long), aes(x=numberoftaxa, y=value)) +
  geom_bar(stat = "identity", aes(fill = variable), colour = "black", position = "dodge", width = 0.75) +
  labs(title = "",
       subtitle = "",
       x = "Number of mediating OTUs",
       y = "Power") +
  scale_fill_manual(values = c("pink", "red", "green", "blue", "gray", "orange")) + 
  facet_grid(sample_size ~ size_type) + # labeller = labeller(include_x = include_x.labs) scale = "free_y"
  theme_Publication(base_size = 20) +
  theme(#legend.spacing.x = unit(0.5, 'cm'),
    legend.text = element_text(margin = margin(r = 15, unit = "pt"))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
p
ggsave("../Figs/S3.power_barplot_ScatteredTaxa.pdf", plot = p, width = 20, height = 15)


##### Simulation: generate empirical FDR curve, N=200, effect size increase from A=0,1e-5,1e-4,1e-3,0.01,0.1,1 B=0.5 ####
# mediating taxa clustered 3/6/15
tmp <- read.csv("../Simulation/IncreaseEffectSize/continuous/RESULT/Type1PowerAllCont.csv", row.names = 1)
tmp <- tmp[tmp$B == 0.5,]
tmp$group <- paste0(tmp$CausalType, ".", tmp$B)
tmp$CausalType[tmp$CausalType == 3] <- 15
tmp$CausalType[tmp$CausalType == 2] <- 6
tmp$CausalType[tmp$CausalType == 1] <- 3
tmp$A[tmp$A == 0] <- 1e-7
p <- ggplot(tmp, aes(x=A, y=efdr.jsmix.perm.BH*100, group=group, color = factor(CausalType))) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(breaks = c("3", "6", "15"),values = c("red", "blue", "green")) +
  #scale_shape_manual(breaks = c("0.5", "1", "2", "4"),values = 1:4) + 
  scale_x_log10(breaks = c(1e-7,1e-5,1e-4,1e-3,0.01,0.1,1),
                labels= c(-Inf,-5:0)) +
  #scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_y_continuous(limits = c(0,5)) +
  theme_Publication(base_size = 20) + 
  theme(legend.position = 'bottom', 
        legend.title = element_text(size=rel(1.2)),
        legend.box = 'vertical') + 
  guides(#shape = guide_legend(title = "Effect size A", byrow = T, nrow = 1),
         color = guide_legend(title = "Number of mediating OTUs", byrow = T, nrow = 1)) +
  #labs(title = "", x = "Effect Size B", y = "Empirical FDR (%)")
  labs(title = "", x = bquote(log[10](A)), y = "Empirical FDR (%)")
p
ggsave("../Figs/S4.lineplot_efdr_fixB_log10.pdf", plot = p, width = 10, height = 10)

# #### Archive: generate empirical pFDR plot, N=200, effect size increase from A=1e-3,0.01,0.1,1 B=0.5####
# tmp <- read.csv("../Simulation/IncreaseEffectSize/continuous/RESULT/Type1PowerAllCont_pFDR.csv", row.names = 1)
# #tmp <- tmp[!(tmp$A %in% c(1e-5, 1e-4)),]
# tmp$group <- paste0(tmp$CausalType, ".", tmp$B)
# tmp$CausalType[tmp$CausalType == 3] <- 15
# tmp$CausalType[tmp$CausalType == 2] <- 6
# tmp$CausalType[tmp$CausalType == 1] <- 3
# p <- ggplot(tmp, aes(x=A, y=Jsmix*100, group=group, color = factor(CausalType))) + 
#   geom_point() + 
#   geom_line() + 
#   scale_color_manual(breaks = c("3", "6", "15"),values = c("red", "blue", "green")) +
#   #scale_shape_manual(breaks = c("0.5", "1", "2", "4"),values = 1:4) + 
#   scale_x_log10(breaks = c(1e-5,1e-4,1e-3, 0.01, 0.1, 1), 
#                 labels= c(-5,-4,-3, -2, -1, 0)) +
#   scale_y_continuous(limits = c(0,11)) +
#   theme_Publication(base_size = 20) + 
#   theme(legend.position = 'bottom', 
#         legend.title = element_text(size=rel(1.2)),
#         legend.box = 'vertical') + 
#   guides(#shape = guide_legend(title = "Effect size A", byrow = T, nrow = 1),
#     color = guide_legend(title = "Number of mediating OTUs", byrow = T, nrow = 1)) +
#   labs(title = "", x = bquote(log[10](A)), y = "Empirical FDR (%)")
# p
# ggsave("../Figs/lineplot_efdr_pFDR.pdf", plot = p, width = 10, height = 10)

# ##### Archive: generate qqplot of null taxa when exists causal taxa (non-null taxa) ####
# # continuous outcome, n = 200, causalType =1
# load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0.5B0.5.Rdata")
# rawp.asym <- rslt_all$rawp.asym.jsmix
# rawp.asym.null <- rawp.asym[which(rslt_all$causalNode == 0)]
# rawp.asym.nonnull <- rawp.asym[which(rslt_all$causalNode == 1)]
# rawp.perm <- rslt_all$rawp.perm.jsmix
# rawp.perm.null <- rawp.perm[which(rslt_all$causalNode == 0)]
# rawp.perm.nonnull <- rawp.perm[which(rslt_all$causalNode == 1)]
# my.pvalue.list.con11.l <-list("PhyloMed.A.null"=na.omit(rawp.asym.null),
#                               #"PhyloMed.A.nonnull"=na.omit(rawp.asym.nonnull),
#                               "PhyloMed.P.null"=na.omit(rawp.perm.null))
#                               #"PhyloMed.P.nonnull"=na.omit(rawp.perm.nonnull))
# postscript("../Figs/qq_nullTaxaCausalType1.eps")  
# .qqunifPlot(my.pvalue.list.con11.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
#             par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
#                                                         col = c("red","green"))))
# dev.off()
# 
# rawp.asym <- rslt_all$rawp.asym.jsmix
# rawp.perm <- rslt_all$rawp.perm.jsmix
# globalp.mat = matrix(NA, nrow = nrow(rslt_all$gp), 6)
# for (i in 1:nrow(rslt_all$gp)) {
#   rawp.asym.rm = na.omit(rawp.asym[i,which(rslt_all$causalNode[i,] == 0)])
#   rawp.perm.rm = na.omit(rawp.perm[i,which(rslt_all$causalNode[i,] == 0)])
#   L = length(rawp.asym.rm)
#   globalp.mat[i,] = c(min(L * rawp.asym.rm/rank(rawp.asym.rm)), 
#                       1 - pchisq(-2 * sum(log(rawp.asym.rm)), df = 2 * L), 
#                       p.hmp(rawp.asym.rm, w = rep(1/L, L), L = L),
#                       min(L * rawp.perm.rm/rank(rawp.perm.rm)), 
#                       1 - pchisq(-2 * sum(log(rawp.perm.rm)), df = 2 * L), 
#                       p.hmp(rawp.perm.rm, w = rep(1/L, L), L = L))
# }
# my.pvalue.list.con11.l <-list("PhyloMed.A"=globalp.mat[,3],
#                               "PhyloMed.P"=globalp.mat[,6])
# postscript("../Figs/qq_nullTaxaCausalType1_globalp.eps")  
# .qqunifPlot(my.pvalue.list.con11.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
#             par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
#                                                         col = c("pink","red"))))
# dev.off()
# 
# # continuous outcome, n =200, causalType =2
# load("../Simulation/ProposedModel_Main/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A0.5B0.5.Rdata")
# rawp.asym <- rslt_all$rawp.asym.jsmix
# rawp.asym.null <- rawp.asym[which(rslt_all$causalNode == 0)]
# rawp.asym.nonnull <- rawp.asym[which(rslt_all$causalNode == 1)]
# rawp.perm <- rslt_all$rawp.perm.jsmix
# rawp.perm.null <- rawp.perm[which(rslt_all$causalNode == 0)]
# rawp.perm.nonnull <- rawp.perm[which(rslt_all$causalNode == 1)]
# my.pvalue.list.con11.l <-list("PhyloMed.A.null"=na.omit(rawp.asym.null),
#                               #"PhyloMed.A.nonnull"=na.omit(rawp.asym.nonnull),
#                               "PhyloMed.P.null"=na.omit(rawp.perm.null))
# #"PhyloMed.P.nonnull"=na.omit(rawp.perm.nonnull))
# postscript("../Figs/qq_nullTaxaCausalType2.eps")  
# .qqunifPlot(my.pvalue.list.con11.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
#             par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
#                                                         col = c("red","green"))))
# dev.off()
# 
# rawp.asym <- rslt_all$rawp.asym.jsmix
# rawp.perm <- rslt_all$rawp.perm.jsmix
# globalp.mat = matrix(NA, nrow = nrow(rslt_all$gp), 6)
# for (i in 1:nrow(rslt_all$gp)) {
#   rawp.asym.rm = na.omit(rawp.asym[i,which(rslt_all$causalNode[i,] == 0)])
#   rawp.perm.rm = na.omit(rawp.perm[i,which(rslt_all$causalNode[i,] == 0)])
#   L = length(rawp.asym.rm)
#   globalp.mat[i,] = c(min(L * rawp.asym.rm/rank(rawp.asym.rm)), 
#                       1 - pchisq(-2 * sum(log(rawp.asym.rm)), df = 2 * L), 
#                       p.hmp(rawp.asym.rm, w = rep(1/L, L), L = L),
#                       min(L * rawp.perm.rm/rank(rawp.perm.rm)), 
#                       1 - pchisq(-2 * sum(log(rawp.perm.rm)), df = 2 * L), 
#                       p.hmp(rawp.perm.rm, w = rep(1/L, L), L = L))
# }
# my.pvalue.list.con11.l <-list("PhyloMed.A"=globalp.mat[,3],
#                               "PhyloMed.P"=globalp.mat[,6])
# postscript("../Figs/qq_nullTaxaCausalType2_globalp.eps")  
# .qqunifPlot(my.pvalue.list.con11.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
#             par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
#                                                         col = c("pink","red"))))
# dev.off()

##### RealData: Sensitivity Analysis ####
load("../Data/Deriveddata/rslt.sens90.rda")
pdf(file = "../Figs/S8.sens90.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plot(Cecal.Sens$NodeID75$sens.out, sens.par = "rho", main = "Mouse cecal study", ylim = c(-2, 2), cex.lab = 1.5, cex.main = 2)
plot(Combo.Sens$NodeID297$sens.out, sens.par = "rho", main = "Human gut study", ylim = c(-2, 2), cex.lab = 1.5, cex.main = 2)
dev.off()
