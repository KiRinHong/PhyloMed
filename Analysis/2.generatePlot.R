rm(list = ls())
setwd("~/Project/QH0020/medtest/PhyloMed/Analysis/")

source("2a.utility.R")
load("../Data/Deriveddata/cecal.full.rda")
load("../Data/Deriveddata/fecal.full.rda")
load("../Data/Deriveddata/rslt.runModel.rda")

library(ggplot2)
library(ggtree)
library(ggpubr)
# box plot show Trt ~ pFat
Trt = data.cecal.full$treatment
Y = data.cecal.full$outcome
dat_tmp = data.frame(Trt =  factor(Trt, levels = c(0,1), labels = c("Control", "Antibiotics")), outcome = Y)

p <- ggplot(dat_tmp, aes(x = Trt, y = outcome, shape = Trt, color = Trt)) +
  # facet_grid(cols = vars(type)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1)) + 
  stat_compare_means(label.x = 1.3, label.y = 30, size = 8) + 
  scale_color_brewer(palette="Set2") + 
  scale_shape_manual(values=c(17,19)) + 
  theme_Publication() + 
  theme(legend.position="none") +
  labs(title = "", x = "Treatment", y = "Body Fat (%)")
p
ggsave("../SuppFigs/boxplot_pFat.pdf", plot = p, width = 10, height = 10)


# cecal full perm tree
pdf("../Figs/Fig4_cecalTreeInfo_o.pdf", width = 15, height = 15)
.plotTreeBoxScatter(data.cecal.full, rslt.cecal.full)
dev.off()
# fecal full perm tree
pdf("../SuppFigs/fecalTreeInfo_o.pdf", width = 15, height = 15)  
.plotTreeBoxScatter(data.fecal.full, rslt.fecal.full)
dev.off()

# cecal full perm qqplot
rawp.list.cecal <- list("PhyloMed"=rslt.cecal.full$node.pval["perm.jsmix",],
                        "JS"=rslt.cecal.full$node.pval["perm.js",],
                        "Sobel"=rslt.cecal.full$node.pval["sobel",])

pdf("../SuppFigs/qqplot_cecal_full.pdf", width = 10, height = 10)
.qqunifPlot(rawp.list.cecal, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
dev.off()
# fecal full perm qqplot
rawp.list.fecal <- list("PhyloMed"=rslt.fecal.full$node.pval["perm.jsmix",],
                        "JS"=rslt.fecal.full$node.pval["perm.js",],
                        "Sobel"=rslt.fecal.full$node.pval["sobel",])
pdf("../SuppFigs/qqplot_fecal_full.pdf", width = 10, height = 10)
.qqunifPlot(rawp.list.fecal, auto.key=list(corner=c(.95,.05),padding.text=4,cex = 1.8), main = "",  pch = c(19,0,2), 
            par.settings = list(superpose.symbol = list(pch = c(19,0,2), cex = 1.5, cex.title = 1.5, col = "red")))
dev.off()

##### Fig2: qqplot continuous outcome #####
# qq-plot, continuous outcome, n.sample = 200
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"],
                            "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"],
                            "CMM"=gp[,"ccmm.norm.tide"])
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample200A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType2Nsample200A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"],
                            "CMM"=gp[,"ccmm.norm.tide"])
# qq-plot, continuous outcome, n.sample = 50
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.con00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType1Nsample50A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/continuous_JC/OUTPUT/shared/Type/allType2Nsample50A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.con01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                            "PhyloMed.P"=gp[,"simes.perm"],
                            "MedTest"=gp[,"dist.jac.MedTest"], 
                            "MODIMA"=gp[,"dist.jac.MODIMA"])

setEPS()# Set postscript arguments
postscript("../Figs/Fig2a_qq_lscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/Fig2b_qq_lscon10.eps")  
.qqunifPlot(my.pvalue.list.con10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/Fig2c_qq_lscon01.eps")  
.qqunifPlot(my.pvalue.list.con01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue","orange"))))
dev.off()
postscript("../Figs/Fig2d_qq_sscon00.eps")  
.qqunifPlot(my.pvalue.list.con00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/Fig2e_qq_sscon10.eps")  
.qqunifPlot(my.pvalue.list.con10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../Figs/Fig2f_qq_sscon01.eps")  
.qqunifPlot(my.pvalue.list.con01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex = 1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()

##### FigS2: qqplot binary outcome #####
# qq-plot, binary outcome, n.sample = 200
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample200A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType2Nsample200A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample200A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType2Nsample200A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.l <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=gp[,"dist.jac.MODIMA"])
# qq-plot, binary outcome, n.sample = 50
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B0.Rdata")
gp <- rslt_all$gp
my.pvalue.list.bin00.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample50A1B0.Rdata")
gp <- rslt_all$gp
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType2Nsample50A1B0.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin10.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=gp[,"dist.jac.MODIMA"])
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType1Nsample50A0B1.Rdata")
gp <- rslt_all$gp
load("../Simulation/binary_JC/OUTPUT/shared/Type/allType2Nsample50A0B1.Rdata")
gp <- rbind(gp, rslt_all$gp)
my.pvalue.list.bin01.s <-list("PhyloMed.A"=gp[,"simes.asym"],
                              "PhyloMed.P"=gp[,"simes.perm"],
                              "MedTest"=gp[,"dist.jac.MedTest"], 
                              "MODIMA"=na.omit(gp[,"dist.jac.MODIMA"])) # 26 replicates produce NA

setEPS()# Set postscript arguments
postscript("../SuppFigs/FigS1a_qq_lsbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../SuppFigs/FigS1b_qq_lsbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../SuppFigs/FigS1c_qq_lsbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.l, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../SuppFigs/FigS1d_qq_ssbin00.eps")  
.qqunifPlot(my.pvalue.list.bin00.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2, 
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../SuppFigs/FigS1e_qq_ssbin10.eps")  
.qqunifPlot(my.pvalue.list.bin10.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()
postscript("../SuppFigs/FigS1f_qq_ssbin01.eps")  
.qqunifPlot(my.pvalue.list.bin01.s, auto.key=list(corner=c(.95,.05),padding.text=4, cex=1.5), main = "",
            par.settings = list(superpose.symbol = list(pch = 19, cex = 1.2,
                                                        col = c("pink","red","green","blue"))))
dev.off()


# power bar plot, continuous outcome and binary outcome
method.sel <- c("dist.jac.MedTest", "dist.jac.MODIMA", "simes.asym", "simes.perm", "ccmm.norm.tide")
power.con.tab <- read.csv("../Simulation/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
power.con.tab <- power.con.tab[order(power.con.tab$CausalType, power.con.tab$N),]
power.con.tab <- power.con.tab[which(power.con.tab$A == 1 & power.con.tab$B == 1), c("CausalType", "N", method.sel)]
rownames(power.con.tab) <-  paste0("con_",c("ss_type1", "ls_type1", "ss_type2", "ls_type2"))
power.bin.tab <- read.csv("../Simulation/binary_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
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
ggsave("../Figs/Fig3_power_barplot.pdf", plot = p, width = 20, height = 15)


