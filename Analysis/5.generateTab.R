rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Analysis/")
#### generate pseudo count table, pseudo count = 0.1, 0.5, 1, A = 0.5, B = 0.5 ####
method.sel <- c("hmp.asym", "hmp.perm")
tab.pseudo <- read.csv("../Simulation/PseudoCount/continuous/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.ori <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.ori <- tab.ori[,c(1:5,19:24)]
tab.ori$CountType <- 3
tab <- rbind(tab.pseudo, tab.ori)
tab <- tab[c(which(tab$A == 0 & tab$B == 0),
             which(tab$A == 0.5 & tab$B == 0),
             which(tab$A == 0 & tab$B == 0.5),
             which(tab$A == 0.5 & tab$B == 0.5)), c("CountType", "CausalType", "N", "A", "B", method.sel)]
tab$CountType[tab$CountType == 1] = 0.1
tab$CountType[tab$CountType == 2] = 1
tab$CountType[tab$CountType == 3] = 0.5
tab <- tab[order(tab$CausalType, tab$N, tab$A, tab$B, tab$CountType),]
write.csv(tab, file = "../Tabs/PseduoCountA5B5.csv", row.names = FALSE)

#### generate Type 1 table, A = 0.5 B = 0.5 ####
method.sel <- c("dist.omn.MedTest", "dist.bonf.MODIMA", "hmp.asym", "hmp.perm", "ldm.med.global", "ccmm.norm.tide")
tab.ori <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.con <- tab.ori[c(which(tab.ori$A == 0 & tab.ori$B == 0),
                     which(tab.ori$A == 0.5 & tab.ori$B == 0),
                     which(tab.ori$A == 0 & tab.ori$B == 0.5)), c("CausalType", "N", "A", "B", method.sel)]
tab.con$OutcomeType <- "continuous"
tab.ori <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
tab.bin <- tab.ori[c(which(tab.ori$A == 0 & tab.ori$B == 0),
                     which(tab.ori$A == 0.5 & tab.ori$B == 0),
                     which(tab.ori$A == 0 & tab.ori$B == 0.5)), c("CausalType", "N", "A", "B", method.sel)]
tab.bin$OutcomeType <- "binary"
tab <- rbind(tab.con, tab.bin)
tab <- tab[order(tab$OutcomeType, tab$N, tab$B, tab$A, tab$CausalType),]
cat(round(tab[11:20,c("hmp.asym", "hmp.perm", "dist.omn.MedTest", "dist.bonf.MODIMA", "ldm.med.global", "ccmm.norm.tide")],4) %>%
      kbl(caption= "",
          format= "latex",
          align = "c") %>%
      kable_paper(full_width = F,latex_options = "scale_down"))
write.csv(tab, file = "../Tabs/Type1A5B5.csv", row.names = FALSE)

#### generate empirical FDR table and discovery rate table ####
tab.con <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.con <- tab.con[which(tab.con$A == 0.5 & tab.con$B == 0.5), 
                   c("CausalType", "N", names(tab.con)[startsWith(names(tab.con), "efdr.sobel")],names(tab.con)[intersect(which(startsWith(names(tab.con), "efdr")), grep("perm", names(tab.con)))])]
tab.con.diffmeds <- tab.con[,c("CausalType","N",names(tab.con)[startsWith(names(tab.con), "efdr.jsmix")])]
tab.bin <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
tab.bin <- tab.bin[which(tab.bin$A == 0.5 & tab.bin$B == 0.5), 
                   c("CausalType", "N", names(tab.bin)[startsWith(names(tab.bin), "efdr.sobel")],names(tab.bin)[intersect(which(startsWith(names(tab.bin), "efdr")), grep("perm", names(tab.bin)))])]
tab.bin.diffmeds <- tab.bin[,c("CausalType","N",names(tab.bin)[startsWith(names(tab.bin), "efdr.jsmix")])]
tab.diffmeds <- rbind(tab.con.diffmeds, tab.bin.diffmeds)
cat(round(tab.diffmeds,3) %>%
      kbl(caption= "",
          format= "latex",
          align = "c") %>%
      kable_paper(full_width = F,latex_options = "scale_down"))

tab.con <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.con <- tab.con[which(tab.con$A == 0.5 & tab.con$B == 0.5), 
                   c("CausalType", "N", names(tab.con)[startsWith(names(tab.con), "dr.sobel")],names(tab.con)[intersect(which(startsWith(names(tab.con), "dr")), grep("perm", names(tab.con)))])]
tab.con.diffmeds <- tab.con[,c("CausalType","N",names(tab.con)[startsWith(names(tab.con), "dr.jsmix")])]
tab.bin <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
tab.bin <- tab.bin[which(tab.bin$A == 0.5 & tab.bin$B == 0.5), 
                   c("CausalType", "N", names(tab.bin)[startsWith(names(tab.bin), "dr.sobel")],names(tab.bin)[intersect(which(startsWith(names(tab.bin), "dr")), grep("perm", names(tab.bin)))])]
tab.bin.diffmeds <- tab.bin[,c("CausalType","N",names(tab.bin)[startsWith(names(tab.bin), "dr.jsmix")])]
tab.diffmeds <- rbind(tab.con.diffmeds, tab.bin.diffmeds)
cat(round(tab.diffmeds,3) %>%
      kbl(caption= "",
          format= "latex",
          align = "c") %>%
      kable_paper(full_width = F,latex_options = "scale_down"))

#### generate FDR table, compare different FDR control method, continuous and binary outcome ####
tab.ori <- read.csv("../Simulation/ProposedModel_Main/continuous_JC/RESULT/Type1PowerAllCont.csv", row.names = 1)
tab.con <- tab.ori[which(tab.ori$A == 0.5 & tab.ori$B == 0.5), 
                   (startsWith(colnames(tab.ori), "efdr.jsmix") + startsWith(colnames(tab.ori), "dr.jsmix")) == 1]
tab.con$OutcomeType <- "continuous"
tab.con.keep <- tab.ori[which(tab.ori$A == 0.5 & tab.ori$B == 0.5), 
                        c("CausalType", "N", "A", "B",
                          "efdr.jsmix.perm.BH","efdr.js.perm.BH", "efdr.sobel.BH",
                          "dr.jsmix.perm.BH","dr.js.perm.BH", "dr.sobel.BH") ]
tab.con.keep$OutcomeType <- "continuous"
tab.ori <- read.csv("../Simulation/ProposedModel_Main/binary_JC/RESULT/Type1PowerAllBin.csv", row.names = 1)
tab.bin <- tab.ori[which(tab.ori$A == 0.5 & tab.ori$B == 0.5), 
                   (startsWith(colnames(tab.ori), "efdr.jsmix") + startsWith(colnames(tab.ori), "dr.jsmix")) == 1]
tab.bin$OutcomeType <- "binary"
tab.bin.keep <- tab.ori[which(tab.ori$A == 0.5 & tab.ori$B == 0.5), 
                        c("CausalType", "N", "A", "B", 
                          "efdr.jsmix.perm.BH","efdr.js.perm.BH", "efdr.sobel.BH",
                          "dr.jsmix.perm.BH","dr.js.perm.BH", "dr.sobel.BH") ]
tab.bin.keep$OutcomeType <- "binary"
tab <- rbind(tab.con, tab.bin)
write.csv(tab, file = "../Tabs/EvalFDR.csv", row.names = FALSE)
tab.keep <- rbind(tab.con.keep, tab.bin.keep)
write.csv(tab.keep, file = "../Tabs/fdrNdrCompA5B5.csv", row.names = FALSE)

